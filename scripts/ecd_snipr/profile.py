"""Deterministic per-candidate molecular profile, boundary-precision analysis,
multichain context and multi-pass loop enumeration.

Everything here is computed from annotations already present in the normalized
UniProt record (see uniprot.normalize). No structure prediction, no scoring:
each result is a named, inspectable record that feeds explicit reason codes in
the screening decision table. Length/cysteine/glycosylation outcomes remain
review warnings, never pass/fail thresholds.
"""

import re
from .common import contains, interval, overlaps

# Informational annotation kinds: parsed and profiled, but never candidate-
# defining. A coordinate problem in one of these is deferred/recorded, never a
# blocking coordinate error (blocking is reserved for candidate-relevant kinds).
INFORMATIONAL_KINDS = {"region", "motif", "glycosylation_site", "site", "binding_site",
                       "active_site", "variant", "mutagenesis", "lipidation"}

# Kinds whose coordinates define or block candidates; prediction-grade evidence
# on these is what prediction_requires_annotation_review is about.
BOUNDARY_DEFINING_KINDS = {"extracellular", "cytoplasmic", "lumenal", "other_topology",
                           "transmembrane", "signal_peptide", "chain", "propeptide",
                           "processed_peptide", "gpi_signal", "gpi_attachment_site"}
PREDICTION_RELEVANT_KINDS = BOUNDARY_DEFINING_KINDS | {"domain", "disulfide"}

# Anchor kinds considered when checking that a candidate boundary agrees with
# the neighboring annotation (e.g. ECD start should abut the signal peptide).
BOUNDARY_ANCHOR_KINDS = ("signal_peptide", "transmembrane", "propeptide", "gpi_signal", "chain")

SEQUON_RE = re.compile(r"(?=(N[^P][ST]))")

# Conservative hetero-oligomer assembly detector for SUBUNIT comment text.
# Matches explicit "heterodimer/heterotrimer/heterotetramer/hetero-oligomer..."
# wording only; "interacts with X" or "part of a complex" phrasing is recorded
# in the raw comment but does NOT fire the conditional code (documented limit).
HETERO_ASSEMBLY_RE = re.compile(r"[Hh]etero[- ]?(?:di|tri|tetra|oligo|multi)mer")


def n_glyco_sequons(fragment, start):
    """N-X-S/T (X != P) sequon positions, 1-based on the reference sequence."""
    return [start + m.start() for m in SEQUON_RE.finditer(fragment)]


def _valid_bounds(feature, seq_len):
    try:
        return interval(feature, seq_len)
    except ValueError:
        return None


def molecular_profile(bounds, fragment, features, seq_len):
    """Per-candidate molecular profile from data already in the record.

    All positions are 1-based inclusive on the reference sequence. Features with
    fuzzy/invalid coordinates are skipped here (they are deferred explicitly by
    the screening partition step, never silently used).
    """
    start, end = bounds
    profile = {"length": len(fragment)}
    cysteines = [start + i for i, aa in enumerate(fragment) if aa == "C"]
    profile["cysteine_count"] = len(cysteines)
    profile["cysteine_odd"] = bool(cysteines and len(cysteines) % 2)
    profile["cysteine_positions"] = cysteines
    disulfides_full, disulfides_partial = [], []
    for f in features:
        if f.get("kind") != "disulfide":
            continue
        b = _valid_bounds(f, seq_len)
        if not b:
            continue
        entry = {"name": f.get("name", ""), "start": b[0], "end": b[1],
                 "evidence_kind": f.get("evidence", {}).get("kind", "")}
        in_start, in_end = start <= b[0] <= end, start <= b[1] <= end
        if in_start and in_end:
            disulfides_full.append(entry)
        elif in_start != in_end:
            disulfides_partial.append(dict(entry, retained_position=b[0] if in_start else b[1]))
    profile["disulfides_fully_contained"] = disulfides_full
    profile["disulfides_partial"] = disulfides_partial
    sequon_positions = n_glyco_sequons(fragment, start)
    profile["n_glyco_sequons"] = {"count": len(sequon_positions), "positions": sequon_positions}
    glyco = []
    for f in features:
        if f.get("kind") != "glycosylation_site":
            continue
        b = _valid_bounds(f, seq_len)
        if not b or not start <= b[0] <= end:
            continue
        desc = f.get("name", "")
        lowered = desc.lower()
        subtype = "N-linked" if lowered.startswith("n-linked") else \
            "O-linked" if lowered.startswith("o-linked") else "other_or_unspecified"
        glyco.append({"position": b[0], "name": desc, "subtype": subtype,
                      "evidence_kind": f.get("evidence", {}).get("kind", "")})
    profile["annotated_glycosylation_in_fragment"] = glyco
    domains_full, domains_cut, domains_outside = [], [], []
    for f in features:
        if f.get("kind") != "domain":
            continue
        b = _valid_bounds(f, seq_len)
        if not b:
            continue
        entry = {"name": f.get("name", ""), "start": b[0], "end": b[1]}
        if contains(bounds, b):
            domains_full.append(entry)
        elif overlaps(bounds, b):
            sides = []
            if b[0] < start:
                sides.append("n_side")
            if b[1] > end:
                sides.append("c_side")
            domains_cut.append(dict(entry, cut_side="+".join(sides) or "interior",
                                    retained_interval=[max(b[0], start), min(b[1], end)]))
        else:
            domains_outside.append(entry)
    profile["domains_fully_contained"] = domains_full
    profile["domains_cut"] = domains_cut
    profile["domains_outside_fragment"] = domains_outside
    for kind, key in (("active_site", "active_sites"), ("binding_site", "binding_sites"),
                      ("site", "sites"), ("lipidation", "lipidation_sites")):
        inside, outside = [], []
        for f in features:
            if f.get("kind") != kind:
                continue
            b = _valid_bounds(f, seq_len)
            if not b:
                continue
            entry = {"name": f.get("name", ""), "start": b[0], "end": b[1],
                     "evidence_kind": f.get("evidence", {}).get("kind", "")}
            (inside if overlaps(bounds, b) else outside).append(entry)
        profile[key + "_in_fragment"] = inside
        profile[key + "_outside_fragment"] = outside
    for kind, key in (("variant", "variants"), ("mutagenesis", "mutagenesis_records")):
        hits = []
        for f in features:
            if f.get("kind") != kind:
                continue
            b = _valid_bounds(f, seq_len)
            if not b or not overlaps(bounds, b):
                continue
            hits.append({"name": f.get("name", ""), "start": b[0], "end": b[1],
                         "feature_id": f.get("feature_id", "")})
        profile[key + "_overlapping"] = hits
    return profile


def _side_analysis(position, direction, features, seq_len, tolerance):
    """Check one candidate boundary against adjacent anchor annotations.

    direction "n": the boundary is `position` = candidate start; an anchor
    feature should end at start-1. direction "c": boundary is candidate end; an
    anchor should start at end+1. A mature chain that starts/ends exactly at
    the candidate boundary is agreement, never a conflict; chains spanning
    across the boundary are not adjacency evidence. Returns agreement status
    plus the evidence: exact / gap_tolerated / gap_beyond_tolerance / overlap /
    missing / chain_terminus.
    """
    best = None
    for f in features:
        if f.get("kind") not in BOUNDARY_ANCHOR_KINDS:
            continue
        b = _valid_bounds(f, seq_len)
        if not b:
            continue
        if f.get("kind") == "chain":
            coincident = (b[0] == position) if direction == "n" else (b[1] == position)
            if not coincident:
                continue  # a chain spanning across the boundary is not adjacency evidence
            gap = 0
        elif direction == "n":
            if b[0] > position:
                continue  # feature lies entirely after the boundary: not an n-side anchor
            gap = position - b[1] - 1  # anchor ends before candidate start
        else:
            if b[1] < position:
                continue  # feature lies entirely before the boundary: not a c-side anchor
            gap = b[0] - position - 1  # anchor starts after candidate end
        candidate = {"kind": f.get("kind"), "name": f.get("name", ""), "start": b[0], "end": b[1],
                     "gap": gap, "evidence_kind": f.get("evidence", {}).get("kind", "")}
        if best is None or abs(gap) < abs(best["gap"]):
            best = candidate
    if direction == "n" and position == 1:
        return {"status": "chain_terminus", "adjacent_feature": None}
    if direction == "c" and position == seq_len:
        return {"status": "chain_terminus", "adjacent_feature": None}
    if best is None:
        return {"status": "missing", "adjacent_feature": None}
    if best["gap"] == 0:
        status = "exact"
    elif best["gap"] > 0:
        status = "gap_tolerated" if best["gap"] <= tolerance else "gap_beyond_tolerance"
    else:
        status = "overlap"
    return {"status": status, "adjacent_feature": best}


def boundary_analysis(bounds, region, features, seq_len, topology, tolerance=5):
    """Which feature justifies each candidate boundary, and does the adjacent
    SIGNAL/TM/PROPEP annotation agree?

    The defining feature is the region the candidate derives from. Adjacent
    disagreement is recorded with explicit reason codes:
    - gap within tolerance -> tm_adjacency_tolerance_used (TM anchor) or
      boundary_adjacency_tolerance_used (other anchor): informational, the
      candidate boundary is used as annotated.
    - overlap or gap beyond tolerance -> boundary_feature_conflict: a named
      verification item (conditional in screening). Note a genuine overlap with
      TM/SP normally already blocks via retained_* checks before this point.
    """
    start, end = bounds
    result = {
        "defining_feature": {"kind": region.get("kind"), "name": region.get("name", ""),
                             "start": region.get("start"), "end": region.get("end"),
                             "boundary_status": region.get("boundary_status", "exact"),
                             "evidence": region.get("evidence", {})},
        "tolerance": tolerance,
        "reason_codes": [],
    }
    # Domain-fragment candidates are defined by the domain itself; adjacency to
    # SP/TM is not the justification and is not forced.
    adjacency_applicable = region.get("kind") in {"extracellular", "chain", "processed_peptide"}
    for direction, position in (("n_side", start), ("c_side", end)):
        if not adjacency_applicable:
            side = {"status": "not_applicable",
                    "note": "boundary defined by the annotated domain/region itself"}
        else:
            side = _side_analysis(position, "n" if direction == "n_side" else "c",
                                  features, seq_len, tolerance)
        status = side["status"]
        anchor = (side.get("adjacent_feature") or {}).get("kind", "")
        code = None
        if status == "gap_tolerated":
            code = "tm_adjacency_tolerance_used" if anchor == "transmembrane" else \
                "boundary_adjacency_tolerance_used"
        elif status in {"gap_beyond_tolerance", "overlap"}:
            code = "boundary_feature_conflict"
        side["reason_code"] = code
        result[direction] = side
    result["reason_codes"] = sorted({result[s]["reason_code"] for s in ("n_side", "c_side")}
                                    - {None})
    return result


def detect_multichain_partners(subunit_comments):
    """Explicit hetero-oligomer assembly sentences from SUBUNIT comments.

    Returns matched sentences verbatim as native context, not proof that an
    isolated antigen fragment needs the partner. Obvious negation is withheld;
    this lexical triage is not a complete natural-language entailment model.
    Only fires on explicit hetero-*mer wording; homo-oligomers, bare
    "interacts with" and "part of a complex" phrasing stay unflagged.
    """
    matches = []
    for comment in subunit_comments or []:
        text = comment.get("text", "") if isinstance(comment, dict) else str(comment)
        for sentence in re.split(r"(?<=[.!?])\s+", text):
            if HETERO_ASSEMBLY_RE.search(sentence) and not re.search(r"\b(?:not|no|neither|never|unable|cannot|fails?\s+to)\b", sentence, re.I):
                cleaned = sentence.strip()
                if cleaned and cleaned not in matches:
                    matches.append(cleaned)
    return matches


def enumerate_external_loops(protein, features, seq_len):
    """Enumerate annotated extracellular loops individually (multi-pass depth).

    Each loop reports its interval, length, the sourced domains it fully
    contains and whether such a domain is eligible as a conditional alternate.
    Loops are never stitched; the enumeration is recorded even when no
    alternate is proposed.
    """
    loops = []
    domains = []
    for f in features:
        if f.get("kind") == "domain":
            b = _valid_bounds(f, seq_len)
            if b:
                domains.append((b, f.get("name", ""), f.get("evidence", {}).get("kind", "")))
    for f in sorted((f for f in features if f.get("kind") == "extracellular"),
                    key=lambda x: x.get("start") or 0):
        b = _valid_bounds(f, seq_len)
        if not b:
            loops.append({"start": f.get("start"), "end": f.get("end"), "length": None,
                          "boundary_status": f.get("boundary_status", "fuzzy"),
                          "domains_fully_contained": [], "eligible_domain_alternate": False,
                          "note": "fuzzy/invalid loop coordinates; deferred to review, not used"})
            continue
        contained = [{"name": name, "start": db[0], "end": db[1], "evidence_kind": kind}
                     for (db, name, kind) in domains if contains(b, db)]
        loops.append({"start": b[0], "end": b[1], "length": b[1] - b[0] + 1,
                      "boundary_status": f.get("boundary_status", "exact"),
                      "domains_fully_contained": contained,
                      "eligible_domain_alternate": bool(contained)})
    return loops
