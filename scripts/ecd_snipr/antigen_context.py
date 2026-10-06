"""Sourced processing units and antigen-repertoire scope, not function prediction."""

from .common import contains, interval, overlaps, sourced


def processed_products(raw, accession, canonical_ids, evidence):
    """Link named mature chains to exact, same-reference comments only.

    Product-scoped localization must never be promoted to the whole precursor.
    Ambiguous names, fuzzy boundaries and foreign isoforms cannot rescue a unit.
    """
    chains = []
    chain_names = []
    for feature in raw.get("features", []):
        if feature.get("type") != "Chain":
            continue
        loc = feature.get("location", {})
        if loc.get("sequence", accession) not in {accession, *canonical_ids}:
            continue
        chain_names.append(feature.get("description", ""))
        chain = {"kind": "chain", "name": feature.get("description", ""),
                 "start": loc.get("start", {}).get("value"), "end": loc.get("end", {}).get("value"),
                 "boundary_status": "exact" if all(loc.get(s, {}).get("modifier", "EXACT") == "EXACT"
                                                     for s in ("start", "end")) else "fuzzy",
                 "evidence": dict(evidence, eco=feature.get("evidences", []))}
        try:
            interval(chain, len(raw["sequence"]["value"]))
        except ValueError:
            continue
        chains.append(chain)
    products = []
    for chain in chains:
        if not chain["name"] or chain_names.count(chain["name"]) != 1:
            continue
        comments = [c for c in raw.get("comments", []) if c.get("commentType") == "SUBCELLULAR LOCATION"
                    and c.get("molecule") == chain["name"]]
        if not comments:
            continue
        products.append(dict(chain, locations=[loc.get("location", {}).get("value", "")
                                              for c in comments for loc in c.get("subcellularLocations", [])],
                             raw_comments=comments, scope="processed_product_not_whole_reference"))
    return products


def _valid(features, length, kind):
    result = []
    for f in features:
        if f.get("kind") != kind or not sourced(f):
            continue
        try:
            interval(f, length)
        except ValueError:
            continue
        result.append(f)
    return result


def _record(f):
    return {k: f.get(k) for k in ("kind", "name", "start", "end", "evidence")}


def assess_antigen_context(protein, bounds):
    """Describe exact sequence coverage without inventing epitope completeness.

    Disjoint chain annotations nominate processing questions, not proven
    simultaneous products. Nested precursor/mature annotations alone do not.
    Repeats are not automatically treated as autonomous domain candidates.
    """
    length = len(protein.get("sequence", ""))
    fs = protein.get("features", [])
    chains = _valid(fs, length, "chain")
    split_pairs = []
    for i, a in enumerate(chains):
        ab = (a["start"], a["end"])
        for b in chains[i + 1:]:
            bb = (b["start"], b["end"])
            if not overlaps(ab, bb) and overlaps(bounds, ab) and overlaps(bounds, bb):
                split_pairs.append({"chains": [_record(a), _record(b)],
                                    "interpretation": "Fragment spans distinct annotated chain segments; processing and co-occurrence require review"})
    intramembrane = _valid(fs, length, "intramembrane")
    excluded = _valid(fs, length, "transmembrane") + intramembrane + _valid(fs, length, "cytoplasmic") + \
        _valid(fs, length, "signal_peptide") + _valid(fs, length, "propeptide") + _valid(fs, length, "gpi_signal")
    products = []
    for p in protein.get("processed_products", []):
        try:
            pb = interval(p, length)
        except ValueError:
            continue
        if not sourced(p) or not any(v.lower().startswith("secreted") or v.lower() in {"extracellular space", "extracellular matrix"}
                                    for v in p.get("locations", [])):
            continue
        if any(overlaps(pb, (f["start"], f["end"])) for f in excluded):
            continue
        state = "sequence_contained" if contains(bounds, pb) else "partial" if overlaps(bounds, pb) else "omitted"
        products.append(dict(_record(p), locations=p.get("locations", []), raw_comments=p.get("raw_comments", []),
                             coverage=state, supported_route="sourced_product_for_separate_review_not_auto_generated"))
    repeats = []
    for f in _valid(fs, length, "repeat"):
        fb = (f["start"], f["end"])
        if overlaps(bounds, fb):
            repeats.append(dict(_record(f), coverage="contained" if contains(bounds, fb) else "cut"))
    flags = []
    if split_pairs:
        flags.append({"code": "processed_chain_segments_spanned", "severity": "review",
                      "detail": "Distinct sourced chain segments overlap the fragment; natural or fusion processing is not established"})
    if any(p["coverage"] != "sequence_contained" for p in products):
        flags.append({"code": "extracellular_processed_product_not_covered", "severity": "review",
                      "detail": "Another exactly annotated extracellular mature product is not fully represented; clarify intended antigen unit"})
    if any(r["coverage"] == "cut" for r in repeats):
        flags.append({"code": "repeat_cut_by_boundary", "severity": "review",
                      "detail": "Fragment cuts a sourced repeat annotation; not proof of misfolding or an autonomous-domain definition"})
    has_epitope = bool(_valid(fs, length, "epitope"))
    context = {
        "processed_chain_pairs_spanned": split_pairs,
        "extracellular_mature_products": products,
        "overlapping_repeats": repeats,
        "reference_intramembrane_segments": [_record(f) for f in intramembrane],
        "full_ecd_label_scope": "complete_annotated_topological_region_not_all_mature_products_or_all_epitopes",
        "epitope_evidence_scope": "mapped_annotations_supplied_not_exhaustive" if has_epitope else "no_mapped_annotations_supplied",
        "external_epitope_database_search": "not_performed",
        "recognition_retention": "not_measured",
        "interpretation": "Sequence/annotation coverage only; a gene can encode several antigenic units. No new sequence or functional label inferred",
    }
    return context, flags
