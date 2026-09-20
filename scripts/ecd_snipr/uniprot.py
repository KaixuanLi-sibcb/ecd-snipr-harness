"""UniProt sequence/features import; preserve ECO and uncertain positions."""

import json
import hashlib
import re
import time
import urllib.error
import urllib.request
from pathlib import Path
from .common import digest, file_hash, now, read_json, write_json
from . import __version__

RETRYABLE_HTTP = {429, 500, 502, 503, 504}
MAX_RETRY_AFTER = 60  # never sleep longer than this on a single Retry-After hint


def retry_delay(exc, attempt):
    """Backoff seconds for a transient failure; honors a bounded Retry-After."""
    if isinstance(exc, urllib.error.HTTPError) and exc.code in RETRYABLE_HTTP:
        try:
            hinted = (exc.headers or {}).get("Retry-After")
            if hinted is not None and str(hinted).strip().isdigit():
                return min(int(str(hinted).strip()), MAX_RETRY_AFTER)
        except (AttributeError, TypeError):
            pass
    return 2 ** attempt


# UniProt feature types mapped into the normalized record. Every entry keeps
# positions, the raw description, exact/fuzzy boundary status and evidence
# (curated vs prediction by ECO code). Candidate-defining kinds are interpreted
# by the engine; informational kinds (region/motif/sites/variants/...) are
# parsed for the molecular profile and never silently dropped.
FEATURE_KINDS = {"Signal": "signal_peptide", "Transmembrane": "transmembrane", "Domain": "domain",
                 "Disulfide bond": "disulfide", "Chain": "chain", "Peptide": "processed_peptide",
                 "Propeptide": "propeptide", "Region": "region", "Motif": "motif",
                 "Glycosylation": "glycosylation_site", "Site": "site", "Binding site": "binding_site",
                 "Active site": "active_site", "Natural variant": "variant", "Mutagenesis": "mutagenesis",
                 "Lipidation": "lipidation"}

# Comment types whose text is machine-readable enough to record verbatim.
# SUBUNIT feeds the deterministic hetero-oligomer scan; FUNCTION/PTM are
# recorded for review only (no automated interpretation).
COMMENT_TEXT_TYPES = {"SUBUNIT": "subunit_comments", "FUNCTION": "function_comments", "PTM": "ptm_comments"}


def normalize(raw):
    accession = raw["primaryAccession"]
    seq = raw["sequence"]["value"]
    version = str(raw.get("entryAudit", {}).get("entryVersion", "unknown"))
    reviewed = raw.get("entryType", "").startswith("UniProtKB reviewed")
    source = {"kind": "curated_annotation" if reviewed else "prediction",
              "source": f"https://www.uniprot.org/uniprotkb/{accession}", "version": version,
              "raw_sha256": digest(raw)}
    p = {"protein_id": accession, "accession": accession, "isoform": raw.get("requested_isoform", ""),
         "isoform_ambiguous": not bool(raw.get("requested_isoform")),
         "gene": (raw.get("genes") or [{}])[0].get("geneName", {}).get("value", ""),
         "taxon_id": raw.get("organism", {}).get("taxonId"), "sequence": seq,
         "reviewed": reviewed, "annotation_score": raw.get("annotationScore"),
         "evidence": source, "features": [], "topology": "unknown", "location": "unknown",
         "notes": ["Canonical fetch does not establish the intended experimental isoform."]}
    # Isoform inventory (recording only, never a screening gate): the main entry
    # document lists annotated isoforms and textual difference features but does
    # NOT carry isoform sequences, so canonical-vs-isoform comparison stays
    # explicitly not evaluated unless the isoform itself was the fetched accession.
    isoforms = []
    events = set()
    for comment in raw.get("comments", []):
        if comment.get("commentType") == "ALTERNATIVE PRODUCTS":
            events.update(comment.get("events", []))
            for iso in comment.get("isoforms", []):
                name = iso.get("name") or {}
                isoforms.append({"name": name.get("value", "") if isinstance(name, dict) else str(name),
                                 "isoform_ids": iso.get("isoformIds", []),
                                 "sequence_status": iso.get("isoformSequenceStatus", "")})
    p["alternative_products"] = ({"isoform_count": len(isoforms), "isoforms": isoforms,
                                  "events": sorted(events)} if isoforms else None)
    displayed = [i for i in isoforms if i["sequence_status"].lower() == "displayed"]
    canonical_ids = {v for i in displayed for v in i["isoform_ids"]}
    p["canonical_isoform_ids"] = sorted(canonical_ids)
    requested = raw.get("requested_isoform", "")
    p["reference_coordinate_status"] = "canonical_entry"
    if requested and requested not in canonical_ids and requested != accession:
        # The main entry JSON is not an isoform-specific coordinate map.
        p["reference_coordinate_status"] = "unverified_noncanonical_request"
        p["isoform_ambiguous"] = True
    canonical_molecules = {accession, *canonical_ids}
    canonical_molecules.update(v for i in displayed for v in (i["name"], "Isoform " + i["name"]))
    principal_name = raw.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
    # A named mature main chain (e.g. LAG3) is not its released shed fragment.
    # Accept only an exact main-name chain ending at the reference terminus.
    for f in raw.get("features", []):
        loc = f.get("location", {})
        if (principal_name and f.get("type") == "Chain" and f.get("description") == principal_name
                and loc.get("sequence", accession) in {accession, *canonical_ids}
                and loc.get("end", {}).get("value") == len(seq)
                and all(loc.get(side, {}).get("modifier", "EXACT") == "EXACT" for side in ("start", "end"))):
            canonical_molecules.add(principal_name)
    p["excluded_annotations"] = []
    p["isoform_differences"] = [
        {"start": (f.get("location", {}).get("start") or {}).get("value"),
         "end": (f.get("location", {}).get("end") or {}).get("value"),
         "description": f.get("description", ""), "feature_id": f.get("featureId", ""),
         "evidence": dict(source, eco=f.get("evidences", []))}
        for f in raw.get("features", []) if f.get("type") == "Alternative sequence"
        and f.get("location", {}).get("sequence", accession) in {accession, *canonical_ids}]
    locations, raw_locations = [], []
    for comment in raw.get("comments", []):
        molecule = comment.get("molecule", "")
        relevant = comment.get("commentType") in {"SUBCELLULAR LOCATION", *COMMENT_TEXT_TYPES}
        applies = not molecule or molecule in canonical_molecules
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for entry in comment.get("subcellularLocations", []):
                value = entry.get("location", {}).get("value", "")
                raw_locations.append(value)
                if applies:
                    locations.append(value)
        if relevant and not applies:
            p["excluded_annotations"].append({"reason": "molecule_context_not_reference",
                                              "molecule": molecule, "raw_comment": comment,
                                              "evidence": source})
            continue
        target = COMMENT_TEXT_TYPES.get(comment.get("commentType"))
        if target:
            texts = [t.get("value", "") for t in comment.get("texts", []) if t.get("value")]
            if texts:
                p.setdefault(target, []).append(
                    {"text": " ".join(texts),
                     "eco": [e.get("evidenceCode", "") for e in comment.get("evidences", [])]})
    p["raw_locations"] = raw_locations
    p["reference_locations"] = locations
    p["has_secreted_reference_annotation"] = any(v.lower() == "secreted" or v.lower().startswith("secreted,") for v in locations)
    if any(v.lower() in {"cell membrane", "cell surface", "plasma membrane"}
           or v.lower().endswith(" cell membrane") for v in locations):
        p["location"] = "plasma_membrane"
    elif p["has_secreted_reference_annotation"]:
        # Secreted precedence over generic organelle membrane: entries annotated
        # e.g. "Secreted" + "Membrane" (secretory-pathway membrane) are secreted
        # proteins, not organelle-membrane residents. Plasma membrane still wins.
        p["location"] = "secreted"
    elif any(v.lower() == "membrane" for v in locations):
        p["location"] = "membrane_unspecified"
    elif any("membrane" in v.lower() for v in locations):
        p["location"] = "other_membrane"
    for f in raw.get("features", []):
        kind = FEATURE_KINDS.get(f.get("type"))
        if f.get("type") == "Lipidation" and "GPI-anchor" in f.get("description", ""):
            kind = "gpi_attachment_site"
        if f.get("type") == "Topological domain":
            desc = f.get("description", "").lower()
            kind = "extracellular" if desc == "extracellular" else "cytoplasmic" if desc == "cytoplasmic" else "lumenal" if "lumen" in desc or "lumin" in desc else "other_topology"
        if not kind:
            continue
        loc = f.get("location", {})
        if loc.get("sequence", accession) not in {accession, *canonical_ids}:
            p["excluded_annotations"].append({"reason": "feature_coordinates_on_other_isoform",
                                              "kind": kind, "raw_feature": f, "evidence": source})
            continue
        s, e = loc.get("start", {}), loc.get("end", {})
        evidence = dict(source, eco=f.get("evidences", []))
        if any(v.get("evidenceCode") in {"ECO:0000255", "ECO:0000256", "ECO:0000259"} for v in f.get("evidences", [])):
            evidence["kind"] = "prediction"
        p["features"].append({"kind": kind, "start": s.get("value"), "end": e.get("value"),
                              "boundary_status": "exact" if s.get("modifier", "EXACT") == e.get("modifier", "EXACT") == "EXACT" else "fuzzy",
                              "name": f.get("description", kind), "feature_id": f.get("featureId", ""),
                              "evidence": evidence, "raw_feature": f})
    tm = [f for f in p["features"] if f["kind"] == "transmembrane"]
    ex = [f for f in p["features"] if f["kind"] == "extracellular"]
    if len(tm) > 1:
        p["topology"] = "multi_pass"
    elif not tm and any(f["kind"] == "gpi_attachment_site" for f in p["features"]):
        p["topology"] = "gpi"
    elif not tm and p["has_secreted_reference_annotation"]:
        p["topology"] = "secreted"
    elif len(tm) == 1 and len(ex) == 1 and all(type(v) is int for v in [tm[0]["start"], tm[0]["end"], ex[0]["start"], ex[0]["end"]]):
        if ex[0]["end"] < tm[0]["start"]:
            p["topology"] = "type_i"
        elif ex[0]["start"] > tm[0]["end"]:
            p["topology"] = "type_ii"
    return p


def fetch(accession, cache, offline=False, refresh=False, opener=None, sleeper=time.sleep):
    if not re.fullmatch(r"[A-Z0-9]{6,10}(?:-[0-9]+)?", accession):
        raise ValueError("Explicit UniProt accession required, no gene-symbol guessing")
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    meta_path = cache / (accession + ".meta.json")
    if not refresh and meta_path.exists():
        try:
            metadata = read_json(meta_path)
            raw_path = cache / metadata.get("raw_file", "absent")
            if raw_path.parent == cache and raw_path.is_file() and metadata.get("raw_response_sha256") == file_hash(raw_path):
                return dict(metadata, status="cached", cache_path=str(raw_path))
        except (ValueError, OSError):
            pass
    if offline:
        return {"status": "missing", "accession": accession, "reason": "offline_cache_absent_or_invalid", "retrieved_at": now()}
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    opener = opener or urllib.request.urlopen
    error = ""
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "ecd-snipr-harness/" + __version__})
            with opener(request, timeout=30) as response:
                payload = response.read()
            raw = json.loads(payload)
            if not raw.get("sequence", {}).get("value"):
                raise ValueError("API response missing sequence")
            if raw.get("primaryAccession") != accession.split("-")[0]:
                raise ValueError("Returned accession differs; explicit ID resolution required")
            payload_hash = hashlib.sha256(payload).hexdigest()
            raw_path = cache / f"{accession}.{payload_hash}.json"
            if raw_path.exists() and file_hash(raw_path) != payload_hash:
                # Preserve an externally corrupted cache object rather than overwrite it.
                raw_path = cache / f"{accession}.{payload_hash}.{time.time_ns()}.json"
            if not raw_path.exists():
                raw_path.write_bytes(payload)
            meta = {"status": "fetched", "accession": accession, "url": url, "retrieved_at": now(), "source_version": str(raw.get("entryAudit", {}).get("entryVersion", "unknown")), "raw_response_sha256": file_hash(raw_path), "raw_file": raw_path.name, "attempts": attempt+1}
            write_json(cache / f"{accession}.{time.time_ns()}.retrieval.json", meta)
            write_json(meta_path, meta)
            return dict(meta, cache_path=str(raw_path))
        except (OSError, ValueError) as exc:
            error = str(exc)
            if isinstance(exc, urllib.error.HTTPError) and exc.code not in RETRYABLE_HTTP:
                break
            if attempt < 2:
                sleeper(retry_delay(exc, attempt))
    return {"status": "error", "accession": accession, "url": url, "retrieved_at": now(), "reason": error}
