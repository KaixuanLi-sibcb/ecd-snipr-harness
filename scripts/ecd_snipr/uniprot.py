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


def normalize(raw):
    accession = raw["primaryAccession"]
    seq = raw["sequence"]["value"]
    version = str(raw.get("entryAudit", {}).get("entryVersion", "unknown"))
    source = {"kind": "curated_annotation" if raw.get("entryType", "").startswith("UniProtKB reviewed") else "prediction",
              "source": f"https://www.uniprot.org/uniprotkb/{accession}", "version": version,
              "raw_sha256": digest(raw)}
    p = {"protein_id": accession, "accession": accession, "isoform": raw.get("requested_isoform", ""),
         "isoform_ambiguous": not bool(raw.get("requested_isoform")),
         "gene": (raw.get("genes") or [{}])[0].get("geneName", {}).get("value", ""),
         "taxon_id": raw.get("organism", {}).get("taxonId"), "sequence": seq,
         "evidence": source, "features": [], "topology": "unknown", "location": "unknown",
         "notes": ["Canonical fetch does not establish the intended experimental isoform."]}
    locations = []
    for comment in raw.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            for entry in comment.get("subcellularLocations", []):
                locations.append(entry.get("location", {}).get("value", ""))
    p["raw_locations"] = locations
    if any(v.lower() in {"cell membrane", "cell surface", "plasma membrane"} for v in locations):
        p["location"] = "plasma_membrane"
    elif any(v.lower() == "secreted" for v in locations):
        # Secreted precedence over generic organelle membrane: entries annotated
        # e.g. "Secreted" + "Membrane" (secretory-pathway membrane) are secreted
        # proteins, not organelle-membrane residents. Plasma membrane still wins.
        p["location"] = "secreted"
    elif any("membrane" in v.lower() for v in locations):
        p["location"] = "other_membrane"
    for f in raw.get("features", []):
        kind = {"Signal": "signal_peptide", "Transmembrane": "transmembrane", "Domain": "domain", "Disulfide bond": "disulfide", "Chain": "chain", "Peptide": "processed_peptide", "Propeptide": "propeptide"}.get(f.get("type"))
        if f.get("type") == "Lipidation" and "GPI-anchor" in f.get("description", ""):
            kind = "gpi_attachment_site"
        if f.get("type") == "Topological domain":
            desc = f.get("description", "").lower()
            kind = "extracellular" if desc == "extracellular" else "cytoplasmic" if desc == "cytoplasmic" else "lumenal" if "lumen" in desc or "lumin" in desc else "other_topology"
        if not kind:
            continue
        loc = f.get("location", {})
        s, e = loc.get("start", {}), loc.get("end", {})
        evidence = dict(source, eco=f.get("evidences", []))
        if any(v.get("evidenceCode") in {"ECO:0000255", "ECO:0000256", "ECO:0000259"} for v in f.get("evidences", [])):
            evidence["kind"] = "prediction"
        p["features"].append({"kind": kind, "start": s.get("value"), "end": e.get("value"),
                              "boundary_status": "exact" if s.get("modifier", "EXACT") == e.get("modifier", "EXACT") == "EXACT" else "fuzzy",
                              "name": f.get("description", kind), "evidence": evidence, "raw_feature": f})
    tm = [f for f in p["features"] if f["kind"] == "transmembrane"]
    ex = [f for f in p["features"] if f["kind"] == "extracellular"]
    if len(tm) > 1:
        p["topology"] = "multi_pass"
    elif not tm and any(f["kind"] == "gpi_attachment_site" for f in p["features"]):
        p["topology"] = "gpi"
    elif not tm and p["location"] == "secreted":
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
            if isinstance(exc, urllib.error.HTTPError) and exc.code not in {429, 500, 502, 503, 504}:
                break
            if attempt < 2:
                sleeper(2 ** attempt)
    return {"status": "error", "accession": accession, "url": url, "retrieved_at": now(), "reason": error}
