"""Batch analysis-set construction from public references or user lists.

Identity rules: explicit accessions are processed directly; gene names are
mapped through a recorded UniProt query, ambiguity is preserved, nothing is
guessed and no input row is dropped. Every entry receives an explicit
disposition (membrane membership / natural cell-surface target / design
scope / extension set / technical failure); nothing is silently deleted.
"""

import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from . import __version__
from .common import digest, file_hash, now, read_json, write_json
from .uniprot import RETRYABLE_HTTP, fetch, normalize, retry_delay

ACCESSION_RE = re.compile(r"[A-Z0-9]{6,10}(?:-[0-9]+)?")
SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
GENE_MAP_FIELDS = "accession,gene_names,protein_name"
# Documented default for a public human membrane-protein reference set.
# Recorded verbatim in the set definition; callers may pass their own query.
DEFAULT_SET_QUERY = '(organism_id:9606) AND (reviewed:true) AND (keyword:"Membrane" OR keyword:"Cell membrane")'


def read_target_list(path):
    """Read a user target list; every row is preserved with its raw value.

    A header line is recognised when a column is named accession or gene.
    Otherwise the first column is treated as the target value. Blank or
    unclassifiable values stay in the list with identity_status unresolved.
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8-sig")
    lines = [ln for ln in text.splitlines()]
    delimiter = "\t" if path.suffix.lower() == ".tsv" or "\t" in (lines[0] if lines else "") else ","
    parsed = list(csv.reader(lines, delimiter=delimiter))
    rows = []
    header = None
    if parsed:
        lowered = [c.strip().lower() for c in parsed[0]]
        if any(c in {"accession", "uniprot", "uniprot_accession", "gene", "gene_name", "gene_symbol"} for c in lowered):
            header = parsed[0]
            parsed = parsed[1:]
    acc_col = gene_col = None
    if header:
        for i, name in enumerate(header):
            key = name.strip().lower()
            if key in {"accession", "uniprot", "uniprot_accession"}:
                acc_col = i
            if key in {"gene", "gene_name", "gene_symbol"}:
                gene_col = i
    for i, fields in enumerate(parsed):
        row_no = i + 2 if header else i + 1
        if not any(f.strip() for f in fields):
            rows.append({"row": row_no, "value": "", "input_type": "blank", "extra": {}})
            continue
        accession = fields[acc_col].strip() if acc_col is not None and acc_col < len(fields) else ""
        gene = fields[gene_col].strip() if gene_col is not None and gene_col < len(fields) else ""
        if not header:
            value = fields[0].strip()
            if ACCESSION_RE.fullmatch(value):
                accession, input_type = value, "accession"
            else:
                gene, input_type = value, "gene_name"
        else:
            value = accession or gene
            input_type = "accession" if accession else "gene_name" if gene else "blank"
            if accession and not ACCESSION_RE.fullmatch(accession):
                rows.append({"row": row_no, "value": accession, "input_type": "invalid_accession", "extra": _extra(header, fields)})
                continue
        extra = _extra(header, fields) if header else {}
        rows.append({"row": row_no, "value": value, "input_type": input_type,
                     "accession": accession, "gene_name": gene, "extra": extra})
    return rows


def _extra(header, fields):
    return {header[i]: fields[i] for i in range(min(len(header), len(fields)))}


def _search_page(query, opener, fields=GENE_MAP_FIELDS, size=25, cursor=None):
    params = {"query": query, "fields": fields, "size": str(size)}
    if cursor:
        params["cursor"] = cursor
    url = SEARCH_URL + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "ecd-snipr-harness/" + __version__})
    opener = opener or urllib.request.urlopen
    with opener(request, timeout=30) as response:
        payload = response.read()
        headers = dict(response.headers.items()) if hasattr(response, "headers") else {}
    link = headers.get("Link", headers.get("link", ""))
    next_url = ""
    match = re.search(r'<([^>]+)>;\s*rel="next"', link)
    if match:
        next_url = match.group(1)
    return json.loads(payload), headers, next_url


def map_gene_name(name, cache, offline=False, refresh=False, opener=None, sleeper=time.sleep):
    """Map a human gene name to UniProt accessions with recorded evidence.

    Exactly one hit resolves the name. Zero or multiple hits are preserved as
    unresolved/ambiguous with the full candidate list; nothing is guessed.
    """
    if not name or not re.fullmatch(r"[A-Za-z0-9._-]+", name):
        return {"status": "unresolved", "gene_name": name, "reason": "gene_name_format_not_supported",
                "candidates": [], "retrieved_at": now()}
    cache = Path(cache) / "gene_map"
    cache.mkdir(parents=True, exist_ok=True)
    key = digest({"gene_name": name.upper()})[:20]
    meta_path = cache / (key + ".json")
    if not refresh and meta_path.exists():
        try:
            record = read_json(meta_path)
            if record.get("gene_name") == name.upper():
                return dict(record, status=record["mapping_status"], cached=True)
        except (ValueError, OSError):
            pass
    if offline:
        return {"status": "missing", "gene_name": name, "reason": "offline_cache_absent", "retrieved_at": now()}
    queries = []
    hits, release = None, "unknown"
    for reviewed_only in (True, False):  # broader query only when the strict one has zero hits
        query = f"gene_exact:{name} AND organism_id:9606" + (" AND reviewed:true" if reviewed_only else "")
        hits, release, error = _query_with_retry(query, opener, sleeper)
        record_q = {"query": query, "retrieved_at": now(), "uniprot_release": release}
        if error:
            queries.append(dict(record_q, error=error))
            return {"status": "error", "gene_name": name, "reason": error, "queries": queries, "retrieved_at": now()}
        queries.append(dict(record_q, hits=len(hits)))
        if hits:
            break
    candidates = [h for h in (hits or []) if h["gene"].upper() == name.upper()] or (hits or [])
    if len(candidates) == 1:
        result = {"mapping_status": "resolved", "accession": candidates[0]["accession"], "candidates": candidates}
    elif candidates:
        result = {"mapping_status": "ambiguous", "candidates": candidates}
    else:
        result = {"mapping_status": "unresolved", "candidates": []}
    record = {"gene_name": name.upper(), "input_gene_name": name, "queries": queries,
              "retrieved_at": now(), "mapping_policy": "gene_exact + organism_id:9606; reviewed first, then all entries; ambiguity preserved"}
    record.update(result)
    write_json(meta_path, record)
    return dict(record, status=record["mapping_status"])


def _fetch_page_with_retry(query, opener, sleeper, size, cursor):
    """One search page; bounded retries on transient codes only, honoring a
    bounded Retry-After hint. Non-retryable HTTP errors raise immediately."""
    for attempt in range(3):
        try:
            return _search_page(query, opener, size=size, cursor=cursor)
        except (OSError, ValueError) as exc:
            if isinstance(exc, urllib.error.HTTPError) and exc.code not in RETRYABLE_HTTP:
                raise
            if attempt == 2:
                raise
            sleeper(retry_delay(exc, attempt))


def _query_with_retry(query, opener, sleeper, size=100, max_pages=20):
    """Return (hits, release, error); pagination is followed so a multi-page
    result can never be silently truncated into a false 'resolved' mapping."""
    hits, release, cursor, pages = [], "unknown", None, 0
    while True:
        try:
            data, headers, next_url = _fetch_page_with_retry(query, opener, sleeper, size, cursor)
        except (OSError, ValueError) as exc:
            return [], "unknown", str(exc)
        release = headers.get("x-uniprot-release", headers.get("X-UniProt-Release", release))
        for entry in data.get("results", []):
            genes = entry.get("genes") or [{}]
            hits.append({"accession": entry.get("primaryAccession", ""),
                         "gene": genes[0].get("geneName", {}).get("value", ""),
                         "reviewed": str(entry.get("entryType", "")).startswith("UniProtKB reviewed")})
        pages += 1
        match = re.search(r"[?&]cursor=([^&]+)", next_url or "")
        if not match:
            return hits, release, ""
        if pages >= max_pages:
            return [], "unknown", "pagination_cap_reached_incomplete_candidate_list"
        cursor = match.group(1)
        sleeper(0.2)


def search_accessions(query, limit=None, offline=False, opener=None, sleeper=time.sleep):
    """List accessions for a public reference-set query; records release and completeness."""
    if offline:
        return {"status": "missing", "reason": "offline_query_not_possible", "query": query, "retrieved_at": now()}
    entries, seen, cursor = [], set(), None
    release, total = "unknown", None
    try:
        while True:
            data, headers, next_url = _fetch_page_with_retry(query, opener, sleeper,
                                                             size=min(500, limit or 500), cursor=cursor)
            release = headers.get("x-uniprot-release", headers.get("X-UniProt-Release", release))
            if total is None:
                raw_total = headers.get("x-total-results", headers.get("X-Total-Results"))
                total = int(raw_total) if raw_total and str(raw_total).isdigit() else None
            for entry in data.get("results", []):
                acc = entry.get("primaryAccession", "")
                if acc and acc not in seen:
                    seen.add(acc)
                    genes = entry.get("genes") or [{}]
                    entries.append({"accession": acc, "gene": genes[0].get("geneName", {}).get("value", ""),
                                    "protein_name": (((entry.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName") or {}).get("value", "")})
                if limit and len(entries) >= limit:
                    entries = entries[:limit]
                    return {"status": "fetched", "query": query, "entries": entries, "retrieved_at": now(),
                            "uniprot_release": release, "total_results": total, "completeness": "truncated_by_limit"}
            if not next_url:
                break
            match = re.search(r"[?&]cursor=([^&]+)", next_url)
            if not match:
                break
            cursor = match.group(1)
            sleeper(0.2)
        completeness = "complete" if total is None or len(entries) >= total else "partial_result_set"
        return {"status": "fetched", "query": query, "entries": entries, "retrieved_at": now(),
                "uniprot_release": release, "total_results": total, "completeness": completeness}
    except (OSError, ValueError) as exc:
        return {"status": "error", "query": query, "reason": str(exc), "retrieved_at": now(),
                "entries": entries, "completeness": "partial_result_set"}


def select_analysis_reference(protein, confirmation="not_performed"):
    """Explicitly choose the analysis reference and record the rationale.

    Selecting the database canonical sequence for screening is not the
    laboratory's experimental isoform confirmation; that confirmation still
    gates downstream fusion assembly. Genuine identity ambiguity or version
    conflicts are preserved by callers, not resolved here. The annotated
    isoform inventory is recorded for reconciliation; sequence-level comparison
    against other isoforms stays explicitly not evaluated because the UniProt
    entry document does not carry isoform sequences.
    """
    p = protein
    if p.get("reference_coordinate_status") == "unverified_noncanonical_request":
        basis = "unresolved_noncanonical"
        rationale = "Requested isoform sequence/coordinates are not established by the canonical entry JSON"
    elif p.get("isoform") and not p.get("isoform_ambiguous"):
        basis, rationale = "explicit_isoform", "Input named an explicit accession/isoform; used as the analysis reference"
    else:
        displayed = sorted({v for i in (p.get("alternative_products") or {}).get("isoforms", [])
                            if i.get("sequence_status", "").lower() == "displayed"
                            for v in i.get("isoform_ids", [])})
        p["isoform"] = displayed[0] if len(displayed) == 1 else p.get("accession", "") + ":canonical"
        p["isoform_ambiguous"] = False
        basis = "database_canonical"
        rationale = ("Screening explicitly selected the database canonical sequence as the analysis reference. "
                     "This is not an experimental isoform confirmation; the laboratory isoform must still be "
                     "confirmed before any fusion assembly.")
    ap = p.get("alternative_products")
    if ap and ap.get("isoform_count", 0) > 1:
        iso_note = (f"UniProt 注释 {ap['isoform_count']} 个 isoform；主条目不含各 isoform 完整序列，"
                    "canonical 与其他 isoform 的序列/边界差异未评估")
    elif ap:
        iso_note = "UniProt 仅注释 canonical isoform（无其他已注释 isoform）"
    else:
        iso_note = "无 ALTERNATIVE PRODUCTS 注释（未注释可变 isoform，不证明不存在）"
    p["reference_selection"] = {"basis": basis, "selected_isoform": p["isoform"], "rationale": rationale,
                                "experimental_isoform_confirmation": confirmation,
                                "selected_sequence_sha256": digest(p.get("sequence", "")),
                                "selected_entry_version": p.get("evidence", {}).get("version"),
                                "annotated_isoform_count": (ap or {}).get("isoform_count"),
                                "isoform_comparison": "not_evaluated",
                                "isoform_comparison_note": iso_note,
                                "note": "Canonical/reference choice gates nothing in screening; experimental isoform confirmation gates fusion assembly."}
    return p


def disposition(protein):
    """Explicit membrane-set / surface-target / scope classification."""
    topology = protein.get("topology", "unknown")
    location = protein.get("location", "unknown")
    tm = [f for f in protein.get("features", []) if f.get("kind") == "transmembrane"]
    external = [f for f in protein.get("features", []) if f.get("kind") == "extracellular"]
    gpi = [f for f in protein.get("features", []) if f.get("kind") == "gpi_attachment_site"]
    membership = "yes" if (tm or gpi or location in {"plasma_membrane", "other_membrane", "membrane_unspecified"}) else \
        "no" if topology in {"secreted", "intracellular"} or location == "secreted" else "undetermined"
    if topology == "secreted" and location == "plasma_membrane":
        surface = "undetermined"  # dual annotated/peripheral presentation is not an integral surface-target proof
    elif topology == "secreted" or location == "secreted":
        surface = "no"
    elif location == "plasma_membrane" and (external or topology in {"type_i", "type_ii", "gpi"}):
        surface = "yes"
    elif location == "plasma_membrane" and topology == "multi_pass":
        surface = "yes"
    elif membership == "no":
        surface = "no"
    else:
        surface = "undetermined"
    reasons = []
    if protein.get("taxon_id") != 9606:
        reasons.append("non_human_record")
    if topology == "intracellular":
        reasons.append("intracellular_only")
    if location == "other_membrane" and not external and topology != "secreted":
        reasons.append("organelle_membrane_not_default_scope")
    if location == "unknown" and topology == "unknown" and not tm:
        reasons.append("membrane_annotation_insufficient")
    extension = "secreted_extension" if (topology == "secreted" or location == "secreted") else "core"
    scope = "out_of_scope" if reasons and any(r in {"non_human_record", "intracellular_only", "organelle_membrane_not_default_scope"} for r in reasons) else "in_scope"
    notes = []
    if topology == "multi_pass":
        notes.append("multi_pass_disposition: no loop stitching; only sourced external domains may be conditional candidates")
    if topology == "gpi":
        notes.append("gpi_disposition: mature boundary and omega residue must be verified")
    if location == "other_membrane":
        notes.append("organelle_disposition: lumenal/extracellular confusion is not resolved automatically")
    if location == "membrane_unspecified":
        notes.append("generic membrane annotation: surface localization unknown; not an organelle exclusion")
    if topology == "unknown":
        notes.append("under_annotated_disposition: retained for screening; likely insufficient_evidence without topology/boundary support")
    return {"membrane_set_membership": membership, "natural_cell_surface_target": surface,
            "scope_status": scope, "scope_reason_codes": reasons, "extension_set": extension,
            "disposition_notes": notes}


def build_set(rows=None, query=None, cache=None, outdir=None, offline=False, refresh=False,
              limit=None, opener=None, sleeper=time.sleep, inclusion_rules=None):
    """Build the analysis set. Single-item failure never halts the batch."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "proteins").mkdir(exist_ok=True)
    query_record = None
    items = []
    if query is not None:
        query_record = search_accessions(query, limit=limit, offline=offline, opener=opener, sleeper=sleeper)
        if query_record["status"] == "error" and not query_record.get("entries"):
            items = []
        else:
            items = [{"row": i + 1, "value": e["accession"], "input_type": "accession",
                      "accession": e["accession"], "gene_name": e.get("gene", ""), "extra": {"protein_name": e.get("protein_name", "")}}
                     for i, e in enumerate(query_record.get("entries", []))]
    else:
        items = rows or []
    entries = []
    for index, item in enumerate(items):
        entry_id = f"entry-{index + 1:05d}"
        entry = {"entry_id": entry_id, "input": {k: item.get(k) for k in ("row", "value", "input_type", "extra") if item.get(k) is not None and item.get(k) != ""},
                 "identity_status": "unresolved", "accession": "", "gene": "",
                 "mapping": None, "fetch": None, "protein_file": "", "protein_sha256": "",
                 "reference_selection": None, "processing_status": "pending", "processing_error": ""}
        try:
            accession = item.get("accession", "")
            if item.get("input_type") == "gene_name" and not accession:
                mapping = map_gene_name(item.get("gene_name", ""), cache or outdir / "cache", offline, refresh, opener, sleeper)
                entry["mapping"] = {k: mapping.get(k) for k in ("mapping_status", "gene_name", "candidates", "queries", "mapping_policy", "reason")}
                if mapping["status"] == "resolved":
                    accession = mapping["accession"]
                elif mapping["status"] in {"ambiguous", "unresolved", "missing", "error"}:
                    entry["identity_status"] = {"ambiguous": "ambiguous", "unresolved": "unresolved"}.get(mapping["status"], "unresolved")
                    entry["processing_status"] = "failed" if mapping["status"] == "error" else "ok"
                    entry["processing_error"] = mapping.get("reason", "") if mapping["status"] == "error" else ""
                    entry["disposition"] = {"scope_status": "unevaluated", "note": "identity not resolved; ambiguity preserved" if mapping["status"] == "ambiguous" else "identity not resolved"}
                    entries.append(entry)
                    continue
            entry["accession"] = accession
            if not accession:
                entry["identity_status"] = "unresolved"
                entry["processing_status"] = "ok"
                entry["disposition"] = {"scope_status": "unevaluated", "note": "no usable accession or gene name in this row"}
                entries.append(entry)
                continue
            status = fetch(accession, cache or outdir / "cache", offline=offline, refresh=refresh, opener=opener, sleeper=sleeper)
            entry["fetch"] = {k: status.get(k) for k in ("status", "accession", "url", "retrieved_at", "source_version", "raw_response_sha256", "reason", "attempts")}
            if status["status"] not in {"fetched", "cached"}:
                entry["identity_status"] = "unresolved"
                entry["processing_status"] = "failed" if status["status"] == "error" else "ok"
                entry["processing_error"] = status.get("reason", "") if status["status"] == "error" else ""
                entry["disposition"] = {"scope_status": "unevaluated", "note": f"fetch_{status['status']}: {status.get('reason', 'not in offline cache')}"}
                entries.append(entry)
                continue
            raw = read_json(status["cache_path"])
            if "-" in accession:
                raw["requested_isoform"] = accession
            protein = select_analysis_reference(normalize(raw))
            protein["protein_id"] = accession
            protein_path = outdir / "proteins" / (accession + ".json")
            write_json(protein_path, protein)
            entry.update(identity_status="resolved", gene=protein.get("gene", "") or item.get("gene_name", ""),
                         protein_file=str(protein_path.relative_to(outdir)), protein_sha256=file_hash(protein_path),
                         reference_selection=protein["reference_selection"],
                         disposition=disposition(protein), processing_status="ok")
        except Exception as exc:  # per-item isolation: one failure never halts the batch
            entry["processing_status"] = "failed"
            entry["processing_error"] = f"{type(exc).__name__}: {exc}"
            entry["disposition"] = {"scope_status": "unevaluated", "note": "technical failure during acquisition"}
        entries.append(entry)
    # Duplicate inputs are recorded and linked, never silently merged or dropped.
    seen = {}
    for entry in entries:
        if entry["processing_status"] == "ok" and entry.get("accession") and entry.get("identity_status") == "resolved":
            if entry["accession"] in seen:
                entry["duplicate_of"] = seen[entry["accession"]]
                entry["disposition"]["duplicate_input"] = True
            else:
                seen[entry["accession"]] = entry["entry_id"]
    # Partial means re-running could add data: fetch/mapping missing or error,
    # any technical failure, or a truncated/partial query result. A gene name
    # with a definitive zero-hit answer is a complete, recorded outcome.
    resolved = sum(e["identity_status"] == "resolved" for e in entries)
    failed = sum(e["processing_status"] == "failed" for e in entries)
    incomplete = failed or any(
        (e.get("fetch") or {}).get("status") in {"missing", "error"} or
        ((e.get("mapping") or {}).get("reason") and (e.get("mapping") or {}).get("mapping_status") != "unresolved")
        for e in entries)
    completeness = "partial" if incomplete else "complete"
    if query_record and query_record.get("completeness") not in {None, "complete"}:
        completeness = "partial"
    if query_record and query_record.get("status") in {"missing", "error"}:
        completeness = "partial"
    definition = {
        "set_id": "set-" + digest({"items": [(i.get("accession"), i.get("gene_name"), i.get("value")) for i in items], "query": query})[:16],
        "created_at": now(), "generator_version": __version__,
        "database": {"name": "UniProtKB", "access": "https://rest.uniprot.org/uniprotkb",
                     "release": (query_record or {}).get("uniprot_release", "per-entry entryVersion recorded per fetch"),
                     "query_record": query_record},
        "query": query, "inclusion_rules": inclusion_rules or (
            [f"UniProt query: {query}"] if query else
            ["User-supplied target list; accessions processed directly",
             "Gene names mapped with gene_exact + organism_id:9606 (reviewed first); ambiguity preserved, no guessing, no row dropped"]),
        "counts": {"input_rows": len(entries), "resolved": resolved,
                   "unresolved": sum(e["identity_status"] != "resolved" for e in entries),
                   "technical_failures": failed,
                   "duplicate_input_rows": sum(bool(e.get("duplicate_of")) for e in entries)},
        "completeness": completeness,
        "entries": entries,
    }
    write_json(outdir / "set_definition.json", definition)
    return definition


def load_set_proteins(set_dir):
    """Yield (entry, protein) pairs from a built set, verifying stored hashes."""
    set_dir = Path(set_dir)
    definition = read_json(set_dir / "set_definition.json")
    proteins = {}
    for entry in definition["entries"]:
        if entry.get("protein_file"):
            path = set_dir / entry["protein_file"]
            if file_hash(path) != entry["protein_sha256"]:
                raise ValueError("Cached protein record checksum mismatch: " + entry["protein_file"])
            proteins[entry["entry_id"]] = read_json(path)
    return definition, proteins
