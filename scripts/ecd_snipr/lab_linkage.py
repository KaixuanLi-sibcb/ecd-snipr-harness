"""Fragment/reference matches are not full-construct or functional equivalence."""

from collections import defaultdict
from pathlib import Path
import re
from .common import digest, file_hash, read_json, sequence, translate, interval, overlaps, contains
from .harness import verify_bundle
from .lab_intake import ACCESSION, ASSAY_FIELDS, field_values, nonblank
from .observations import ENDPOINTS, normalize_observation


def identity_hint(value):
    """A checkbox/status or composite note is never a shared construct identifier."""
    text = str(value).strip()
    if text.lower() in {"y", "n", "yes", "no", "na", "n/a", "x", "pending", "true", "false", "none"}:
        return None
    return text if re.fullmatch(r"[A-Za-z0-9_.-]+", text) else None


def load_references(set_dir, wanted):
    if not set_dir:
        return {}, [], None
    root = Path(set_dir).resolve()
    definition = read_json(root / "set_definition.json")
    references, errors, conflicted = {}, [], set()
    for entry in definition["entries"]:
        accession = entry.get("accession")
        if accession not in wanted:
            continue
        if accession in conflicted:
            continue
        try:
            path = (root / entry["protein_file"]).resolve()
            if not path.is_relative_to(root) or file_hash(path) != entry["protein_sha256"]:
                raise ValueError("Reference file outside set or checksum mismatch")
            protein = read_json(path)
            if protein["accession"] != accession:
                raise ValueError("Reference identity mismatch")
            if accession in references and digest(references[accession]) != digest(protein):
                references.pop(accession)
                conflicted.add(accession)
                raise ValueError("Conflicting reference versions for same accession; none selected")
            references[accession] = protein
        except (OSError, ValueError, KeyError, TypeError) as exc:
            errors.append({"accession": accession, "status": "reference_error", "reason": str(exc)})
    return references, errors, {"sha256": file_hash(root / "set_definition.json"),
                                "set_id": definition["set_id"], "completeness": definition["completeness"]}


def fragment_records(rows, references):
    result = []
    for row in rows:
        aas = [f for f in row["fields"] if f["field"] in {"ecd_aa", "ecd_unspecified"} and nonblank(f["raw_value"])]
        nts = [f for f in row["fields"] if f["field"] == "ecd_nt" and nonblank(f["raw_value"])]
        for source in aas or nts:
            record = {"fragment_record_id": digest([row["row_id"], source["cell_id"], "fragment"]),
                      "row_id": row["row_id"], "sheet": row["sheet"], "row": row["row"],
                      "source_cell": source["cell"], "source_cell_id": source["cell_id"],
                      "agid": field_values(row, "agid"), "receiver_plasmid": field_values(row, "receiver_plasmid"),
                      "receiver_label": field_values(row, "receiver_label"), "gene": field_values(row, "gene"),
                      "accession": field_values(row, "accession"), "sequence": "", "sequence_digest": "",
                      "sequence_status": "unresolved", "reference_status": "not_checked",
                      "matches": [], "issues": [], "full_construct_identity": "not_established",
                      "experimental_isoform_confirmation": "not_performed"}
            result.append(record)
            if row["layout_status"] != "mapped":
                record["issues"].append("layout_review_required_before_sequence_interpretation")
                continue
            if source.get("formula_present") or source.get("raw_formula"):
                record["issues"].append("formula_cached_sequence_not_recalculated")
            try:
                if source["field"] == "ecd_nt":
                    seq = sequence(translate(str(source["raw_value"])))
                    record["issues"].append("translated_cds_not_independently_confirmed_protein")
                else:
                    seq = sequence(str(source["raw_value"]))
                    if source["field"] == "ecd_unspecified":
                        if set(seq) <= set("ACGT"):
                            raise ValueError("Untyped ECD alphabet cannot distinguish DNA from protein")
                        record["issues"].append("amino_acid_kind_inferred_from_untyped_ECD_alphabet")
                record.update(sequence=seq, sequence_digest=digest(seq), sequence_status="valid_fragment_string")
            except ValueError as exc:
                record["sequence_status"] = "sequence_error_or_type_unresolved"
                record["issues"].append(str(exc))
                continue
            if aas and nts:
                if len(nts) != 1:
                    record["issues"].append("multiple_CDS_cells_no_automatic_pairing")
                else:
                    try:
                        if nts[0].get("formula_present") or nts[0].get("raw_formula"):
                            record["issues"].append("formula_CDS_cache_unconfirmed")
                        if translate(str(nts[0]["raw_value"])) != seq:
                            record["issues"].append("ECD_nt_translation_vs_aa_mismatch")
                    except ValueError as exc:
                        record["issues"].append("CDS_check:" + str(exc))
            ids = sorted(set(str(v).strip() for v in record["accession"]))
            if len(ids) != 1 or not ACCESSION.fullmatch(ids[0]):
                record["reference_status"] = "identity_unresolved_no_gene_guess"
                continue
            protein = references.get(ids[0])
            if not protein:
                record["reference_status"] = "reference_not_supplied_or_unavailable"
                continue
            try:
                full = sequence(protein["sequence"])
                positions, offset = [], 0
                while (pos := full.find(seq, offset)) >= 0:
                    positions.append((pos + 1, pos + len(seq)))
                    offset = pos + 1
                for bounds in positions:
                    match = {"reference": protein["protein_id"], "isoform": protein.get("isoform"),
                             "reference_sequence_digest": digest(full), "reference_evidence": protein.get("evidence", {}),
                             "start": bounds[0], "end": bounds[1], "coordinate_system": "1-based-inclusive",
                             "topology": protein.get("topology"), "retained_exclusions": [], "cut_domains": [],
                             "omitted_extracellular_domains": [], "epitope_status": "not_evaluated", "annotation_issues": []}
                    external = []
                    for f in protein.get("features", []):
                        if f.get("kind") == "extracellular":
                            try:
                                external.append(interval(f, len(full)))
                            except ValueError:
                                pass
                    for feature in protein.get("features", []):
                        try:
                            fb = interval(feature, len(full))
                        except ValueError:
                            match["annotation_issues"].append("invalid_feature:" + str(feature.get("kind")))
                            continue
                        kind = feature["kind"]
                        if kind in {"transmembrane", "intramembrane", "signal_peptide", "gpi_signal", "cytoplasmic", "propeptide"} and overlaps(bounds, fb):
                            match["retained_exclusions"].append({"kind": kind, "start": fb[0], "end": fb[1]})
                        if kind == "domain" and not contains(bounds, fb):
                            if overlaps(bounds, fb):
                                match["cut_domains"].append(feature)
                            elif any(contains(e, fb) for e in external):
                                match["omitted_extracellular_domains"].append(feature)
                    record["matches"].append(match)
                record["reference_status"] = "unique_exact_fragment_match" if len(positions) == 1 else "ambiguous_multiple_positions" if positions else "no_exact_match_not_proof_of_wrong_construct"
            except (ValueError, KeyError, TypeError) as exc:
                record["reference_status"] = "reference_error"
                record["issues"].append(str(exc))
    return result


def identity_conflicts(fragments):
    groups, conflicts = defaultdict(list), []
    for fragment in fragments:
        for field in ("agid", "receiver_plasmid", "receiver_label"):
            for label in fragment[field]:
                key = identity_hint(label)
                if key:
                    groups[(field, key)].append(fragment)
    for (field, label), records in groups.items():
        seqs = {r["sequence_digest"] for r in records if r["sequence_digest"]}
        accessions = {str(a).strip() for r in records for a in r["accession"]}
        if len(seqs) > 1 or len(accessions) > 1:
            conflicts.append({"type": "shared_label_different_fragment_or_reference", "label_field": field,
                              "label": label, "rows": [r["row_id"] for r in records],
                              "unique_fragment_strings": len(seqs), "accessions": sorted(accessions),
                              "interpretation": "possible intended variants or identity conflict; never auto-merge"})
    return conflicts


def proposed_links(rows, fragments):
    """Suggest review edges only; same gene/fragment never establishes a receptor."""
    index = defaultdict(list)
    fields = ("agid", "gene", "accession", "receiver_plasmid", "receiver_label")
    for f in fragments:
        if f["sequence"]:
            for field in fields:
                for value in f[field]:
                    key = identity_hint(value)
                    if key:
                        index[(field, key)].append(f)
    edges, readiness = [], []
    for row in rows:
        assays = [c for c in row["fields"] if c["field"] in ASSAY_FIELDS and nonblank(c["raw_value"])]
        if not assays:
            continue
        matches = {}
        if row["layout_status"] == "mapped":
            for field in fields:
                for value in field_values(row, field):
                    for f in index.get((field, identity_hint(value)), []):
                        key = f["fragment_record_id"]
                        item = matches.setdefault(key, {"assay_row_id": row["row_id"], "fragment_record_id": key,
                                                        "fragment_row_id": f["row_id"], "fragment_sequence_digest": f["sequence_digest"],
                                                        "matching_fields": [], "status": "proposed_not_confirmed_construct_link"})
                        if field not in item["matching_fields"]:
                            item["matching_fields"].append(field)
        edges.extend(matches.values())
        unique_seqs = {m["fragment_sequence_digest"] for m in matches.values()}
        readiness.append({"row_id": row["row_id"], "sheet": row["sheet"], "row": row["row"],
                          **{field: field_values(row, field) for field in fields}, "layout_status": row["layout_status"],
                          "nonblank_fields": [c["field"] for c in assays],
                          "candidate_fragment_records": len(matches), "candidate_unique_fragment_strings": len(unique_seqs),
                          "link_readiness": "layout_review" if row["layout_status"] != "mapped" else
                          "multiple_fragment_variants_review" if len(unique_seqs) > 1 else
                          "fragment_links_proposed_not_receiver_equivalence" if matches else "no_identity_link_found",
                          "actual_receiver_identity": "requires_explicit_construct_audit",
                          "assay_semantics": "see_per_cell_semantic_status"})
    return edges, readiness


def assay_records(rows, semantics=None, constructs=None):
    semantics, constructs = semantics or {}, constructs or {}
    reviews = semantics.get("reviews", [])
    if not isinstance(reviews, list) or any(not isinstance(r, dict) or not isinstance(r.get("cell_id"), str)
            or not r["cell_id"] or r.get("decision") not in {"pending", "confirmed"} for r in reviews):
        raise ValueError("Semantic reviews require explicit cell IDs and pending/confirmed decisions")
    if len({r.get("cell_id") for r in reviews}) != len(reviews):
        raise ValueError("Duplicate semantic review cell selectors")
    review_map = {r["cell_id"]: r for r in reviews}
    results, template = [], []
    for row in rows:
        for cell in row["fields"]:
            field = cell["field"]
            if field not in ASSAY_FIELDS:
                continue
            reported = {"self_activation": "basal_activation", "surface_mfi": "surface_expression"}.get(field, "unconfirmed")
            kind = "project_status_not_function" if field == "snipr_assay_status" else "downstream_serology_not_receiver_endpoint" if field == "serology_ec50" else "assay_readout"
            o = {"observation_id": cell["cell_id"], "row_id": row["row_id"], "field": field,
                 "reported_endpoint": reported, "record_kind": kind, "raw_value": cell["raw_value"],
                 "raw_formula": cell["raw_formula"], "raw_header": cell["raw_header"],
                 "formula_present": cell.get("formula_present", False),
                 "formula_attributes": cell.get("formula_attributes", {}),
                 "shared_formula_anchors": cell.get("shared_formula_anchors", []),
                 "raw_record_digest": digest(cell), "endpoint": reported, "unit": "unconfirmed",
                 "construct_id": "", "context": {}, "semantic_status": "unconfirmed",
                 "source_ref": {"workbook_sha256": cell["workbook_sha256"], "sheet": cell["sheet"],
                                "cell": cell["cell"], "row": cell["row"]}}
            review = review_map.get(cell["cell_id"])
            if review and review.get("decision") == "confirmed":
                if review.get("raw_record_digest") != digest(cell):
                    raise ValueError("Semantic review raw record hash mismatch")
                if not all(isinstance(review.get(k), str) and review[k].strip() for k in ("reviewer", "reviewed_at", "rationale")):
                    raise ValueError("Confirmed semantics need reviewer, date and rationale")
                if kind != "assay_readout" or row["layout_status"] != "mapped":
                    raise ValueError("Project status/serology/suspect layout cannot be promoted to receiver endpoint")
                if (cell.get("formula_present") or cell["raw_formula"]) and review.get("formula_cache_reviewed") is not True:
                    raise ValueError("Formula cache must be explicitly reviewed")
                definition = review.get("definition", {})
                allowed = {"endpoint", "unit", "denominator", "statistic", "gate_path", "channel", "construct_id", "context"}
                if not isinstance(definition, dict) or set(definition) - allowed or definition.get("endpoint") not in ENDPOINTS:
                    raise ValueError("Invalid semantic definition fields/endpoint")
                if not isinstance(definition.get("context", {}), dict):
                    raise ValueError("Assay context must be an object")
                o.update(definition, semantic_status="confirmed_definition_not_function", semantic_review=review)
            normalized = normalize_observation(o, constructs)
            if o["semantic_status"].startswith("confirmed") and o["construct_id"] in constructs:
                audited = constructs[o["construct_id"]]
                local_fragments = fragment_records([row], {})
                for f in local_fragments:
                    if any(i in f["issues"] for i in ("formula_cached_sequence_not_recalculated", "formula_CDS_cache_unconfirmed")) and review.get("sequence_formula_cache_reviewed") is not True:
                        normalized["issues"].append("sequence_formula_cache_requires_review")
                    if not f["sequence"] or "ECD_nt_translation_vs_aa_mismatch" in f["issues"]:
                        normalized["issues"].append("source_fragment_requires_review")
                    elif f["sequence"] != audited.get("antigen_sequence"):
                        normalized["issues"].append("reviewed_construct_conflicts_with_row_fragment")
            if o["semantic_status"] == "unconfirmed":
                normalized["issues"].append("semantic_definition_unconfirmed")
            if row["layout_status"] != "mapped":
                normalized["issues"].append("layout_requires_review")
            if normalized["issues"] and normalized["eligibility"] != "missing":
                normalized["eligibility"] = "unresolved"
            if normalized["eligibility"] != "eligible_measurement_not_success_label":
                # Partial numeric conversion is not an approved measurement.
                normalized["normalized_value"] = None
                normalized.pop("normalized_unit", None)
            results.append(normalized)
            if nonblank(cell["raw_value"]) and kind == "assay_readout":
                template.append({"cell_id": cell["cell_id"], "raw_record_digest": digest(cell),
                                 "source_ref": o["source_ref"], "field": field, "raw_value": cell["raw_value"],
                                 "decision": "pending", "reviewer": "", "reviewed_at": "", "rationale": "",
                                 "definition": {"endpoint": reported, "unit": "", "denominator": "", "statistic": "",
                                                "channel": "", "gate_path": "", "construct_id": "", "context": {}}})
    unknown = set(review_map) - {o["observation_id"] for o in results}
    if unknown:
        raise ValueError("Semantic review selects cells absent from this intake")
    return results, template


def audited_constructs(run_dir):
    if not run_dir:
        return {}, None
    root = Path(run_dir)
    if not verify_bundle(root):
        raise ValueError("Construct audit bundle integrity failed")
    records = read_json(root / "construct_audit.json")
    if len({c["construct_id"] for c in records}) != len(records):
        raise ValueError("Duplicate construct IDs in audited run")
    for record in records:
        if record.get("audit_status") == "sequence_consistent_not_functionally_validated":
            try:
                if sequence(record.get("antigen_sequence", "")) not in sequence(record.get("fusion_sequence", "")):
                    raise ValueError("Antigen fragment absent from supplied full fusion")
            except ValueError as exc:
                record.update(audit_status="requires_review", audit_issues=record.get("audit_issues", []) + [str(exc)])
    return {c["construct_id"]: c for c in records}, file_hash(root / "manifest.json")


def reviewed_assay_conflicts(observations):
    groups = defaultdict(list)
    for o in observations:
        if o["eligibility"] == "eligible_measurement_not_success_label":
            key = digest([o.get(k) for k in ("construct_id", "endpoint", "context", "unit", "statistic", "channel", "gate_path", "denominator")])
            groups[key].append(o)
    return [{"type": "duplicate_confirmed_context_different_values", "context_key": key,
             "observations": [o["observation_id"] for o in group], "resolution": "pending_no_pooling"}
            for key, group in groups.items() if len({str(o["raw_value"]) for o in group}) > 1]
