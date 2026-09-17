"""Strict structural checks before any candidate generation."""

from .common import unique
from .design import CONFIGURATION
from .criteria import CONSIDERATIONS, FORMS, SCOPES, STATUSES

EVIDENCE_FIELDS = "record_version candidate_id gene_symbol species uniprot_accession ensembl_gene_id evidence_type evidence_value normalized_value evidence_status source_name source_authority_level source_url source_version retrieved_at query adapter_name adapter_version raw_response_sha256 confidence conflict_status failure_mode license_note".split()
EVIDENCE_STATUSES = {"confirmed_live", "confirmed_local_xlsx", "confirmed_fixture", "inferred_live", "inferred_local_xlsx", "missing", "not_applicable", "conflict", "error"}


def evidence_v2(record):
    if set(record) != set(EVIDENCE_FIELDS) or any(not isinstance(v, str) for v in record.values()):
        raise ValueError("EvidenceRecord v2 requires exactly the existing string fields")
    if record["record_version"] != "2.0" or record["evidence_status"] not in EVIDENCE_STATUSES:
        raise ValueError("Invalid EvidenceRecord v2 version/status")
    if record["source_authority_level"] not in {"local_lab", "fixture", "curated_database", "public_database", "computed", "unknown"}:
        raise ValueError("Invalid evidence authority")
    if record["conflict_status"] not in {"none", "value_conflict", "source_conflict", "version_conflict", "missing_conflict", "not_checked"}:
        raise ValueError("Invalid conflict status")
    if record["evidence_status"] in {"missing", "error"} and (record["normalized_value"] or record["evidence_value"]):
        raise ValueError("Missing/error evidence must not contain asserted values")


def validate_project(project):
    if not isinstance(project, dict) or project.get("schema_version") != "1.0":
        raise ValueError("Project schema_version must be 1.0")
    if project.get("configuration") != CONFIGURATION:
        raise ValueError("This skill supports Antibody-Sender -> Antigen-Receiver only")
    if not isinstance(project.get("proteins"), list):
        raise ValueError("proteins list required")
    if not project["proteins"]:
        raise ValueError("At least one normalized protein record required")
    for field in ("proteins", "existing_constructs", "observations", "evidence_records_v2", "reviews"):
        if field in project and (not isinstance(project[field], list) or any(not isinstance(v, dict) for v in project[field])):
            raise ValueError(field + " must contain objects")
    for field, key in (("proteins", "protein_id"), ("existing_constructs", "construct_id"), ("observations", "observation_id")):
        if any(not isinstance(v.get(key), str) or not v[key] for v in project.get(field, [])):
            raise ValueError(key + " required")
        unique(project.get(field, []), key)
    for p in project["proteins"]:
        if not isinstance(p.get("evidence", {}), dict):
            raise ValueError("Protein evidence must be an object")
        if "isoform_ambiguous" in p and type(p["isoform_ambiguous"]) is not bool:
            raise ValueError("isoform_ambiguous must be a boolean")
        if p.get("candidate_policy", "full_ecd_only") not in {"full_ecd_and_domains", "full_ecd_only"}:
            raise ValueError("Unknown candidate_policy")
        if p.get("topology", "unknown") not in {"type_i", "type_ii", "gpi", "multi_pass", "secreted", "intracellular", "unknown"}:
            raise ValueError("Unknown topology category")
        for field in ("features", "candidate_regions", "considerations"):
            if not isinstance(p.get(field, []), list) or any(not isinstance(f, dict) for f in p.get(field, [])):
                raise ValueError(field + " must be a list of feature objects")
        if any(not f.get("kind") for f in p.get("features", [])):
            raise ValueError("Every feature needs a kind")
        if any(not isinstance(f.get("evidence", {}), dict) for f in p.get("features", []) + p.get("candidate_regions", [])):
            raise ValueError("Feature evidence must be an object")
        for region in p.get("candidate_regions", []):
            if region.get("antigen_form_type", "domain_fragment") not in FORMS:
                raise ValueError("Unknown antigen_form_type")
            if not isinstance(region.get("release_evidence", {}), dict):
                raise ValueError("release_evidence must be an evidence object")
            if not isinstance(region.get("release_context", ""), str):
                raise ValueError("release_context must be a description string")
        for record in p.get("considerations", []):
            if record.get("kind") not in CONSIDERATIONS or record.get("status") not in STATUSES or record.get("scope") not in SCOPES:
                raise ValueError("Invalid contextual consideration kind/status/scope")
            if not isinstance(record.get("evidence", {}), dict) or not isinstance(record.get("context", {}), dict):
                raise ValueError("Consideration evidence/context must be objects")
            if not isinstance(record.get("reason", ""), str):
                raise ValueError("Consideration reason must be a string")
    scaffold = project.get("scaffold")
    if scaffold is not None:
        if not isinstance(scaffold, dict) or not isinstance(scaffold.get("modules", []), list):
            raise ValueError("Scaffold must be a module configuration object or null")
        if any(not isinstance(m, dict) for m in scaffold.get("modules", [])):
            raise ValueError("Scaffold modules must be objects")
    for o in project.get("observations", []):
        if not isinstance(o.get("context", {}), dict):
            raise ValueError("Observation context must be an object")
    rules = project.get("rules", {})
    for name in ("short_length", "long_length"):
        if name in rules and (type(rules[name]) is not int or rules[name] <= 0):
            raise ValueError("Length thresholds must be positive integers")
    if rules.get("short_length", 50) >= rules.get("long_length", 500):
        raise ValueError("short_length must be smaller than long_length")
    if not 0 < rules.get("cysteine_fraction", 0.08) <= 1:
        raise ValueError("cysteine_fraction must be in (0,1]")
    for r in project.get("evidence_records_v2", []):
        evidence_v2(r)
    if project.get("reviews"):
        if any(not r.get("review_key") for r in project["reviews"]):
            raise ValueError("Every review must have a scoped review_key")
        unique(project["reviews"], "review_key")


def validate_artifacts(candidates, observations):
    unique(candidates, "candidate_id")
    for c in candidates:
        if c["sequence"] and len(c["sequence"]) != c["end"] - c["start"] + 1:
            raise ValueError("Candidate interval/sequence invariant failed")
        if c["assembly_status"] != "assembled_unvalidated" and c["fusion_sequence"]:
            raise ValueError("Unapproved or blocked fusion exported")
    for o in observations:
        if o["eligibility"] == "eligible_measurement_not_success_label" and o["issues"]:
            raise ValueError("Assay eligibility inconsistent")
