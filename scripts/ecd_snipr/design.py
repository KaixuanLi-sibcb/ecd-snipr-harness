"""Sourced ECD candidates and scaffold assembly with scoped human review."""

import re
from .common import contains, digest, interval, overlaps, sequence, sourced, translate
from .criteria import assess_context
from .profile import (INFORMATIONAL_KINDS, PREDICTION_RELEVANT_KINDS, boundary_analysis,
                      detect_multichain_partners, molecular_profile, n_glyco_sequons)

CONFIGURATION = "antibody_sender_antigen_receiver"
DEFAULT_RULES = {"short_length": 50, "long_length": 500, "cysteine_fraction": 0.08,
                 # Review heuristics, never pass/fail thresholds: dense glycosylation
                 # fires at >= min sequons AND <= max residues per sequon; boundary
                 # adjacency tolerance is the unannotated gap we still accept as
                 # agreement-with-tolerance (recorded, not silent).
                 "dense_glyco_min_sequons": 4, "dense_glyco_max_residues_per_sequon": 50,
                 "boundary_gap_tolerance": 5}


def flag(code, severity="review", detail=""):
    return {"code": code, "severity": severity, "detail": detail}


def scaffold_check(scaffold):
    if not scaffold:
        return [flag("scaffold_missing", "block", "Candidate fragments only; no fusion sequence")]
    errors = []
    if scaffold.get("record_type") == "public_reference_architecture" or scaffold.get("status") == "reference_only":
        errors.append(flag("reference_architecture_not_lab_scaffold", "block", "Continue fragment analysis; public reference cannot authorize fusion assembly"))
    if scaffold.get("configuration") != CONFIGURATION:
        errors.append(flag("wrong_configuration", "block"))
    for k in ("scaffold_id", "version", "source", "host_cell", "reporter"):
        if scaffold.get(k) in {None, "", "unknown", "pending"}:
            errors.append(flag("scaffold_metadata_missing", "block", k))
    if scaffold.get("status") not in {"verified_local", "synthetic_fixture"}:
        errors.append(flag("scaffold_not_verified", "block"))
    modules = scaffold.get("modules", [])
    roles = [m.get("role") for m in modules]
    required = ("signal_peptide", "antigen_slot", "transmembrane", "transcription_factor")
    if any(roles.count(r) != 1 for r in required):
        errors.append(flag("scaffold_required_module_count", "block"))
    elif not roles.index("signal_peptide") == 0 < roles.index("antigen_slot") < roles.index("transmembrane") < roles.index("transcription_factor"):
        errors.append(flag("scaffold_module_order", "block"))
    if scaffold.get("signal_policy") != "scaffold_signal_only":
        errors.append(flag("signal_policy_unresolved", "block"))
    if scaffold.get("junctions_verified") is not True:
        errors.append(flag("junctions_unverified", "block"))
    for module in modules:
        if module.get("role") == "antigen_slot":
            if module.get("sequence"):
                errors.append(flag("antigen_slot_must_be_empty", "block"))
            continue
        try:
            sequence(module.get("sequence", ""))
        except ValueError as exc:
            errors.append(flag("scaffold_sequence_invalid", "block", str(exc)))
        if not module.get("name") or not module.get("source"):
            errors.append(flag("scaffold_module_source_missing", "block"))
    return errors


def propose(protein, scaffold=None, rules=None, reviews=()):
    rules = dict(DEFAULT_RULES, **(rules or {}))
    base_flags = []
    try:
        seq = sequence(protein.get("sequence", ""))
    except ValueError as exc:
        return [], [flag("reference_sequence_invalid", "block", str(exc))]
    if protein.get("taxon_id") != 9606:
        base_flags.append(flag("not_human", "block"))
    if not protein.get("accession") or not protein.get("isoform") or protein.get("isoform_ambiguous"):
        base_flags.append(flag("isoform_unresolved", "block"))
    if not sourced(protein):
        base_flags.append(flag("reference_provenance_missing", "block"))
    topology = protein.get("topology", "unknown")
    if protein.get("location") != "plasma_membrane" and topology != "secreted":
        base_flags.append(flag("plasma_membrane_not_confirmed", "review", protein.get("location", "unknown")))
    if topology == "secreted":
        base_flags.append(flag("secreted_antigen_scope_extension", detail="Tethered antigen candidate, not a native cell-surface target"))
        if protein.get("location") != "secreted":
            base_flags.append(flag("secreted_location_unconfirmed", "block"))
    if topology == "type_ii":
        base_flags.append(flag("type_ii_attachment_orientation_change", detail="Native N-terminal attachment becomes C-terminal attachment in this receiver"))
    if topology == "gpi":
        base_flags.append(flag("gpi_anchor_replaced", detail="Verify mature cleavage boundary and membrane presentation"))
    if topology == "unknown":
        base_flags.append(flag("topology_unknown", "block"))
    if topology == "intracellular":
        return [], [flag("outside_extracellular_antigen_scope", "block", "Scope decision, not proof that all artificial display is impossible")]
    features = []
    for f in protein.get("features", []):
        try:
            interval(f, len(seq))
            if not sourced(f):
                raise ValueError("Feature lacks source/version/evidence kind")
            features.append(f)
        except ValueError as exc:
            if f.get("kind") in INFORMATIONAL_KINDS:
                # Informational annotations (variant/mutagenesis/site/...) never
                # block; a coordinate problem there is recorded, not fatal. In
                # batch screening these are deferred earlier by partition_features.
                base_flags.append(flag("informational_annotation_deferred", "review",
                                       f"{f.get('kind')}:{f.get('name', '')} — {exc}"))
            else:
                base_flags.append(flag("feature_annotation_invalid", "block", str(exc)))
    extracellular = [f for f in features if f["kind"] == "extracellular"]
    tm = [f for f in features if f["kind"] == "transmembrane"]
    chains = [f for f in features if f["kind"] == "chain"]
    omega = [f for f in features if f["kind"] == "gpi_attachment_site"]
    if topology == "secreted" and tm:
        base_flags.append(flag("secreted_topology_conflict", "block"))
    if topology in {"type_i", "type_ii"}:
        if len(tm) != 1:
            base_flags.append(flag("single_pass_topology_incomplete_or_conflicting", "block"))
        elif any((f["end"] >= tm[0]["start"] if topology == "type_i" else f["start"] <= tm[0]["end"]) for f in extracellular):
            base_flags.append(flag("topology_orientation_conflict", "block"))
    if topology == "multi_pass" and len(tm) < 2:
        base_flags.append(flag("multipass_topology_incomplete", "block"))
    if topology == "gpi" and (len(omega) != 1 or omega[0]["start"] != omega[0]["end"]):
        base_flags.append(flag("gpi_processing_boundary_unconfirmed"))
    if topology == "type_i" and not any(f["kind"] == "signal_peptide" for f in features):
        base_flags.append(flag("signal_peptide_annotation_missing"))
    if len(extracellular) > 1 or topology == "multi_pass":
        base_flags.append(flag("discontinuous_or_multipass", "review", "No concatenation; individual sourced domains only, native conformation not guaranteed"))
    regions = []
    if len(extracellular) == 1 and topology not in {"multi_pass", "secreted"}:
        regions.append(dict(extracellular[0], antigen_form_type="mature_gpi" if topology == "gpi" else "full_ecd", rationale="Complete annotated continuous extracellular region", origin="baseline"))
    if topology == "secreted":
        for chain in chains:
            regions.append(dict(chain, antigen_form_type="mature_secreted", rationale="Annotated processed chain; autonomy and tethered presentation not established", origin="mature_chain"))
        if not chains:
            base_flags.append(flag("mature_chain_boundary_missing"))
        elif len(chains) > 1:
            base_flags.append(flag("multiple_processed_chains_dependency_unknown"))
    if topology == "gpi" and not extracellular and len(omega) == 1:
        # Keep the omega residue; the removed C-terminal signal begins at omega+1.
        for chain in chains:
            if chain["end"] == omega[0]["start"] == omega[0]["end"]:
                extracellular.append(chain)
                regions.append(dict(chain, antigen_form_type="mature_gpi", rationale="Annotated mature chain ending at the sourced omega residue (included)", origin="mature_chain"))
    if protein.get("candidate_policy", "full_ecd_only") == "full_ecd_and_domains":
        for domain in (f for f in features if f["kind"] == "domain"):
            if any(contains(interval(e, len(seq)), interval(domain, len(seq))) for e in extracellular):
                regions.append(dict(domain, antigen_form_type="domain_fragment", rationale="Complete sourced domain; independent folding and epitope coverage require review", origin="domain_alternative"))
    for region in protein.get("candidate_regions", []):
        regions.append(dict(region, origin="proposed"))
    if not regions:
        return [], base_flags + [flag("no_supported_continuous_candidate", "block", "Needs curated domain candidate or a different presentation approach")]
    candidates = []
    seen = set()
    for region in regions:
        risks = list(base_flags)
        try:
            start, end = interval(region, len(seq))
        except ValueError as exc:
            candidates.append(_invalid(protein, region, risks + [flag("candidate_boundary_invalid", "block", str(exc))]))
            continue
        form = region.get("antigen_form_type", "domain_fragment" if region["origin"] == "proposed" else "full_ecd")
        if (start, end, form) in seen:
            continue
        seen.add((start, end, form))
        bounds = (start, end)
        fragment = seq[start-1:end]
        if region["origin"] == "domain_alternative":
            risks.append(flag("domain_fragment_independence_unproven"))
        if not sourced(region) or not region.get("rationale"):
            risks.append(flag("candidate_support_missing", "block"))
        permitted = chains if topology == "secreted" else extracellular
        if not any(contains(interval(f, len(seq)), bounds) for f in permitted):
            risks.append(flag("candidate_not_within_extracellular_annotation", "block"))
        if form == "mature_secreted" and (topology != "secreted" or not any(interval(c, len(seq)) == bounds for c in chains)):
            risks.append(flag("mature_secreted_form_not_supported", "block"))
        if form == "mature_gpi" and (topology != "gpi" or len(omega) != 1 or end != omega[0]["start"] or omega[0]["start"] != omega[0]["end"]):
            risks.append(flag("mature_gpi_boundary_not_supported", "block"))
        if topology == "gpi" and len(omega) == 1 and end > omega[0]["start"]:
            risks.append(flag("retained_gpi_c_terminal_signal", "block"))
        if form == "shed_product":
            release = region.get("release_evidence", {})
            if not sourced({"evidence": release}) or release.get("kind") == "prediction":
                risks.append(flag("shed_release_evidence_missing", "block", "Motif or guessed cleavage is not a documented released form"))
            if not region.get("release_context"):
                risks.append(flag("shed_release_context_missing", "block"))
            risks.append(flag("shed_product_tethering_unvalidated"))
        if form == "full_ecd" and not any(interval(e, len(seq)) == bounds for e in extracellular):
            risks.append(flag("full_ecd_label_boundary_mismatch", "block"))
        covered, cut, lost, epitopes = [], [], [], []
        epitope_annotations = [f for f in features if f["kind"] == "epitope"]
        for f in features:
            b = interval(f, len(seq))
            k = f["kind"]
            name = f.get("name", f"{k}:{b[0]}-{b[1]}")
            if k in {"transmembrane", "signal_peptide", "gpi_signal", "cytoplasmic", "propeptide"} and overlaps(bounds, b):
                risks.append(flag("retained_" + k, "block", name))
            if k == "domain":
                if contains(bounds, b):
                    covered.append(name)
                elif overlaps(bounds, b):
                    cut.append(name)
                    risks.append(flag("domain_cut", detail=name))
                elif any(overlaps(b, interval(e, len(seq))) for e in permitted):
                    lost.append(name)
                    risks.append(flag("domain_omitted_epitope_scope_changed", detail=name))
            if k == "disulfide" and (start <= b[0] <= end) != (start <= b[1] <= end):
                risks.append(flag("disulfide_partner_removed", detail=name))
            if k == "epitope":
                state = "sequence_retained_not_binding_proven" if contains(bounds, b) else "partially_removed" if overlaps(bounds, b) else "removed"
                epitopes.append({"name": name, "state": state, "antibody_id": f.get("antibody_id", "unknown")})
                if state != "sequence_retained_not_binding_proven":
                    risks.append(flag("known_epitope_loss", detail=name))
        if not epitope_annotations:
            risks.append(flag("epitope_coverage_unknown", detail="No mapped epitopes; full ECD is not proof of preserved recognition"))
        if not covered and not cut:
            risks.append(flag("domain_annotation_missing"))
        if len(fragment) < rules["short_length"]:
            risks.append(flag("short_fragment_geometry_review"))
        if len(fragment) > rules["long_length"]:
            risks.append(flag("long_fragment_geometry_review"))
        if fragment.count("C") / len(fragment) >= rules["cysteine_fraction"]:
            risks.append(flag("cysteine_rich_folding_review"))
        sequons = n_glyco_sequons(fragment, start)
        if sequons:
            risks.append(flag("potential_n_glycosylation", detail="Sequon only; occupancy and effect unknown"))
        if len(sequons) >= rules["dense_glyco_min_sequons"] and \
                len(fragment) / len(sequons) <= rules["dense_glyco_max_residues_per_sequon"]:
            risks.append(flag("dense_glycosylation", "review",
                              f"{len(sequons)} sequons in {len(fragment)} aa (>= {rules['dense_glyco_min_sequons']} and "
                              f"<= {rules['dense_glyco_max_residues_per_sequon']} aa/sequon); heuristic review item, not a threshold"))
        cysteine_count = fragment.count("C")
        if cysteine_count and cysteine_count % 2:
            risks.append(flag("free_thiol_odd_cysteine", "review",
                              f"Odd cysteine count ({cysteine_count}): potential unpaired thiol; review folding/interchain-bond context"))
        multichain = detect_multichain_partners(protein.get("subunit_comments"))
        if multichain:
            risks.append(flag("multichain_partner_required", "review",
                              "SUBUNIT hetero-oligomer annotation: " + " | ".join(multichain)))
        criteria_review, context_risks = assess_context(protein, bounds)
        risks.extend(context_risks)
        missing_criteria = [r["kind"] for r in criteria_review if r["status"] == "missing"]
        incomplete_criteria = any(r["status"] == "missing" or any(v["effective_status"] in {"missing", "unresolved"} for v in r["records"]) for r in criteria_review)
        if missing_criteria:
            risks.append(flag("contextual_risk_evidence_missing", detail=",".join(missing_criteria)))
        # Prediction-evidence review is scoped to candidate-defining features:
        # the candidate's own region plus boundary-defining/domain/disulfide
        # annotations. A predicted glycosylation site or variant record no
        # longer raises this boundary-evidence flag.
        if any(f["evidence"]["kind"] == "prediction" for f in [region] + features
               if sourced(f) and (f is region or f.get("kind") in PREDICTION_RELEVANT_KINDS)):
            risks.append(flag("prediction_requires_annotation_review"))
        profile = molecular_profile(bounds, fragment, features, len(seq))
        boundary = boundary_analysis(bounds, region, features, len(seq), topology,
                                     tolerance=rules["boundary_gap_tolerance"])
        for side in ("n_side", "c_side"):
            code = boundary[side].get("reason_code")
            if not code:
                continue
            adjacent = (boundary[side].get("adjacent_feature") or {})
            # Severity stays "review" here; screening maps
            # boundary_feature_conflict to a conditional reason and the
            # tolerance codes to informational reason codes.
            risks.append(flag(code, "review",
                              f"{side}: defining {boundary['defining_feature']['kind']} "
                              f"{boundary['defining_feature'].get('start')}-{boundary['defining_feature'].get('end')} vs adjacent "
                              f"{adjacent.get('kind', 'none')} {adjacent.get('start', '?')}-{adjacent.get('end', '?')} "
                              f"({boundary[side]['status']})"))
        candidate = {
            "protein_id": protein["protein_id"], "gene": protein.get("gene", ""),
            "accession": protein["accession"], "isoform": protein.get("isoform", ""),
            "reference_sha256": digest(seq), "topology": topology,
            "antigen_form_type": form, "criteria_review": criteria_review,
            "risk_review_status": "incomplete" if incomplete_criteria else "records_present_requires_context_review",
            "release_evidence": region.get("release_evidence", {}), "release_context": region.get("release_context", ""),
            "start": start, "end": end, "coordinate_system": "1-based-inclusive",
            "length": len(fragment), "sequence": fragment, "origin": region["origin"],
            "rationale": region.get("rationale", ""), "boundary_evidence": region.get("evidence", {}),
            "boundary_analysis": boundary, "molecular_profile": profile,
            "multichain_partners": multichain,
            "domains_retained": covered, "domains_cut": cut, "domains_omitted": lost,
            "epitope_review": epitopes, "potential_glycosylation_sites": sequons,
            "risks": risks, "scaffold_issues": scaffold_check(scaffold),
            "scaffold_id": (scaffold or {}).get("scaffold_id", ""),
            "scaffold_version": (scaffold or {}).get("version", ""),
            "functional_status": "not_experimentally_validated",
        }
        candidate["candidate_id"] = "cand-" + digest([protein["protein_id"], seq, start, end, form])[:16]
        candidate["review_key"] = digest({"protein": protein, "region": region, "scaffold": scaffold, "rules": rules})
        review = next((r for r in reviews if r.get("review_key") == candidate["review_key"]), None)
        required_ack = {r["code"] for r in risks if r["severity"] == "review"}
        approved = bool(review and review.get("decision") == "approve" and review.get("reviewer") and review.get("rationale") and review.get("reviewed_at") and required_ack <= set(review.get("acknowledged_risks", [])))
        blocked = any(r["severity"] == "block" for r in risks)
        candidate["design_status"] = "blocked" if blocked else "approved_for_assembly" if approved else "needs_review"
        candidate["assembly_status"] = "blocked" if candidate["scaffold_issues"] or blocked else "awaiting_review" if not approved else "assembled_unvalidated"
        candidate["fusion_sequence"] = ""
        candidate["connection_map"] = []
        if candidate["assembly_status"] == "assembled_unvalidated":
            offset = 1
            for module in scaffold["modules"]:
                aa = fragment if module["role"] == "antigen_slot" else sequence(module["sequence"])
                candidate["fusion_sequence"] += aa
                candidate["connection_map"].append({"role": module["role"], "name": module.get("name", "antigen"), "start": offset, "end": offset + len(aa) - 1})
                offset += len(aa)
        candidate["synthetic_only"] = protein.get("evidence", {}).get("kind") == "synthetic_fixture" or (scaffold or {}).get("status") == "synthetic_fixture"
        candidates.append(candidate)
    return candidates, base_flags


def _invalid(protein, region, risks):
    return {"candidate_id": "invalid-" + digest([protein["protein_id"], region])[:16],
            "protein_id": protein["protein_id"], "start": region.get("start"), "end": region.get("end"),
            "risks": risks, "design_status": "blocked", "assembly_status": "blocked",
            "sequence": "", "fusion_sequence": "", "review_key": "", "scaffold_issues": []}


def audit_construct(record, proteins):
    issues = []
    p = proteins.get(record.get("protein_id"))
    if not p:
        return dict(record, audit_status="blocked", audit_issues=["protein_unresolved"])
    try:
        seq = sequence(record.get("antigen_sequence", ""))
        bounds = interval(record, len(sequence(p["sequence"])))
        if sequence(p["sequence"])[bounds[0]-1:bounds[1]] != seq:
            issues.append("actual_sequence_differs_from_reference_interval")
        if record.get("antigen_cds") and translate(record["antigen_cds"]) != seq:
            issues.append("antigen_cds_translation_mismatch")
        if record.get("fusion_cds") and translate(record["fusion_cds"]) != sequence(record.get("fusion_sequence", "")):
            issues.append("fusion_cds_translation_mismatch")
    except ValueError as exc:
        issues.append(str(exc))
    if not record.get("scaffold_id") or not record.get("scaffold_version") or not record.get("fusion_sequence"):
        issues.append("complete_construct_identity_missing")
    if not record.get("source_ref"):
        issues.append("construct_source_missing")
    if record.get("fusion_sequence"):
        try:
            sequence(record["fusion_sequence"])
        except ValueError as exc:
            issues.append(str(exc))
    return dict(record, audit_status="requires_review" if issues else "sequence_consistent_not_functionally_validated", audit_issues=issues)
