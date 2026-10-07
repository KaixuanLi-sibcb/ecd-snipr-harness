"""Name-blind, uncalibrated mechanism-priority triage for antigen receivers.

Positive annotations and retained coordinates are distinguished from native
whole-protein context, sequence heuristics and measured receptor function.
"""
import re
import hashlib
import json
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .common import digest, file_hash, interval, now, read_json, sequence, sourced, tsv, write_json
from .sequence_tools import PARAMETERS, descriptors, sequence_hash, validate_prediction
from .annotation_scope import association_subject, ligand_subject, statements
from .core_evidence import check_core, export_core

POLICY = {"version": "2.0", "sequence_parameters": PARAMETERS,
          "comment_subject_policy": "explicit_reference_subject_or_unresolved",
          "core_evidence_policy": "exact_fragment_domain_and_topology_crosscheck_v1",
          "tiers": ["elevated_review_priority", "context_dependent_review", "sequence_alert_only", "no_specific_signal_detected", "not_assessed"],
          "identity_features": False, "lab_outcomes_used": False, "probability": False,
          "negative_annotation_means_safe": False, "sequence_alert_can_establish_high_risk": False}
POLICY_SHA256 = digest(POLICY)
REFERENCE_READ_TIMEOUT = 5
LIMIT = "Mechanism-priority triage is an uncalibrated hypothesis, not a prediction of SNIPR self-activation or functional success."
SELF_ASSOC = re.compile(r'\b(?:homodimer\w*|homotrimer\w*|homotetramer\w*|homo[- ]?oligomer\w*|self[- ]associat\w*|homophilic\w*|self[- ]ligand)\b', re.I)
SHEDDING = re.compile(r'\b(?:shedd\w*|shed|ectodomain cleavage|cleaved by ADAM\d*)\b', re.I)
ENDOGENOUS_BINDING = re.compile(r'\b(?:sialic acid|sialylated|glycan|immunoglobulin|IgG|IgM|Fc region)\b', re.I)
BINDING = re.compile(r'\b(?:bind\w*|recogniz\w*|receptor for|ligand)\b', re.I)
NEGATION = re.compile(r'\b(?:not|no|cannot|unable|lack\w*|inhibit\w*|prevent\w*)\b', re.I)


def read_reference(path, expected_hash):
    # A file-provider read can block indefinitely despite the path existing.
    # A bounded child isolates that failure without unkillable reader threads.
    try:
        process = subprocess.run([sys.executable, '-c',
            'from pathlib import Path; import sys; sys.stdout.buffer.write(Path(sys.argv[1]).read_bytes())', str(path)],
            capture_output=True, timeout=REFERENCE_READ_TIMEOUT)
    except subprocess.TimeoutExpired as exc:
        raise ValueError('reference_read_timeout') from exc
    if process.returncode:
        raise ValueError('reference_read_error')
    if hashlib.sha256(process.stdout).hexdigest() != expected_hash:
        raise ValueError('Reference path/hash conflict')
    return json.loads(process.stdout)


def assess_risk(candidate, protein=None, predictions=(), use_biopython=False, interpro=None):
    """Do not inspect gene/accession or outcome labels in the decision logic."""
    result = {"policy": POLICY, "policy_sha256": POLICY_SHA256, "claim_limit": LIMIT,
        "signals": [], "excluded_or_unresolved_evidence": [], "tool_evidence": [],
        "functional_confidence": "not_established", "functional_risk": "undetermined_not_low_risk",
        "success_probability": None, "source_scope": "annotation_and_sequence_only_no_lab_labels",
        "missing_information": ["actual_receiver_backbone_and_junction", "host_and_expression_regime",
            "matched_unstimulated_controls", "validated_endpoint_acceptance_criteria"]}
    if not candidate.get('sequence') or candidate.get('design_status') == 'blocked':
        return dict(result, review_priority='not_assessed', assessment_status='blocked_design_no_functional_transfer', reason_codes=['candidate_not_usable'])
    try:
        seq = sequence(candidate['sequence'])
    except ValueError as exc:
        return dict(result, review_priority='not_assessed', assessment_status='invalid_sequence', reason_codes=[str(exc)])
    tools = descriptors(seq, use_biopython)
    result['sequence_analysis'] = tools
    result['core_evidence'] = check_core(candidate, protein, interpro, predictions)
    full = candidate.get('antigen_form_type') in {'full_ecd','mature_gpi','mature_secreted'}
    validated_reference = False
    if protein is not None:
        try:
            reference = sequence(protein['sequence'])
            start, end = interval(candidate, len(reference))
            if reference[start-1:end] != seq or candidate.get('reference_sha256') != digest(reference):
                raise ValueError('Candidate is not exact fragment of this reference')
            validated_reference = sourced(protein)
            if not validated_reference:
                raise ValueError('Reference source/version/kind missing')
        except (ValueError, KeyError, TypeError) as exc:
            return dict(result, review_priority='not_assessed', assessment_status='reference_conflict', reason_codes=[str(exc)])
    if not validated_reference:
        result['missing_information'].append('source_bound_native_annotations_not_supplied')

    def signal(code, tier, endpoints, evidence, note, action):
        result['signals'].append({"code": code, "priority": tier, "endpoints_to_check": endpoints,
            "evidence": evidence, "interpretation": note, "next_action": action,
            "receiver_outcome_inference": "unvalidated"})

    if validated_reference:
        for index, f in enumerate(protein.get('features', [])):
            try:
                a,b = interval(f,len(reference))
            except ValueError:
                continue
            if not (start <= a <= b <= end) or not sourced(f):
                continue
            name = f.get('name','')
            ev = dict(f.get('evidence',{}), reference_interval=[a,b], feature_index=index, text=name,
                      applicability='annotated_feature_retained_in_fragment')
            if f.get('kind') == 'disulfide' and re.search(r'\binterchain\b', name, re.I):
                signal('retained_interchain_disulfide', 'elevated_review_priority', ['surface_expression','basal_activation'], ev,
                    'Native covalent chain association is retained; it does not establish fusion clustering or self-activation.',
                    'Determine partner and mature-chain identity; compare surface expression and matched baseline without changing epitopes blindly.')
            if f.get('kind') in {'site','region'} and re.search(r'\b(?:cleavage|cleaved|shedd\w*)\b',name,re.I):
                if NEGATION.search(name):
                    result['excluded_or_unresolved_evidence'].append(dict(ev, reason='negated_or_ambiguous_feature_text'))
                else:
                    signal('retained_annotated_processing_site', 'elevated_review_priority', ['surface_expression','basal_activation','recognition_retention'],ev,
                        'A native processing annotation intersects the proposed fragment, not proof the new receiver is cleaved.',
                        'Check protease and junction context; compare retained receptor/fragment and unstimulated reporter.')
        for field in ('ptm_comments','subunit_comments','function_comments'):
            for index, comment in enumerate(protein.get(field, [])):
                clauses = ((s, part) for s in statements(comment) for part in re.split(r'[.;]\s+', s.get('text','')))
                for statement, clause in clauses:
                    if not (SELF_ASSOC.search(clause) or SHEDDING.search(clause) or (ENDOGENOUS_BINDING.search(clause) and BINDING.search(clause))):
                        continue
                    ev = dict(protein['evidence'], field=field, comment_index=index, text=clause,
                              eco=statement.get('eco',[]), statement_evidences=statement.get('evidences',[]),
                              applicability='native_reference_context_site_not_localized', fragment_is_full_topological_interval=full)
                    if NEGATION.search(clause):
                        result['excluded_or_unresolved_evidence'].append(dict(ev, reason='negated_or_mixed_clause_requires_manual_interpretation'))
                        continue
                    if SELF_ASSOC.search(clause):
                        subject = association_subject(clause)
                        if subject != 'reference_subject':
                            result['excluded_or_unresolved_evidence'].append(dict(ev, reason=subject))
                            continue
                        intracellular = re.search(r'\b(?:cytoplasm\w*|intracellular|transmembrane|kinase domain)\b',clause,re.I)
                        if intracellular and not re.search(r'\bextracellular\b',clause,re.I):
                            result['excluded_or_unresolved_evidence'].append(dict(ev, reason='native_association_not_localized_to_retained_ecd'))
                        else:
                            explicit_external = bool(re.search(r'\b(?:homophilic\w*|self[- ]ligand)\b|\bhomo\w*.{0,60}(?:via|through|mediated by).{0,40}extracellular|extracellular domain.{0,30}mediates.{0,30}homo',clause,re.I))
                            tier = 'elevated_review_priority' if full and explicit_external else 'context_dependent_review'
                            signal('native_self_association_context', tier, ['basal_activation','recognition_retention'],ev,
                                'Native self-association/adhesion is a receptor-context question, not demonstrated autonomous signaling.',
                                'Test receiver-only and nonbinding-sender baseline with matched receptor abundance; establish relevant interface and orientation.')
                    if SHEDDING.search(clause) and field == 'ptm_comments':
                        signal('native_shedding_context', 'context_dependent_review', ['surface_expression','basal_activation','recognition_retention'],ev,
                            'Native shedding is reported, but the cleavage position/protease applicability to this fragment is not established by this text.',
                            'Localize the processing site and assess retention plus actual receiver junction; measure intact surface receptor and baseline.')
                    if ENDOGENOUS_BINDING.search(clause) and BINDING.search(clause):
                        if ligand_subject(clause) != 'reference_subject':
                            result['excluded_or_unresolved_evidence'].append(dict(ev, reason='ligand_binding_subject_unresolved'))
                            continue
                        signal('endogenous_or_medium_ligand_context', 'context_dependent_review', ['basal_activation','recognition_retention'],ev,
                            'Native glycan/immunoglobulin binding may complicate cell/medium controls; ligand availability and force coupling are unknown.',
                            'Check actual host, sender, tags and medium; include nonbinding/ligand-context controls, not a gene-specific blacklist.')

    for risk in candidate.get('risks', []):
        if risk.get('code') in {'domain_cut','disulfide_partner_removed','processed_chain_segments_spanned'}:
            signal('structural_unit_disruption', 'elevated_review_priority', ['surface_expression','recognition_retention'],
                   {'kind':'candidate_annotation_check','record':risk,'source':candidate.get('boundary_evidence',{})},
                   'Specific annotated structural/processing unit is interrupted; generic motif counts are not used.',
                   'Review an intact-unit alternative without assuming better functional response or preserved epitope repertoire.')
    if tools['hydrophobic_segments']:
        signal('local_hydrophobic_segment', 'sequence_alert_only', ['surface_expression'],
               {'kind':'computed_descriptor','algorithm':tools['algorithm'],'parameters':PARAMETERS,'segments':tools['hydrophobic_segments']},
               'Hydropathy alone cannot distinguish buried cores, transmembrane helices or aggregation.',
               'Check annotation/structure and optional DeepTMHMM on the exact sequence; do not automatically delete the segment.')
    if tools['stp_enriched_segments']:
        signal('stp_enriched_stalk_geometry', 'sequence_alert_only', ['recognition_retention','induced_response'],
               {'kind':'computed_descriptor','segments':tools['stp_enriched_segments'],'parameters':PARAMETERS},
               'A serine/threonine/proline-enriched interval suggests a stalk/compositional context, not measured O-glycosylation or self-activation.',
               'Check annotated O-glycosylation/mucin regions and actual attachment geometry; retain unknown epitope effects.')
    if tools['low_complexity_segments']:
        signal('low_complexity_sequence_context', 'sequence_alert_only', ['recognition_retention','induced_response'],
               {'kind':'computed_descriptor','segments':tools['low_complexity_segments'],'parameters':PARAMETERS},
               'Window entropy detects compositional bias, not disorder, condensate formation or SNIPR activation.',
               'Use a local disorder predictor or structural/domain evidence if this interval affects the connection; do not truncate merely for complexity.')
    for check in result['core_evidence']['checks']:
        if check['code'] in {'domain_boundary_cut','retained_native_exclusion_region'}:
            signal('core_' + check['code'], 'context_dependent_review', ['surface_expression','recognition_retention'],check,
                   'Exact fragment/domain or exclusion-region conflict needs review; boundaries are not silently replaced.',
                   'Reconcile source/sequence/version and preserve epitope tradeoffs before assembly.')
    for prediction in predictions:
        try:
            record = validate_prediction(prediction,candidate,protein)
        except (ValueError,KeyError,TypeError) as exc:
            result['excluded_or_unresolved_evidence'].append({'reason':'rejected_tool_evidence','detail':str(exc),'record':prediction})
            continue
        result['tool_evidence'].append(record)
        if record['software'] in {'DeepTMHMM','DeepTMHMM2','SignalP'} and any(r['kind'] in {'transmembrane','signal_peptide','inside','cytoplasmic','reentrant','interfacial','transit_peptide'} for r in record['regions']):
            signal('predicted_membrane_or_signal_conflict', 'context_dependent_review', ['surface_expression'],record,
                   'Predicted regions in the proposed fragment conflict with a soluble extracellular route; prediction is not experimental evidence.',
                   'Reconcile reference topology and sequence before assembly; do not silently replace curated boundaries.')
        elif record.get('scores') and record['software'] in {'IUPred2A','IUPred3','AIUPred'}:
            result.setdefault('disorder_descriptors',[]).append({'software':record['software'],'version':record['version'],
                'fraction_above_0_5':sum(v>0.5 for v in record['scores'])/len(seq),
                'interpretation':'prediction_only_not_functional_risk_threshold'})
    tiers = [s['priority'] for s in result['signals']]
    priority = next((p for p in POLICY['tiers'][:3] if p in tiers),'no_specific_signal_detected')
    return dict(result, review_priority=priority,
                assessment_status='annotation_and_sequence_evaluated' if validated_reference else 'sequence_only_annotation_coverage_missing',
                reason_codes=sorted({s['code'] for s in result['signals']}),
                rule_context_limit='Native associations and processing are not transferred as measured receiver outcomes.')


def export_risk(outdir,candidates):
    outdir = Path(outdir)
    rows, signals = [], []
    usable = [c for c in candidates if c.get('sequence') and c.get('design_status') != 'blocked']
    for c in candidates:
        r = c.get('receiver_risk') or assess_risk(c)
        c['receiver_risk'] = r
        metadata = {k:c.get(k) for k in ('candidate_id','protein_id','gene','accession','isoform','start','end','is_primary','screening_recommendation','design_status')}
        rows.append(dict(metadata, review_priority=r['review_priority'], reason_codes=r['reason_codes'],
            functional_confidence=r['functional_confidence'], assessment_status=r['assessment_status'],
            missing_information=r['missing_information'], signals=r['signals']))
        signals.extend(dict(candidate_id=c['candidate_id'], **s) for s in r['signals'])
    tsv(outdir/'receiver_risk.tsv',rows,[*metadata,'review_priority','reason_codes','functional_confidence','assessment_status','missing_information','signals'] if rows else ['candidate_id','review_priority'])
    tsv(outdir/'receiver_risk_signals.tsv',signals,['candidate_id','code','priority','endpoints_to_check','evidence','interpretation','next_action','receiver_outcome_inference'])
    priorities = dict(Counter(c['receiver_risk']['review_priority'] for c in usable))
    primary = [c for c in usable if c.get('is_primary')]
    cross = defaultdict(Counter)
    for c in usable:
        cross[str(c.get('screening_recommendation','unclassified'))][c['receiver_risk']['review_priority']] += 1
    summary = {'policy':POLICY,'policy_sha256':POLICY_SHA256,'candidate_records':len(candidates),
        'usable_candidate_denominator':len(usable),'primary_candidate_denominator':len(primary),
        'priorities':priorities,'primary_priorities':dict(Counter(c['receiver_risk']['review_priority'] for c in primary)),
        'design_by_priority':{k:dict(v) for k,v in cross.items()},
        'signal_candidate_counts':dict(Counter(code for c in usable for code in c['receiver_risk']['reason_codes'])),
        'annotation_assessment_states':dict(Counter(c['receiver_risk']['assessment_status'] for c in usable)),
        'claim_limit':LIMIT,'biological_accuracy':'not_estimated','lab_outcome_labels_used':False}
    summary['core_evidence'] = export_core(outdir,candidates)
    write_json(outdir/'receiver_risk_summary.json',summary)
    return summary


def audit_risk(run_dir, reference_set, outdir, tool_evidence=None, resume=False, interpro_cache=None, interpro_live=False):
    from . import __version__
    from .harness import engine_hash, verify_bundle
    source, root, setdir = Path(run_dir).resolve(), Path(outdir).resolve(), Path(reference_set).resolve()
    if not verify_bundle(source):
        raise ValueError('Source bundle integrity failed')
    if root == source or source in root.parents or root == setdir or setdir in root.parents:
        raise ValueError('Risk output must be outside frozen source run and reference set')
    definition = read_json(setdir/'set_definition.json')
    source_hash, set_hash = file_hash(source/'manifest.json'), file_hash(setdir/'set_definition.json')
    predictions = read_json(tool_evidence) if tool_evidence else []
    if not isinstance(predictions,list):
        raise ValueError('Tool evidence must be a list of normalized records')
    candidates = read_json(source/'candidates.json')
    ids = {c['candidate_id'] for c in candidates}
    accessions = {c.get('accession') or c.get('protein_id') for c in candidates}
    if any((p.get('accession') not in accessions if p.get('input_scope') == 'full_reference' else p.get('candidate_id') not in ids) for p in predictions):
        raise ValueError('Unknown candidate/reference in tool evidence; no silently ignored records')
    for p in predictions:
        raw_path = p.get('raw_output_path')
        if not raw_path or file_hash(raw_path) != p.get('raw_output_sha256'):
            raise ValueError('Local predictor raw output file/hash required; provenance not verified')
        if p.get('raw_input_path') and file_hash(p['raw_input_path']) != p.get('raw_input_sha256'):
            raise ValueError('Predictor input FASTA changed after import')
    inputs = {'source_manifest_sha256':source_hash,'set_definition_sha256':set_hash,'engine_sha256':engine_hash(),
              'policy_sha256':POLICY_SHA256,'tool_evidence_sha256':file_hash(tool_evidence) if tool_evidence else None}
    key = digest(inputs)
    dest = root/key[:20]
    entries = {e.get('accession'):e for e in definition['entries'] if e.get('protein_file')}
    proteins, errors = {}, []
    def load(pid):
        entry = entries.get(pid)
        try:
            if not entry:
                raise ValueError('Reference missing from supplied set')
            path = (setdir/entry['protein_file']).resolve()
            if not path.is_relative_to(setdir):
                raise ValueError('Reference path/hash conflict')
            return pid,read_reference(path,entry['protein_sha256']),None
        except (ValueError,OSError,KeyError) as exc:
            return pid,None,{'accession':pid,'reason':str(exc)}
    reference_ids = list(dict.fromkeys(c.get('accession') or c.get('protein_id') for c in candidates))
    # Bound both process count and each file read; order remains deterministic.
    with ThreadPoolExecutor(max_workers=4) as pool:
        for pid,p,error in pool.map(load,reference_ids):
            proteins[pid] = p
            if error:
                errors.append(error)
    interpro_records = {}
    if interpro_live and not interpro_cache:
        raise ValueError('InterPro live retrieval needs an explicit cache directory')
    if interpro_cache:
        from .interpro import fetch_domains
        for pid,p in proteins.items():
            if p is None:
                continue
            # Only public human canonical annotations leave the host, never sequences.
            if p.get('taxon_id') != 9606 or p.get('evidence',{}).get('kind') == 'synthetic_fixture':
                interpro_records[pid] = {'status':'missing','reason':'not_a_public_human_reference'}
            else:
                interpro_records[pid] = fetch_domains(pid,interpro_cache,offline=not interpro_live)
    inputs['interpro_records_sha256'] = digest({pid:{k:v for k,v in r.items() if k not in {'snapshot','retrieved_at','status','raw_snapshot_path'}} for pid,r in interpro_records.items()})
    key = digest(inputs)
    dest = root/key[:20]
    for c in candidates:
        pid = c.get('accession') or c.get('protein_id')
        try:
            selected = [p for p in predictions if (p.get('accession') == pid if p.get('input_scope') == 'full_reference' else p.get('candidate_id') == c['candidate_id'])]
            c['receiver_risk'] = assess_risk(c,proteins[pid],selected,interpro=interpro_records.get(pid))
        except (ValueError,KeyError,TypeError,AttributeError) as exc:
            errors.append({'candidate_id':c['candidate_id'],'reason':str(exc)})
            c['receiver_risk'] = dict(assess_risk(c),review_priority='not_assessed',assessment_status='technical_assessment_error',reason_codes=['processing_error'])
    if file_hash(source/'manifest.json') != source_hash or not verify_bundle(source) or file_hash(setdir/'set_definition.json') != set_hash:
        raise ValueError('Source changed during risk audit')
    with ThreadPoolExecutor(max_workers=4) as pool:
        for pid,p,error in pool.map(load,[pid for pid,p in proteins.items() if p is not None]):
            if error:
                raise ValueError('Reference changed or became unavailable during audit: '+pid)
    for record in interpro_records.values():
        if record.get('raw_snapshot_path') and file_hash(record['raw_snapshot_path']) != record['raw_response_sha256']:
            raise ValueError('InterPro snapshot changed during audit')
    for prediction in predictions:
        if file_hash(prediction['raw_output_path']) != prediction['raw_output_sha256']:
            raise ValueError('Predictor output changed during audit')
        if prediction.get('raw_input_path') and file_hash(prediction['raw_input_path']) != prediction['raw_input_sha256']:
            raise ValueError('Predictor input changed during audit')
    if dest.exists():
        if resume and not errors and verify_bundle(dest):
            return {'run_dir':str(dest),'execution':'verified_cache_hit','summary':read_json(dest/'receiver_risk_summary.json')}
        raise FileExistsError('Risk bundle exists/corrupt or reference changed; preserved. Choose a new output root')
    dest.mkdir(parents=True)
    summary = export_risk(dest,candidates)
    enrichment_errors = [dict(accession=pid,reason=r.get('reason','source_unavailable')) for pid,r in interpro_records.items() if r['status'] == 'error']
    enrichment_conflicts = [c['candidate_id'] for c in candidates if c.get('receiver_risk',{}).get('core_evidence',{}).get('coverage',{}).get('interpro') == 'conflict']
    tool_conflicts = [c['candidate_id'] for c in candidates if 'rejected' in c.get('receiver_risk',{}).get('core_evidence',{}).get('coverage',{}).values()]
    summary.update(reference_errors=errors, enrichment_errors=enrichment_errors, enrichment_conflicts=enrichment_conflicts,
        tool_evidence_conflicts=tool_conflicts,
        completeness='partial' if errors or enrichment_errors or enrichment_conflicts or tool_conflicts or definition.get('completeness') != 'complete' else 'complete',
        enrichment_coverage='partial' if any(r['status'] not in {'cached','fetched'} for r in interpro_records.values()) else 'complete' if interpro_records else 'not_run')
    write_json(dest/'receiver_risk_summary.json',summary)
    write_json(dest/'candidates.json',candidates)
    write_json(dest/'core_source_records.json',{pid:{k:v for k,v in r.items() if k != 'snapshot'} for pid,r in interpro_records.items()})
    write_json(dest/'risk_policy_snapshot.json',{'policy':POLICY,'implementation_sha256':inputs['engine_sha256'],
        'freeze_basis':'rules applied before importing laboratory outcome workbook','prior_cases_seen_during_development':True,
        'retrospective_evaluation_is_not_held_out':True})
    write_json(dest/'source_link.json',dict(inputs,source_run=str(source),reference_set=str(setdir),
        design_rules_changed=False,lab_outcomes_read=False,reference_errors=errors))
    artifacts = {p.name:file_hash(p) for p in sorted(dest.iterdir()) if p.is_file()}
    write_json(dest/'manifest.json',dict(inputs,state='complete',analysis_completeness=summary['completeness'],version=__version__,
        run_key=key,finished_at=now(),artifacts=artifacts))
    return {'run_dir':str(dest),'execution':'computed','summary':summary}


def run_sequence_tools(run_dir,outdir,use_biopython=False):
    from .harness import engine_hash,verify_bundle
    source,root = Path(run_dir).resolve(),Path(outdir).resolve()
    if not verify_bundle(source):
        raise ValueError('Source bundle integrity failed')
    if root == source or source in root.parents:
        raise ValueError('Sequence-tool output must be outside frozen source')
    root.mkdir(parents=True,exist_ok=False)
    records = []
    for c in read_json(source/'candidates.json'):
        if c.get('sequence') and c.get('design_status') != 'blocked':
            records.append({'candidate_id':c['candidate_id'],'descriptors':descriptors(c['sequence'],use_biopython)})
    write_json(root/'sequence_descriptors.json',records)
    summary = {'candidates':len(records),'biopython_states':dict(Counter(r['descriptors']['biopython']['status'] for r in records)),
        'prediction_states':{name:'not_run' for name in ('IUPred','DeepTMHMM','SignalP')},
        'source_manifest_sha256':file_hash(source/'manifest.json'),'engine_sha256':engine_hash(),'claim_limit':LIMIT}
    write_json(root/'sequence_tools_summary.json',summary)
    if not verify_bundle(source):
        raise ValueError('Source changed during sequence analysis')
    write_json(root/'manifest.json',dict(summary,state='complete',artifacts={p.name:file_hash(p) for p in root.iterdir() if p.is_file()}))
    return {'run_dir':str(root),'summary':summary}
