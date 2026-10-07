import copy
import json
import tempfile
import unittest
import urllib.error
from pathlib import Path
from unittest.mock import patch

from test_design import feature
from test_receiver_risk import case
from test_screening import Response, make_protein, make_set
from ecd_snipr.annotation_scope import association_subject
from ecd_snipr.common import digest, file_hash, read_json, write_json
from ecd_snipr.core_evidence import check_core
from ecd_snipr.interpro import fetch_domains, normalized_domains
from ecd_snipr.predictor_import import import_prediction
from ecd_snipr.receiver_risk import assess_risk, audit_risk
from ecd_snipr.sequence_tools import sequence_hash, validate_prediction
from ecd_snipr.screening import run_screening
from ecd_snipr.harness import verify_bundle


def entry(accession, start=20, end=40, kind='domain', database='pfam', fragments=None):
    return {'metadata': {'accession':'PF_TEST', 'source_database':database, 'name':'Synthetic domain', 'type':kind},
            'proteins':[{'accession':accession.lower(), 'protein_length':100,
                'entry_protein_locations':[{'fragments':fragments or [{'start':start,'end':end,'dc-status':'CONTINUOUS'}]}]}]}


def interpro_record(p, entries=()):
    snapshot = {'accession':p['accession'], 'sequence':p['sequence'], 'sequence_sha256':sequence_hash(p['sequence']),
                'version':'synthetic-1', 'completeness':'complete', 'pages':[{'results':list(entries)}]}
    return {'status':'cached', 'version':'synthetic-1', 'raw_response_sha256':digest(snapshot), 'snapshot':snapshot}


def prediction(p, **kwargs):
    return dict({'software':'DeepTMHMM', 'version':'synthetic', 'parameters':{}, 'source':'synthetic fixture',
                 'accession':p['accession'], 'sequence_sha256':sequence_hash(p['sequence']), 'input_scope':'full_reference',
                 'coordinate_system':'1-based-inclusive', 'raw_output_sha256':'a'*64, 'regions':[]}, **kwargs)


class ScopeTests(unittest.TestCase):
    def test_partner_homodimer_not_target_association(self):
        clauses = ['Interacts with PARTNER_A homodimer.',
                   'In turn, the hexamer interacts with PARTNER_A/PARTNER_B homodimer to form the complex.',
                   'Forms a complex with a homodimer of ANOTHER_PROTEIN.',
                   'Is part of a complex containing a homodimer of ANOTHER_PROTEIN.']
        for text in clauses:
            p,c=case(); p['subunit_comments']=[{'text':text}]
            r=assess_risk(c,p)
            self.assertNotIn('native_self_association_context',r['reason_codes'])
            self.assertTrue(r['excluded_or_unresolved_evidence'])

    def test_reference_subjects_preserved(self):
        for text in ('Homodimer.', 'Forms a homodimer.', 'The receptor is a homotrimer.', 'Mediates homophilic binding.'):
            self.assertEqual(association_subject(text), 'reference_subject')

    def test_unresolved_subject_not_positive(self):
        p,c=case(); p['function_comments']=[{'text':'Binding to LIGAND leads to homodimerization.'}]
        r=assess_risk(c,p)
        self.assertFalse(r['signals'])
        self.assertEqual(r['excluded_or_unresolved_evidence'][0]['reason'],'subject_unresolved')

    def test_partner_ligand_binding_not_target_binding(self):
        p,c=case(); p['function_comments']=[{'text':'Interacts with a receptor that binds IgG.'}]
        self.assertNotIn('endogenous_or_medium_ligand_context',assess_risk(c,p)['reason_codes'])

    def test_per_statement_citations_kept(self):
        from ecd_snipr.uniprot import normalize
        from test_screening import Raw
        raw=Raw.entry()
        raw['comments'].append({'commentType':'SUBUNIT','texts':[{'value':'Homodimer.','evidences':[{'evidenceCode':'ECO:0000269','source':'PubMed','id':'12345'}]}]})
        p=normalize(raw)
        self.assertEqual(p['subunit_comments'][0]['statements'][0]['evidences'][0]['id'],'12345')
        p,c=case(); p['subunit_comments']=[{'statements':[{'text':'Homodimer.','evidences':[{'source':'PubMed','id':'12345'}]}]}]
        self.assertEqual(assess_risk(c,p)['signals'][0]['evidence']['statement_evidences'][0]['id'],'12345')


class CoreTests(unittest.TestCase):
    def test_no_optional_evidence_is_not_failure(self):
        p,c=case(); r=check_core(c,p)
        self.assertEqual(r['coverage']['interpro'],'not_run')
        self.assertNotIn('missing',r['reason_codes'])
        self.assertEqual(r['status'],'checks_completed_no_conflict')

    def test_input_not_mutated(self):
        p,c=case(); old=copy.deepcopy((p,c)); check_core(c,p)
        self.assertEqual(old,(p,c))

    def test_exact_sequence_mismatch_not_checked(self):
        p,c=case(); c['sequence']='A'*len(c['sequence'])
        self.assertEqual(check_core(c,p)['status'],'reference_conflict')

    def test_retained_cut_and_omitted_domains(self):
        p,c=case(); r=interpro_record(p,[entry(p['accession'],20,40),entry(p['accession'],5,30),entry(p['accession'],80,95)])
        a=check_core(c,p,r)
        observed=[d['retention'] for d in a['domains'] if d['database']=='pfam']
        self.assertEqual(observed,['retained','cut','omitted'])
        self.assertIn('domain_boundary_cut',a['reason_codes'])
        self.assertEqual(c['design_status'],'needs_review')

    def test_family_not_autonomous_domain(self):
        p,c=case(); r=interpro_record(p,[entry(p['accession'],1,100,kind='family')])
        self.assertFalse([d for d in check_core(c,p,r)['domains'] if d['database']=='pfam'])

    def test_discontinuous_domain_not_stitched(self):
        p,c=case(); r=interpro_record(p,[entry(p['accession'],fragments=[{'start':20,'end':30,'dc-status':'N_TERMINAL_DISC'},{'start':80,'end':90,'dc-status':'C_TERMINAL_DISC'}])])
        d=check_core(c,p,r)['domains'][-1]
        self.assertEqual(d['retention'],'cut')
        self.assertTrue(d['discontinuous'])
        self.assertEqual(len(d['fragments']),2)

    def test_interpro_sequence_version_conflict_visible(self):
        p,c=case(); r=interpro_record(p)
        r['snapshot']['sequence']='A'*100
        a=check_core(c,p,r)
        self.assertEqual(a['coverage']['interpro'],'conflict')
        self.assertTrue(a['rejected_evidence'])

    def test_interpro_noncanonical_coordinates_not_transferred(self):
        p,c=case(); p['reference_coordinate_status']='unverified_noncanonical_request'
        self.assertEqual(check_core(c,p,interpro_record(p))['coverage']['interpro'],'conflict')

    def test_duplicate_domain_hits_are_not_votes(self):
        p,c=case(); r=interpro_record(p,[entry(p['accession']),entry(p['accession'],database='interpro')])
        a=check_core(c,p,r)
        self.assertEqual([x for x in a['checks'] if x['code']=='domain_annotations_available'][0]['independent_vote_count'],'not_used')

    def test_retained_native_tm_stays_conflict(self):
        p,c=case(); p['features'].append(feature('transmembrane',30,40))
        self.assertIn('retained_native_exclusion_region',check_core(c,p)['reason_codes'])

    def test_outside_native_tm_not_conflict(self):
        p,c=case(); self.assertNotIn('retained_native_exclusion_region',check_core(c,p)['reason_codes'])

    def test_invalid_secondary_annotation_not_discarded(self):
        p,c=case(); p['features'].append(feature('domain',0,30))
        self.assertEqual(check_core(c,p)['status'],'checks_partial_evidence_unresolved')


class ProjectionTests(unittest.TestCase):
    def test_native_sp_and_tm_outside_fragment_not_conflict(self):
        p,c=case(); r=prediction(p,regions=[{'kind':'signal_peptide','start':1,'end':c['start']-1},{'kind':'transmembrane','start':c['end']+1,'end':90}])
        b=validate_prediction(r,c,p)
        self.assertEqual(b['regions'],[])
        self.assertEqual(len(b['not_retained_regions']),2)
        self.assertNotIn('predicted_membrane_or_signal_conflict',assess_risk(c,p,[r])['reason_codes'])

    def test_tm_partial_overlap_correct_coordinates(self):
        p,c=case(); r=prediction(p,regions=[{'kind':'transmembrane','start':5,'end':15}])
        b=validate_prediction(r,c,p)
        self.assertEqual((b['regions'][0]['start'],b['regions'][0]['end']),(1,5))
        self.assertEqual(b['regions'][0]['reference_interval'],[5,15])
        self.assertIn('predicted_membrane_or_signal_conflict',assess_risk(c,p,[r])['reason_codes'])

    def test_outside_not_claimed_cell_surface(self):
        p,c=case(); r=prediction(p,regions=[{'kind':'outside','start':11,'end':70}])
        a=assess_risk(c,p,[r])
        self.assertNotIn('predicted_membrane_or_signal_conflict',a['reason_codes'])
        self.assertIn('organelle lumen',a['tool_evidence'][0]['interpretation_limit'])

    def test_full_sequence_hash_and_reference_required(self):
        p,c=case(); r=prediction(p)
        with self.assertRaises(ValueError): validate_prediction(r,c)
        r['sequence_sha256']='b'*64
        with self.assertRaises(ValueError): validate_prediction(r,c,p)

    def test_disorder_scores_project_without_function_threshold(self):
        p,c=case(); r=prediction(p,software='IUPred3',scores=[0.1]*(c['start']-1)+[0.8]*len(c['sequence'])+[0.1]*(100-c['end']))
        r.pop('regions')
        b=validate_prediction(r,c,p)
        self.assertEqual(b['scores'],[0.8]*len(c['sequence']))
        a=assess_risk(c,p,[r]); self.assertEqual(a['disorder_descriptors'][0]['fraction_above_0_5'],1)
        self.assertEqual(a['review_priority'],'no_specific_signal_detected')


class ImportTests(unittest.TestCase):
    def setup_files(self, root):
        p,_=case(); p['accession']='P00001'
        fasta=root/'input.fasta'; fasta.write_text('>P00001\n'+p['sequence']+'\n')
        return p,fasta

    def test_iupred_import_exact_residues(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); p,fasta=self.setup_files(root); raw=root/'iupred.txt'
            raw.write_text('# synthetic\n'+''.join(f'{i} {aa} 0.8\n' for i,aa in enumerate(p['sequence'],1)))
            r=import_prediction('IUPred3',raw,fasta,p,'synthetic',{'mode':'long'})
            self.assertEqual(len(r['scores']),100)
            self.assertEqual(r['raw_input_sha256'],file_hash(fasta))
            raw.write_text('1 A 0.8\n')
            with self.assertRaises(ValueError): import_prediction('IUPred3',raw,fasta,p,'synthetic',{})

    def test_deeptmhmm_import_sequence_and_labels(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); p,fasta=self.setup_files(root); raw=root/'prediction.3line'
            raw.write_text('>P00001\n'+p['sequence']+'\n'+'S'*10+'O'*60+'M'*20+'I'*10+'\n')
            r=import_prediction('DeepTMHMM',raw,fasta,p,'synthetic',{})
            self.assertEqual([x['kind'] for x in r['regions']],['signal_peptide','outside','transmembrane','inside'])
            raw.write_text('>P00001\n'+p['sequence']+'\n'+'Z'*100+'\n')
            with self.assertRaises(ValueError): import_prediction('DeepTMHMM',raw,fasta,p,'synthetic',{})

    def test_signalp_positive_gff3_only(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); p,fasta=self.setup_files(root); raw=root/'signalp.gff3'
            raw.write_text('P00001\tSignalP\tsignal_peptide\t1\t10\t.\t+\t.\tID=sp\n')
            r=import_prediction('SignalP',raw,fasta,p,'synthetic',{'organism':'eukarya'})
            self.assertEqual(r['regions'][0]['end'],10)
            raw.write_text('# no explicit per-sequence result\n')
            with self.assertRaises(ValueError): import_prediction('SignalP',raw,fasta,p,'synthetic',{})

    def test_input_fasta_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); p,fasta=self.setup_files(root); fasta.write_text('>P00001\nAAAA\n')
            with self.assertRaises(ValueError): import_prediction('IUPred3',root/'absent',fasta,p,'synthetic',{})


class FetchTests(unittest.TestCase):
    def setup_opener(self, accession='P00001', first_entry=None):
        p=make_protein(accession)
        protein={'metadata':{'accession':accession,'sequence':p['sequence'],'length':100,'source_organism':{'taxId':9606}}}
        page={'count':1,'next':None,'results':[first_entry or entry(accession)]}
        def opener(request, **kwargs):
            return Response(protein if '/api/protein/' in request.full_url else page, {'InterPro-Version':'synthetic-1'})
        return p,opener

    def test_fetch_cache_and_offline_versioned(self):
        p,opener=self.setup_opener()
        with tempfile.TemporaryDirectory() as tmp:
            r=fetch_domains(p['accession'],tmp,offline=False,opener=opener)
            self.assertEqual(r['status'],'fetched')
            self.assertEqual(len(normalized_domains(r,p)),1)
            cached=fetch_domains(p['accession'],tmp)
            self.assertEqual(cached['status'],'cached')
            self.assertEqual(r['raw_response_sha256'],cached['raw_response_sha256'])

    def test_offline_missing_is_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(fetch_domains('P00001',tmp)['status'],'missing')

    def test_string_taxon_is_equivalent_and_raw_bytes_preserved(self):
        p,opener=self.setup_opener()
        def string_taxon(request, **kwargs):
            r=opener(request,**kwargs)
            payload=json.loads(r.read())
            if '/api/protein/' in request.full_url:
                payload['metadata']['source_organism']['taxId']='9606'
            return Response(payload,r.headers)
        with tempfile.TemporaryDirectory() as tmp:
            r=fetch_domains('P00001',tmp,offline=False,opener=string_taxon)
            self.assertEqual(r['status'],'fetched')
            for retrieval in r['snapshot']['retrievals']:
                self.assertEqual(file_hash(Path(tmp)/retrieval['raw_response_file']),retrieval['response_sha256'])

    def test_corrupt_original_response_prevents_cached_checked_state(self):
        p,opener=self.setup_opener()
        with tempfile.TemporaryDirectory() as tmp:
            r=fetch_domains('P00001',tmp,offline=False,opener=opener)
            raw=Path(tmp)/r['snapshot']['retrievals'][0]['raw_response_file']
            raw.write_text('corrupt')
            self.assertEqual(fetch_domains('P00001',tmp)['status'],'missing')

    def test_zero_domain_matches_explicitly_complete_not_negative_function(self):
        p,opener=self.setup_opener()
        def empty(request, **kwargs):
            if '/api/entry/' in request.full_url:
                return Response({'count':0,'next':None,'results':[]},{'InterPro-Version':'synthetic-1'})
            return opener(request,**kwargs)
        with tempfile.TemporaryDirectory() as tmp:
            r=fetch_domains('P00001',tmp,offline=False,opener=empty)
            self.assertEqual(r['status'],'fetched')
            self.assertEqual(normalized_domains(r,p),[])

    def test_malformed_source_is_error_not_batch_crash(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=fetch_domains('P00001',tmp,offline=False,opener=lambda *a,**k:Response({'metadata':None},{'InterPro-Version':'synthetic-1'}))
            self.assertEqual(r['status'],'error')

    def test_cache_corruption_preserved(self):
        p,opener=self.setup_opener()
        with tempfile.TemporaryDirectory() as tmp:
            r=fetch_domains('P00001',tmp,offline=False,opener=opener)
            path=Path(tmp)/r['raw_file']; path.write_text('corrupt')
            self.assertEqual(fetch_domains('P00001',tmp)['status'],'missing')
            new=fetch_domains('P00001',tmp,offline=False,opener=opener)
            self.assertEqual(new['status'],'fetched')
            self.assertEqual(path.read_text(),'corrupt')
            self.assertNotEqual(new['raw_file'],r['raw_file'])

    def test_query_error_and_retry_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch('ecd_snipr.interpro.time.sleep') as sleep:
                r=fetch_domains('P00001',tmp,offline=False,opener=lambda *a,**k: (_ for _ in ()).throw(OSError('unavailable')),sleeper=sleep)
            self.assertEqual(r['status'],'error')
            self.assertEqual(sleep.call_count,2)

    def test_cross_release_pages_error(self):
        p,opener=self.setup_opener()
        def mismatched(request, **kwargs):
            r=opener(request,**kwargs)
            if '/api/entry/' in request.full_url: r.headers['InterPro-Version']='different'
            return r
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIn('release_changed',fetch_domains('P00001',tmp,offline=False,opener=mismatched)['reason'])

    def test_pagination_count_no_silent_partial(self):
        p,opener=self.setup_opener()
        def truncated(request, **kwargs):
            if '/api/entry/' in request.full_url:
                return Response({'count':2,'next':None,'results':[entry('P00001')]},{'InterPro-Version':'synthetic-1'})
            return opener(request,**kwargs)
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(fetch_domains('P00001',tmp,offline=False,opener=truncated)['status'],'error')

    def test_external_next_url_not_followed(self):
        p,opener=self.setup_opener()
        def bad(request, **kwargs):
            if '/api/entry/' in request.full_url:
                return Response({'count':1,'next':'https://external.invalid/','results':[entry('P00001')]},{'InterPro-Version':'synthetic-1'})
            return opener(request,**kwargs)
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIn('pagination_url',fetch_domains('P00001',tmp,offline=False,opener=bad)['reason'])


class IntegrationTests(unittest.TestCase):
    def test_screen_exports_core_without_scaffold(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); make_set(root/'set',[make_protein('P00001')])
            result=run_screening(root/'set',root/'run',pilot=True)
            dest=Path(result['run_dir'])
            self.assertTrue(verify_bundle(dest))
            self.assertTrue((dest/'candidate_core_evidence.tsv').exists())
            self.assertEqual(read_json(dest/'candidates.json')[0]['screening_recommendation'],'standard_candidate')

    def test_overlay_enrichment_no_design_mutation_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); p,opener=self.setup_opener_for_public()
            make_set(root/'set',[p]); run=run_screening(root/'set',root/'run',pilot=True)['run_dir']
            before=file_hash(Path(run)/'manifest.json')
            fetch_domains('P00001',root/'cache',offline=False,opener=opener)
            a=audit_risk(run,root/'set',root/'overlay',interpro_cache=root/'cache')
            self.assertEqual(a['summary']['core_evidence']['coverage_states']['interpro'],{'evaluated':1})
            self.assertTrue(verify_bundle(a['run_dir']))
            self.assertEqual(file_hash(Path(run)/'manifest.json'),before)
            self.assertEqual(audit_risk(run,root/'set',root/'overlay',resume=True,interpro_cache=root/'cache')['execution'],'verified_cache_hit')

    def setup_opener_for_public(self):
        p,opener=FetchTests().setup_opener()
        p['taxon_id']=9606
        p['evidence']['kind']='curated_annotation'
        return p,opener

    def test_live_failure_partial_not_negative(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); p,_=self.setup_opener_for_public(); make_set(root/'set',[p])
            run=run_screening(root/'set',root/'run',pilot=True)['run_dir']
            with patch('ecd_snipr.interpro.fetch_domains',return_value={'status':'error','reason':'query_failed'}):
                a=audit_risk(run,root/'set',root/'overlay',interpro_cache=root/'cache',interpro_live=True)
            self.assertEqual(a['summary']['completeness'],'partial')
            self.assertEqual(a['summary']['core_evidence']['coverage_states']['interpro'],{'error':1})

    def test_generated_core_outputs_private(self):
        from manage_skill import forbidden
        for name in ('candidate_core_evidence.tsv','candidate_domain_coverage.tsv','core_evidence_summary.json','core_source_records.json'):
            self.assertTrue(forbidden(Path('examples')/name))

    def test_independent_core_verifier(self):
        from verify_core_evidence import verify_core
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); make_set(root/'set',[make_protein('P00001')])
            result=run_screening(root/'set',root/'run',pilot=True)
            self.assertTrue(verify_core(result['run_dir'],root/'set')['passed'])

    def test_independent_verifier_catches_resealed_derived_error(self):
        from verify_core_evidence import verify_core
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); make_set(root/'set',[make_protein('P00001')])
            result=run_screening(root/'set',root/'run',pilot=True); run=Path(result['run_dir'])
            candidates=read_json(run/'candidates.json')
            candidates[0]['receiver_risk']['core_evidence']['domains'][0]['retention']='omitted'
            write_json(run/'candidates.json',candidates)
            manifest=read_json(run/'manifest.json'); manifest['artifacts']['candidates.json']=file_hash(run/'candidates.json')
            write_json(run/'manifest.json',manifest)
            self.assertTrue(verify_bundle(run))
            self.assertFalse(verify_core(run,root/'set')['passed'])
