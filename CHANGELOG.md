# Changelog

## 0.3.0 - 2026-09-16

- Lab-workbook cross-check fix (same version, regression test added): UniProt entries annotated "Secreted" plus a generic organelle membrane location (e.g. REN "Secreted"+"Membrane", TFPI "Secreted"+"Microsome membrane") were misclassified as other_membrane and ruled out_of_scope. `normalize` now gives Secreted precedence over generic membrane; plasma membrane still wins over both; raw_locations are preserved. Found via validation_v031 three-way comparison, not tuned to force agreement. 124 tests green.

- New default main line: user target list or public UniProt query -> identity resolution and explicit analysis-reference selection -> batch fetch/cache/normalization -> sourced antigen fragments -> per-target screening recommendation -> set-level summary. Real scaffold assembly and experiment linkage remain optional downstream branches in `run` and are unchanged.
- New `acquisition` module: `build-set` builds a membrane-protein analysis set; explicit accessions are processed directly, gene names are mapped via recorded `gene_exact + organism_id:9606` queries (reviewed first, then all entries); ambiguity is preserved, no guessing, no input row dropped, duplicates linked not merged. Every entry gets an explicit disposition distinguishing membrane-set membership, natural cell-surface target, design scope and secreted extension; single-item failure never halts the batch; incomplete acquisition is marked partial.
- Explicit analysis-reference selection: screening may select the database canonical sequence with a recorded rationale (`reference_selection`); this is not laboratory isoform confirmation, which still gates fusion assembly. Genuine identity ambiguity, version conflicts and non-canonical coordinate problems are preserved.
- New independent `screening_recommendation` layer (`screening.py`), separate from design_status / assembly_status / functional_status. Transparent decision table: standard_candidate / conditional_candidate / no_standard_route (route limitation) / insufficient_evidence (evidence state, not a biological failure); out-of-scope via scope_status, technical failures via processing_status. No composite scores. Every recommendation carries reason_codes, rationale, evidence, missing info and alternates with reasons. Standard candidates require positive sequence/topology/boundary evidence; missing epitope or contextual-risk literature is reported as not evaluated and never auto-downgrades a candidate with reliable boundaries. Prediction-evidence boundaries in reviewed entries stay review-level risks; in unreviewed entries they make the candidate conditional. Fuzzy secondary annotations (shed-form or alternative splice chains, fuzzy domains) are explicitly deferred and recorded; sole essential annotations still block.
- Sourced lab experience rules can be imported (`--lab-rules`); they annotate or downgrade with provenance and never upgrade. When absent, every record states 尚未纳入.
- Real design content per recommendation: reference identity, isoform, candidate ID, antigen form, exact 1-based coordinates, length, sequence, rationale, main risks, missing info and alternates, plus junction notes. Without a verified scaffold only candidate-fragment FASTA is exported; no guessed fusion sequence.
- New summary deliverables: `protein_screening.tsv` (one row per analysis object), `screening.json`, extended `candidate_plan.tsv` / `candidates.json` / `candidate_fragments.fasta`, `summary.json` (set definition, completeness, separate gene/reference/candidate counts, class counts, out-of-scope and technical-failure counts, all ratios with explicit denominators), `PI_SUMMARY.md` (conclusions first), and a dependency-free SVG figure `screening_overview.svg` with `screening_overview_data.json` source data, labelled 候选设计覆盖/初筛建议 (pilot runs labelled PILOT), never experimental success rates.
- New CLI commands `build-set` and `screen` chain set-building -> normalization -> screening -> summary with cache/version/resume and per-item failure isolation; no per-protein hand-written JSON. Exit code 3 marks a partial (incomplete) set.
- tests/test_screening.py adds 30 tests (123 total): screening without scaffold/review/experiments, canonical-vs-isoform distinction, missing optional risk info, core sequence errors still blocked, topology/mature-boundary/domain-alternate handling, multi-reference/multi-candidate/duplicate/out-of-scope counting, per-item failure isolation, partial-run honesty, lab-rule provenance, and non-regression of the legacy assembly/endpoint branch.
- Live public pilot: 12 UniProt entries across type-I/II, GPI, multi-pass GPCR, multichain/secreted extension (logs and review in the validation output directory). Software tests demonstrate contract execution, not biological accuracy.

## 0.2.0 - 2026-09-16

- Independently reviewed screening criteria; separate annotation readiness, assembly authorization and four experimental endpoints. No no-evidence-to-low-risk conversion.
- Prefer full antigen forms by default; domain truncation now explicit opt-in, preserving epitope-scope review.
- Add annotated mature secreted/GPI forms, propeptide exclusion and sourced shed-product proposals. Retain the GPI omega residue and distinguish multichain dependence.
- Add source/context records for processing, native shedding, junction cleavage, oligomerization, aggregation and culture interference; export `criteria_review.tsv` and contextual claims. These are not new functional predictors.
- Reassess patent-family dependence and public synNotch descriptions, keeping actual Ag-SNIPR equivalence unconfirmed.
- Candidate IDs now include antigen form; old run bundles remain valid historical artifacts, but new rules/code invalidate prior approvals and caches.

## 0.1.0 - 2026-09-16

- New independent Antibody-Sender -> Antigen-Receiver skill, leaving the older therapeutic-antigen prioritization workflow intact.
- Deterministic source-bound ECD/domain candidates, epitope/domain/disulfide risk review and scoped human-reviewed scaffold assembly.
- Read-only workbook inventory and explicit mapping; four context-specific experimental endpoints without inferred success labels.
- Content-addressed immutable runs, checksums, retry/cache states, evidence/construct provenance and conflict preservation.
- Synthetic offline regression suite, clean package/install allowlist and documented biological limitations.
- Source-checked patent/Addgene reference catalog kept separate from actual lab scaffold; sender PDGFR-LC sequence IDs explicitly excluded from receiver assembly.
- Public-reference input does not block ECD candidate analysis and cannot authorize fusion export.
