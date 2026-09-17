# Data contracts

## Analysis set (batch screening, v0.3.0)

`build-set` produces a set directory: `set_definition.json` plus one JSON per resolved reference under `proteins/`. The definition records the database name/release, the exact query (or the user list source), inclusion rules, creation time, completeness (`complete`/`partial`) and one entry per input row:

- `input`: original row, raw value, input type (accession / gene_name / blank / invalid_accession) and any extra columns. Rows are never dropped.
- `identity_status`: resolved / ambiguous / unresolved. Gene-name mapping keeps the recorded queries, the hit list and the mapping policy; ambiguity is preserved and never guessed away.
- `fetch`: status, URL, retrieval time, source entry version and raw-response SHA256 for accession fetches.
- `protein_file` / `protein_sha256`: normalized reference record, hash-verified on load.
- `reference_selection`: basis (`database_canonical` or `explicit_isoform`), selected isoform, rationale and `experimental_isoform_confirmation`. Canonical selection is an analysis choice, not laboratory isoform confirmation.
- `disposition`: membrane_set_membership, natural_cell_surface_target, scope_status with reason codes, extension_set (core / secreted_extension) and explicit notes for multi-pass, GPI, organelle and under-annotated entries.
- `processing_status` / `processing_error`: a single-item failure is recorded and never halts the batch.
- `duplicate_of`: duplicate inputs are linked to the first occurrence, not merged or dropped.

## Screening record (v0.3.0)

One record per analysis object in `screening.json` / `protein_screening.tsv`: entry/protein identity, reference selection, scope_status (+ reason codes), extension set, processing_status, screening_recommendation (one of the four classes or null for out-of-scope/technical failure), reason_codes, brief rationale, evidence (boundary/topology/reference/interval), missing_info, not_evaluated, candidates_evaluated, primary_candidate_id, alternates with `reason_not_primary`, deferred_annotations, lab_rules application, and functional_status (always `not_experimentally_validated`).

Lab rules file: a JSON object with a `rules` list. Each rule needs `rule_id`, a non-empty `match` on gene/accession/topology, `effect` in `annotate` / `downgrade_to_conditional` / `downgrade_to_no_standard_route`, `rationale`, and sourced `evidence` (kind local_experiment or curated_annotation, source, version). Rules annotate or downgrade only.

## Project

The schema in [project.schema.json](../schemas/project.schema.json) documents the serialized interface; `contracts.validate_project` and downstream coordinate/scaffold/assay checks implement runtime validation. The runtime is not a general-purpose JSON Schema engine. No external dependency is required.

`schema_version=1.0`, `configuration=antibody_sender_antigen_receiver`, `proteins`, optional `scaffold`, `existing_constructs`, `observations`, `reviews`, `rules`, `evidence_records_v2`.

Each protein needs a stable local `protein_id`, explicit accession/isoform, `taxon_id=9606`, actual amino-acid sequence, topology/location and evidence. A missing/ambiguous isoform blocks sequence export for that candidate. Do not label an arbitrary canonical sequence as the intended isoform.

Sequence normalization strips whitespace and uppercases only. The engine accepts the 20 standard amino acids; unknown residues, stops, gaps and nonstandard residues require explicit manual resolution, not substitution. Source sequences must not be silently altered to a structure-paper variant.

## Annotation and candidate coordinates

- Every coordinate is **1-based inclusive** on that protein's exact stored sequence; extraction is `seq[start-1:end]`.
- Feature fields: `kind`, integer `start/end`, optional `name`, `boundary_status` (`exact` or `fuzzy`), `evidence`.
- Supported interpretation: `signal_peptide`, `transmembrane`, `extracellular`, `cytoplasmic`, `lumenal`, `domain`, `disulfide`, `epitope`, `gpi_signal`, `gpi_attachment_site`, `chain`, `propeptide`. `processed_peptide` is preserved but not automatically considered an independent fold. Unknown annotations may be retained but are not automatically interpreted.
- Evidence object: `kind` in `curated_annotation / prediction / local_experiment / synthetic_fixture`, `source`, `version`; optional `raw_sha256`, ECO records, URL, retrieval time. Recorded annotations are not automatically experimental confirmations.
- Disulfide start/end refer to the two cysteines; crossing a candidate boundary raises a risk. Epitope start/end conservatively span the mapped epitope; sequence containment does not imply preserved conformation. Discontinuous/conformational epitope maps should retain their primary source; explicit spatial compatibility prediction is not implemented.
- `candidate_regions` adds exact sourced intervals with `rationale`. No unsourced LLM interval can be assembled. `candidate_policy=full_ecd_only` is default; `full_ecd_and_domains` explicitly enables domain proposals. Explicit candidates remain supported; an unspecified form is a `domain_fragment`, not a presumed full ECD.
- Overlapping/alternative/fuzzy topology annotations are not arbitrarily merged. A complete domain is a candidate alternative, not a demonstrated autonomous fold; domain interfaces may span neighboring domains.

Antigen-form and contextual-consideration fields are specified in [screening criteria](screening-criteria.md). Version 0.2 includes antigen form in candidate identity; do not match old and new candidates only by their truncated hash. Match reference hash, isoform, coordinates, sequence and form explicitly. Historical run bundles are preserved; changed code or rules invalidate assembly approval.
- Thresholds `short_length=50`, `long_length=500`, `cysteine_fraction=0.08` are configurable **review heuristics**, not empirically fitted SNIPR thresholds. No probability or ranking score uses them.

## Actual receiver scaffold

Required fields: `scaffold_id`, `version`, `source`, `host_cell`, `reporter`, `configuration`, `status=verified_local`, `signal_policy=scaffold_signal_only`, `junctions_verified=true`, ordered `modules`.

Each module has `name`, `role`, exact amino-acid `sequence`, and `source`. Exactly one `antigen_slot` has no sequence. Exactly one N-terminal `signal_peptide` precedes it; exactly one `transmembrane` and one `transcription_factor` follow it. Supply actual linker, regulatory extracellular segment, juxtamembrane or other modules in their verified order; do not add them from literature defaults. The configuration is a caller assertion of reviewed lab material, not something the engine independently measures.

`synthetic_fixture` is permitted solely for synthetic software tests and marks resulting sequences accordingly. **No laboratory scaffold template sequence is bundled.** Do not change synthetic status to verified_local and treat it as real.

The antigen ectodomain being inserted is not synonymous with the regulatory extracellular module called ECD in SNIPR literature. Different modules remain independently named in `connection_map` with coordinates on the assembled fusion. Missing/unsupported configurations produce fragment proposals only, never guessed fusions.

## Existing constructs

`construct_id`, `protein_id`, `start`, `end`, `antigen_sequence`, `scaffold_id`, `scaffold_version`, `fusion_sequence`, `source_ref`; optional `antigen_cds`, `fusion_cds`.

Sequence/CDS mismatches remain review issues. A terminal stop is allowed only in CDS translation, not amino-acid strings. No CDS/codon optimization is generated. Multiple constructs of one gene remain separate; duplicated construct IDs are rejected. Assembly approval is not experimental validation.

Exact antigen matches are listed as `antigen_only_not_equivalent_receptor` unless the full fusion and scaffold identity also match. Actual observations are linked by construct ID; they are not inherited by a newly truncated or differently connected receptor.

## Observations

See [assay interpretation](assay-interpretation.md). Preserve `raw_value`, `source_ref`, `endpoint`, `unit`, gate/channel/statistic, and the entire assay context. The project file may retain original fields as metadata. Missing values must not become zeros.

## Reviews

`review_key`, `candidate_id`, `decision=pending/approve/reject`, `reviewer`, `reviewed_at` (ISO8601), `rationale`, `acknowledged_risks`. `required_risks` in templates is explanatory. Approval requires every review-level risk code to be explicitly acknowledged. Duplicate review keys are rejected; edit the input review list to select the current decision while retaining old immutable runs. This is a traceable review record, not a cryptographic signature or access-control system.

## EvidenceRecord v2 interoperability

Imported legacy v2 records are validated and exported without modification; extra properties are rejected. Their original statuses (including missing, conflict and error) remain untouched. `evidence_links.json` references each original record by canonical JSON SHA256. A candidate/gene identifier by itself does not bind an experiment to a construct, sender or endpoint.

New computations use `claims.jsonl` with `evidence_class=deterministic_rule`, `status=computed`. **They are intentionally not coerced into legacy confirmed_fixture/confirmed_live statuses.** The old v2 enum has no honest generic local-computation status. New claims carry native context and checksums; legacy context linkage remains unresolved unless a later explicit mapping supplies it. This is backward-compatible consumption, not a silent v2 schema extension.

Hashes labelled `reference_sha256`, `claim_id`, `review_key`, and run keys use canonical JSON serialization (including string quotes); source-file and artifact hashes use raw bytes. The algorithm and encoding are stable and recorded in code. Do not compare hashes from different conventions as if identical.
