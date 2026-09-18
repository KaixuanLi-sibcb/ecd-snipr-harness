# Screening criteria: evidence before gates

## Purpose

Design review for **Antibody-Sender -> Antigen-Receiver**, not therapeutic target ranking or a universal receptor-function classifier. The candidate unit is `(reference sequence, isoform, antigen form, exact interval)`. Functional observations additionally require actual construct, sender/antibody and experimental conditions.

## Adopted, conditional and deferred

| Criterion | Treatment | Reason / limitation |
|---|---|---|
| Human sequence identity, accession, isoform and provenance | Required | Canonical is a reference choice, not proof of the experimental isoform |
| Reviewed annotation | Preferred, not an absolute biological gate | A documented noncanonical isoform or exact local sequence must not be discarded solely by this preference |
| Topology and compartment | Required | Intracellular-only records are outside the default antigen scope; unknown topology is missing evidence, not intracellular classification |
| Continuous full ECD | Default candidate | Coordinates must be supported on the intended sequence; SP/TM-derived boundary guesses are not silently promoted to annotations |
| Mature secreted antigen | Conditional scope extension | Requires secreted context, no TM conflict and an annotated chain; not evidence of native membrane accessibility |
| Multiple processed chains | Review each chain, preserve dependencies | An annotated chain need not be independently stable or represent all conformational epitopes |
| GPI mature antigen | Conditional candidate | Include the omega residue; remove the downstream GPI signal. Require sourced attachment and mature boundaries |
| Known shed product | Explicit sourced proposal only | Exact boundary, release evidence and release context required; a cleavage motif is insufficient |
| Type II topology | Orientation review | Native N-terminal tether becomes C-terminal tether in this supported receiver layout |
| Multi-pass | No whole-ECD synthesis / no loop concatenation | An explicitly sourced autonomous external domain may be reviewed; original membrane-dependent epitopes may be lost |
| Domain boundaries and inter-domain contacts | Required review | Full annotated domain does not prove autonomous folding |
| Known / unknown epitope coverage | Required review | No mapped epitope does not mean no epitope. Preserve full antigen baseline and document changes in screenable repertoire |
| Native processing / autocleavage / shedding | Sourced contextual flag | Native processing cannot be transferred directly into receiver background predictions |
| Junction cleavage | Separate contextual review | Requires actual fusion context; no unvalidated motif-to-self-activation classifier |
| Oligomerization vs aggregation | Separate records | Native assembly can support structure or recognition; it is not synonymous with pathological aggregation |
| Short / long / cysteine-rich / odd-cysteine / glycosylated | Configurable review heuristics | No universal optimum, success probability or automatic rejection; glycans/disulfides can be essential; odd cysteine parity and sequon density are review warnings, never thresholds |
| Culture / B-cell interference | Conditional context review | Need actual co-culture system, interacting partner and applicable evidence; gene names alone do not establish interference |
| Four functional endpoints | Experimental records only | Surface expression, recognition retention, basal activation and induced response remain independent |
| Disease value / tissue specificity / known drugs | Deferred from receiver engineering | Potentially useful to prioritize applications, not to decide whether an antigen form can be evaluated |
| Human-mouse identity | Deferred from receiver engineering | A project-specific preference is not a generic constructability criterion |
| Global red/yellow/green success | Not adopted | Keep design, assembly and experimental status separate. If colors are added to a report, they must label a named readiness dimension |

## Implemented input contract

Supported antigen forms: `full_ecd`, `mature_secreted`, `mature_gpi`, `shed_product`, `domain_fragment`.

`chain` is an annotated processed chain, not automatically a secreted antigen. Only explicit `topology=secreted`, `location=secreted` without TM conflict enables the secreted route. A GPI chain must end at the sourced `gpi_attachment_site`, a single inclusive coordinate. The GPI signal is omega+1 onward; uncertainty does not justify shortening by one residue. These are candidate rules, not a recommendation to preserve native GPI anchoring in the receiver.

For a shed candidate supply `antigen_form_type=shed_product`, exact bounds, `rationale`, coordinate `evidence`, separate `release_evidence` (not a motif-only prediction), and `release_context`. Current automated support requires containment in an annotated external region (or processed chain in a secreted reference). Conflicts require annotation resolution, not guessed bounds. Protease identity may be recorded in context; it is not always known even when release is established.

Each optional `considerations` record has:

- `kind`: `processing`, `native_shedding`, `junction_cleavage`, `oligomerization`, `aggregation`, `culture_interference`.
- `status`: `reported`, `not_observed_in_context`, `missing`, `not_applicable`.
- `scope`: `native_protein`, `receiver_construct`, `assay_context`.
- `reason`, `evidence` with source/version/kind, `context`; optional exact `start/end`.
- Receiver/assay scope requires `context.construct_id` and `context.condition_id`. Context should retain host, sender, antibody, medium and batch where applicable. The software checks record structure, not whether a paper truly supports the claim; source interpretation remains a review task.

Every candidate has six contextual review rows. No record yields `missing/not_assessed`, never low risk. Contradictory or differently scoped observations are preserved separately (`mixed_context_records`), not pooled. A reported site absent from the candidate remains a contextual record, not proof the fusion cannot be cleaved elsewhere. Native shedding and self-cleavage are not automatically receiver basal activation.

`candidate_plan.tsv` separates `antigen_form_type`, `design_status`, `assembly_status`, `risk_review_status`. Details appear in `criteria_review.tsv` and native contextual claims. No automatic green/safe label. Legacy EvidenceRecord v2 remains losslessly preserved, not relabeled with incompatible statuses.

## Screening recommendation layer (v0.3.0)

`screening_recommendation` is independent of design/assembly/functional status and is assigned by a transparent decision table — no uncalibrated composite scores.

| Class | Requirement | Explicitly not |
|---|---|---|
| `standard_candidate` | In scope; a usable candidate with positive sequence, topology and exact sourced boundary evidence; no conditional issue on the primary candidate | Proof of expression, recognition or activation |
| `conditional_candidate` | A sourced candidate exists but carries a named issue: type-II orientation change, GPI mature-boundary check, multichain dependency (multiple processed chains or SUBUNIT-recorded hetero-oligomer partnership), secreted-extension presentation, multi-pass external-domain choice, a domain cut by the candidate boundary, a boundary/adjacent-annotation conflict, unreviewed-entry prediction boundaries, unconfirmed plasma-membrane annotation, or a sourced lab rule | A failure label |
| `no_standard_route` | In scope with understood topology, but no supported route under current rules (e.g. multi-pass without a sourced external domain) | "This protein can never be used" — it is a route limitation |
| `insufficient_evidence` | Identity, topology or key boundary evidence missing or conflicting | A biological failure |

Decision order: technical failure (`processing_status`) and out-of-scope (`scope_status`) are recorded first and are not recommendations; then insufficient-evidence codes; then usable-candidate selection. The primary candidate is chosen by an explicit deterministic rule — form priority (full ECD > mature GPI > mature secreted > shed > domain fragment), fewest conditional issues, fewest review risks, then candidate ID — and alternates are kept with their reason for not being primary.

Isoform inventory is recorded, not guessed. The analysis reference records `annotated_isoform_count` and `isoform_comparison: not_evaluated`: the UniProt entry document lists annotated isoforms and textual alternative-sequence differences but does not carry isoform sequences, so a sequence-level canonical-vs-isoform comparison cannot be performed in batch screening and is never silently assumed. An annotated alternative-sequence difference overlapping the candidate interval is surfaced under `evidence.isoform_differences_within_candidate` plus a `missing_info` note — as a record only. It never downgrades, blocks or reclassifies: a textual difference feature is not the experimental isoform, and the validated policy for exactly this situation (e.g. the LAG3 shed-form annotation) is record-and-confirm, not auto-conditional. A genuine conditional flag would require the actual reviewed isoform sequence in hand, which is the assembly-gate confirmation, not a screening input.

Positive evidence is required for a standard candidate; "no risk found" is not evidence. Conversely, missing epitope maps or missing contextual-risk literature never auto-downgrade a candidate with reliable boundaries: the aspects are listed under `missing_info` / `not_evaluated`. Prediction-evidence boundaries in a reviewed entry remain a review-level risk; in an unreviewed entry they make the recommendation conditional. Fuzzy or unsourced *secondary* annotations (for example a shed-form chain annotation on a type-I receptor, or alternative splice chains next to an exact main chain) are deferred with an explicit `deferred_annotations` record; a fuzzy *sole essential* annotation (the only signal peptide, the only chain of a secreted protein) still blocks. The strict assembly path in `run` keeps every fuzzy feature blocking.

Laboratory experience rules (`--lab-rules`) must carry rule_id, match (gene/accession/topology), effect (`annotate`, `downgrade_to_conditional`, `downgrade_to_no_standard_route`), rationale and sourced evidence. They annotate or downgrade only — never upgrade. Without an imported rule set every record states 尚未纳入.

## Deep screening layer (v0.4.0)

The screen now reads all candidate-relevant annotation already present in the UniProt entry and records it explicitly. Nothing here is a composite score; every outcome is a named code with named evidence.

### Fuller feature extraction

`normalize` parses SIGNAL, PROPEP, CHAIN, PEPTIDE, TOPO_DOM, TRANSMEM, DOMAIN, REGION, MOTIF, DISULFID, CARBOHYD, SITE, BINDING, ACT_SITE, VARIANT, MUTAGEN and LIPIDATION (non-GPI lipidations kept separate from the GPI omega site), each with positions, description, feature ID, exact/fuzzy boundary status and curated-vs-prediction evidence (ECO codes 0000255/0000256/0000259 mark prediction). The entry records `reviewed` status, `annotation_score`, and the machine-readable SUBUNIT/FUNCTION/PTM comment texts. Informational kinds (region/motif/glycosylation/site/binding/active_site/variant/mutagenesis/lipidation) never produce blocking coordinate errors: fuzzy or invalid ones are deferred with an explicit record in screening, and remain a review note in the strict assembly path. Candidate-relevant kinds keep blocking semantics.

### Per-candidate molecular profile

Each candidate carries `molecular_profile`: length; cysteine count, positions and odd parity; disulfide bonds fully contained vs crossing the boundary (with the retained cysteine position); N-glyco sequon scan (N-X-S/T, X≠P) with positions; annotated CARBOHYD sites inside the fragment with N/O subtype and evidence kind; domains fully contained, cut by the boundary (identity, cut side, retained interval) or outside; active/binding/SITE/lipidation sites inside vs outside; natural-variant and mutagenesis records overlapping the fragment. Variant/mutagenesis overlaps are records only — they never change a class.

### Boundary precision analysis

`boundary_analysis` names the defining feature for each candidate (e.g. TOPO_DOM extracellular 20–291 or the mature chain) and checks the adjacent anchor annotation on each side (SIGNAL/TM/PROPEP/GPI-signal/chain; a mature chain coincident with the boundary counts as agreement, never as conflict). Outcomes per side: `exact`, `gap_tolerated` (unannotated gap ≤ `boundary_gap_tolerance`, default 5), `gap_beyond_tolerance`, `overlap`, `missing`, `chain_terminus`, or `not_applicable` (domain-fragment boundaries are defined by the domain itself). Domain-fragment candidates are not forced to justify against SP/TM.

### New reason codes and their severity mapping

| Code | Level | Evidence used | Rationale for the mapping |
|---|---|---|---|
| `domain_cut_by_boundary` | Conditional | DOMAIN feature overlapped but not contained by the candidate interval | Cutting an annotated domain changes the fold/epitope scope; the domain identity and cut side are named in the rationale (e.g. CD59: the UPAR/Ly6 annotation extends past the omega residue into the GPI signal, so the mature-boundary cut is correct but must be verified) |
| `multichain_partner_required` | Conditional | SUBUNIT comment containing explicit hetero-oligomer wording (`heterodimer/heterotrimer/heterotetramer/hetero-oligomer/heteromultimer`); matched sentences quoted verbatim | A documented heteromeric assembly means ECD surfaces may engage a different chain; this is the "multichain dependency" conditional category. Conservative by design: homo-oligomers, bare "interacts with" and "part of a complex" phrasing are recorded in the raw comment but do not fire the code |
| `boundary_feature_conflict` | Conditional | Adjacent anchor annotation overlapping the boundary or separated by a gap beyond tolerance | A boundary disagreement is a named verification item. A genuine TM/SP overlap with the candidate still blocks via the retained-* checks before this point |
| `tm_adjacency_tolerance_used` / `boundary_adjacency_tolerance_used` | Informational reason code | Adjacent anchor within the tolerated unannotated gap | Documents that tolerance was applied; never changes the class |
| `free_thiol_odd_cysteine` | Warning (review) | Odd cysteine count in the fragment | Potential unpaired thiol; cysteine context is a review heuristic, not a pass/fail threshold |
| `dense_glycosylation` | Warning (review) | ≥ `dense_glyco_min_sequons` (4) sequons AND ≤ `dense_glyco_max_residues_per_sequon` (50) residues per sequon | Glycan density may matter for recognition/expression; heuristic only, occupancy is not evaluated |

`prediction_requires_annotation_review` is now scoped to the candidate's own region plus boundary-defining, domain and disulfide features: a prediction-grade glycosylation site or variant record no longer raises a boundary-evidence flag. Reviewed-entry handling (risk stays review-level) and unreviewed-entry handling (conditional) are unchanged.

### Multi-pass loop enumeration

For multi-pass references every annotated extracellular loop is enumerated individually in `multipass_loops`: interval, length, fully contained sourced domains, eligibility as a conditional alternate, and links to any proposed alternate candidate. The enumeration is recorded even when no alternate exists (e.g. CXCR4: four loops of 38/11/21/21 residues, none containing a sourced domain, so no alternate is eligible and the route stays no_standard_route). Loops are never stitched.

### Checks applied vs not evaluated

`PI_SUMMARY.md` carries a "checks applied" section listing exactly what the deterministic screen evaluated and what remains not evaluated (epitope literature, contextual-risk literature unless records are supplied, isoform sequence-level comparison, glycan occupancy, structure/folding/aggregation/cleavage prediction, the four experimental endpoints, laboratory isoform confirmation). `summary.json` adds `reason_code_tallies` with explicit denominators; per-code tallies need not sum to class counts because one reference can carry several codes.

## Scientific sources and boundaries

- [UniProt sequence processing](https://www.uniprot.org/help/sequence_processing): precursor processing may include SP and propeptide removal; absent processing annotation is not proof of absent processing.
- [UniProt GPI transamidase annotation](https://www.uniprot.org/uniprotkb/Q92643/entry) and [GPI anchor definition](https://www.uniprot.org/locations/SL-9902): attachment is to the mature terminal omega residue, not to a residue preceding it. This supports coordinate semantics, not SNIPR function.
- [Ming et al., 2022](https://pubmed.ncbi.nlm.nih.gov/35761082/): LAG3 has mapped ligand/antibody interfaces and oligomeric context; truncation may change the recognition repertoire. It is a case-specific precedent, not a universal truncation predictor.
- [Zhu et al., 2022](https://doi.org/10.1016/j.cell.2022.03.023): modular receptor research distinguishes recognition and receptor regulatory modules. Its ECD terminology must not be equated with every target protein's complete extracellular antigen.
- [Public-reference review](public-reference-architecture.md): source descriptions are not independent experimental validation and cannot substitute for the lab scaffold.

## What remains unimplemented

No protease-site, aggregation, glycan-occupancy, B-cell-interference or SNIPR-function predictor is introduced. Contextual source records currently need explicit curation. No automatic domain-interface/epitope literature mining and no optimization of mutations, linkers or assay conditions. Review thresholds are not fitted to the laboratory's data. Screening recommendations are a deterministic first pass; targeted literature review and agent explanation are optional enhancements, not required steps.
