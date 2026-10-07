---
name: ecd-snipr-harness
description: Design and review sourced human membrane-protein extracellular fragments for a laboratory SNIPR Antigen-Receiver paired with an Antibody-Sender. Use for batch ECD candidate planning, per-target screening recommendations over a defined protein set, scaffold-aware fusion review, epitope-retention checks, and context-specific experimental evidence tracking. Not for antibody design, antigen-sender qDis engineering, therapeutic target ranking, or unsupported SNIPR success prediction.
---

## Current integrated workflow (v0.10.2)

Read [the complete workflow and standards](references/workflow-and-standards.zh-CN.md)
before a batch run. Keep design advice, core-evidence checks, mechanism-review
priority and four experimental endpoints independent. For scheduled local tools,
follow [the portable HPC guide](references/hpc-execution.md); require exact native
CLI concordance before adopting the DeepTMHMM2 compatibility worker. Never package
institutional records, generated results, configurations, licenses or model assets.

# ECD-SNIPR Harness

## Scope and routing
Target configuration is **Antibody-Sender -> Antigen-Receiver**. The antigen is on the SNIPR receiving cell, not on the LC-bearing sending scaffold. Treat literature on antigen-sender qDis as background only.

The default main line is batch screening: **a user target list or a public human UniProt set -> identity resolution and explicit analysis-reference selection -> batch fetch/cache/normalization -> sourced antigen fragments -> a screening recommendation per target -> set-level summary.** Real scaffold assembly and experimental-endpoint linkage are optional downstream branches (`run`); their absence never blocks screening, and no scaffold, linker, review or experiment label may be fabricated to keep a run going.

Read [workflow](references/workflow.md) for execution and task-state rules; [data contracts](references/data-contracts.md) before mapping inputs or assembling sequences; [screening criteria](references/screening-criteria.md) for the recommendation decision table; [assay interpretation](references/assay-interpretation.md) for experimental data; [evidence and methods](references/evidence-and-methods.md) when explaining predictions or citations. CLI examples are in [README](README.md).

When the actual scaffold is not yet supplied, read the [source-checked public reference descriptions](references/public-reference-architecture.md). CN114437232B / WO2022095916A1 belong to one patent family, not independent replication; Addgene #79127 is a synNotch reference, not the current Ag-SNIPR. Never use sender SEQ ID NO.5/6 as the receiver scaffold. Continue ECD analysis while keeping assembly pending.

## Required behavior
1. Build the analysis set honestly. Record the data version, query, inclusion rules and fetch completeness. Accessions are processed directly; gene names are mapped with recorded evidence, ambiguity preserved, no guessing, no row dropped. Multi-pass, GPI, organelle and under-annotated entries get explicit dispositions — nothing is silently deleted. Distinguish membrane-set membership, natural cell-surface target, and design-evaluation scope; secreted proteins are a separately counted extension set.
2. Screening may explicitly select the database canonical/reference sequence as the analysis reference with a recorded rationale. This is not experimental isoform confirmation, which still gates final fusion assembly. The annotated isoform inventory is recorded with the selection and sequence-level comparison against other isoforms is explicitly `not_evaluated` (the UniProt entry document carries no isoform sequences); annotated differences overlapping a candidate are recorded as evidence only, never an automatic downgrade. Genuine identity ambiguity, version conflicts and non-canonical coordinate problems are preserved, never smoothed over.
3. Give each in-scope target one screening_recommendation from the transparent decision table: `standard_candidate` / `conditional_candidate` / `no_standard_route` (route limitation, not "unusable forever") / `insufficient_evidence` (missing core information, not a biological failure). A standard candidate needs positive sequence, topology and boundary evidence — never a mere "no risk found". Missing epitope or contextual-risk literature is reported as not evaluated, not used to downgrade a candidate with reliable boundaries.
4. Every recommendation carries actual design content: reference identity, isoform, candidate ID, antigen form, exact 1-based coordinates, length, candidate sequence, selection rationale, main risks, missing info and alternates with reasons. Type-II orientation change, GPI mature-boundary check and multi-pass external-domain conditionality stay distinct. No loop stitching, no random window-sliding, no default mutation optimization, no mechanical truncation, no forced uniform candidate counts.
5. Coordinate errors and retained native TM/SP/cytoplasmic tail remain blocking. Length, cysteine and glycosylation motifs are warnings, not pass/fail thresholds. Invalid topology/exclusion annotations cannot be discarded merely because another feature of the same kind is valid. Secondary annotations may be explicitly deferred, without claiming their impact is absent. Molecular profiles and adjacent-boundary checks are deterministic annotation summaries, not folding or function predictions. SUBUNIT heteromer wording is native context requiring review, not proof that an ECD requires the partner.
6. Use the **actual laboratory scaffold** with explicit module sequences/order/version, signal policy, host, reporter and verified junctions for the assembly branch. Without it, deliver candidate-fragment FASTA plus junction notes and export no fusion. Require a named, authorized human review for assembly; the agent never self-approves.
7. Keep four endpoints separate: surface expression, recognition retention, basal activation, induced response. Associate each observation with a specific construct, sender/antibody, batch, replicate and condition. Blank is missing; Y/N is not a functional label.
8. Lead reports with research conclusions and coverage with explicit denominators — not approval counts. Label outputs 候选设计覆盖/初筛建议, never experimental success rates; label small runs pilot.

## Execution
Read [local predictor execution](references/local-predictors-v0100.md) before using
`configure-predictors`, `predict-local` or `validate-local-tools`. Install tools/models
outside the skill. Actual execution and version/hash verification are separate from
an importer; DeepTMHMM2 is not legacy DeepTMHMM and uses CPU on Apple Silicon.
Requested unavailable tools make the tool run partial. Never submit private sequences
to hosted services, relabel empty output as negative, or call disorder/topology
agreement SNIPR functional accuracy. Preserve frozen candidates and four endpoints.

Read [core evidence v0.9.0](references/core-evidence-v090.md). Preserve the candidate
engine and design classes; independently cross-check exact fragments and
structural units. `core-evidence` supports optional public-accession-only InterPro
retrieval with exact reference-sequence equality, cache/version/hash checks and
visible conflicts. `import-prediction` reads existing local outputs plus original
FASTA; it never executes a predictor or submits a private sequence. Full-reference
coordinates require explicit projection, not silent transfer. Partner association
is not reference self-association; unresolved subjects stay unresolved. Family
classification is not an autonomous domain, correlated hits are not independent
votes, and no conflict detected never establishes low functional risk.

Read [receiver risk v0.8.0](references/receiver-risk-v080.md). Every new `screen` and
`run` includes an independent mechanism-priority layer, with source and scope checks.
Use `receiver-risk` plus the exact reference set for immutable historical overlays.
Use `sequence-tools` for local descriptors and optional Biopython; normalized local
IUPred/DeepTMHMM/SignalP predictions require fragment identity, hash, coordinates,
version, parameters and actual raw output file. Never call an unavailable tool checked.
Do not turn native shedding/dimerization or generic glycosylation into measured
SNIPR self-activation. Missing context never establishes low risk. Do not tune rules
or sequence cutoffs against a supplied failure list. Report whole-cohort flag
prevalence alongside retrospective matches; absent eligible successes means no
accuracy/specificity estimate. This is risk-priority triage, not a validated predictor.

Read the [v0.7.0 receiver-function policy](references/receiver-function-v070.md).
Always present fragment design and receiver function as separate axes. `standard_candidate`
does not imply low basal activation; `conditional_candidate` is not a failure prediction.
Export the four endpoint evidence states, scope and next actions; unknown risk remains
`undetermined_not_low_risk`, even when no warning is found. Never add gene-specific rules
from a failure list or report generic review flags as successful failure predictions.
Native protein biology cannot establish function of a new receiver. Never transfer outcomes
between antigen-only matches or scaffold versions. Use `receiver-audit --run-dir RUN --outdir NEW`
to overlay a verified historical bundle without changing its design results.

Read [antigen-unit scope](references/antigen-units-v051.md) before interpreting a fragment as representative of a gene. Review product-specific omissions, spans across annotated processed chains and cut repeats. A full topological ECD is not a complete mature-product/epitope repertoire. Missing epitope mappings mean none were supplied, not that an external search was negative. v0.5.1 keeps the candidate generator intact; unsupported segmented mature-product routes remain explicit gaps, not guessed concatenations.

For current selection and evidence semantics read [methodology v0.5.0](references/methodology-v050.md). Keep design class separate from annotation support; never upgrade feature evidence using entry reviewed status. Explain primary/alternate repertoire tradeoffs, route gaps and receiver review questions. The pending stratified review queue is a purposive audit aid, not a biological accuracy sample or assembly approval.

Use `python3 scripts/ecd_snipr_cli.py --help`. Core and offline tests require Python >=3.10, no external packages. Always use a private output directory outside the installed skill.

- `screen --list targets.tsv --cache CACHE --outdir OUT [--pilot] [--resume]` chains build-set -> normalize -> screening -> summary. `build-set --query` builds from a public UniProt query; `screen --set SETDIR` screens a previously built set.
- `run --resume` (legacy project branch) reuses only hash-verified bundles and preserves corrupted/failed attempts. Incomplete acquisition or processing is marked `partial` and reported as such, never as a full run (CLI exit code 3).
- Use `lab-evidence` for workbook recovery and reconciliation; read the [private evidence contract](references/lab-evidence-v060.md). Raw intake needs no actual scaffold. Preserve duplicate headers/unknown cells; quarantine suspect layouts and count source rows separately from constructs. Proposed gene/fragment links never confirm receiver identity. Cell-specific semantics plus audited constructs gate quantitative interpretation only. `inspect` and `map` remain low-level helpers.
- Use `normalize-uniprot` for local public annotations, and optional `fetch-uniprot` only for explicitly public accessions. Never upload private workbook contents, scaffold sequences or unpublished constructs to prediction services without specific authorization.
- Lab experience rules import via `--lab-rules` with provenance; when absent, outputs state 尚未纳入.

Before reporting completion run `make all-checks-offline` in the source package and `verify-run` for produced bundles. DNA synthesis orders, publishing and model training are separate actions, not implied by this skill.

For v0.4.2 reference handling and migration read [audit corrections](references/audit-v042.md): use the annotated Displayed isoform or an explicit `accession:canonical` marker, never assume `-1`. Other-isoform features and unmatched molecule-scoped comments remain recorded but are not transferred. Generic `Membrane` means unspecified location, not organelle-only. Candidate-level recommendations must not inherit the primary reference's positive class. Rebuild normalized inputs from preserved raw caches after parser changes; old hashes validate history, not current interpretation.

## Partial completion and stopping points
Missing scaffold stops fusion export only. Missing/ambiguous isoform or invalid coordinates block that candidate, not unrelated proteins. A single item's fetch, mapping or processing failure — including a corrupted cached protein record at screening time — is recorded as that entry's failure and marks the run partial; it never halts the batch. Uninterpretable assay fields stay unresolved. Conflicts remain visible. If a predictor or API is not run, report that fact. No trained ECD-SNIPR function predictor is included.
