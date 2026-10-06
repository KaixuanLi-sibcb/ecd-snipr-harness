# Receiver mechanism-priority triage v0.8.0

v0.8.1 / policy 1.1 recognizes the general `self-ligand` synonym as native
self-association context, including negation and fragment-scope checks. The gap was
noticed during a developmental audit after the first frozen cohort evaluation;
that original evaluation is preserved. The corrected output must not be presented
as independent blind or held-out validation. No gene list or numeric threshold was added.

## Purpose and decision contract

Separate fragment design, mechanism-priority triage and measured receiver function.
The target is Antibody-Sender -> Antigen-Receiver. Candidate generation, primary
selection, hard coordinate/TM/SP/tail checks and four endpoint measurement contracts
are unchanged. A routine fragment can warrant elevated receiver review without
being relabeled a failed design. No success probability is produced.

Positive signals need exact reference-sequence correspondence and source/version
provenance. Reference identities and laboratory outcomes are not decision features.
The algorithm retains the evidence type; curated entries are not upgraded to
experimental support. Inherited comment summaries retain global entry provenance
and available ECO codes; they do not supply per-text PMID or localized interaction
sites when the normalized reference lacks those details.

| Priority | Transparent trigger | What remains unproved |
|---|---|---|
| Elevated review | Retained annotated processing site, interchain disulfide, interrupted annotated structural unit, or explicit extracellular/homophilic association for full topological interval | Cleavage, folding, clustering, low baseline and inducibility in the new receiver |
| Context-dependent | Native unlocalized shedding/homo-association; native glycan/immunoglobulin binding; exact-fragment topology prediction conflict | Site retention, ligand availability, culture effect, protease and force coupling |
| Sequence only | Local hydropathy, S/T/P enrichment or low window entropy | TM identity, aggregation, disorder, glycan occupancy, SNIPR activity |
| No specific signal | None detected by the methods actually run | Absence of risk, general compatibility or negative evidence |
| Not assessed | Unusable candidate, sequence/source mismatch or reference conflict | Any functional conclusion |

The earliest applicable priority wins, not a sum of risk counts. Native homodimer
wording without an ECD interface remains context-dependent. Explicit intracellular
or transmembrane association is excluded from ECD transfer. Negated or mixed text
is retained unresolved, not a positive signal. Sites outside the fragment cannot
trigger retained-site risk. Generic missing literature, length, odd cysteine and
potential N-glycosylation are not elevated-risk triggers. Regex extraction is a
conservative triage aid; mixed statements, by-similarity assignments and unlocalized
interactions need review. Native homophilic binding is not proof that homotypic
receiver contacts activate the regulatory core.

Four endpoints remain independent: surface expression, recognition retention,
basal activation and induced response. Unknown backbone/host/expression/controls
remain missing. Tags, medium and sender matter for potential glycan/Fc interactions.
Receiver confidence is not established by these signals. Changing or removing an
interaction/processing motif may destroy epitopes; no mutations or automatic
truncations are proposed by this risk layer.

## Algorithms and integrations actually available

| Method | Execution | Interpretation and limits |
|---|---|---|
| Kyte-Doolittle | Local sliding mean, 19 residues, alert >=1.6; average GRAVY also exported | Hydrophobicity descriptor; no TM/aggregation classification or functional calibration |
| Shannon window entropy | Local 12-residue windows <=2.2 bits | Composition bias; not the SEG algorithm, disorder, LC-driven condensation or SNIPR prediction |
| S/T/P composition | Local 30-residue windows >=0.65 fraction | Stalk/glycosylation geometry question, not glycan occupancy or failure cutoff |
| Biopython ProtParam 1.86 | Optional pinned extra, local CPU; MW, pI, charge and independent KD profile cross-check | Physicochemical properties, not fold or receptor-function prediction. No instability-index pass/fail rule |
| IUPred2A / IUPred3 / AIUPred | Import normalized local results, one score per exact fragment residue | Disorder prediction only. No automatic invocation, model distribution or functional threshold |
| DeepTMHMM / SignalP | Import normalized local predicted TM/SP regions | Topology conflict prompts review, never silent replacement of annotation |

Window cutoffs are explicit uncalibrated descriptor alerts, not learned SNIPR
thresholds. A sequence alert alone cannot become elevated review priority. Short
sequences have no forced window. No GPU or online sequence upload is required.
Biopython uses the Biopython License Agreement/BSD-3-Clause dual license; external
predictors retain their own terms. The IUPred site describes academic availability
and commercial enquiry requirements. Do not redistribute their models as MIT assets.

Normalized predictor JSON is a list; each record requires `candidate_id`,
`sequence_sha256` (SHA256 of plain uppercase fragment bytes), `software`, `version`,
`parameters`, `source`, `input_scope=candidate_fragment`,
`coordinate_system=1-based-inclusive`, `raw_output_path` and `raw_output_sha256`.
The CLI checks the actual local raw file hash. Disorder tools need finite 0..1
`scores` of exact fragment length; topology tools need `regions` with `kind`,
`start`, `end`. Coordinates from a full reference must not be passed as fragment
coordinates. No raw-text format is guessed. Unknown candidate IDs abort the import;
invalid record content is retained as rejected evidence. The hash verifies provenance,
not the biological correctness of the predictor or hand-authored normalization.

## Harness and output

`screen`/`run` export `receiver_risk.tsv`, `receiver_risk_signals.tsv`,
`receiver_risk_summary.json` and `candidates.json.receiver_risk` by default.
`protein_screening.tsv` includes the primary candidate's separate review priority.
Neither original screening recommendation nor sequence changes.

`receiver-risk --run-dir RUN --reference-set SET --outdir NEW [--tool-evidence FILE]`
reads a verified historical bundle and exact hash-verified reference records. It
creates a new content-addressed bundle, binds rule/engine/input hashes, retains
single-reference errors and marks analysis partial. Resume rechecks reference hashes.
Each reference read uses a bounded child process (five-second timeout) so a blocked
file-provider/cache read cannot hang the batch. Timeouts remain explicit missing
annotation coverage; they cannot be converted into low risk. Reference bytes are
hash-checked again before export. This bounds I/O, not biological computation.
At most four readers run concurrently; acquisition order and candidate decisions
remain deterministic. No additional biological rules are introduced by I/O concurrency.
Original bundles are checked before/after and never edited. Missing records do not
drop candidates. `sequence-tools` writes a separate immutable local descriptor
bundle; unavailable optional tools are explicit, not silently simulated.

```bash
python3 scripts/ecd_snipr_cli.py receiver-risk --run-dir RUN --reference-set SET --outdir NEW --resume
python3 scripts/ecd_snipr_cli.py sequence-tools --run-dir RUN --outdir NEW_TOOLS --biopython
python3 scripts/ecd_snipr_cli.py verify-run --run-dir NEW_BUNDLE
make receiver-risk-fixture
make all-checks-offline
make install-user
```

Use new output directories; workbooks, sequences, overlays and tool outputs are
private by default and excluded from release assets. Installation is allowlisted,
backs up the old skill and does not push to GitHub.

## Validation and anti-leakage

Freeze the policy plus implementation hash before importing outcome labels; run it
unchanged across all eligible public candidates, not just observed failures. Report
flag prevalence, annotation coverage and candidate-versus-reference denominators.
Private failures are developmental retrospective checks because examples were
already seen, not blinded/held-out tests. Do not alter thresholds after inspecting
the matches. Generic alerts and an unknown state cannot count as self-activation
prediction hits. Compare elevated and context-only tiers separately.

Without eligible independent successes, matched actual constructs, assay definitions
and predefined labels there is no accuracy, specificity, PPV or calibrated sensitivity
estimate. Antibody-screening outcome and SNIPR Assay Y/N are not receiver-functional
labels. True validation needs frozen rules and prospective protein/family-grouped
constructs, multiple backbone/expression conditions and matched controls.

## Primary sources and scope

- [Zhu et al., Cell 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9108009/): synthetic receptor modules and proteolysis regulate signaling stringency. Its regulatory ECD is not synonymous with this project's antigen fragment; it supplies motivation, not antigen-receiver validation.
- [Enhanced synthetic Notch, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7069970/): expression-dependent ligand-independent activation warrants matched expression/control context, not a universal sequence cutoff.
- [UniProt disulfide documentation](https://www.uniprot.org/help/disulfid): distinguishes intra/interchain and predicted/experimentally supported annotation; chain association is not functional failure.
- [Biopython ProtParam documentation](https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html) and [source/license](https://github.com/biopython/biopython): available deterministic sequence properties; no SNIPR outcome model.
- [NCBI masking guidance](https://www.ncbi.nlm.nih.gov/books/NBK569845/): SEG is a specific low-complexity masking algorithm. The built-in Shannon descriptor is deliberately not claimed to be SEG.
- [IUPred official download and terms](https://iupred.elte.hu/download): optional future/local disorder evidence, not run by default.

Sources checked 2026-10-06. No source establishes a general antigen-ECD self-activation classifier.
