# ECD-SNIPR Harness

An evidence-bounded, LLM-assisted workflow for antigen-fragment design and receiver-context review in the **Antibody-Sender -> Antigen-Receiver** configuration.

[![CI](https://github.com/KaixuanLi-sibcb/ecd-snipr-harness/actions/workflows/ci.yml/badge.svg)](https://github.com/KaixuanLi-sibcb/ecd-snipr-harness/actions/workflows/ci.yml)
[![Version](https://img.shields.io/badge/version-0.10.2-208B83.svg)](CHANGELOG.md)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![License](https://img.shields.io/badge/workflow-MIT-blue.svg)](LICENSE)

**Public release: workflow code, documentation, synthetic tests and public-method examples only.** Institutional source workbooks, actual constructs, assay labels, retrospective comparisons and generated analysis outputs are not published. External predictors and model weights are not redistributed.

## Overview

Given a defined human membrane-protein set or target list, the workflow proposes sourced antigen fragments for further evaluation at the antigen-recognition position of a SNIPR receiver. It reports exact sequences, reference coordinates, alternatives, rationale, conflicts and unresolved questions. It does not design antibodies or sender PDGFR-LC constructs.

The first pass is deterministic. An agent can orchestrate tools, investigate difficult cases and explain evidence; an LLM does not invent sequences, choose unchecked boundaries or approve constructs. This is not a trained ECD-SNIPR compatibility classifier.

| Independent layer | Question answered | Not established |
|---|---|---|
| Fragment design | Is there a supported candidate under current routes? | Expression or receptor function |
| Core evidence | Does the exact fragment agree with sequence, domain and predicted topology evidence? | Correct folding or preserved epitopes |
| Mechanism-priority review | Which specific signals warrant receiver-context checks? | Failure prediction or calibrated risk |
| Functional evidence | Which construct-matched endpoint measurements exist? | Success without defined acceptance criteria |

Read the [complete Chinese workflow and standards](references/workflow-and-standards.zh-CN.md). Formal technical contracts remain in [workflow](references/workflow.md), [screening criteria](references/screening-criteria.md) and [data contracts](references/data-contracts.md).

## Core capabilities

- Build a versioned UniProt set; preserve ambiguity, duplicate relationships, scope exclusions and per-entry failures.
- Select a canonical/reference or requested isoform without asserting laboratory isoform confirmation.
- Prefer full mature antigen forms; retain sourced alternatives and domain/epitope tradeoffs.
- Check coordinates, exact slicing, native SP/TM/intramembrane/cytoplasmic overlap, domain cuts and disulfide crossings.
- Cross-check optional public-accession InterPro/Pfam annotations and exact full-reference local predictions.
- Execute IUPred2A, DeepTMHMM2 and licensed SignalP 6 locally, with optional portable scheduled HPC workers.
- Export transparent design recommendations, separate review priorities, candidate FASTA and reconciled coverage.
- Optionally assemble a real human-reviewed scaffold and reconcile private measurements into four endpoint states.

## Workflow

~~~mermaid
flowchart TD
    A[Defined public set or user list] --> B[Identity resolution and reference selection]
    B --> C[Versioned fetch, cache and normalization]
    C --> D[Scope, topology and sourced mature boundaries]
    D --> E[Candidate generation and four-class design advice]
    E --> F[Exact fragment and structural-unit checks]
    P[Optional local SignalP / DeepTMHMM2 / IUPred2A] --> F
    I[Optional public InterPro/Pfam annotations] --> F
    F --> G[Independent mechanism-priority review]
    G --> H[Per-protein results, candidate plans, FASTA and coverage]
    H --> V[Independent raw / sequence / artifact verification]
    H -. Real scaffold and named human review .-> S[Optional full fusion assembly]
    S -. Matched constructs and defined assays .-> X[Expression, recognition, basal and induced endpoints]
~~~

Absent scaffold, review or experiments do not block candidate delivery. They do block unsupported fusion export or functional claims. Frozen designs are never silently edited by a prediction overlay.

## Screening standards

| Recommendation | Interpretation |
|---|---|
| `standard_candidate` | Positive sequence/topology/boundary evidence supports a routine proposal |
| `conditional_candidate` | A sourced proposal exists with a named orientation, structural-unit or native-context issue |
| `no_standard_route` | No supported proposal under current routes; not an assertion the protein can never be used |
| `insufficient_evidence` | Core identity/topology/boundaries are insufficient or conflicting; not biological failure |

Scope and processing status are separate. Set membership, natural cell-surface localization and design scope are distinct. Secreted extensions use a separate denominator.

| Topology | Current route |
|---|---|
| Type I | Sourced continuous mature ECD; remove native SP/TM/tail |
| Type II | Sourced external region plus attachment-orientation review |
| GPI | Sourced mature boundary including omega residue; remove downstream GPI signal |
| Multi-pass | No loop concatenation; only sourced domain alternatives in one annotated external interval |
| Processed / multi-chain | Preserve product scope, dependencies and omitted repertoire |
| Organelle / uncertain location | Lumen-facing topology is not natural cell-surface proof |

Coordinates are **1-based inclusive on the exact selected reference**. Invalid coordinates, sequence mismatch and retained annotated native exclusions remain blocking. Length, cysteine, glycans, hydropathy and disorder are descriptors/review signals, not success thresholds. Primary selection uses explicit structural-integrity, mapped-epitope, antigen-form and stable-ID tradeoffs, not a weighted score.

Review tiers are elevated review, context-dependent review, sequence-only alert, no specific signal and not assessed. **No detected signal is not low functional risk.** Annotated and predicted conflicts retain different evidence types; neither establishes self-activation.

## Quick start

Python >= 3.10; core/offline tests use the standard library only.

~~~bash
git clone https://github.com/KaixuanLi-sibcb/ecd-snipr-harness.git
cd ecd-snipr-harness
make all-checks-offline

# Offline synthetic end-to-end run
python3 scripts/ecd_snipr_cli.py smoke --outdir outputs/demo

# Public pilot; network or populated cache required
python3 scripts/ecd_snipr_cli.py screen --list examples/public_core_pilot.tsv \
  --cache cache --outdir outputs/public-pilot --pilot
~~~

Use the actual immutable bundle path printed by the command, not a guessed ID. `verify-run --run-dir BUNDLE` checks every exported artifact.

## Defined membrane-reference set

~~~bash
python3 scripts/ecd_snipr_cli.py build-set \
  --query '(organism_id:9606) AND (reviewed:true) AND (keyword:"Membrane" OR keyword:"Cell membrane")' \
  --cache cache --outdir sets/human-membrane
python3 scripts/ecd_snipr_cli.py screen --set sets/human-membrane \
  --outdir outputs/human-membrane --resume
~~~

This defines membrane-related reviewed human entries, not all genes/isoforms or only plasma-membrane targets. Release, query, reference policy and completeness are recorded. Limited queries and unresolved/fetch/processing gaps stay partial. Input rows, genes, references and fragments are different units. Installation and tests do not launch a full-set task.

## Local tools and independent cross-checks

~~~bash
python3 scripts/ecd_snipr_cli.py configure-predictors \
  --iupred-python /outside/repo/iupred-env/bin/python \
  --dtm-python /outside/repo/dtm-env/bin/python \
  --dtm-models /outside/repo/dtm-models \
  --signalp-python /outside/repo/signalp-env/bin/python \
  --signalp-models /outside/repo/signalp-models \
  --output /outside/repo/predictor-config.json

python3 scripts/ecd_snipr_cli.py validate-local-tools \
  --run-dir SCREENING_BUNDLE --reference-set sets/human-membrane \
  --config /outside/repo/predictor-config.json \
  --outdir outputs/tool-validation --timeout 3600 --resume
~~~

See [local setup](references/local-predictors-v0100.md) and [portable HPC execution](references/hpc-execution.md). Importing predictions is not running a model. DeepTMHMM2 is not legacy DeepTMHMM. SignalP requires an official licensed package supplied by the user. Executables, code/data, models, original FASTA and raw outputs are hash-bound. Unsupported residues are retained without substitution.

Optional public-accession InterPro enrichment does not submit sequences:

~~~bash
python3 scripts/ecd_snipr_cli.py core-evidence --run-dir SCREENING_BUNDLE \
  --reference-set sets/human-membrane --outdir outputs/core-review \
  --interpro-cache cache/interpro --interpro-live
python3 scripts/verify_core_evidence.py --run-dir CORE_BUNDLE \
  --reference-set sets/human-membrane
~~~

Omit `--interpro-live` for cache-only use. This retrieves existing matches, not InterProScan execution. The independent verifier re-derives sequence/domain retention; it is a deterministic second pass, not independent biological validation.

## Outputs

| Artifact | Purpose |
|---|---|
| `protein_screening.tsv` | Per-object scope, topology, class, primary proposal and rationale |
| `candidate_plan.tsv`, `candidates.json`, `candidate_fragments.fasta` | Exact primary/alternate fragment plans and sequences |
| `candidate_comparison.tsv`, `antigen_context.tsv` | Structural, epitope and mature-product tradeoffs |
| `candidate_core_evidence.tsv`, `candidate_domain_coverage.tsv` | Fragment cross-checks and source coverage |
| `receiver_risk.tsv`, `receiver_risk_signals.tsv` | Independent source-bounded review priorities |
| `receiver_function.tsv`, `receiver_review_plan.tsv` | Four endpoint states and outstanding checks |
| `summary.json`, `PI_SUMMARY.md`, figure/source data | Coverage with explicit denominators, never success rate |
| Manifests, state/events and source records | Reproducibility, partial/error states and hashes |

Optional [private evidence recovery](references/lab-evidence-v060.md) preserves sheet/row/cell hashes. An experiment-performed flag is not a success label; positive antibody-screening records are not complete receiver-function validation. Inputs and recovered outputs stay private.

## Validation and limitations

Offline checks cover CLI/contracts, synthetic smoke, unit tests, scheduled-worker contracts, privacy and allowlisted packaging. CI runs on Python 3.10 and 3.13 without private workbooks, live APIs or vendor models. Release receipts are recorded in the [changelog](CHANGELOG.md).

Software tests, annotation agreement and vendor-CLI concordance are **not SNIPR biological accuracy**. A valid benchmark needs actual construct/scaffold/sender/condition matching and predefined expression, recognition, basal and induced labels. Previously inspected development cases are not held out. No accuracy, sensitivity, specificity or calibrated probability is claimed.

Epitope preservation, glycan occupancy, actual-fusion geometry/processing and experimental acceptance remain unestablished. HPA/GTEx, Open Targets/ChEMBL, therapeutic ranking and human-mouse preference are outside this receiver-engineering workflow. InterPro is optional; no new AlphaFold/ESM-derived SNIPR model or full-fusion predictor is claimed.

## Installation, development and governance

~~~bash
make all-checks-offline
make install-user
# Default destination: ~/.agents/skills/ecd-snipr-harness
~~~

Installation validates a staged allowlisted copy and backs up the previous skill. Analysis outputs are not installed. Source data, caches, predictor environments and licensed assets must stay outside the repository. Automated checks supplement, not replace, content/history review.

~~~text
scripts/ecd_snipr/   deterministic core and evidence adapters
scripts/hpc/         scheduled workers, reconciliation, CLI concordance
references/         workflow, criteria, contracts and deployment
schemas/            machine-readable contracts
tests/ examples/    synthetic regressions and public-method examples
agents/ SKILL.md    bounded agent instructions
~~~

## Next milestone

Freeze rules/versions, independently review representative designs and collect construct-matched four-endpoint measurements under predefined acceptance criteria. Evaluate by protein/family or prospectively. Targeted structural/epitope analyses should answer specific unresolved questions, not optimize a generic score against a known failure list.

## Methodological origin and license

The workflow formalizes antigen-fragment engineering and experimental-review practices developed through internal discussions in the Meng Lab, Shanghai Institute of Biochemistry and Cell Biology, Chinese Academy of Sciences. It is a reproducible scaffold, not a claim of biological validation.

Workflow code is [MIT licensed](LICENSE). External software, data and weights retain their own terms; this license does not permit their redistribution.
