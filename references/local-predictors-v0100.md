# Local predictor execution and evaluation

## Purpose and limits

This layer checks signal-peptide, transmembrane/topological and disorder evidence
on the exact **full analysis-reference sequence**, then projects the results onto
frozen ECD candidates. It does not redesign antibodies, guess the receiver scaffold,
automatically trim a candidate, or predict SNIPR basal/induced activation.

Native full-reference predictions and the actual fusion are different inputs.
Removing an annotated native signal peptide does not remove the need for the real
fusion's signal peptide. A topology conflict is a review question, not measured
failure. Disorder can be functional or context-dependent and is not an automatic
self-activation warning. Reference-sequence annotation agreement is not an
independent experimental benchmark and must not be called accuracy.

## Installations outside the skill

Core/offline CI has no predictor dependency. Use separate environments outside the
repository. Do not package third-party software, licensed weights or private inputs.

| Tool | Setup and actual role | Constraints |
|---|---|---|
| IUPred2 | Authors' [packaged implementation](https://github.com/saezlab/iupred), `iupred iupred2`; CPU energy-based disorder | CC-BY-NC-ND-4.0 as supplied by authors; retain unmodified original code. Installed package version is recorded, not mislabeled as IUPred3. |
| DeepTMHMM2 | [Author repository](https://github.com/fteufel/DeepTMHMM2), `dtm2`; full-reference topology with `--simplify-io` | Explicitly version 2, not legacy DeepTMHMM. CPU on Apple Silicon; author warns of incorrect MPS embeddings. ESM2 650M and ten checkpoints needed. Pyproject declares MIT; do not redistribute assets in this skill. |
| SignalP 6 | [DTU academic package](https://services.healthtech.dtu.dk/services/SignalP-6.0/), `signalp6`; eukarya/fast local analysis | Academic email/license download plus model weights required. Package absence is not deployment success. Official instructions require Python <=3.10, torch<2 and numpy<2 for 6.0h; follow the received package's instructions. Slow-sequential is a later boundary-review option. |

Example local installation (pin the retrieved Git commits in your setup log):

```bash
python3 -m venv /path/outside/repo/iupred-env
/path/outside/repo/iupred-env/bin/python -m pip install /path/to/unmodified/iupred-source
python3 -m venv /path/outside/repo/deeptmhmm2-env
/path/outside/repo/deeptmhmm2-env/bin/python -m pip install /path/to/unmodified/DeepTMHMM2-source
```

DeepTMHMM2 checkpoints are available in the author repository. Complete the ESM2
download **during setup**, not during batch analysis:

```bash
TORCH_HOME=/path/to/DeepTMHMM2/checkpoints/torch \
  /path/outside/repo/deeptmhmm2-env/bin/python -c \
  "from esm import pretrained; pretrained.load_model_and_alphabet('esm2_t33_650M_UR50D')"
```

This retrieves public model assets; no sequence is submitted. Configuration
requires both ESM2 files under `<model_dir>/torch/hub/checkpoints` and hashes all
model files. CPU execution is conservative for local use; CPU throughput and long
sequence memory demand must be measured on a pilot, not extrapolated as completion.

For SignalP, request the DTU fast package with your academic email; after receiving
it, install into a compatible separate environment and point to its `models`
directory. This workflow never fills an academic-license form with invented user
details. Validate the adapter against your actual received package before declaring
installation accepted. Positive GFF3 and explicit summary SP/OTHER formats are supported; blank
GFF3 alone is not a checked no-SP prediction.

## Configure and run

```bash
python3 scripts/ecd_snipr_cli.py configure-predictors \
  --iupred-python /path/to/iupred-env/bin/python \
  --dtm-python /path/to/deeptmhmm2-env/bin/python \
  --dtm-models /path/to/DeepTMHMM2/checkpoints \
  --output /private/output/predictor-config.json
# Add --signalp-python and --signalp-models only after real installation.

python3 scripts/ecd_snipr_cli.py validate-local-tools \
  --run-dir /path/to/verified/screening/run \
  --reference-set /path/to/analysis-set \
  --config /private/output/predictor-config.json \
  --outdir /private/output/tool-validation \
  --interpro-cache /path/to/public/interpro-cache \
  --timeout 3600 --resume
```

Default requested tools include SignalP: missing it gives `partial`/exit 3, not
success. `--tools IUPred2A DeepTMHMM2` explicitly runs only those two; completion
then applies only to that declared scope. `predict-local` alone covers all resolved
references, including references without ECD candidates. `--limit N` is a pilot
subset and stays partial against the original set. No sequence is truncated to fit
a model or memory limit. A failed item remains `error`; it cannot become a negative.

IUPred uses a local worker with one package import and per-reference error isolation.
DeepTMHMM2 and SignalP use timeout-bounded per-reference subprocesses with raw
stdout/stderr and exact command records. This is robust but repeatedly loads deep
models. The [optional scheduled path](hpc-execution.md) adds model-loaded GPU shards
and SignalP CPU batches, with its own pilot/concordance requirements. No universal
throughput is claimed. Completed immutable batches may be reused only
after raw-file, reference, code and model checks; interrupted attempts are preserved
and need a new output root. This is completed-bundle reuse, not unfinished-item resume.

Outputs include `predictor_status.tsv`, `tool_evidence.json`, `predictor_summary.json`,
raw per-reference input/output/logs, the original manifests and a derived core/risk
overlay. Predictions without any candidate reference are separately counted, not
silently passed into candidate assessment. `verify-predictions --run-dir RUN` checks
nested outputs and original FASTA files; `verify-run` checks the derived overlay.

## What can be evaluated now?

1. **Software correctness:** actual execution, exact reference identity, one score
   or topology label per residue, coordinate projection, hashes and error isolation.
2. **Annotation consistency:** differences from UniProt and InterPro, with explicit
   missing annotations and prediction-only evidence. Not biological accuracy.
3. **Retrospective lab comparison:** only after matching the actual fragment,
   receptor/scaffold, sender, host, batch and endpoint. Positive antibody clones and
   `SNIPR Assay=Y` do not establish successful receiver function. Failure descriptions
   are not interchangeable endpoint labels. Previously examined cases are development
   examples, not held-out data.
4. **Functional benchmark:** freeze rules before testing; collect independently
   labeled surface-expression, binding-retention, basal and induced-response outcomes
   under specified conditions. Split by protein/homolog family; report confusion
   matrices, sensitivity/specificity and uncertainty per endpoint. No valid label set
   means no accuracy number, even if every program executed successfully.
