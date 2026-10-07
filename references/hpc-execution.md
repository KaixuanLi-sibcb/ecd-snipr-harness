# Portable local and scheduled predictor execution

This guide publishes workflow helpers only. It contains no institutional paths,
accounts, job receipts, reference cohorts, private sequences or licensed software.
Workers preserve the candidate engine: predictions are evidence overlays, not
automatic replacement boundaries or a receptor-function model.

## 1. Prerequisites and scope

Use Python >=3.10 for the workflow, separate predictor environments and model
directories **outside the repository**. Follow [local installation guidance](local-predictors-v0100.md).
The core/CI needs none of these tools. Request SignalP directly from
[DTU](https://services.healthtech.dtu.dk/services/SignalP-6.0/) and comply with the
received academic/commercial terms. No download token, license, binary or weight
belongs in Git. IUPred authors' terms and model licensing are separate from MIT
workflow code.

The public worker has two execution modes:

| Method | Entry | Behavior |
|---|---|---|
| Conservative local path | `validate-local-tools` | Vendor CLI for deep tools, CPU/per-reference subprocesses, complete-bundle reuse |
| Scheduled path | `scripts/hpc/` | IUPred model loaded once, DeepTMHMM2 model-loaded per shard, SignalP CPU batches; raw reimport and candidate filtering |

On institutional HPC, login/submit hosts may run only lightweight preparation,
submission and status checks. Use `qsub` for inference, environment compilation,
large verification and overlays on allocated compute nodes. Queue names, parallel
environments, GPU/memory resource syntax, modules and CUDA versions are site-specific;
obtain them from the administrator. This repository deliberately does not guess them.

## 2. Pin external implementations

The compatibility worker was developed against the following author-source commits:

- [DeepTMHMM2](https://github.com/fteufel/DeepTMHMM2/tree/b05e27adf50738405e1de0b4a6b7072bb145fd3a):
  `b05e27adf50738405e1de0b4a6b7072bb145fd3a`.
- [Packaged IUPred2A](https://github.com/saezlab/iupred/tree/655d42c28c571cc42e12f93588bab4dd727a48cf):
  `655d42c28c571cc42e12f93588bab4dd727a48cf`.
- SignalP: official fast academic package supplied by the user; record actual
  distribution/module version (`signalp6` / `signalp`) and model hashes, not a
  version inferred from the archive filename.

DeepTMHMM2 is not legacy DeepTMHMM. Complete the ten checkpoints and ESM2 assets
during setup; no implicit downloads during analysis. Do not use MPS. CPU and CUDA
have different resource/performance requirements; CUDA unavailable is an error,
not a silent CPU fallback. New upstream releases must pass a fresh concordance
pilot before adopting this API compatibility path.

At the pinned DeepTMHMM2 commit, the Python API supplied a one-element batch to
scalar postprocessing. `dtm2_api_compat.py` reproduces the author CLI's single-batch
unpacking and calls vendor inference/ensemble/postprocessing functions. It does
not edit vendor models or change label rules. This compatibility path has synthetic
contract tests; **run `check_dtm_cli.py` on your installed tool and representative
public pilot before full inference**. Header, exact reference and every topology
residue must agree. Concordance is engineering validation, not biological accuracy.

Store source commits, `pip freeze`, platform/CUDA versions and received package
metadata privately with each installation. `configure-predictors` separately
hashes configured software, executables and model files; a configured status is
not an executed-tool validation result.

## 3. Prepare a declared reference scope

First run the documented `build-set` and `screen`; keep their actual returned
paths. The jobs manifest is generated outside the reference set and never overwrites
an existing file. Duplicate inputs are linked, unresolved source rows preserved,
bad references isolated and unsupported residues retained without substitution.

```bash
python3 scripts/hpc/prepare_jobs.py --reference-set /local/analysis-set \
  --output /local/prediction-work/prediction_jobs.json

# Separate small execution/concordance pilot; remains partial against full scope
python3 scripts/hpc/prepare_jobs.py --reference-set /local/analysis-set \
  --output /local/prediction-work/pilot_jobs.json --pilot-limit 4
```

The strict 20-AA workflow contract excludes U or other unsupported characters
explicitly. This is not a claim that every vendor tool intrinsically rejects them.
Do not substitute U with C or truncate long sequences to force completion.

Configure local installations with the exact full command in the README. Prefer
absolute shared paths for the set, config, environments, code and outputs. A
configuration created on another machine must not be reused with different assets.

## 4. Scheduler worker example

`examples/hpc_worker.sh` requires an SGE job allocation and explicit paths; it
does not contain a lab configuration. Supply the appropriate **site-approved CPU
or GPU queue and resources** when invoking qsub. The following exports are examples,
not institutional settings:

```bash
export WORKFLOW_HOME=/shared/code/ecd-snipr-harness
export CORE_PYTHON=/shared/envs/core/bin/python
export JOB_MANIFEST=/shared/work/pilot_jobs.json
export PREDICTOR_CONFIG=/shared/work/predictor-config.json
export PREDICTION_OUT=/shared/work/pilot-results
export TOOL=DeepTMHMM2 DEVICE=cuda SHARDS=2
# Set QUEUE, PARALLEL_ENV and site-specific GPU/memory resources first.
qsub -cwd -q "$QUEUE" -pe "$PARALLEL_ENV" 4 -t 1-2 \
  -v WORKFLOW_HOME,CORE_PYTHON,JOB_MANIFEST,PREDICTOR_CONFIG,PREDICTION_OUT,TOOL,DEVICE,SHARDS \
  examples/hpc_worker.sh
```

For IUPred2A set `TOOL=IUPred2A`, `DEVICE=cpu`, select a CPU queue and matching
shard array. For SignalP set `TOOL=SignalP`, select a CPU queue and submit **one**
non-array job; `SHARDS=1`, `BATCH_SIZE=512` if desired. Never run duplicate SignalP
jobs against one output root. Logs and raw predictions are local/private outputs.
Change to the full jobs manifest only after the pilot and CLI concordance pass.

Raw worker CLI options are inspectable with `--help`. If a scheduler is not SGE,
write its allocation wrapper around the same Python workers; do not pretend the
SGE wrapper supports another scheduler. Reserve enough resources for long sequences.
Timeout/errors preserve failed attempts and require a new output root to retry;
do not overwrite an accepted result.

## 5. Require native CLI concordance for DeepTMHMM2

Run this on an allocated compute node with a new output directory and the same
device/assets as the adapted pilot:

```bash
python3 scripts/hpc/check_dtm_cli.py --jobs /shared/work/pilot_jobs.json \
  --config /shared/work/predictor-config.json \
  --adapted-root /shared/work/pilot-results/DeepTMHMM2 \
  --outdir /shared/work/new-cli-concordance --device cuda
```

This executes the unmodified vendor CLI. A mismatch stops acceptance; do not
discount it as an annotation conflict or proceed to a full run.

## 6. Reconcile and overlay, independently of scheduler status

Run these on an allocated CPU node. A finished job or zero exit code alone does
not establish complete coverage. Reconciliation requires the frozen jobs contract,
source hashes, raw FASTA/output and parser reproduction for every requested tool.

```bash
python3 scripts/hpc/reconcile.py --jobs /shared/work/prediction_jobs.json \
  --run /shared/work/full-results --outdir /shared/work/new-accepted \
  --tools IUPred2A DeepTMHMM2 SignalP
python3 scripts/ecd_snipr_cli.py verify-predictions --run-dir /shared/work/new-accepted
python3 scripts/hpc/finalize.py --prediction-dir /shared/work/new-accepted \
  --screening-run /shared/work/FROZEN_SCREENING_BUNDLE \
  --overlay-outdir /shared/work/new-candidate-review
```

`finalize.py` filters predictions to explicit candidate references and counts the
remaining records; it never silently ignores unknown candidate IDs. The reference
set, original candidate sequences and recommendations remain unchanged. Verify
the returned overlay using `verify-run` and `verify_core_evidence.py`.

`reconcile.py` exits 0 when every *eligible* requested prediction reimports,
otherwise 1; exclusions, unresolved source rows, pilot scope or partial source
still keep `completeness=partial`. `finalize.py` exits 3 for a partial predictor
scope or partial overlay. These semantics distinguish eligible completion from full-scope completion.
Keep partial outputs with their explicit denominators; do not call them full success.

## 7. Offline release checks

```bash
make hpc-fixture
make all-checks-offline
```

These test synthetic raw-result reconciliation, tamper detection, duplicate/missing
reference accounting, partial scopes, unsupported sequences and immutable overlays.
They do not execute vendor models in CI. Actual vendor installation, native CLI
concordance, throughput and hardware suitability require a separate local pilot.

Never publish work manifests, per-reference outputs, predictor configs, environment
receipts, internal host/account metadata, experimental labels or vendor packages.
