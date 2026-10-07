#!/usr/bin/env bash
set -euo pipefail

# SGE example: compute only after a scheduler allocation, never on a login host.
: "${JOB_ID:?Submit this script with qsub; do not execute it on a login host}"
: "${WORKFLOW_HOME:?Absolute workflow checkout path required}"
: "${CORE_PYTHON:?Core Python executable required}"
: "${JOB_MANIFEST:?Frozen jobs JSON required}"
: "${PREDICTOR_CONFIG:?Hashed predictor configuration required}"
: "${PREDICTION_OUT:?Private/local output root required}"
: "${TOOL:?Choose IUPred2A, DeepTMHMM2 or SignalP}"

export PYTHONPATH="$WORKFLOW_HOME/scripts${PYTHONPATH:+:$PYTHONPATH}"
cd "$WORKFLOW_HOME"
if [[ "$TOOL" == SignalP ]]; then
  exec "$CORE_PYTHON" scripts/hpc/run_signalp.py \
    --jobs "$JOB_MANIFEST" --config "$PREDICTOR_CONFIG" \
    --outdir "$PREDICTION_OUT" --threads "${NSLOTS:-1}" \
    --batch-size "${BATCH_SIZE:-512}"
fi
case "$TOOL" in IUPred2A|DeepTMHMM2) ;; *) echo 'Unsupported tool' >&2; exit 2;; esac
TOOL_PYTHON="$("$CORE_PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["tools"][sys.argv[2]]["python"])' "$PREDICTOR_CONFIG" "$TOOL")"
SHARD="${SGE_TASK_ID:-1}"
if [[ "$SHARD" == undefined ]]; then SHARD=1; fi
exec "$TOOL_PYTHON" scripts/hpc/run_tool_v2.py \
  --jobs "$JOB_MANIFEST" --config "$PREDICTOR_CONFIG" \
  --outdir "$PREDICTION_OUT" --tool "$TOOL" \
  --shard "$SHARD" --shards "${SHARDS:-1}" --device "${DEVICE:-cpu}"
