# ECD–SNIPR Harness

**Evidence-bounded antigen-fragment design for SNIPR receivers — at the scale of the human membrane proteome.**

`ecd-snipr-harness` turns "a new membrane target arrived; which fragment should we use, and why?" into a batch-executable, fully traceable pipeline. For every target it delivers a concrete candidate fragment (coordinates + sequence + rationale + open issues), a screening recommendation, and a reconciled summary across the whole set.

Working configuration: **Antibody-Sender → Antigen-Receiver** — the antibody sits on the sender cell; the antigen fragment occupies the antigen-recognition position of the SNIPR receiver. The tool designs the *fragment*, not the antibody and not the sender construct.

## What it does

- Builds a defined analysis set from a user target list or a public UniProt query (identity resolution with recorded evidence; ambiguity preserved, never guessed, no row dropped).
- Fetches, caches and normalizes UniProt annotations (topology, localization, signal peptide / propeptide / transmembrane / topological domains / domains / processing products), with content-addressed caching, resume, and per-item failure isolation.
- Selects the analysis reference explicitly (database canonical/reference sequence, rationale recorded — distinct from experimental isoform confirmation, which is still required before any fusion assembly).
- Generates source-bounded candidate antigen fragments by topology: type-I/II full ECDs, mature GPI forms, sourced shed forms, and sourced domain alternates for multi-pass proteins (extracellular loops are never stitched). No random window sliding, no default mutation "optimization", no mechanical truncation.
- Classifies each reference with a transparent decision table — `standard_candidate` / `conditional_candidate` / `no_standard_route` / `insufficient_evidence` — with `reason_codes`, evidence, and missing information attached. Out-of-scope objects (`scope_status`) and technical failures (`processing_status`) are tracked separately, never disguised as biological conclusions.
- Emits the research deliverables: `protein_screening.tsv`, `candidate_plan.tsv`, `candidates.json`, `candidate_fragments.fasta`, `summary.json` (all ratios carry explicit denominators), `PI_SUMMARY.md`, and a coverage figure with its source data.

## What it does not do

- It does **not** predict experimental success. Screening coverage ≠ functional validation; surface expression, antibody recognition, basal activity and induced response remain separate experiments.
- It does **not** fabricate fusion sequences: without a real, human-reviewed SNIPR scaffold, only candidate fragments and junction notes are delivered.
- It does **not** run per-protein LLM reasoning, structure prediction (AlphaFold/ESM) or paid APIs in the batch pass; the first pass is fully deterministic.
- Blocking checks (coordinate errors, retained native TM/SP/cytoplasmic tail) stay hard failures; length/cysteine/glycosylation-motif notes are warnings, not thresholds.
- Optional lab-experience rules (`--lab-rules`) may annotate or downgrade, never upgrade; when absent, outputs state "not yet incorporated".

## Requirements

Python ≥ 3.10, standard library only. No dependencies, no GPU, no network needed for the offline test suite.

## Quickstart

```bash
# Offline synthetic smoke test (no network)
python3 scripts/ecd_snipr_cli.py smoke --outdir /tmp/ecd_smoke

# Screen a target list (accessions or gene names; every row preserved)
python3 scripts/ecd_snipr_cli.py screen --list examples/target_list_example.tsv \
    --cache cache --outdir runs/pilot --pilot

# Build a defined set from a public UniProt query, then screen it
python3 scripts/ecd_snipr_cli.py build-set \
    --query '(organism_id:9606) AND (reviewed:true) AND (keyword:"Membrane" OR keyword:"Cell membrane")' \
    --cache cache --outdir sets/human_membrane
python3 scripts/ecd_snipr_cli.py screen --set sets/human_membrane --outdir runs/human_membrane
```

Incomplete retrievals are marked `partial` and never reported as complete. `--resume` continues interrupted batches from the cache.

## Tests

```bash
python3 -m unittest discover -s tests -p 'test_*.py'
```

124 checks (v0.3.0), including: screening completes with no scaffold/review/experiments; public-reference selection vs experimental-isoform confirmation stay distinct; missing optional risk literature does not blanket-downgrade; core sequence errors remain blocked; topology/mature-boundary/domain-alternate handling; multi-reference/multi-candidate/duplicate/out-of-scope counting; single-item failure isolation; legacy scaffold assembly and four-endpoint linkage unregressed.

## Repository layout

```
scripts/ecd_snipr/     core package (acquisition · screening · design · reporting · provenance …)
scripts/ecd_snipr_cli.py    command-line entry point
schemas/               input/claim JSON schemas
references/            workflow, data contracts, screening criteria, validation policy
examples/              synthetic fixtures and an example public target list
tests/                 offline unittest suite
agents/openai.yaml     agent harness declaration
SKILL.md               skill-level task definition and boundaries
```

## Provenance and privacy

Runs are content-addressed and checksummed (`verify-run` re-checks every artifact). Raw workbooks, private construct sequences and experimental results are private research data: they live outside the package, and the packaging allowlist plus `manage_skill.py` checks exclude them from any release. This repository contains code, schemas, documentation and synthetic/public examples only.

## Versioning

See [CHANGELOG.md](CHANGELOG.md). Current release: **0.3.0** (batch set-building, independent screening-recommendation layer, reconciled summary outputs; candidate engine and validity checks preserved from 0.2.x).

## License

MIT — see [LICENSE](LICENSE).
