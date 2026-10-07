# Core evidence cross-checks - v0.9.0

## Purpose

Preserve sourced candidate generation, coordinates, primary selection and four
design classes. Strengthen the separate question: does this *exact fragment*
retain an annotated structural unit, and does the evidence actually concern it?
No failure workbook, gene blacklist, fitted cutoff or success probability is used.

## Default checks

`screen` and `run` automatically export exact reference-fragment correspondence,
UniProt domain/repeat retention and native SP/TM/intramembrane/cytoplasmic overlap.
Coordinate and assembly blocking checks remain in the original engine. These new
records never override an old approval or silently rewrite a supported candidate.

Statement attribution is conservative: `Homodimer` / `Forms a homodimer` and
explicit reference-subject homophilic descriptions are native context; `interacts
with a partner homodimer` is excluded, and unrecognized subjects stay unresolved.
This is a restricted grammar, not semantic understanding of arbitrary prose.
Negation, intracellular scope and fragment applicability still apply. New UniProt
normalizations retain text-specific ECO/PubMed records. Historical normalized
records lacking those citations are not retroactively enriched or called checked.

## InterPro/Pfam integration

The [official InterPro API](https://interpro-documentation.readthedocs.io/en/latest/api.html)
is queried by **public human canonical accession only**. No sequence is submitted.
The workflow downloads the reference sequence and paginated matched entries,
records the `InterPro-Version` header, retrieval times and content hashes, and
requires exact equality with the chosen UniProt analysis sequence. A release or
sequence mismatch is a visible conflict, not an automatic coordinate transfer.

This retrieves existing domain predictions; it does **not** execute InterProScan
or a new Pfam HMM search. No heavy dependency or GPU is required. Database family
classification is not treated as an autonomous domain boundary. Domain/repeat/
superfamily hits remain distinct; repeats/superfamilies are not autonomous folds.
Discontinuous match fragments stay separate, never stitched. Member and integrated
entries can describe the same unit and are not independent votes. A boundary cut
prompts review, not an assertion that the receiver will fail. Retention does not
prove correct folding, a preserved conformational epitope or good geometry.

```bash
python3 scripts/ecd_snipr_cli.py core-evidence --run-dir /local/verified/run \
  --reference-set /local/analysis-set --outdir /local/new/core-evidence \
  --interpro-cache /local/interpro-cache --interpro-live
# Omit --interpro-live to use checksum-verified local cache only.
python3 scripts/ecd_snipr_cli.py verify-run --run-dir /local/new/core-evidence/BUNDLE
python3 scripts/verify_core_evidence.py --run-dir /local/new/core-evidence/BUNDLE \
  --reference-set /local/analysis-set
```

Network failure => `error`; unavailable offline cache => `missing`; not requested
=> `not_run`. No hits after complete retrieval are source coverage, not a negative
SNIPR result. Fetch/pagination errors or sequence conflicts mark enrichment partial.
Optional missing tools do not block candidate delivery.

## Local external-tool outputs

The [licensed-tool guidance](https://interproscan6.readthedocs.io/stable/licensed_applications/)
distinguishes the workflow license from third-party models. No proprietary predictor
weights are installed, redistributed or auto-downloaded. Never submit private
fusion sequences to a web service by default.

| Tool | Supported local import | Interpretation |
|---|---|---|
| SignalP | Positive single-reference GFF3 regions + exact original input FASTA | Native precursor signal region; not a secretion verdict for an isolated ECD or unknown fusion |
| DeepTMHMM | Single-reference 3-line sequence/label output + original FASTA | Predicted SP/TM/inside/outside; outside can be organelle lumen, not natural cell-surface proof |
| IUPred2A/IUPred3/AIUPred | Ordered residue/score table + original FASTA | Per-residue disorder prediction; no SNIPR risk threshold |

`import-prediction` reads actual existing outputs, not executes a model. The FASTA
accession/sequence, output identity/residue order, ranges, version, parameters and
both raw-file hashes are checked. Import declares `full_reference` coordinates.
Software version and invocation parameters are caller-supplied provenance assertions;
the importer does not verify a vendor executable or manufacture proof it was run.
Exact slicing projects results into fragment coordinates. Native SP/TM outside a
fragment cannot become a retained-fragment conflict. Disorder scores are sliced,
not re-run on the truncated construct; truncation-dependent changes remain unknown.
Existing hash-bound `candidate_fragment` normalized imports remain supported.
Actual full-fusion prediction is **not implemented** and cannot be faked as a
reference run. SignalP empty GFF alone is not a checked no-SP result; unsupported
formats, DeepTMHMM2 labels and multiple sequences are rejected, not guessed.

```bash
python3 scripts/ecd_snipr_cli.py import-prediction --software DeepTMHMM \
  --input /local/actual/prediction.3line --input-fasta /local/actual/input.fasta \
  --reference /local/analysis-set/proteins/ACCESSION.json --version ACTUAL_VERSION \
  --parameters '{"mode":"actual invocation setting"}' --output /local/new/predictions.json
python3 scripts/ecd_snipr_cli.py core-evidence --run-dir /local/verified/run \
  --reference-set /local/analysis-set --tool-evidence /local/new/predictions.json \
  --outdir /local/new/with-tools
```

## Artifacts and traceability

- `candidate_core_evidence.tsv`: sequence/interval checks, source coverage, conflicts
  and rejected evidence, one row per candidate including blocked records.
- `candidate_domain_coverage.tsv`: exact retained/cut/omitted annotation instances
  with original discontinuous intervals, database/type/version/hash.
- `core_evidence_summary.json`: count denominators and method coverage states.
- `core_source_records.json`: optional InterPro snapshot hashes/statuses; no claim
  that a live lookup succeeded when it did not.
- Existing `candidates.json` embeds `receiver_risk.core_evidence`; all new files are
  bound by the immutable run manifest and excluded from packages/Git.

No legacy EvidenceRecord v2 is modified. New calculations keep native computed/
prediction context, following the [existing contract](data-contracts.md).

## Acceptance, deployment and next boundary

```bash
make all-checks-offline
python3 scripts/run_core_pilot.py --outdir /local/new/public-pilot --live
make install-user
```

The public pilot uses four public accessions spanning type I, type II, GPI and
multi-pass. Unsupported multi-pass entries remain in the denominator. Outputs are
local and the pilot is not a full proteome run or a biological validation set.
Installed validation/smoke must be rerun after installation. The install helper
backs up the previous installation before replacement.

`verify_core_evidence.py` independently re-derives exact reference slicing and
domain retention, checks raw InterPro response hashes and rejects unsupported
functional claims. It can detect an incorrect derived label even if a manifest
has been resealed. It is a deterministic second pass, not an independent biological
benchmark or a separate LLM agent.

Structure/interface mapping, PISA, aggregation prediction and actual receiver
geometry remain future targeted analyses; no unchecked PDB residue numbering or
generic AlphaFold confidence becomes receiver compatibility. The most important
experimental gap remains matched actual constructs and four independently defined
endpoints. Developmental case audits are never held-out accuracy.
