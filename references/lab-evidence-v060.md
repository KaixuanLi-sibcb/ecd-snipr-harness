# Private laboratory evidence workflow

## Purpose

Public screening proposes sourced antigen fragments. This optional branch recovers
what the laboratory has recorded, compares fragments to supplied references, and
identifies what prevents an observation from being interpreted. Missing imported
fields are not missing experiments. Candidate generation rules remain unchanged.

This is a cell-level inventory, not a whole-XLSX reserialization. Drawings, comments,
external links and formatting semantics are not interpreted. Preserve the original
workbook and checksum. Formulas and caches are retained, not recalculated. Shared
formulas preserve their index, attributes and anchor(s); a follower without formula
text is not a literal value. Exact cached sequence matches are tentative, not proof
that a formula has been recalculated against current inputs.

## Flow

1. Hash and read without modifying the source. Preserve exact sheet/row/cell,
   original values, formulas, types, styles and merge ranges. Never fill merges down.
2. Detect repeated header sections with a versioned explicit dictionary. Normalize
   Unicode header syntax, recognize SNIPR and historical SNPIR spellings, preserve
   every duplicate column. Unknown columns remain available in raw outputs.
3. Quarantine displaced/composite identity layouts. Header detection is provisional,
   not universal table understanding. Explicit source-hash-bound region maps can
   resolve them, with a rationale. No positional guess for a headerless table.
4. Validate amino-acid strings and available paired CDS translations. Untyped ECD
   strings can be marked amino-acid-like by alphabet; DNA-compatible ones remain
   ambiguous. Never strip tags or edit a sequence to force a match.
5. Compare against exact, hash-verified references supplied through `--reference-set`.
   Retain every exact occurrence and ambiguity. Record native SP/TM/tail overlaps,
   cut domains and omitted annotated extracellular domains. A mismatch can reflect
   alternative isoforms, tags or engineering; it does not establish laboratory error.
6. Propose cross-row links using exact named identity fields. Y/N/checkmarks in a
   plasmid column are status markers, not shared plasmid identifiers; composite
   labels/notes are not split into guessed identities. Same gene, AgID or
   fragment never proves the same full receiver. Preserve duplicates and variants.
7. Keep raw self-activation and surface-expression headings as reported endpoints.
   BFP/myc/mRuby remains endpoint-unconfirmed. Y/N project status and downstream
   serology cannot become receiver-success labels. Layout-quarantined cells are
   excluded from interpretable-readout counts, but retained for review.
8. With human-confirmed cell-specific semantics and an audited actual construct,
   apply the existing four-endpoint quantitative contract. No automatic pooling,
   success label, ratio between unmatched conditions or functional probability.
9. Export an immutable, hash-verified private bundle and a finite confirmation queue.
   The pipeline makes no external requests and changes no public rankings.

## CLI and agent deployment

```bash
make all-checks-offline
python3 scripts/ecd_snipr_cli.py lab-evidence \
  --input /absolute/private/source.xlsx \
  --outdir /absolute/private/lab-evidence-output \
  --reference-set /absolute/local/public-analysis-set

# Only after responsible staff confirm definitions and actual construct links:
python3 scripts/ecd_snipr_cli.py lab-evidence \
  --input /absolute/private/source.xlsx \
  --outdir /absolute/private/lab-evidence-output \
  --reference-set /absolute/local/public-analysis-set \
  --semantics /absolute/private/confirmed-definitions.json \
  --construct-run /absolute/private/verified-construct-run

python3 scripts/ecd_snipr_cli.py verify-run --run-dir /absolute/private/returned-run-id
```

Use the returned directory, not an assumed run ID. `--resume` verifies source/code
identities and checksums. Corrupt outputs are never replaced; choose a new root for
a failed attempt. Corrupted references are isolated and mark execution `partial`
(exit 3); malformed source/contract fails the invocation. Interpretation can remain
pending when technical intake is complete. Missing scaffold does not block intake.

Install only allowlisted assets with the existing backup-preserving installer.
Run installed validate, smoke, lab-evidence fixture and tests; outputs must be
outside the installation. Never upload these local results or private definitions.

## Input contract

Optional `--layout`: `source_sha256`, `regions`; each region has exact `sheet`,
integer `start_row`/`end_row`, `columns` (Excel letters to dictionary field names)
and `rationale`. Regions cannot overlap and replace the entire map for their rows.
Unspecified columns remain unmapped. Layout confirmation is not assay confirmation.

`semantic_review_template.json` is generated pending. Review a copy outside the
immutable run. Top level: `schema_version=1.0`, matching `source_sha256`, `reviews`.
Confirmed items require `cell_id`, `raw_record_digest`, `reviewer`, `reviewed_at`,
`rationale`, and `definition`: endpoint, unit, denominator, statistic, channel,
gate_path, construct_id and full context from [assay interpretation](assay-interpretation.md).
Measurement formula caches require `formula_cache_reviewed=true`; formula-derived
sequence/CDS cells additionally require `sequence_formula_cache_reviewed=true`.
Wildcards,
unknown selectors and raw-value overwrites are prohibited. Review is an auditable
assertion, not an authenticated signature or proof of reviewer authority.

`--construct-run` must be an intact `run` bundle with `construct_audit.json`. No
missing backbone, linker, approval or sequence is invented. Fragment consistency
does not prove actual receptor function or calibrated measurement validity.

## Artifacts and count units

| Output | Meaning |
|---|---|
| raw_inventory.json / source_cells.tsv | Original cells, duplicate headers, formulas and provenance |
| row_records.json / header_sections.json | Provisional mapping and layout issues |
| fragment_reference_audit.tsv / fragments.json | Sequence, translation and reference-coordinate audit |
| construct_link_candidates.tsv | Proposed links only; no implicit experiment transfer |
| row_readiness.tsv | Per-source-row linkage/layout work remaining |
| assay_evidence.tsv / observations.json | Raw readouts/status and semantic eligibility |
| assay_definition_groups.tsv | Field discussion aid, not pooled experiment groups |
| semantic_review_template.json | Pending cell-specific definitions and construct selection |
| evidence_records_v2.jsonl / source_claims.jsonl | Source presence, not biological validity |
| conflicts.json / reference_errors.json | Non-overwritten conflicts and technical errors |
| summary.json / LAB_EVIDENCE_REPORT.md | Reconciled results and confirmation queue |

Count source cells, rows, fragment strings, audited full constructs and interpretable
measurements separately. Cross-sheet duplicate records are not independent
replicates. Missingness reasons overlap. None of these counts is a model sample-size
claim. No experimental success rate is exported.

EvidenceRecord v2 `confirmed_local_xlsx` confirms only a raw XLSX cell. Its normalized
value stays empty and confidence explicitly means source presence only. Source
SHA and cell query preserve provenance. CSV/TSV use native `source_claims.jsonl`
without falsely claiming XLSX origin. Computed fragment matches stay in native
outputs, not v2 confirmed biological claims. Digests of strings use canonical JSON;
file checksums use raw bytes, consistent with the existing contract.

## Limits and next steps

No FCS parsing, gate inference, unit inference, threshold calibration, new isoform
search, antibody-function transfer, model training or epitope validation is included.
First confirm a complete actual receiver and the dictionary/conditions for a small
set of observations. Then count independent constructs, endpoint coverage and batch
confounding before considering retrospective modeling or prospective validation.
