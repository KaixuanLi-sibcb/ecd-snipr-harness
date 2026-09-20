# v0.5.0 Public Reference-Set Execution

This is an aggregate execution record for public annotations, not a biological benchmark. No private source workbook, construct sequence or experimental label was used in this run. Detailed local outputs and caches are intentionally not committed.

## Scope and snapshot

- Configuration: Antibody-Sender -> Antigen-Receiver; antigen-fragment proposals only.
- UniProtKB release: `2026_03`, release date `2026-09-02`.
- Live query completed: `2026-09-20T05:08:24.390394+00:00`, without a limit or pilot flag.
- Query: `(organism_id:9606) AND (reviewed:true) AND (keyword:"Membrane" OR keyword:"Cell membrane")`.
- Retrieved: 7,783 unique accessions and 7,746 unique gene labels; all rows accounted for, no duplicate input, no technical failure.
- Entry annotation cache retrieval: `2026-09-18T09:00:20.179631+00:00` through `2026-09-18T09:02:35.702570+00:00`. Raw content hashes were checked and all entries were normalized again with v0.5.0. This was not a fresh download of every entry.
- Reference policy: Displayed canonical isoform ID, or an explicit canonical marker when no unique ID is available; not laboratory isoform confirmation.

Complete means complete **within that query and snapshot**. Unreviewed predictions, all noncanonical isoform sequences and entries outside the keyword query are not included. Membrane-set membership, native surface accessibility and current design scope remain distinct.

## Reconciled results

| Recommendation | Core references | Secreted extension | All evaluated |
|---|---:|---:|---:|
| Standard candidate | 566 | 0 | 566 |
| Conditional candidate | 942 | 301 | 1,243 |
| No supported route | 1,595 | 2 | 1,597 |
| Insufficient core evidence | 2,246 | 43 | 2,289 |
| **Total** | **5,349** | **346** | **5,695** |

The other 2,088 references are outside the declared design scope and remain recorded. There are 1,809 references with a primary candidate, 2,129 candidate records and 2,105 unblocked fragments exported as FASTA. The 24 blocked candidate records are retained but not exported as usable sequence. No complete receptor fusion was exported.

## What was checked

The local offline validation, smoke, privacy and package suite passed with 200 unit tests. Independent full-run acceptance passed 23 checks covering query/row accounting, normalized-reference hashes, candidate ID uniqueness, coordinate-to-sequence equality, foreign-isoform exclusion, primary/alternate blocking consistency, exact FASTA membership, pending review labels and preservation of historical bundles. These checks test software/output consistency, not the scientific truth of every annotation or recommendation. GitHub CI runs the offline suite only, not this live query or any private-data task.

| Provenance key | Value |
|---|---|
| Engine SHA256 | `6198ea04b19203921336f859f4db7ac01ec960f02ce57a4cbaa38b6a3612fb17` |
| Local immutable run ID | `9413cdd70f88e61d78f1` |
| Run manifest SHA256 | `1d140b0db351d2334d08fb3c148bbde46c15ae446523dbaf4c9e1a26fa33d5d6` |
| Summary SHA256 | `3338fe6b8b795bcd6ee69481e094fa2f34948abd176294c4f1e9d5eacce01a3f` |

These hashes identify retained local artifacts; the table does not imply that those full artifacts are hosted in Git. Network timestamps and changed database contents may produce a different run ID on a later execution.

## Remaining limitations

- All 1,809 primary candidates have unknown mapped-epitope status in this run. An absence of known epitope loss does not establish preserved recognition.
- 31 primary candidates have annotated domain/disulfide boundary concerns. Existing sourced fragments still require structural and epitope-scope review.
- The 93-row stratified review queue remains pending. It is purposive sampling, not a probability sample for overall accuracy estimation.
- Surface expression, recognition retention, basal activation and induced response are unmeasured. Actual receiver backbone and junction configuration were not supplied; no fusion sequence was invented.
- `SLC24A1` (`O60721`) has one explanatory sentence in `reason_codes`, copied through a shared list with `missing_info`. This is a nonblocking field-format issue; its `insufficient_evidence` recommendation and absence of a usable candidate are unchanged. Consumers should not assume all v0.5.0 reason-code values are machine tokens. The original run was not silently edited.
- Nonstandard residues can be outside the current sequence validator's supported alphabet; a software block does not establish an erroneous natural protein sequence.
- Annotation-driven scope and route decisions are still reviewable. A zero technical-failure count does not mean zero scientific limitations.

## Reproduction entry point

```bash
make all-checks-offline
python3 scripts/ecd_snipr_cli.py build-set --query default \
  --cache cache --outdir sets/human_membrane_v050
python3 scripts/ecd_snipr_cli.py screen --set sets/human_membrane_v050 \
  --outdir runs/human_membrane_v050 --review-per-stratum 2 --resume
python3 scripts/ecd_snipr_cli.py verify-run --run-dir runs/human_membrane_v050/RUN_ID
```

Use a new output directory and the actual returned `RUN_ID`. Exact snapshot replay requires the matching retained public cache; a current live query may differ. Query and annotation service: [UniProt REST](https://rest.uniprot.org/); reference-selection policy: [UniProt canonical and isoforms](https://www.uniprot.org/help/canonical_and_isoforms).
