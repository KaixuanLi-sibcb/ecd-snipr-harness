# v0.4.2 Reference and Recommendation Audit

This is a targeted correction, not a new design engine or a biological performance benchmark.

## Corrected contracts

1. Select the UniProt **Displayed** isoform ID instead of appending `-1`. If no unique ID is annotated, use `accession:canonical` as a reference marker. A noncanonical request without a matching sequence/coordinate map remains unresolved.
2. Retain foreign-isoform features and unmatched molecule-specific comments in `excluded_annotations`, but do not apply them to the canonical sequence. Exact principal-name mature-chain comments ending at the reference terminus are recognized separately from shed fragments. Broader processed-form/isoform-name reconciliation is still a review task.
3. Generic **Membrane** becomes `membrane_unspecified`, not organelle-only. Named cell-membrane sublocations are recognized. A documented secreted form with peripheral cell-surface annotation remains a separate secreted-extension route; integral surface targeting is not inferred.
4. Invalid core topology/exclusion features remain blocking even when another annotation of the same kind is valid. Secondary feature deferral says only that the annotation was not used, not that it has no effect.
5. Each candidate receives its own recommendation/reasons; the protein-level primary recommendation is stored separately. Blocked alternatives cannot inherit a positive class.
6. A selected public reference with `experimental_isoform_confirmation != performed` blocks fusion export even after a candidate review. This does not block screening or fragment export. Confirmation is an asserted, review-bound record, not a measured result or authentication mechanism.
7. Native heteromer text triggers a context-review warning, not a claim of necessary ECD partner dependence. The lexical filter withholds obvious negation but does not replace contextual literature review.
8. The package validator checks that manifest, Python package and pyproject versions agree.

## What the audit does not establish

Passing tests validates these software contracts. It does not establish surface expression, epitope preservation, background activation or induced response. A query-complete dataset is not an exhaustive human membrane-proteome/isoform census. Hash verification checks stored content, not the scientific correctness of its classification. The native subunit, glycosylation, cysteine and adjacency checks remain computational annotation summaries.

## Migration

Preserve old raw caches and immutable runs. Build a **new** normalized set from the raw cache using `build-set --list ... --offline --outdir NEW_SET`, then run `screen --set NEW_SET --outdir NEW_OUTPUT`. Screening an old normalized set does not apply parser corrections. Do not replace old summaries with new counts or reuse old assembly approvals. Obtain approval before a full rerun or publication.

## Sources

- [UniProt canonical and isoform policy](https://www.uniprot.org/help/canonical_and_isoforms).
- [UniProt generic membrane location](https://www.uniprot.org/locations/SL-0162).
- [UniProtKB user manual](https://web.expasy.org/docs/userman.html): sequence-specific annotation and feature locations.

Regression tests use synthetic entries only. Public cached pilot results are local validation artifacts, not packaged institutional data or a functional accuracy estimate.
