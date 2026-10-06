# Antigen-unit scope and research-goal acceptance

## What is and is not achieved

The software can enumerate a declared public reference set, retain unresolved/out-of-scope records, and produce sourced sequence fragments plus design-review questions. This meets the computational triage goal. It does not yet establish a validated receptor-construct library, complete epitope representation, or reliable SNIPR expression/activation predictions.

Sequence coordinates, reference version, processing form, annotation evidence, actual receiver and assay context are distinct layers. No numerical success probability can fill missing experimental evidence.

## Biological units

A gene/reference may produce more than one mature chain. The `full_ecd` label refers to a complete annotated topological interval, not every mature extracellular product or all potential antibody epitopes. Product-specific localization must not be inherited by the whole precursor. UniProt [Chain documentation](https://www.uniprot.org/help/chain) and [feature manual](https://web.expasy.org/docs/userman.html) distinguish processed chains, topology and repeated sequence regions.

v0.5.1 links only uniquely named, exact, same-reference Chain annotations to exact matching SUBCELLULAR LOCATION molecule names. No fuzzy/foreign-isoform/ambiguous-name mapping is used. Other source annotations remain preserved rather than silently transferred. An exactly localized secreted chain is a product for separate review, not an automatically authorized fusion.

## Transparent checks

| Check | Evidence required | Consequence | Not inferred |
|---|---|---|---|
| processed_chain_segments_spanned | Two disjoint, sourced chain intervals both overlap a candidate | Conditional review of processing and co-occurrence | Fusion cleavage, simultaneous mature products or obligatory partner dependence |
| extracellular_processed_product_not_covered | Exact product-scoped secreted localization, valid chain boundaries, no annotated TM/intramembrane/SP/cytoplasmic/propeptide/GPI-signal overlap, not fully inside fragment | Conditional review of the intended antigen unit | Every gene epitope is lost or a replacement sequence is safe |
| repeat_cut_by_boundary | Candidate partly overlaps an exact sourced Repeat feature | Conditional boundary review; included in annotated-integrity selection axis | Autonomous folding unit or inevitable misfolding |
| segmented_single_pass_needs_product_specific_design | Multiple annotated extracellular intervals lie on one side of a single TM, but current routes give no fragment | Explicit route limitation, not unknown orientation | Permission to join intervals or reconstruct processing |
| intramembrane_segmented_route_requires_review | One TM plus separately annotated membrane-embedded sequence interrupts extracellular intervals, with no supported candidate | Explicit embedded-segment route gap | A second full TM span or cleavage between mature products |
| retained_intramembrane | Soluble fragment intersects an exactly annotated membrane-embedded segment | Block that fragment; fuzzy exclusion coordinates remain essential unresolved evidence | All other external domains are unusable |

Nested precursor/chain annotations alone do not trigger the disjoint-product rule. Repeats may be short sequence repetitions or domain-like units; no repeat alone generates an alternate candidate. Invalid repeat coordinates remain deferred informational annotations, not a substitute for exact structural evidence.

## Repertoire and experiments

`epitope_evidence_scope=no_mapped_annotations_supplied` is an input-coverage statement. `external_epitope_database_search=not_performed` prevents a negative-search claim. Even mapped sequence retention is not binding preservation. Screening with an unknown antibody repertoire generally favors retaining mature antigen diversity, but a mechanistic tradeoff cannot be settled by a generic length threshold.

The SNIPR study separates LBD from the receptor's regulatory ECD, TMD and JMD; module effects on expression and signaling do not validate arbitrary antigen domains at the receiving recognition position. [Zhu et al., Cell 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9108009/). Its design rules cannot be imported as experimentally calibrated antigen-fragment thresholds.

## Implementation and outputs

- `antigen_context.py`: read-only interpretation of source annotations; no sequence generation.
- `uniprot.py`: product-scoped source links, repeat and intramembrane import, multiple-external-interval orientation. Intramembrane evidence is retained separately, never counted as an additional full membrane span.
- `design.py`: add named review flags before human assembly approval, so new issues cannot bypass review.
- `screening.py`: conditional reasons remain independent of assembly/functional status.
- `antigen_context.tsv`: one row per candidate, detailed JSON also in candidates.json.
- `processed_product_review.tsv`: reference/product-scoped records, including cases where no candidate exists; no new validated-design count.
- Summary counts processing/repertoire and repeat flags on primary candidates with an explicit primary denominator. Review strata preserve these cases separately.

The run manifest checksums all new outputs. Schemas remain backward-compatible additions; rebuild normalized sets from preserved raw sources when applying parser changes. Do not overwrite old bundles or reuse old review approvals.

## Next evidence needed

1. Independently review source boundaries and biological-unit selection across ordinary, conditional and unsupported routes; annotation-derived labels are not independent ground truth.
2. Supply the actual versioned receiver and connection configuration. Keep candidate-only delivery available before this.
3. Associate expression, recognition, background and induction measurements with construct, sender/antibody, batch, gating, denominator and controls. No blank/Y-N label inference.
4. Compare policies on held-out proteins/families and then prospectively. Rule-consistency tests, independent annotation review and prospective biological validation are different forms of evidence.

No private laboratory evidence is included in this specification. No model or adaptive rule is claimed to have been calibrated using institutional results.
