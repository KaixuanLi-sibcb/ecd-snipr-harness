# Four endpoint contract

| Endpoint | Question | What it does not prove |
|---|---|---|
| `surface_expression` | Is this actual receptor detectable on the cell surface? | Correct antigen conformation, binding or signaling |
| `recognition_retention` | Does the intended antibody recognize this fragment in this receptor context? | Low background or induced transcriptional response |
| `basal_activation` | What is the reporter output under a defined no-binding condition? | Stimulated functionality |
| `induced_response` | What is output under a defined binding-sender condition? | Specificity or retained epitope breadth without appropriate controls |

## Raw fields

- An explicitly defined self-activation field may be assigned to `basal_activation`; missing stimulus/control details still prevent it becoming an eligible quantitative observation.
- BFP+/myc+/mRuby3+ is not a self-explanatory denominator. Record which events are selected, which are reporter-positive, gate hierarchy, channel, thresholding method and stimulus. Preserve ambiguous text.
- `SNIPR Assay (Y/N)` is a raw project-status field until its definition is confirmed. Do not map Y to success or N to failure. A later antibody-screening result is not an ECD-receiver expression outcome.
- Blanks, `x`, NA and pending remain missing; measured numeric zero is retained as zero only with interpretable units/context.
- Formulas retain both formula text and cached value; this reader does not recalculate them or infer number formats. Percent/fraction must be explicitly declared, not guessed from numeric magnitude.

## Minimum eligible measurement

`observation_id`, known and sequence-audited `construct_id`, one of the four endpoint names, numeric finite `raw_value`, `source_ref`, `channel`, `gate_path`, declared `unit`, and complete context:

`batch_id`, `host_cell`, `sender_id`, `antibody_id`, `stimulus`, `condition_id`, `timepoint`, `replicate_id`.

Use explicit `none` or `not_applicable` for truly absent entities. Do not silently fill them. Additional context should include receiver/sender abundance, ratio, tag/reporter definitions, detection reagent, acquisition/settings and sample-processing details wherever available. Add these to `context`; they remain part of the immutable record and conflict key.

Fractions/percentages require `denominator` and are normalized to [0,1] only within valid declared bounds. `fluorescence_au` requires `statistic=mean/median/geometric_mean`; the abbreviation MFI alone is insufficient. Unknown units produce `unresolved`.

Basal requires `stimulus=none/nonbinding_control`. Induced and recognition measurements require `stimulus=binding_sender` in this specific platform; soluble antibody assays should be preserved as raw supporting records until an explicit alternative assay contract is added, not silently treated as equivalent sender stimulation.

`eligible_measurement_not_success_label` means the software contract is met, not a validated biological label. The engine computes neither fold-induction across unmatched samples nor binary success classes. Comparisons must explicitly match relevant conditions, report controls and uncertainty, and account for replicate/batch structure.

## Conflict behavior

Different antibodies, senders, batches, replicates or conditions are different records, not automatically contradictory. Multiple different nonblank values for the exact same construct/endpoint/context/unit/gate/statistic are flagged as a potential duplicate-context conflict. They remain visible; a human determines whether this is a true data error or incompletely described repeat. No value is overwritten.

## Modeling gate

First count distinct full-sequence constructs, constructs with eligible endpoint measurements, class distributions if outcome definitions are prospectively agreed, and batch confounding. Raw row count and unique fragment strings are not independent training-example counts.

Begin with interpretable baselines and explicit uncertainty. Keep all constructs from the same protein together; stronger validation holds out homologous families. Sender, backbone and batch dependence must be assessed, not hidden in random row splitting. Retrospective fit is not prospective validation. This release does not train models, infer success labels, or provide a minimum-number shortcut for model readiness.
