# Receiver-function evidence policy v0.7.0

## What this change solves

Two questions must remain separate: can the reference support a fragment proposal,
and what has been established for that fragment in a specific receiver and assay?
The old design decision table is preserved. The new layer makes the second question
a first-class output rather than a footnote. It does not solve functional prediction
by renaming failed examples, adding gene lists or calling all candidates high risk.

The unit is candidate fragment x complete receptor/backbone version x assay context.
Use Antibody-Sender -> Antigen-Receiver. Do not substitute the antigen-sender LC system.

## Evidence decision table

| Available evidence | Endpoint state | Allowed conclusion |
|---|---|---|
| No applicable endpoint measurement | not_measured | Function risk undetermined, not low |
| Antigen-only match, different backbone, invalid audit or unresolved measurement | related_or_unresolved_only | Related evidence retained, no outcome transfer |
| Eligible measurement for the exact complete receiver | measured_in_recorded_context | A measurement exists in that context, not automatic success |
| Different eligible values for one complete measurement-context key | conflicting_measurements | Keep all values and sources; no silent choice |

All four endpoints are independent: surface expression, recognition retention,
basal activation and induced response. A zero is a numeric observation, not missing;
it is also not universal safety. Unknown gates/units/stimuli remain unresolved.
The existing input normalizer still gates measurement eligibility. No functional
acceptance thresholds or binary labels are added in this release.

## Mechanism review, not retrospective explanation

Six axes organize existing evidence: folding/trafficking; epitope repertoire;
attachment geometry; processing/junctions; assembly interactions; expression/assay
context. Each names its input signals and next check. Positive sourced context
records stay scoped to native protein, receiver or assay; missing, unresolved,
negative-in-context and out-of-fragment records remain visible separately.

Length, cysteine and glycosylation heuristics are unchanged. They cannot become a
failure score. Native shedding is not proof of cleavage in the fusion. A native
complex is not proof of ECD partner dependence or pathological receptor clustering.
No signal means not established, never a validated negative.

The actual backbone, junction/regulatory-core sequence, host, expression regime,
reporter and stimulation controls matter. Where absent, the software lists these
questions rather than inventing their values. Quantitative matched-expression and
control comparisons are not yet computed. The receiver architecture is not inferred
from a protein's UniProt sequence.

## Literature basis and limits

- [Zhu et al., Cell 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9108009/): receptor-module composition, cell context and reporter response must be distinguished; expression alone does not establish signaling. The paper's regulatory ECD/hinge is not synonymous with the target-antigen fragment used here. These results justify separate evidence axes, not a predictor for arbitrary antigen receivers.
- [Engineering of an enhanced synthetic Notch receptor, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7069970/): ligand-independent activation was studied with receptor-expression dependence. This motivates recording expression and control context; it supplies neither a universal background threshold nor validation for the laboratory's receiver configuration.

Sources checked 2026-10-06. Neither paper is evidence that a particular untested
antigen fragment will self-activate in a different scaffold or host.

## Anti-leakage policy and tests

Private failure cases may reveal missing measurement dimensions, but must not be
used to relabel the historical design output or selected to claim accuracy.
There is no gene lookup, automatic failure-table loading, learned cutoff, risk-count
score, model training or prevalence estimate. Synthetic tests check name invariance,
missingness, changed scaffolds, fragment-only matches, context conflicts and all four
endpoints. These are software tests, not biological accuracy tests.

For a future predictor: define endpoints and thresholds independently, include
both failures and successes, freeze data/rules, group by protein/homologous family
and actual construct, separate batches, and evaluate untouched or prospective data.
Cases already used in development are not held-out test cases.

## Execution and immutable migration

New `screen` and `run` bundles automatically contain `receiver_function.tsv`,
`receiver_mechanism_review.tsv`, `receiver_function_summary.json` and
`RECEIVER_FUNCTION.md`. Existing TSV/JSON/FASTA deliverables are preserved.

```bash
python3 scripts/ecd_snipr_cli.py screen --set /path/to/set --outdir /path/to/new-output --resume
python3 scripts/ecd_snipr_cli.py receiver-audit --run-dir /path/to/old-verified-run \
  --outdir /path/to/private/receiver-overlay --resume
python3 scripts/ecd_snipr_cli.py verify-run --run-dir /path/to/overlay-bundle
make all-checks-offline
```

The overlay does not rerun candidate generation, read a failure workbook or change
any original sequence, primary selection or classification. It binds the original
manifest hash and new engine hash. It refuses output inside the original bundle,
checks source integrity before and after reading, and preserves corrupt output for
inspection instead of overwriting. Choose a new output root to retry a corrupt overlay.

## Deployment and agent handoff

```bash
make install-user
python3 "$HOME/.agents/skills/ecd-snipr-harness/scripts/manage_skill.py" validate
python3 "$HOME/.agents/skills/ecd-snipr-harness/scripts/ecd_snipr_cli.py" smoke \
  --outdir /path/to/private/installed-smoke
```

The installer validates an allowlisted staging copy and backs up an existing install
outside skill discovery. Do not copy the working tree, private outputs or workbooks.
Installation does not publish or push anything. Preserve old packages and run bundles.

Suggested task prompt:

> Use ecd-snipr-harness. Preserve the original run and private data. Present fragment
> design and four receiver-function evidence states separately. Audit applicable
> constructs and contexts, retain unknowns and conflicts, and list next measurements.
> Do not infer low background from standard_candidate or count generic warnings as
> failure predictions. Do not train or tune rules from my failure list.
