# Validation and acceptance

## Offline acceptance

Run `make all-checks-offline`. The suite checks:

- Coordinate slicing roundtrip, range errors, fuzzy coordinates, unknown sequence residues.
- Full-ECD versus sourced-domain alternatives; known and unknown epitope coverage.
- Default full-form-only versus explicitly enabled truncations; mature secreted chains and multichain dependence; sourced shed release versus motif-only proposals.
- GPI omega inclusion, downstream-signal rejection, propeptide retention, uncertain mature boundaries and explicit scope exclusions.
- Six contextual criteria: absent/unresolved evidence never becomes low risk; no native-shedding-to-basal-activation inference, no gene-name culture-interference label, no pooling contradictory contexts.
- Domain truncation, disulfide partner removal, native signal/TM retention.
- Multi-pass loop non-concatenation, type-II attachment flag, GPI review and lumen exclusion.
- No scaffold/no fusion, module count/order, explicit signal policy, scoped risk acknowledgement.
- Review invalidation after scaffold, reference/rule/code change; synthetic fusion assembly branch.
- CDS translation/internal stops, actual sequence mismatches, construct identities.
- Percent/fraction semantics, MFI statistic, missing values versus zero, Y/N non-label handling.
- Sender/stimulus/context boundaries and potential duplicate-context conflicts.
- Source-cell mapping, XLSX formula/merge preservation, exact worksheet names.
- EvidenceRecord v2 lossless pass-through and invalid extra-field rejection.
- Immutable output, complete-bundle checksum validation, corruption preservation and bounded network retries with mocked responses.
- Batch acquisition: every list row preserved (blank/invalid included), gene-name mapping evidence and ambiguity, duplicate linking, per-item failure isolation, partial-set marking, offline cache reuse, query release/limit recording.
- Screening layer: completes with no scaffold/review/experiments; canonical analysis reference vs experimental isoform confirmation kept distinct; missing optional risk literature never forces undeterminable; core sequence errors still blocked; type-II/GPI/multi-pass/secreted/multichain handling; multi-candidate primary/alternate rules; duplicate/out-of-scope/extension counting; partial runs never masquerade as complete; lab-rule provenance and no-upgrade; fuzzy secondary annotation deferral versus sole-essential blocking; reviewed versus unreviewed prediction-boundary handling; legacy scaffold assembly and four-endpoint linkage not regressed.

The fixture is deliberately synthetic; it does not assess any real antigen or laboratory receptor. A main smoke run emits 3 candidates and zero fusions because no actual scaffold is supplied. Empty fusion/evidence-v2 files may be legitimate; verify schemas, counts and state, not a blanket nonempty-file rule.

## Public pilot

After synthetic tests pass, a small live pilot fetches roughly a dozen public human UniProt entries spanning type-I, type-II, GPI, multi-pass (GPCR), multichain and secreted-extension classes, screens them deterministically (cache and outputs outside the package), and is manually reviewed for identity mapping, class assignment and coordinate correctness. The pilot is labelled PILOT in its figure and PI_SUMMARY. It demonstrates the software contract on real annotation shapes; it is not evidence of biological accuracy, and it does not start a full-scale run.

## Installation checks

Install an allowlisted clean package. Run installed `manage_skill.py validate`, installed CLI smoke and installed unit tests; write outputs outside the installation. Verify compatibility symlink target and package archive membership.

## Real-data intake regression

When private workbooks are available, `inspect` can test source preservation and mapping without running new protein screening. A successful import proves file parsing, not biological labeling. Keep raw input unchanged and compare before/after hashes. Exact definitions of ambiguous experimental columns remain pending until supplied by the laboratory.

## Not covered by software tests

Actual host, SNIPR module sequence, junction choices, receptor surface expression, native epitope conformation, basal activity and induced response are not validated by synthetic tests. Predictor performance and full-proteome coverage are not benchmarked. Future retrospective and prospective evaluation must remain distinguishable.

## Privacy boundaries

`privacy-check` inspects tracked paths when Git is present and the package allowlist in all cases. No public repository is created or changed by these commands. Generated outputs, raw workbooks, snapshots and cache are excluded. Review all source/documentation diffs before any future public release; naming rules do not detect every possible disclosure of unpublished biology.
