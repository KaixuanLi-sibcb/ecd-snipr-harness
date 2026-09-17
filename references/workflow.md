# Workflow and operational boundaries

## Unit of analysis

**Specific antigen fragment + specific receiver scaffold + specific antibody/sender + assay condition.** A gene name is not this unit. Gene, isoform, fragment, full fusion, sender, batch and observation identities remain distinct.

The human membrane proteome is a reference universe, not an implicit denominator. Any set built from a public query records the database release, taxon 9606, reviewed/unreviewed policy, the exact query, plasma-membrane versus all-membrane scope, canonical/isoform policy and the deduplication key in `set_definition.json`. Query truncation (`--limit`) and fetch gaps mark the set partial; a partial set is never reported as a full run.

## Execution phases

0. **Set building** (`build-set` / the chained first step of `screen`): user target list or public UniProt query. Accessions are processed directly; gene names go through recorded `gene_exact + organism_id:9606` mapping (reviewed entries first, then all), ambiguity preserved, no guessing, no row dropped, duplicates linked. Every entry gets an explicit disposition: membrane-set membership, natural cell-surface target, design scope, secreted extension, technical failure. Fetch/mapping caches are content-addressed and resumable; per-entry UniProt entry versions are recorded.
1. **Analysis-reference selection**: screening explicitly selects the database canonical sequence (or an explicitly requested accession-isoform) as the analysis reference and records the rationale. This is distinct from laboratory isoform confirmation, which still gates final fusion assembly. Genuine identity ambiguity, version conflicts and non-canonical coordinate problems are preserved.
2. **Inventory**: use `inspect`; hash the source bytes and preserve exact sheet/row/cell references, including formulas and cached values. XLSX formatting is not interpreted as units; formulas are not recalculated. Merge ranges are recorded, never silently filled down.
3. **Semantic mapping**: use `map` with explicit columns. Resolve aliases and intended isoforms outside the sequence-design engine. An agent can research ambiguity but must retain alternatives and sources. Do not turn raw mapped rows directly into outcome labels.
4. **Reference annotation**: import local UniProt JSON or public accession fetch results. Coordinate evidence must apply to the exact sequence/isoform. Import prediction outputs only with tool/model/database version and source. Preserve conflicting predictions, do not overwrite curated records. In batch screening, fuzzy/unsourced secondary annotations are deferred with an explicit record; sole essential annotations still block, and the strict assembly path keeps every fuzzy feature blocking.
5. **Candidate generation**: prefer a continuous full antigen form, preserving input order. Full-domain alternatives require explicit opt-in (`auto` opts in for multi-pass only); processed secreted/GPI forms and documented shed forms follow the [criteria policy](screening-criteria.md). A user-sourced experimental fragment is a proposal, not proof of performance. There is no random sliding-window generation or mutation optimization.
6. **Screening recommendation**: one transparent decision table per in-scope reference (standard / conditional / no standard route / insufficient evidence), with reason codes, rationale, evidence, missing info and alternates. See the [criteria policy](screening-criteria.md). Optional sourced lab rules annotate or downgrade only.
7. **Risk review**: validate coordinates and sequence, detect forbidden native segments, list retained/cut/omitted domains and known epitopes, disulfide boundary crossings, cysteine and potential glycosylation risks. Separately record processing, shedding, junction cleavage, oligomerization, aggregation and culture-interference sources/contexts, including missingness. Length thresholds prompt review, not hard biological cutoffs. No imported human-mouse ranking weight is used here.
8. **Scaffold review** (optional downstream branch, `run`): require actual lab module sequences/order and metadata. Confirm that the antigen is in the receiving receptor's antigen slot. Native signal peptide/TM/tail/GPI-signal overlap blocks assembly. This implementation supports one explicitly supplied N-terminal scaffold signal, one antigen slot and one receiver TM; unsupported architectures stop rather than being improvised.
9. **Human decision** (optional downstream branch): pending review template -> named reviewer approves with rationale and risk acknowledgements -> rerun. No self-approval by the agent. Block-level errors cannot be overridden by a review. Annotation changes create a new review key; source history remains in prior run bundles.
10. **Experiment linkage** (optional downstream branch): audit actual antigen/fusion sequences and optional CDS translations; map observations only to explicit constructs. Exact antigen sequence correspondence alone is not full receptor equivalence. Keep all four endpoints and contexts separate.
11. **Export and verify**: protein_screening.tsv, screening.json, candidate tables/FASTA, summary.json with explicit denominators, PI_SUMMARY.md, coverage figure with source data, state, manifest hashes. `verify-run` must pass. Count computationally supported candidates separately from expression-tested, recognition-tested, baseline-tested and response-tested constructs. Reports lead with research conclusions, not approval counts.

## Harness behavior

`run_key = SHA256(canonical project JSON + engine source checksum)`. Inputs, rules, scaffold, reviews, observations and imported evidence all contribute. Review keys additionally bind the candidate proposal and annotation snapshot; engine checksum is included through the effective rule set.

Outputs are immutable per run. `--resume` only accepts a complete bundle whose recorded files match their SHA256 values. A damaged/incomplete bundle is retained and a numbered new attempt is created. Exceptions create an error event/state and return a nonzero exit code. Missing biological inputs yield documented blocked/needs-review records, not a task crash or fabricated sequence.

Runtime is sequential and CPU-only. Candidate slicing and table export are lightweight; versioned API fetch is optional and separate. Three attempts with bounded delay cover transient errors; non-retryable HTTP errors stop. Cache reuse is explicit provenance, not a claim of a fresh live query. Use `--refresh` only when a newer public annotation is wanted. No continuous monitor or remote job scheduler is installed.

## Protein classes

| Class | Default design route | Caveat |
|---|---|---|
| Type I | Sourced mature N-terminal continuous ECD; optional sourced domains | Remove native signal/TM/tail through verified boundaries, not guessed shifts |
| Type II | Sourced C-terminal extracellular portion | Receiver fusion usually changes the side of attachment; orientation/geometry needs review |
| GPI-anchored | Sourced mature extracellular portion, including omega residue | GPI processing boundary and terminal attachment must be confirmed; downstream GPI signal must not remain |
| Secreted (scope extension) | Sourced mature chain, not SP-to-end guessing | Exclude annotated propeptides; multichain and tethered-presentation dependencies require review |
| Shed form | Explicit sourced release boundaries and release context | Native release does not prove suitability as a tethered receiver antigen |
| Multi-pass | Only a sourced extracellular domain contained within a single annotated external interval | No loop concatenation or promise of native conformational epitopes; alternate presentation may be necessary |
| Organelle membrane | Do not label lumen as external cell surface | Any alternative engineering route requires a separately justified scope |

## What may be autonomous

The agent may inspect, normalize, run deterministic validation, retrieve public annotations, prepare unresolved records and explain conflicts. It may not guess the lab's scaffold/host/reporter, fill unknown experimental conditions, silently approve candidates, upload private sequences, order synthesis, or declare a functional construct library solely from software tests.
