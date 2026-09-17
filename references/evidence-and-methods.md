# Evidence and method selection

## Architectural evidence boundaries

SNIPR's ligand-binding domain and its regulatory extracellular module are distinct. Modular replacement in the original work does not establish that every human antigen ectodomain functions as a receiver recognition module. Treat expression, background and induced output separately. [Zhu et al., Cell 2022](https://www.bu.edu/khalillab/_papers/zhu-cell2022.pdf).

SNAP-synNotch uses a different adaptor architecture. Its expression, labeling and signaling measurements illustrate why surface presence is not sufficient, but are not a direct validation dataset for this lab's antigen-receiver design. The article is from **2023**, not 2020. [Ruffo et al., Nature Communications](https://pubmed.ncbi.nlm.nih.gov/37160880/).

LAG3 illustrates an important truncation risk: individual domain names do not guarantee autonomous folds, and removing regions can change antibody recognition scope. In the reported structure, D3 and D4 form a connected unit; engineered structural variants must not silently replace the native sequence. This supports mandatory review, not a universal LAG3/SNIPR design rule. [Ming et al., Nature Immunology 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC10191176/).

## Available versus planned tools

| Resource | Useful evidence | Current integration | Boundary |
|---|---|---|---|
| [UniProt](https://www.uniprot.org/help/canonical_and_isoforms) | Sequence, topology, domains, ECO and isoforms | Local JSON import; optional public accession fetch with cache/retry | Intended experimental isoform is not selected automatically; reviewed entry is not proof every feature is experimentally established |
| [SignalP 6](https://services.healthtech.dtu.dk/services/SignalP-6.0/) | Signal-peptide/cutting-site prediction | Not executed by this release | Predicted targeting is not measured receptor expression |
| [DeepTMHMM](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/) | Membrane topology prediction | Not executed by this release | Topology does not establish plasma-membrane rather than organelle localization |
| [InterProScan](https://github.com/ebi-pf-team/interproscan6) | Domain-family annotations | Explicit feature input only | Domains may depend on neighbors or oligomers; do not slice them blindly |
| [AlphaFold DB](https://alphafold.ebi.ac.uk/faq) | Structure/region-confidence evidence | Agent-assisted review only | Low confidence is not a deletion instruction; no validated receptor-background prediction |
| [IUPred3](https://iupred3.elte.hu/) | Sequence disorder evidence | Not executed by this release | Disorder may include functional linkers or epitopes |
| [ESM](https://github.com/facebookresearch/esm) | Protein-sequence representations | Future features only | No ECD-SNIPR success predictor is bundled or claimed |

Prediction adapters require pinned model/database versions, coordinate conventions, licensing review for the actual deployment, and local privacy-compatible execution. No third-party models or databases are redistributed in this package. Academic availability is not automatically a commercial redistribution license.

## Present heuristics

Length, cysteine content, potential N-X-S/T sequons and topology are interpretable review descriptors. They do not estimate surface expression, epitope conformation, background or induced response. The framework records this lack of experimental validation explicitly. It does not optimize mutations or replace lost epitopes with generated sequence.

## Evidence ladder

1. Reference annotation with source/version and applicable sequence.
2. Prediction with method/version and coordinate mapping.
3. Deterministic rule-derived candidate and risk ledger.
4. Actual construct and context-matched experimental measurements.
5. Prospective validation across representative and held-out targets.

Do not promote a claim along this ladder because a file exists, software tests pass, or a candidate looks plausible. Missing source data is missing evidence, not negative evidence. This is a selected methodological reference set, not a claim of an exhaustive updated literature search or a universal validated classifier.
