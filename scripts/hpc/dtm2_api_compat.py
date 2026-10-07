"""Use the author's CLI batch-unpacking logic without modifying vendor code.

Pinned upstream API passes a one-element batch into scalar postprocessing.
Keep model inference, ensemble voting and label conversion in vendor functions.
"""
def unpack_single_batch(topologies, membrane_probabilities, marginals, length):
    if len(topologies) != 1 or not isinstance(topologies[0], str):
        raise ValueError('Expected exactly one vendor topology string in batch')
    topology = topologies[0]
    if len(topology) != length:
        raise ValueError('Vendor topology length does not match full reference')
    if len(membrane_probabilities) != 1 or (marginals and len(marginals) != 1):
        raise ValueError('Unexpected vendor batch dimensions')
    return topology, membrane_probabilities[0], marginals[0] if marginals else None


class ReferencePredictor:
    def __init__(self, model_dir, device):
        from deeptmhmm2_predictor import DeepTMHMM2
        self.vendor = DeepTMHMM2(model_dir=model_dir, device=device)
        self.device = device

    def predict(self, sequences, simplify_io=True, progress_bar=False):
        if len(sequences) != 1 or simplify_io is not True:
            raise ValueError('One explicit full reference and simplified I/O required')
        import numpy as np
        import torch
        from deeptmhmm2_predictor.api import Prediction, PredictionCollection
        from deeptmhmm2_predictor.constants import DISPLAY_PROT_TYPE_MAPPING, MEMBRANE_TYPE_NAMES
        from deeptmhmm2_predictor.utils.esm_embed import get_embedding
        from deeptmhmm2_predictor.utils.postprocessing import (
            aggregate_memtype_predictions, aggregate_topology_predictions,
            infer_structural_type, suppress_implausible_predictions,
            translate_ab_topology, get_named_regions)
        accession, sequence = next(iter(sequences.items()))
        embedding = get_embedding(sequence, self.vendor._esm_model,
                                  self.vendor._esm_alphabet, self.device).unsqueeze(0)
        mask = torch.ones(1, len(sequence), device=self.device)
        all_topologies, all_probabilities, all_marginals = [], [], []
        # Match the official CLI precision and unpacking, with batch size one.
        with torch.amp.autocast(self.device):
            for index in range(len(self.vendor._models.topology_models)):
                t, p, m, _ = self.vendor._models.predict_single_model(
                    index, embedding, mask, compute_marginals=False)
                topology, probabilities, marginals = unpack_single_batch(t, p, m, len(sequence))
                all_topologies.append(topology)
                all_probabilities.append(probabilities)
                all_marginals.append(marginals)
        raw, marginals = aggregate_topology_predictions(all_topologies, all_marginals)
        indices, probabilities = aggregate_memtype_predictions(all_probabilities)
        raw = suppress_implausible_predictions(raw)
        type_code = infer_structural_type(raw)
        main = max(indices, key=lambda i: probabilities[i])
        membrane_names = [MEMBRANE_TYPE_NAMES[i] for i in indices]
        if type_code in {'I', 'S', 'T'}:
            membrane_names = []
            probabilities = np.zeros_like(probabilities)
        return PredictionCollection([Prediction(
            id=accession, sequence=sequence, type=DISPLAY_PROT_TYPE_MAPPING[type_code],
            membrane_types=membrane_names, membrane_types_probabilities=probabilities,
            topology_string=translate_ab_topology(raw, main, simplify_io=True),
            segments=get_named_regions(raw, main, simplify_io=True),
            _topology_ab=raw, _marginals=marginals, _main_mem_idx=main, _simplify_io=True)])
