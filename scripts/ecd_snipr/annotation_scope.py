"""Conservative statement attribution, not a general natural-language parser."""
import re


def association_subject(clause):
    """Only promote explicit/implicit reference-subject templates.

    Interacting with a partner homodimer is not being a homodimer. Unknown
    grammar remains unresolved rather than guessing the biological subject.
    """
    text = clause.strip()
    implicit = r'^(?:homo(?:dimer|trimer|tetramer|[- ]?oligomer)\w*|self[- ]ligand\b)'
    predicate = r'^(?:(?:it|the (?:protein|receptor|ectodomain|extracellular domain))\s+)?(?:(?:forms?|is|exists? as)\s+(?:(?:a|an|the)\s+)?homo(?:dimer|trimer|tetramer|[- ]?oligomer)\w*|self[- ]associates?\b|mediates? homophilic\b|binds? (?:to )?itself\b)'
    if re.search(implicit, text, re.I) or re.search(predicate, text, re.I):
        return 'reference_subject'
    if re.search(r'\b(?:interacts? with|associates? with|binds? to|complex (?:with|of)|consists? of|composed of)\b', text, re.I):
        return 'partner_or_complex_subject'
    return 'subject_unresolved'


def ligand_subject(clause):
    if re.search(r'^(?:(?:it|the (?:protein|receptor))\s+)?(?:binds?|recogniz\w*|receptor for|acts? as (?:a )?receptor for)\b', clause.strip(), re.I):
        return 'reference_subject'
    return 'subject_unresolved'


def statements(comment):
    """Preserve per-text citations when normalized UniProt supplies them."""
    return comment.get('statements') or [{'text': comment.get('text', ''),
                                         'evidences': comment.get('evidences', []),
                                         'eco': comment.get('eco', [])}]
