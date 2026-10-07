"""Accession-only InterPro domain retrieval with sequence-bound snapshots.

No sequence submission or InterProScan execution. Cached canonical annotations
are usable only after exact sequence comparison with the chosen reference.
"""
import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path
from urllib.parse import urlparse

from .common import file_hash, now, read_json, sequence, write_json
from .sequence_tools import sequence_hash
from .uniprot import RETRYABLE_HTTP, retry_delay

BASE = 'https://www.ebi.ac.uk/interpro/api/'


def fetch_domains(accession, cache, offline=True, refresh=False, opener=None, sleeper=time.sleep):
    if not re.fullmatch(r'[A-Z0-9]{6,10}', accession):
        return {'status': 'missing', 'reason': 'canonical_accession_required_no_isoform_transfer', 'accession': accession}
    root = Path(cache)
    root.mkdir(parents=True, exist_ok=True)
    meta = root / (accession + '.meta.json')
    if meta.exists() and not refresh:
        try:
            record = read_json(meta)
            path = root / record['raw_file']
            payload = path.read_bytes()
            if path.parent != root or hashlib.sha256(payload).hexdigest() != record['raw_response_sha256']:
                raise ValueError('cache_hash_mismatch')
            snapshot = json.loads(payload)
            for retrieval in snapshot.get('retrievals', []):
                name = retrieval.get('raw_response_file')
                if not name or Path(name).name != name or file_hash(root/name) != retrieval['response_sha256']:
                    raise ValueError('original_response_cache_hash_mismatch')
            return dict(record, status='cached', snapshot=snapshot, raw_snapshot_path=str(path.resolve()))
        except (OSError, ValueError, KeyError, TypeError):
            pass
    if offline:
        return {'status': 'missing', 'accession': accession, 'reason': 'offline_cache_absent_or_invalid'}
    opener = opener or urllib.request.urlopen
    pages, versions, retrievals = [], set(), []

    def get(url):
        parsed = urlparse(url)
        if parsed.scheme != 'https' or parsed.netloc != 'www.ebi.ac.uk' or not parsed.path.startswith('/interpro/api/'):
            raise ValueError('unexpected_interpro_pagination_url')
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'ecd-snipr-harness'})
                with opener(request, timeout=30) as response:
                    payload = response.read()
                    headers = response.headers
                    version = headers.get('InterPro-Version', '')
                    http_status = getattr(response, 'status', 200)
                payload_hash = hashlib.sha256(payload).hexdigest()
                raw_path = root / ('response.' + payload_hash + '.json')
                if raw_path.exists() and file_hash(raw_path) != payload_hash:
                    raw_path = root / ('response.' + payload_hash + '.' + str(time.time_ns()) + '.json')
                if not raw_path.exists():
                    raw_path.write_bytes(payload)
                if not version:
                    raise ValueError('interpro_version_header_missing')
                versions.add(version)
                if len(versions) != 1:
                    raise ValueError('interpro_release_changed_during_fetch')
                retrievals.append({'url': url, 'response_sha256': payload_hash,
                                   'raw_response_file': raw_path.name, 'retrieved_at': now(), 'version': version})
                if http_status == 204:
                    return {'count':0, 'results':[], 'next':None}
                return json.loads(payload)
            except (OSError, ValueError) as exc:
                if isinstance(exc, urllib.error.HTTPError) and exc.code == 204:
                    return {'count': 0, 'results': [], 'next': None}
                if isinstance(exc, ValueError) or (isinstance(exc, urllib.error.HTTPError) and exc.code not in RETRYABLE_HTTP) or attempt == 2:
                    raise
                sleeper(retry_delay(exc, attempt))

    try:
        protein = get(BASE + 'protein/uniprot/' + accession + '/?extra_fields=sequence')
        metadata = protein['metadata']
        if metadata.get('accession', '').upper() != accession or str(metadata.get('source_organism', {}).get('taxId')) != '9606':
            raise ValueError('interpro_accession_or_species_mismatch')
        ref = sequence(metadata['sequence'])
        if metadata.get('length') != len(ref):
            raise ValueError('interpro_length_mismatch')
        url = BASE + 'entry/all/protein/uniprot/' + accession + '/?page_size=200'
        seen = set()
        total = None
        while url:
            if url in seen or len(pages) >= 100:
                raise ValueError('interpro_pagination_cycle_or_limit')
            seen.add(url)
            page = get(url)
            if not isinstance(page.get('results'), list) or type(page.get('count')) is not int:
                raise ValueError('malformed_interpro_page')
            total = page['count'] if total is None else total
            if page['count'] != total:
                raise ValueError('interpro_count_changed')
            pages.append(page)
            url = page.get('next')
        if sum(len(p['results']) for p in pages) != total:
            raise ValueError('incomplete_interpro_pages')
        snapshot = {'accession': accession, 'sequence': ref, 'sequence_sha256': sequence_hash(ref),
                    'version': next(iter(versions)), 'protein': protein, 'pages': pages,
                    'retrievals': retrievals, 'completeness': 'complete', 'retrieved_at': now()}
        payload = (json.dumps(snapshot, sort_keys=True, ensure_ascii=False) + '\n').encode()
        checksum = hashlib.sha256(payload).hexdigest()
        path = root / (accession + '.' + checksum + '.json')
        if path.exists() and file_hash(path) != checksum:
            path = root / (accession + '.' + checksum + '.' + str(time.time_ns()) + '.json')
        if not path.exists():
            path.write_bytes(payload)
        record = {'status': 'fetched', 'accession': accession, 'raw_file': path.name,
                  'raw_response_sha256': checksum, 'hash_scope':'snapshot_including_original_response_hashes',
                  'version': snapshot['version'], 'retrieved_at': snapshot['retrieved_at']}
        write_json(root / (accession + '.' + str(time.time_ns()) + '.retrieval.json'), record)
        write_json(meta, record)
        return dict(record, snapshot=snapshot, raw_snapshot_path=str(path.resolve()))
    except (OSError, ValueError, KeyError, TypeError, AttributeError, StopIteration) as exc:
        return {'status': 'error', 'accession': accession, 'reason': str(exc), 'retrievals': retrievals}


def normalized_domains(record, protein):
    snapshot = record['snapshot']
    ref = sequence(protein['sequence'])
    accession = protein.get('accession')
    if snapshot.get('accession') != accession or sequence(snapshot['sequence']) != ref or snapshot.get('sequence_sha256') != sequence_hash(ref):
        raise ValueError('interpro_reference_sequence_conflict')
    if protein.get('reference_coordinate_status') == 'unverified_noncanonical_request':
        raise ValueError('noncanonical_reference_coordinates_not_verified')
    if snapshot.get('completeness') != 'complete' or not snapshot.get('version'):
        raise ValueError('interpro_snapshot_incomplete_or_unversioned')
    result = []
    for page in snapshot['pages']:
        for entry in page['results']:
            meta = entry['metadata']
            if meta.get('type') not in {'domain', 'repeat', 'homologous_superfamily'}:
                continue  # Family assignment is not an autonomous domain boundary.
            for match in entry.get('proteins', []):
                if match.get('accession', '').upper() != accession or match.get('protein_length') != len(ref):
                    raise ValueError('interpro_match_identity_conflict')
                for location in match.get('entry_protein_locations', []):
                    fragments = location.get('fragments', [])
                    if not fragments:
                        raise ValueError('interpro_empty_domain_location')
                    for f in fragments:
                        if type(f.get('start')) is not int or type(f.get('end')) is not int or not 1 <= f['start'] <= f['end'] <= len(ref):
                            raise ValueError('interpro_domain_coordinate_conflict')
                    result.append({'kind': meta['type'], 'name': meta.get('name', ''),
                        'entry_accession': meta['accession'], 'database': meta['source_database'],
                        'integrated': meta.get('integrated'), 'fragments': fragments,
                        'start': min(f['start'] for f in fragments), 'end': max(f['end'] for f in fragments),
                        'representative': location.get('representative'), 'model': location.get('model'),
                        'evidence': {'kind': 'prediction', 'source': BASE + 'entry/' + meta['source_database'] + '/' + meta['accession'] + '/',
                                     'version': snapshot['version'], 'raw_response_sha256': record['raw_response_sha256'],
                                     'reference_sequence_sha256': sequence_hash(ref)}})
    return result
