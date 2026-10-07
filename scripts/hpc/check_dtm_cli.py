#!/usr/bin/env python3
"""Require residue-wise equality between compatibility adapter and vendor CLI."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ecd_snipr.common import file_hash, now, write_json


def read_records(path):
    lines = [s.strip() for s in Path(path).read_text().splitlines() if s.strip()]
    if len(lines) % 3:
        raise ValueError('Incomplete vendor three-line records')
    records = {}
    for index in range(0, len(lines), 3):
        header, sequence, topology = lines[index:index+3]
        accession = header[1:].split()[0]
        if not header.startswith('>') or accession in records or len(sequence) != len(topology):
            raise ValueError('Vendor output identity/length conflict')
        records[accession] = (header, sequence, topology)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--adapted-root', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--device', choices=['cpu','cuda'], default='cpu')
    parser.add_argument('--timeout', type=int, default=1800)
    args = parser.parse_args()
    jobs = json.loads(args.jobs.read_text())['jobs']
    if not jobs or args.timeout<=0:
        parser.error('Nonempty pilot and positive timeout required')
    folder = args.outdir.resolve()
    folder.mkdir(parents=True, exist_ok=False)
    fasta = folder/'input.fasta'
    fasta.write_text(''.join('>'+j['accession']+'\n'+j['sequence']+'\n' for j in jobs))
    tool = json.loads(args.config.read_text())['tools']['DeepTMHMM2']
    command = [tool['executable'], str(fasta), str(folder/'vendor'), '--model-dir', tool['model_dir'],
               '--device', args.device, '--simplify-io', '--batch-size', '1']
    with (folder/'stdout.log').open('w') as out, (folder/'stderr.log').open('w') as err:
        subprocess.run(command, stdout=out, stderr=err, check=True, timeout=args.timeout)
    actual = read_records(folder/'vendor/predicted_topologies.3line')
    if set(actual) != {j['accession'] for j in jobs}:
        raise ValueError('CLI coverage mismatch')
    for accession in actual:
        adapted = read_records(args.adapted_root/accession/'predicted_topologies.3line')
        if adapted.get(accession) != actual[accession]:
            raise ValueError('Adapter differs from official CLI for '+accession)
    result = dict(checked_at=now(), total=len(jobs), exact_matches=len(actual),
                  input_sha256=file_hash(fasta), command=command,
                  output_sha256=file_hash(folder/'vendor/predicted_topologies.3line'),
                  validation='Exact header, reference sequence and residue topology equality; not biological accuracy')
    write_json(folder/'concordance.json', result)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
