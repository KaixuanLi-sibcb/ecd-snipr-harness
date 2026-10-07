#!/usr/bin/env python3
"""Local IUPred2 API worker; one import, isolated errors and original score files."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', required=True)
    parser.add_argument('--mode', choices=['long', 'short'], default='long')
    args = parser.parse_args()
    from iupred.iupred2a import iupred
    for job in json.loads(Path(args.jobs).read_text()):
        try:
            scores, _ = iupred(job['sequence'], mode=args.mode)
            if len(scores) != len(job['sequence']):
                raise ValueError('Incomplete per-residue output')
            path = Path(job['output'])
            if path.exists():
                raise FileExistsError(path)
            lines = ['# IUPred2 API output: '+job['accession'], '# pos aa iupred2']
            lines += [f'{i}\t{aa}\t{float(score):.17g}' for i, (aa, score) in
                      enumerate(zip(job['sequence'], scores), 1)]
            path.write_text('\n'.join(lines)+'\n')
            result = {'accession':job['accession'], 'status':'computed'}
        except Exception as exc:
            result = {'accession':job['accession'], 'status':'error',
                      'reason':type(exc).__name__+': '+str(exc)}
        print(json.dumps(result), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
