#!/usr/bin/env python3
"""Report solved counts and paired iteration metrics from selector_sweep JSON.

Any instance with a time limit/error is censored from iteration comparisons.
PAR2 assigns unsolved runs twice their fixed iteration budget. SGM uses +1000.
"""
import argparse
import json
import math
from pathlib import Path


def gm(values):
    return math.exp(sum(math.log(max(v, 1)) for v in values) / len(values)) if values else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path, nargs='+')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rows = {}
    for path in args.results:
        for result in json.loads(path.read_text()):
            rows.setdefault(result['name'], {})[result['config']] = result
    assert all(set(r) == {'A', 'L', 'S'} for r in rows.values()), 'each instance must have A,L,S'
    solved = lambda r: r['status'] == 'OPTIMAL' and r.get('kkt_pass') is True
    censored = [n for n, r in rows.items() if any(v['status'] not in ('OPTIMAL', 'ITERATION_LIMIT') for v in r.values())]
    valid = [n for n in rows if n not in censored]
    common = [n for n in valid if all(solved(r) for r in rows[n].values())]
    summary = {'instances': len(rows), 'solved': {c: sum(solved(r[c]) for r in rows.values()) for c in ('A', 'L', 'S')},
               'censored': censored, 'common_solved': len(common), 'metrics': {}}
    for config in ('A', 'L', 'S'):
        costs = [rows[n][config]['iterations'] if solved(rows[n][config]) else 2 * rows[n][config]['iteration_budget'] for n in valid]
        summary['metrics'][config] = {'common_solved_iteration_GM': gm([rows[n][config]['iterations'] for n in common]),
                                      'PAR2_SGM': gm([v + 1000 for v in costs]) - 1000 if costs else None}
    # Ratios may be below one; do not floor them as iteration counts are floored.
    def ratio_gm(values):
        return math.exp(sum(math.log(v) for v in values) / len(values)) if values else None
    summary['selector_over_A_common_GM'] = ratio_gm([max(rows[n]['S']['iterations'], 1) / max(rows[n]['A']['iterations'], 1) for n in common])
    summary['selector_over_A_shifted_PAR2_GM'] = ratio_gm([(1000 + (rows[n]['S']['iterations'] if solved(rows[n]['S']) else 2 * rows[n]['S']['iteration_budget'])) /
                                                        (1000 + (rows[n]['A']['iterations'] if solved(rows[n]['A']) else 2 * rows[n]['A']['iteration_budget'])) for n in valid])
    summary['worst_regressions'] = sorted([{'name': n, 'A': rows[n]['A']['iterations'], 'S': rows[n]['S']['iterations'],
                                          'ratio': max(rows[n]['S']['iterations'], 1) / max(rows[n]['A']['iterations'], 1)} for n in common], key=lambda r: r['ratio'], reverse=True)[:10]
    encoded = json.dumps(summary, indent=2) + '\n'
    print(encoded, end='')
    if args.output:
        args.output.write_text(encoded)


if __name__ == '__main__':
    main()
