#!/usr/bin/env python3
"""Report actual policy runs. Censor time limits in CPU iteration comparisons."""
import argparse
import collections
import json
import math
from pathlib import Path


def solved(r):
    return r['status'] == 'OPTIMAL' and r['kkt_pass']


def report(rows):
    by = collections.defaultdict(dict)
    for r in rows:
        by[r['name']][r['config']] = r
    policies = sorted(set(r['config'] for r in rows))
    if any(set(by[n]) != set(policies) for n in by):
        raise ValueError('incomplete policy comparison')
    common = [n for n in by if all(solved(by[n][c]) for c in policies)]
    result = dict(count=len(by), backend=rows[0]['backend'], common_names=common,
                  policies={}, switches=[dict(name=r['name'], chosen=r['chosen'],
                      discarded_iterations=r['discarded_iterations']) for r in rows
                      if r['config'] == 'M' and r['discarded_iterations']])
    for c in policies:
        selected = [by[n][c] for n in by]
        entry = dict(solved=sum(solved(r) for r in selected),
            observed_seconds=sum(r['solver_wall_sec'] for r in selected),
            par2_seconds=sum(r['solver_wall_sec'] if solved(r) else 2*r['time_limit'] for r in selected),
            common_solved_seconds=sum(by[n][c]['solver_wall_sec'] for n in common),
            feature_seconds=sum(r['feature_wall_sec'] for r in selected))
        if common and rows[0]['backend'] == 'cpu':
            entry['common_iteration_sum'] = sum(by[n][c]['iterations'] for n in common)
            entry['common_iteration_gm'] = math.exp(sum(math.log(max(1, by[n][c]['iterations'])) for n in common) / len(common))
        if rows[0]['backend'] == 'cpu':
            for key in ('observed_seconds', 'par2_seconds', 'common_solved_seconds'):
                entry.pop(key)
        result['policies'][c] = entry
    if 'M' in policies and 'P' in policies:
        pair = [n for n in by if solved(by[n]['P']) and solved(by[n]['M'])]
        result['MP_common_names'] = pair
        result['worst_time_regressions'] = sorted([dict(name=n,
            ratio=by[n]['M']['solver_wall_sec']/max(1e-12,by[n]['P']['solver_wall_sec']),
            P_seconds=by[n]['P']['solver_wall_sec'], M_seconds=by[n]['M']['solver_wall_sec'])
            for n in pair], key=lambda r:r['ratio'], reverse=True)[:5]
        if rows[0]['backend'] == 'cpu':
            result.pop('worst_time_regressions')
        result['worst_iteration_regressions'] = sorted([dict(name=n,
            ratio=by[n]['M']['iterations']/max(1,by[n]['P']['iterations']),
            P=by[n]['P']['iterations'], M=by[n]['M']['iterations'])
            for n in pair], key=lambda r:r['ratio'], reverse=True)[:5]
        # A CPU run stopped by wall time is censored, including its pairing.
        eligible = [n for n in by if not any(by[n][c]['status']=='TIME_LIMIT' for c in policies)]
        result['cpu_iteration_eligible_names'] = eligible
        result['cpu_iteration_censored_names'] = sorted(set(by)-set(eligible))
        result['cpu_par2_iterations'] = {c:sum(by[n][c]['iterations'] if solved(by[n][c])
             else 2*by[n][c]['iteration_limit'] for n in eligible) for c in policies}
        if rows[0]['backend'] != 'cpu':
            for key in ('worst_iteration_regressions', 'cpu_iteration_eligible_names',
                        'cpu_iteration_censored_names', 'cpu_par2_iterations'):
                result.pop(key)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('results', type=Path)
    p.add_argument('--output', type=Path)
    a=p.parse_args()
    r=report(json.loads(a.results.read_text()))
    if a.output:
        a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))


if __name__ == '__main__':
    main()
