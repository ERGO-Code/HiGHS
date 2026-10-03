#!/usr/bin/env python3
"""Train a depth-two, asymmetric-loss tree from verified CPU iteration labels.

Time saving is predicted from iteration savings and the live default's measured
seconds/iteration. CPU parallel wall times and CUDA censored runs are not labels.
A 3x penalty on regression and leaf support/family guards favor continuing P.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import numpy as np


def family(name):
    # Conservative: merge numbered sizes, B/filtered aliases and common prefixes.
    s = name.lower().replace('-filtered', '')
    if s in ('distillation', 'distillation_dual'):
        return 'distillation'
    for prefix in ('neos', 'neso', 'sc', 'sctap', 'qap', 'nug', 'pds', 'osa',
                   'fome', 'cre', 'ken', 'ship', 'savsched', 'psched', 'watson'):
        if s.startswith(prefix):
            return prefix
    s = re.sub(r'\d+', '', s).strip('-_')
    if s.endswith('b') and name.lower().endswith('b'):
        s = s[:-1]
    return s


def labels(data):
    rows = {}
    for file in ('cpu_all.csv', 'cpu_long.csv', 'cpu_patched.csv'):
        for row in csv.DictReader((data / file).open()):
            if row['config'] in ('P', 'L', 'A'):
                rows[row['instance'], row['config']] = row
    result = {}
    for name in (data / 'instances.txt').read_text().split():
        triple = [rows.get((name, c), {}) for c in ('P', 'L', 'A')]
        if all(r.get('status') == 'OPTIMAL' and r.get('kkt_pass') == 'True' for r in triple):
            result[name] = [int(r['iterations']) for r in triple]
    return result


def fit(X, costs, groups, depth=2, min_leaf=8):
    nodes = []
    def leaf(ids):
        baseline = costs[ids, 0].sum()
        delta = costs[ids] - costs[ids, 0, None]
        loss = np.where(delta > 0, 3 * delta, delta).sum(axis=0)
        action = int(np.argmin(loss))
        ratio = 1 + loss[action] / max(1, baseline)
        regression_fraction = float(np.mean(delta[:, action] > 0))
        # Thresholds are fixed before examining validation or external results.
        if ratio >= .75 or regression_fraction > .20 or len(set(groups[ids])) < 3:
            action, ratio = 0, 1.
        return dict(action=('P', 'L', 'A')[action], ratio=float(ratio),
                    baseline_iterations=float(np.mean(costs[ids, 0])),
                    support=len(ids), families=len(set(groups[ids])),
                    regression_fraction=regression_fraction,
                    loss=float(loss[action]))
    def grow(ids, level):
        index = len(nodes)
        record = leaf(ids)
        nodes.append(record)
        if level == depth or len(ids) < 2 * min_leaf:
            return index
        best, split = record['loss'], None
        for j in range(X.shape[1]):
            values = np.unique(X[ids, j])
            if len(values) < 2:
                continue
            cuts = np.unique(np.quantile(values, np.linspace(.1, .9, 9)))
            for cut in cuts:
                left, right = ids[X[ids, j] <= cut], ids[X[ids, j] > cut]
                if min(len(left), len(right)) < min_leaf:
                    continue
                score = leaf(left)['loss'] + leaf(right)['loss']
                if score < best - .01:
                    best, split = score, (j, float(cut), left, right)
        if split is None:
            return index
        j, cut, left, right = split
        nodes[index] = dict(feature_index=j, threshold=cut,
                            left=grow(left, level+1), right=grow(right, level+1))
        return index
    grow(np.arange(len(X)), 0)
    return nodes


def predict(nodes, x):
    node = nodes[0]
    while 'feature_index' in node:
        node = nodes[node['left'] if x[node['feature_index']] <= node['threshold'] else node['right']]
    return node


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--features', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    study = Path(__file__).resolve().parents[1]
    data = study / 'data'
    truth = labels(data)
    external = (data / 'heldout-instances.txt').read_text().split() + (data / '2026-10-02-large20/instances.txt').read_text().split()
    external_families = {family(n) for n in external}
    rows = [r for r in json.loads(a.features.read_text()) if r['name'] in truth and 'primal_200' in r['features']]
    excluded = [r['name'] for r in rows if family(r['name']) in external_families]
    rows = [r for r in rows if family(r['name']) not in external_families]
    keys = sorted(set.intersection(*(set(r['features']) for r in rows)) - {'seconds_per_iteration'})
    # Family split is determined from names alone, independent of label outcomes.
    validation_families = {family(r['name']) for r in rows if int(hashlib.sha256(family(r['name']).encode()).hexdigest()[:8], 16) % 4 == 0}
    train = [r for r in rows if family(r['name']) not in validation_families]
    valid = [r for r in rows if family(r['name']) in validation_families]
    X = np.array([[r['features'][k] for k in keys] for r in train])
    groups = np.array([family(r['name']) for r in train])
    C = np.array([truth[r['name']] for r in train], dtype=float)
    C[:, 1:] += 200  # Default pilot is discarded on cold cross-scaling switches.
    tree = fit(X, C, groups)
    evaluation = []
    for r in valid:
        leaf = predict(tree, np.array([r['features'][k] for k in keys]))
        action = leaf['action']
        counts = truth[r['name']]
        selected = counts[('P', 'L', 'A').index(action)] + (200 if action != 'P' else 0)
        evaluation.append(dict(name=r['name'], family=family(r['name']), chosen=action,
                               P=counts[0], L=counts[1], A=counts[2], selected=selected,
                               ratio=selected/max(1, counts[0])))
    a.output.mkdir(parents=True, exist_ok=True)
    text = []
    for i, n in enumerate(tree):
        if 'feature_index' in n:
            n['feature'] = keys[n['feature_index']]
            text.append(f'split {i} {n["feature"]} {n["threshold"]:.17g} {n["left"]} {n["right"]}')
        else:
            text.append(f'leaf {i} {n["action"]} {n["ratio"]:.17g} {n["baseline_iterations"]:.17g} {n["support"]}')
    (a.output / 'model.txt').write_text('\n'.join(text) + '\n')
    report = dict(schema='one-pilot-tree/1', default='P', depth=2,
                  regression_penalty=3, minimum_leaf=8, minimum_families=3,
                  minimum_predicted_saving=.25, maximum_training_regression_fraction=.2,
                  label_kind='verified CPU iterations, plus 200 on switch; no timeout labels',
                  runtime_gate='savings > 1.5*(elapsed pilot + feature cost + estimated fresh setup)',
                  features=keys, tree=tree,
                  label_source_sha256={f: hashlib.sha256((data/f).read_bytes()).hexdigest()
                      for f in ('cpu_all.csv', 'cpu_long.csv', 'cpu_patched.csv')},
                  training_labels={r['name']: dict(zip(('P','L','A'), truth[r['name']])) for r in train},
                  validation_labels={r['name']: dict(zip(('P','L','A'), truth[r['name']])) for r in valid},
                  training_names=[r['name'] for r in train],
                  validation_names=[r['name'] for r in valid],
                  family_map={r['name']: family(r['name']) for r in rows},
                  external_families=sorted(external_families), excluded_external_family_names=excluded,
                  excluded_early_solves=[r['name'] for r in json.loads(a.features.read_text()) if 'primal_200' not in r['features']],
                  validation=evaluation)
    if evaluation:
        report['validation_summary'] = dict(count=len(evaluation),
            baseline_iteration_sum=sum(r['P'] for r in evaluation),
            selected_iteration_sum=sum(r['selected'] for r in evaluation),
            switches=sum(r['chosen'] != 'P' for r in evaluation),
            worst_ratio=max(r['ratio'] for r in evaluation))
    (a.output / 'training.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(train=len(train), families=len(set(groups)),
                         validation=len(valid), tree=tree,
                         summary=report.get('validation_summary')), indent=2))


if __name__ == '__main__':
    main()
