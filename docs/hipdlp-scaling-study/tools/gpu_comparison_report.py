#!/usr/bin/env python3
"""Archive compact GPU evidence and render comparable wall-time tables."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--expected-runs', type=int, default=30)
    args = p.parse_args()
    records = json.loads((args.root/'results.json').read_text())
    chosen = {(r['name'],r['solver']): r for r in records if r.get('selected')}
    assert len(chosen) == sum(bool(r.get('selected')) for r in records), 'Duplicate selected run'
    assert len(chosen) == args.expected_runs, 'Comparison is incomplete'
    names = list(dict.fromkeys(r['name'] for r in records))
    solvers = ['rewrite','hpr','cupdlpx','cupdlpx-no-gm']
    metadata = json.loads((args.root/'metadata.json').read_text())
    for path, expected in metadata['binaries'].items():
        assert sha(Path(path)) == expected, path
    args.evidence.mkdir(parents=True, exist_ok=True)
    for name in ('results.json','metadata.json','manifest.json',
                 'hpr-instrumentation.patch','cupdlpx-instrumentation.patch'):
        shutil.copyfile(args.root/name, args.evidence/name)
    for file in (args.root/'runs').glob('*/*/*'):
        relative = file.relative_to(args.root)
        if file.suffix == '.json' or file.name.endswith('_summary.txt'):
            target = args.evidence/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(file, target)
        elif file.suffix == '.log':
            target = args.evidence/(str(relative)+'.gz')
            target.parent.mkdir(parents=True, exist_ok=True)
            with file.open('rb') as source, target.open('wb') as out:
                with gzip.GzipFile(filename='', mode='wb', fileobj=out, mtime=0) as zipped:
                    shutil.copyfileobj(source,zipped)
    for run in (args.root/'runs').glob('*/*'):
        hashes = {p.name: sha(p) for p in sorted(run.glob('*.txt'))
                  if not p.name.endswith('_summary.txt')}
        if hashes:
            (args.evidence/run.relative_to(args.root)/'solution-sha256.json').write_text(
                json.dumps(hashes, indent=2)+'\n')
    evidence = 'hipdlp-scaling-study/data/' + args.evidence.name
    lines = ['# Expanded GPU comparison after preprocessing optimization', '',
             'Date: 2026-10-04 (Europe/London). GPU: RTX 4070 Ti SUPER, 16 GB.',
             'CPU: Intel Core i7-14700. CUDA 12.8.61; driver 580.173.02.',
             '`OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` for every solver.', '',
             'The rewrite source is `'+metadata['commits']['rewrite']+'`. HPR-LP uses',
             'the official C implementation at `'+metadata['commits']['hpr']+'`;',
             'cuPDLPx is `'+metadata['commits']['cupdlpx']+'`.', '',
             '## Measurement and acceptance', '',
             'All runs use GPU device 0 with presolve disabled and run serially.',
             'Rewrite uses auto scaling and PID. HPR-C requests its default CR',
             '(with its upstream guard). cuPDLPx is measured both with its default',
             '12-sweep geometric stage and with `--geo_mean_iter 0`.', '',
             'Every exported primal, row-dual and reduced-cost vector has 17 digits.',
             'A persistent, read-only host verifier loads each LP once and checks',
             '`hipdlp-original-kkt/2` at `1e-4` between solver runs. Its work is outside',
             'all benchmark timers. Rewrite also retains its production KKT screen',
             'and host acceptance check inside the solver.', '',
             'For external solvers the predeclared native tolerance grid is `1e-4`,',
             '`1e-6`, `1e-8`. The first run on this grid passing common KKT is selected;',
             'a time limit or missing export ends the search. This does not identify',
             'the earliest possible KKT-passing iterate or the fastest intermediate',
             'tolerance. It can require extra accuracy from an external solver.',
             'Each attempt has a 300-second native budget and an outer safety timeout.',
             'Native budgets retain each implementation\'s scope; cuPDLPx uses CPU',
             '`clock()` internally. The reported comparison timers all use wall time.', '',
             'A table entry is one selected attempt. Tolerance calibration and failed',
             'attempts are retained in the evidence and are excluded from that entry.',
             'These are single runs with a fixed order; small timing differences do',
             'not establish a statistically reliable ranking.', '',
             'Six added models come from `/srv/mps_da/*.mps.gz`. Canonicalization',
             'removed 0 entries and added 0 bound-only columns on each of these six.',
             'Exact HiGHS input audits pass for matrices, costs, bounds, names, sense',
             'and offset. NESO uses the previously audited canonical inputs described',
             'in [the original comparison](hipdlp-scaling-study/NESO_GPU_COMPARISON.md).',
             'The failed default-geometric NESO runs are retained from that study;',
             'the current NESO runs retest rewrite, HPR-C and cuPDLPx without the',
             'geometric stage using the new common timers.', '',
             '## End-to-end process wall time', '',
             'Seconds from process launch to exit, including program startup, input',
             'reading, GPU preparation, solve, full-precision export and cleanup.',
             '`FAIL` means the selected run did not pass common KKT.', '']

    def table(field):
        result = ['| Model | Rewrite | HPR-LP-C | cuPDLPx default | cuPDLPx no geometric |',
                  '|---|---:|---:|---:|---:|']
        for name in names:
            cells=[]
            for solver in solvers:
                r=chosen.get((name,solver))
                if r is None:
                    cells.append('prior FAIL')
                elif not r['verification']['kkt_pass']:
                    cells.append('FAIL' + (' ('+format(r[field],'.3f')+')' if field in r else ''))
                else:
                    cells.append(format(r[field],'.3f') + (' †' if r['native']['status'] != 'OPTIMAL' else ''))
            result.append('| '+name+' | '+' | '.join(cells)+' |')
        return result

    lines += table('process_wall_sec')
    lines += ['', '## Model sizes', '',
              '| Model | Rows | Columns | Nonzeros |',
              '|---|---:|---:|---:|']
    for name in names:
        d = chosen[(name,'rewrite')]['dimensions']
        lines.append(f"| {name} | {d['rows']:,} | {d['cols']:,} | {d['nnz']:,} |")
    lines += ['', '## Solver-call wall time', '',
              'Seconds around the call receiving the already parsed LP and returning',
              'the solution. Includes each solver\'s model preparation, scaling, GPU',
              'initialization, iterations and postsolve; excludes text parsing and',
              'solution export. HPR-C and cuPDLPx have passive CLI timer/export patches.', '']
    lines += table('solver_wall_sec')
    if any(r['verification']['kkt_pass'] and r['native']['status'] != 'OPTIMAL' for r in chosen.values()):
        lines += ['', '† Common KKT passed at a native limit rather than native OPTIMAL.']
    lines += ['', '## Input reading and rewrite preprocessing', '',
              'All values are seconds. Preprocessing is already included in the',
              'rewrite solver-call time above. Parsing is outside that call.', '',
              '| Model | Rewrite read | HPR-C read | cuPDLPx default read | cuPDLPx no geometric read | Rewrite preprocessing |',
              '|---|---:|---:|---:|---:|---:|']
    for name in names:
        values = [format(chosen[(name,s)]['read_wall_sec'],'.3f')
                  if (name,s) in chosen and 'read_wall_sec' in chosen[(name,s)] else '—' for s in solvers]
        values.append(format(chosen[(name,'rewrite')]['native']['preprocessing_seconds'],'.3f'))
        lines.append('| '+name+' | '+' | '.join(values)+' |')
    lines += ['',
              '## Selected native tolerances and iterations', '',
              '| Model | Solver | Native tolerance | Iterations | Native status | Common KKT |',
              '|---|---|---:|---:|---|---|']
    for name in names:
        for solver in solvers:
            r=chosen.get((name,solver))
            if r:
                lines.append('| '+ ' | '.join((name,solver,r['native_tolerance'],str(r['native'].get('iterations','—')),
                             r['native']['status'],'PASS' if r['verification']['kkt_pass'] else 'FAIL'))+' |')
    lines += ['', '## Interpretation', '',
              'The rewrite has the shortest solver-call time on `25fv47`,',
              '`pypsa-eur-240t5n`, `nug20`, and `mod2` in this sample. Including input',
              'reading changes the ordering on `nug20`. `pds-20` is close at the',
              'solver-call level. `a2864` still exposes a substantial rewrite',
              'preprocessing/solve cost and slower input reading.', '',
              'The NESO preprocessing optimization improves rewrite time but does',
              'not close its gap to HPR-C on these two models. See the controlled',
              '[before/after measurements](hipdlp-preprocessing.md), where iteration',
              'counts and exported solution vectors are unchanged.', '',
              'The current cuPDLPx no-geometric NESO runs already pass common KKT',
              'at native `1e-6`; the previous study used `1e-8` for its tight',
              'ablation. In particular, the shorter 2005 time here must not be',
              'attributed to a cuPDLPx code optimization. Both native tolerance',
              'and the reported timing scope differ from that older table.', '',
              'cuPDLPx geometric scaling helps `pypsa-eur-240t5n` and hurts `mod2`',
              'in this experiment; disabling it globally is not supported by these',
              'results. Its min/max geometric scaling is distinct from CR.', '',
              'These rankings are for reaching the common KKT acceptance contract',
              'under the stated tolerance grid. They do not isolate kernel speed',
              'or measure all solvers at identical internal stopping conditions.', '',
              '## Reproduction and evidence', '',
              'The benchmark tools are under `docs/hipdlp-scaling-study/tools/`:',
              '`gpu_comparison_rewrite.cpp`, `gpu_comparison_verify.cpp`,',
              '`gpu_comparison_run.py`, and `gpu_comparison_report.py`.',
              'The input normalization and audit reuse `neso_canonicalize.py` and',
              '`neso_input_audit.cpp`. Numerical solver sources are unchanged by',
              'this expanded comparison. External CLI patches, binary hashes, input',
              'and solution hashes, exact argv, all calibration attempts, KKT metrics and compressed',
              'logs are retained in [the evidence directory]('+evidence+'/).', '',
              'Build the two adapters against the same Release CUDA library. For example:', '',
              '```sh',
              'c++ -O3 -std=c++11 -Ihighs -Ibuild-scaling-cuda \\',
              '  docs/hipdlp-scaling-study/tools/gpu_comparison_rewrite.cpp \\',
              '  -Lbuild-scaling-cuda/lib -Wl,-rpath,"$PWD/build-scaling-cuda/lib" \\',
              '  -lhighs -o /tmp/hipdlp-expanded-gpu/rewrite',
              'c++ -O3 -std=c++11 -Ihighs -Ibuild-scaling-cuda \\',
              '  docs/hipdlp-scaling-study/tools/gpu_comparison_verify.cpp \\',
              '  -Lbuild-scaling-cuda/lib -Wl,-rpath,"$PWD/build-scaling-cuda/lib" \\',
              '  -lhighs -o /tmp/hipdlp-expanded-gpu/verify-server',
              '```', '',
              'External build options are recorded in the original comparison.',
              'Apply the archived passive patches to the pinned external sources.',
              'Place audited canonical inputs in the `models/` subdirectory using',
              'the filenames recorded in each command JSON. Keep all solver binaries',
              'and the shared library fixed for the whole serial run.', '',
              '```sh',
              'python3 docs/hipdlp-scaling-study/tools/gpu_comparison_run.py \\',
              '  --root /tmp/hipdlp-expanded-gpu --external /tmp/neso-gpu-comparison \\',
              '  --names 25fv47,pds-20,pypsa-eur-240t5n,nug20,mod2,a2864,neso-2245,neso-2005 \\',
              '  --limit 300',
              '```', '']
    args.report.write_text('\n'.join(lines))
    hashes={str(p.relative_to(args.evidence)):sha(p) for p in sorted(args.evidence.rglob('*'))
            if p.is_file() and p.name != 'sha256.json'}
    (args.evidence/'sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    print('Archived',len(records),'attempts;',len(chosen),'selected runs.')


if __name__ == '__main__':
    main()
