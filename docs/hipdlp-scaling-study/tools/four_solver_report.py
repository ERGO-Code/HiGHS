#!/usr/bin/env python3
"""Archive and tabulate all 80 attempts without tolerance selection."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

SOLVERS = ['legacy','rewrite','hpr','cupdlpx']
LABELS = ['Legacy','Rewrite','HPR-LP-C','cuPDLPx']


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    args = p.parse_args()
    root, archive = args.root, args.archive
    data = json.loads((root/'results.json').read_text())
    names = json.loads((root/'names.json').read_text())
    assert len(data) == len(names)*4 == 80
    assert len({(r['name'],r['solver']) for r in data}) == 80
    meta = json.loads((root/'metadata.json').read_text())
    for path, expected in meta['artifacts'].items():
        assert sha(Path(path)) == expected, path
    assert sha(root/'build/lib/libhighs.so') == meta['library_sha256']
    manifest = json.loads((root/'manifest.json').read_text())
    assert len(manifest)==20 and all(r['audit']['equal_after_highs_default_filter'] for r in manifest)
    inputs = {r['name']:r for r in manifest}
    for r in manifest:
        assert sha(Path(r['target'])) == r['target_sha256']
    archive.mkdir(parents=True,exist_ok=True)
    for f in ['results.json','manifest.json','metadata.json','names.json','prepare.log','run.log']:
        shutil.copy2(root/f,archive/f)
    for source,name in [(root/'hpr-passive-cli.patch','hpr-passive-cli.patch'),
                        (root/'cupdlpx-passive-cli.patch','cupdlpx-passive-cli.patch'),
                        (root/'build/CMakeCache.txt','hipdlp-CMakeCache.txt')]:
        with gzip.open(archive/(name+'.gz'),'wb') as out:
            out.write(source.read_bytes())
        # Preserve raw evidence byte-for-byte inside gzip, including whitespace.
        if (archive/name).exists(): (archive/name).unlink()
    shutil.copy2(root/'prepare.py',archive/'executed-input-preparation.py')
    shutil.copy2(root/'metadata.py',archive/'executed-metadata-capture.py')
    for f in ['configure','build']:
        source=Path('/tmp/hipdlp-four20-'+f+'.log')
        with gzip.open(archive/(f+'.log.gz'),'wb') as out:
            out.write(source.read_bytes())
    rows=[]
    by={(r['name'],r['solver']):r for r in data}
    for r in data:
        assert r['native_tolerance']==1e-4 and r['presolve'] and r['scaling']=='Ruiz10+PC(alpha=1)'
        assert r['sha256'] == inputs[r['name']]['target_sha256']
        assert r['dimensions'] == {k:inputs[r['name']]['audit'][k] for k in ['rows','cols','nnz']}
        assert r['executable_sha256']==sha(Path(r['command'][0]))
        source=root/'runs'/r['name']/r['solver']
        dest=archive/'runs'/r['name']/r['solver'];dest.mkdir(parents=True,exist_ok=True)
        for f in source.iterdir():
            if f.name=='run.log':
                with gzip.open(dest/'run.log.gz','wb') as out: out.write(f.read_bytes())
            elif f.suffix=='.json' or f.name.endswith('_summary.txt'):
                shutil.copy2(f,dest/f.name)
        for export in r['exports'].values():
            assert sha(Path(export['path']))==export['sha256']
        log=(source/'run.log').read_text()
        n=r['native']; v=r['verification']
        if n['status'] == 'OPTIMAL':
            assert r['exit_code'] == 0
            if r['solver'] in ['legacy','rewrite']:
                assert 'Presolve status:' in log
                if r['solver']=='legacy' and n.get('iterations',0)>0:
                    assert 'HiPDLP first order PDLP solver on a GPU' in log
            elif r['solver']=='hpr':
                assert re.search(r'Presolver:\s+PSLP',log)
                assert re.search(r'Curtis-Reid:\s+Disabled',log)
                assert re.search(r'Bounds/Cost:\s+Disabled',log)
            else:
                assert 'Running presolver (PSLP' in log
                assert 'Ruiz scaling (10 iterations)' in log
                assert 'Pock-Chambolle scaling (alpha=1.0000)' in log
                assert re.search(r'geo_mean_iter\s*:\s*0',log)
                assert re.search(r'bound_obj_rescaling\s*:\s*off',log)
        row=dict(model=r['name'],solver=r['solver'],native_status=n['status'],native_tolerance=r['native_tolerance'],
                 common_kkt_pass=v['kkt_pass'],iterations=n.get('iterations'),process_wall_sec=r['process_wall_sec'],
                 read_wall_sec=r.get('read_wall_sec'),solver_wall_sec=r.get('solver_wall_sec'),
                 read_plus_solve_wall_sec=r.get('read_plus_solve_wall_sec'),**r['dimensions'])
        for k in ['presolve_seconds','algorithm_seconds','postsolve_seconds','reduced_rows','reduced_cols','reduced_nnz',
                  'rewrite_preprocessing_seconds','backend_setup_seconds','operator_norm_seconds','iterations_seconds']:
            row[k]=n.get(k)
        for k in ['objective','primal','dual','gap','complementarity','stationarity']:row[k]=v.get(k)
        row['cupdlpx_presolve_reported_sec'] = None
        row['cupdlpx_scaling_reported_sec'] = None
        if r['solver']=='hpr':
            match=re.search(r'Presolved \(PSLP\) Model: nRow = (\d+), nCol = (\d+), nnz = (\d+)',log)
            if match: row.update(zip(['reduced_rows','reduced_cols','reduced_nnz'],map(int,match.groups())))
        if r['solver']=='cupdlpx':
            match=re.search(r'reduced problem\s*:\s*(\d+) rows, (\d+) columns, (\d+) nonzeros',log)
            if match: row.update(zip(['reduced_rows','reduced_cols','reduced_nnz'],map(int,match.groups())))
            summaries=list(source.glob('*_summary.txt'))
            if summaries:
                summary=dict(line.split(': ',1) for line in summaries[0].read_text().splitlines() if ': ' in line)
                for label,key in [('Presolve Time (sec)','cupdlpx_presolve_reported_sec'),('Precondition time (sec)','cupdlpx_scaling_reported_sec')]:
                    if label in summary: row[key]=float(summary[label])
        rows.append(row)
    with (archive/'results.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    metrics={}
    all_native=[name for name in names if all(by[name,s]['native']['status']=='OPTIMAL' for s in SOLVERS)]
    all_kkt=[name for name in names if all(by[name,s]['verification']['kkt_pass'] for s in SOLVERS)]
    for s in SOLVERS:
        rs=[by[name,s] for name in names]
        metrics[s]=dict(native_optimal=sum(r['native']['status']=='OPTIMAL' for r in rs),
                        common_kkt_pass=sum(r['verification']['kkt_pass'] for r in rs),
                        native_time_limit=sum('TIME_LIMIT' in r['native']['status'] for r in rs),
                        total_process_seconds=sum(r['process_wall_sec'] for r in rs))
        if all_native:
            metrics[s]['shifted_geomean_process_seconds_common_native']=math.exp(sum(math.log1p(by[name,s]['process_wall_sec']) for name in all_native)/len(all_native))-1
            metrics[s]['shifted_geomean_solver_seconds_common_native']=math.exp(sum(math.log1p(by[name,s]['solver_wall_sec']) for name in all_native)/len(all_native))-1
    (archive/'summary.json').write_text(json.dumps(dict(solvers=metrics,all_native_optimal_models=all_native,all_common_kkt_models=all_kkt),indent=2)+'\n')
    lines=['# 20 个 LP 的四求解器 GPU 比较：Ruiz + PC、1e-4、presolve on','',
           '2026-10-04。RTX 4070 Ti SUPER 16 GB；Intel i7-14700；CUDA 12.8.61。每个模型、每个求解器运行一次，共 80 次，固定顺序串行执行。', '',
           '## 比较口径','',
           '- 四家原生停止容差都固定为 `1e-4`，没有按结果挑选更紧容差。Ruiz 10 轮 + PC（alpha=1）；关闭 CR、几何平均、L2 和额外 bounds/cost 缩放。Legacy/rewrite 使用 PID，其他启发式保留各家默认。',
           '- Legacy/rewrite 都来自同一 HiPDLP 提交，均启用 CUDA；legacy 日志明确显示 GPU。`hipdlp_scaling=legacy` 只是 rewrite 的 Ruiz+PC 模式，legacy 列实际调用旧求解器。',
           '- 两个 HiPDLP 路径实际执行 `Highs::presolve()`，在 reduced LP 上求解，再执行无 basis 的 `postsolve()`，没有 simplex cleanup。HPR 使用 PSLP 0.0.8；cuPDLPx 使用 PSLP 0.0.11。Presolve 算法并不相同，缩放作用的矩阵也可能不同。',
           '- 主时间采用外部墙钟。读取不含 presolve；求解调用含 presolve、模型准备、scaling、GPU 初始化、迭代及 postsolve。进程总时间另含启动、导出及清理。统一 KKT 检查在计时外进行。',
           '- 原生限时 300 秒，保留各家内部计时范围。HiPDLP 的 presolve 和后续求解分别受限；cuPDLPx 内部使用 CPU clock。外部墙钟兜底 780 秒。超时值不能解释为达标求解时间。',
           '- 环境固定 `CUDA_VISIBLE_DEVICES=0, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1`；HiGHS threads=1。没有同时运行其他求解器或输入准备任务。单次结果不能支持接近时间的统计排名。', '',
           '**精度标签：**所有导出的 x/y/z 都用 `hipdlp-original-kkt/2` 在原始规范化模型上检查 `1e-4`。`†` 表示原生 OPTIMAL，但未通过该共同检查；`T` 表示原生超时；`E` 表示其他失败。原生 `1e-4` 参数相同，不代表 postsolve 后共同指标也相同；带标签的时间不能称为共同精度下的胜出。','',
           '## 汇总','', '| 求解器 | 原生 OPTIMAL / 20 | 原始模型共同 KKT 通过 / 20 | 超时 / 20 |', '|---|---:|---:|---:|']
    for s,label in zip(SOLVERS,LABELS):
        m=metrics[s];lines.append(f"| {label} | {m['native_optimal']} | {m['common_kkt_pass']} | {m['native_time_limit']} |")
    lines.extend(['','## 相同原生容差下的时间概览','',
                  f'四家均原生 OPTIMAL 的共同子集有 {len(all_native)} 个模型；四家均通过共同 KKT 的交集有 {len(all_kkt)} 个。下表是前者的时间描述，不是共同 KKT 达标的性能排名。', '',
                  '| 求解器 | 总时间 SGM-1（秒） | 求解调用 SGM-1（秒） | 20 次进程时间之和（秒） |',
                  '|---|---:|---:|---:|'])
    for s,label in zip(SOLVERS,LABELS):
        m=metrics[s]
        values=[m.get('shifted_geomean_process_seconds_common_native'),m.get('shifted_geomean_solver_seconds_common_native'),m['total_process_seconds']]
        lines.append('| '+label+' | '+' | '.join('—' if value is None else f'{value:.3f}' for value in values)+' |')
    lines.extend(['','SGM-1 = `exp(mean(log(1+t))) - 1`，以 1 秒平移，避免极短用时过度影响几何平均。后面的逐项表保留全部结果。'])
    faster=sum(by[name,'rewrite']['process_wall_sec']<by[name,'legacy']['process_wall_sec'] for name in names)
    a=by['a2864','rewrite']; n=by['neso-2005','rewrite']; old=by['neso-2005','legacy']
    lines.extend(['','## 这轮结果说明什么','',
                  f"- 按上述配置的单次总时间，rewrite 在 {faster}/20 个算例上低于 legacy；NESO-2005、bnl2、pilot87 上较长的求解抵消了这些收益。原生状态与共同 KKT 验收结果必须一起看。",
                  f"- a2864 的 rewrite 总时间 {a['process_wall_sec']:.3f} 秒：读取 {a['read_wall_sec']:.3f} 秒，HiGHS presolve {a['native']['presolve_seconds']:.3f} 秒，rewrite 预处理 {a['native']['rewrite_preprocessing_seconds']:.3f} 秒。这里读取和 presolve 是主要成本；即使省去全部 rewrite 预处理，总时间的改善也不足 1%。",
                  f"- NESO-2005 的 rewrite 使用 {n['native']['iterations']:,} 步，legacy 使用 {old['native']['iterations']:,} 步。Rewrite 的迭代阶段 {n['native']['iterations_seconds']:.3f} 秒、预处理 {n['native']['rewrite_preprocessing_seconds']:.3f} 秒；这次 rewrite 较慢需要同时检查停止条件和迭代路径，不能归因于 scaling 一项。",
                  '- 本测试传给 rewrite 的 OriginalModel 是 presolve 后的 LP，因此其 GPU 筛查/host 验证保证的是 reduced LP 的共同 KKT。Postsolve 后重新检查原始模型，只有 10/20 通过；例如 NESO-2005 的共同 relative dual residual 为约 0.591。要比较真正相同的最终精度，需要统一 postsolve 后的验收，并记录达到该条件的时间。本轮保持原生容差 1e-4，并原样验证返回向量。'])
    def table(title,key):
        lines.extend(['',f'## {title}','', '| 算例 | Legacy | Rewrite | HPR-LP-C | cuPDLPx |','|---|---:|---:|---:|---:|'])
        for name in names:
            cells=[]
            for s in SOLVERS:
                r=by[name,s];status=r['native']['status']
                flag=' T' if 'TIME_LIMIT' in status else ' E' if status!='OPTIMAL' else ' †' if not r['verification']['kkt_pass'] else ''
                value=r.get(key)
                cells.append(('—' if value is None else f'{value:.3f}')+flag)
            lines.append('| '+name+' | '+' | '.join(cells)+' |')
    table('进程总时间（秒）','process_wall_sec')
    table('求解调用时间，含 presolve / postsolve（秒）','solver_wall_sec')
    table('读取时间（秒，标签仍表示该次求解的最终状态）','read_wall_sec')
    lines.extend(['','## HiPDLP 时间拆分（秒）','',
                  '| 算例 | Legacy presolve | Legacy 后续求解 | Rewrite presolve | Rewrite 预处理 | Rewrite 后续求解 |',
                  '|---|---:|---:|---:|---:|---:|---:|'])
    for name in names:
        a,b=by[name,'legacy']['native'],by[name,'rewrite']['native']
        vals=[a.get('presolve_seconds'),a.get('algorithm_seconds'),b.get('presolve_seconds'),b.get('rewrite_preprocessing_seconds'),b.get('algorithm_seconds')]
        lines.append('| '+name+' | '+' | '.join('—' if v is None else f'{v:.3f}' for v in vals)+' |')
    lines.extend(['','Rewrite 预处理已包含在“后续求解”中，不能重复相加。Postsolve、GPU 初始化、范数估计和迭代的可用细分见 CSV；外部原生计时保留在日志/summary，不混充统一墙钟。','',
                  '## 输入规模与审计','',
                  '| 算例 | 行 | 列 | 非零元 | 删除 <=1e-9 的矩阵项 | 补充声明的 bound-only 列 |','|---|---:|---:|---:|---:|---:|'])
    for r in manifest:
        a=r['audit'];lines.append(f"| {r['name']} | {a['rows']} | {a['cols']} | {a['nnz']} | {r['removed']} | {r['bound_only_columns']} |")
    lines.extend(['','输入来自 `/srv/mps_da/*.mps.gz`，统一解压后计时读取 MPS；不计一次性解压/规范化。规范化显式删除绝对值 <=1e-9 的矩阵项，并补齐只出现在 BOUNDS 的变量声明。20 个模型均通过与 HiGHS 默认读入结果逐项比较的审计（矩阵、目标、上下界、名称、目标方向和偏移）。这属于明确的输入政策，NESO 与 stocfor3 不是未经修改的原始文本。', '',
                  '## 版本及复现','',
                  f"- HiPDLP: `{meta['hipdlp_commit']}`。Release，`CUPDLP_GPU=ON`、`HIGHS_HIPDLP_CUDA=ON`，独立构建目录 `{root/'build'}`。",
                  f"- HPR-LP-C: `{meta['hpr_commit']}`；cuPDLPx: `{meta['cupdlpx_commit']}`。仅 CLI 的被动墙钟计时和 17 位解导出补丁，数值求解源码未修改。",
                  '- 运行工具：`four_solver_prepare.py`（输入）、`four_solver_presolve.cpp`（HiPDLP 适配器）、`four_solver_run.py`（80 次串行运行）、`gpu_comparison_verify.cpp`（共同验证）、`four_solver_report.py`（归档）。',
                  '- 证据目录包含所有逐次命令、原生结果、统一验证结果、CSV、完整压缩日志、输入/二进制/导出向量 SHA-256、构建选项和被动 CLI 补丁。大型原始模型与解向量留在 `/tmp/hipdlp-four20`，不提交到 Git。', '',
                  '全部命令和环境见 `metadata.json` 及每次 `command.json`。直接重新运行：', '',
                  '```sh',
                  'python docs/hipdlp-scaling-study/tools/four_solver_run.py --root /tmp/hipdlp-four20 --external /tmp/neso-gpu-comparison --limit 300',
                  '```', '',
                  '已有结果会被跳过；独立重测需使用新 root，准备对应的 models、names.json、hipdlp-presolve 和 verify-server。编译适配器时应链接该次 GPU 构建的 libhighs 并设置相应 rpath。'])
    args.report.write_text('\n'.join(lines)+'\n')
    index={str(f.relative_to(archive)):sha(f) for f in sorted(archive.rglob('*')) if f.is_file() and f.name!='sha256.json'}
    (archive/'sha256.json').write_text(json.dumps(index,indent=2)+'\n')
    print(json.dumps(metrics,indent=2));print('common native:',len(all_native),'common KKT:',len(all_kkt))


if __name__ == '__main__':
    main()
