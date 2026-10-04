import hashlib,json,os,pathlib,subprocess,time
root=pathlib.Path('/tmp/hipdlp-cr-guard-study')
models=pathlib.Path('/tmp/hipdlp-expanded-gpu/models')
old=json.loads(pathlib.Path('/tmp/hipdlp-expanded-gpu/results.json').read_text())
records=[]
env=dict(os.environ,CUDA_VISIBLE_DEVICES='0',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
for name in ('25fv47','pds-20','pypsa-eur-240t5n','nug20','mod2','a2864','neso-2245','neso-2005'):
    print('START',name,flush=True)
    out=root/'runs'/name
    out.mkdir(parents=True,exist_ok=True)
    model=models/(name+'-canonical.mps')
    record=dict(name=name,baseline=next(r for r in old if r['name']==name and r['solver']=='rewrite'),commands=[])
    for phase,cmd in [('preprocess',[str(root/'new/preprocess'),str(model),'auto','1']),
                      ('solve',[str(root/'new/profile'),str(model),'auto',str(out/'solution.bin'),'120'])]:
        record['commands'].append(cmd)
        with (out/(phase+'.log')).open('w') as log:
            run=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=300)
        record[phase+'_exit']=run.returncode
        lines=(out/(phase+'.log')).read_text().splitlines()
        parsed=[json.loads(line) for line in lines if line.startswith('{')]
        if parsed: record[phase]=parsed[-1]
        else: raise RuntimeError((name,phase,lines[-10:]))
    solution=out/'solution.bin'
    record['solution_sha256']=hashlib.sha256(solution.read_bytes()).hexdigest()
    prior=pathlib.Path('/tmp/neso-preprocessing-study/optimized')/(name+'.solution')
    if prior.exists():
        record['same_solution_as_previous_auto']=record['solution_sha256']==hashlib.sha256(prior.read_bytes()).hexdigest()
    (out/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    records.append(record)
    (root/'results.json').write_text(json.dumps(records,indent=2)+'\n')
    print('DONE',name,'CR rejected',record['preprocess']['curtis_reid_rejected'],'solve',round(record['solve']['total_wall_seconds'],3),'KKT',record['solve']['kkt_pass'],'iterations',record['solve']['iterations'],flush=True)
