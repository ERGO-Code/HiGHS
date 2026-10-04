import gzip,json,sys,shutil,subprocess
from pathlib import Path
sys.path.insert(0,'/home/yzhou/Documents/Github/hipdlp2/docs/hipdlp-scaling-study/tools')
from neso_canonicalize import canonicalize
root=Path('/tmp/hipdlp-four20')
for d in ['models','original','runs']: (root/d).mkdir(exist_ok=True)
old=json.loads(Path('/tmp/hipdlp-expanded-gpu/manifest.json').read_text())
records=[]
names='25fv47,pds-20,pypsa-eur-240t5n,nug20,mod2,a2864,neso-2245,neso-2005,steiner_t,savsched1,watson_1,graph40-40,fome13,neos1,pds-100,cre-b,bnl2,ship12l,pilot87,stocfor3'.split(',')
(root/'names.json').write_text(json.dumps(names,indent=2)+'\n')
for name in names:
    target=root/'models'/(name+'-canonical.mps')
    existing=next((x for x in old if Path(x['target']).name==target.name),None)
    if existing:
        if not target.exists(): target.symlink_to(existing['target'])
        record=dict(existing,target=str(target),reused_audit=True)
    else:
        source=root/'original'/(name+'.mps')
        with gzip.open('/srv/mps_da/'+name+'.mps.gz','rb') as f,source.open('wb') as out: shutil.copyfileobj(f,out)
        record=canonicalize(source,target,1e-9)
        audit=subprocess.run(['/tmp/hipdlp-expanded-gpu/input-audit',str(source),str(target)],capture_output=True,text=True,check=True)
        record['audit']=json.loads(audit.stdout)
    records.append(dict(record,name=name))
    (root/'manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    print(name,record['audit'],flush=True)
