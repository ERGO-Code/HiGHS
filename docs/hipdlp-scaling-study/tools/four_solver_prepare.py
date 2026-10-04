#!/usr/bin/env python3
"""Prepare the predeclared 20-model LP suite and audit exact HiGHS ingestion."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
from neso_canonicalize import canonicalize

NAMES = '25fv47,pds-20,pypsa-eur-240t5n,nug20,mod2,a2864,neso-2245,neso-2005,steiner_t,savsched1,watson_1,graph40-40,fome13,neos1,pds-100,cre-b,bnl2,ship12l,pilot87,stocfor3'.split(',')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, default=Path('/srv/mps_da'))
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--audit', type=Path, required=True)
    args = p.parse_args()
    for sub in ['models','original','runs']:
        (args.root/sub).mkdir(parents=True, exist_ok=True)
    (args.root/'names.json').write_text(json.dumps(NAMES,indent=2)+'\n')
    records = []
    for name in NAMES:
        source = args.root/'original'/(name+'.mps')
        target = args.root/'models'/(name+'-canonical.mps')
        with gzip.open(args.source/(name+'.mps.gz'),'rb') as f, source.open('wb') as out:
            shutil.copyfileobj(f,out)
        record = canonicalize(source,target,1e-9)
        audit = subprocess.run([str(args.audit),str(source),str(target)],capture_output=True,text=True,check=True)
        record.update(name=name,audit=json.loads(audit.stdout))
        assert record['audit']['equal_after_highs_default_filter']
        records.append(record)
        (args.root/'manifest.json').write_text(json.dumps(records,indent=2)+'\n')
        print(name,record['audit'],flush=True)


if __name__ == '__main__':
    main()
