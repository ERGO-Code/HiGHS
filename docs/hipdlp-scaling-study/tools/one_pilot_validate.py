#!/usr/bin/env python3
"""Integration checks for useful default work and cold-switch budget accounting."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--driver', type=Path, required=True)
    p.add_argument('--instance', type=Path, required=True)
    p.add_argument('--backend', choices=['cpu','cuda'], default='cpu')
    a=p.parse_args()
    def run(config, limit, model):
        command=[str(a.driver.resolve()),str(a.instance.resolve()),'adlittle',
                 a.backend, config,'60',str(limit),str(model)]
        process=subprocess.run(command,text=True,capture_output=True,check=True)
        lines=[s[11:] for s in process.stdout.splitlines() if s.startswith('ONE_RESULT ')]
        assert len(lines)==1, process.stdout
        return json.loads(lines[0])
    with tempfile.TemporaryDirectory(prefix='hipdlp-onepilot-check-') as directory:
        model=Path(directory)/'model.txt'
        model.write_text('leaf 0 P 1 1000000 12\n')
        baseline=run('P',100000,'-')
        observed=run('M',100000,model)
        assert baseline['status']==observed['status']=='OPTIMAL'
        assert baseline['iterations']==observed['iterations']
        assert observed['discarded_iterations']==0
        for candidate in ('L','A'):
            # Artificial leaf deliberately forces a cold switch; this tests
            # solver orchestration, not predictive model quality.
            model.write_text(f'leaf 0 {candidate} 0.001 100000000 12\n')
            selected=run('M',100000,model)
            standalone=run(candidate,100000,'-')
            assert selected['status']==standalone['status']=='OPTIMAL'
            assert selected['chosen']==candidate
            assert selected['iterations']==standalone['iterations']+200
            assert selected['discarded_iterations']==200
            limited=run('M',441,model)
            assert limited['iterations']==441
            assert limited['discarded_iterations']==200
            short=run('M',200,model)
            assert short['iterations']==200
            assert short['discarded_iterations']==0
            assert short['chosen']=='P'
        model.write_text('split 0 unknown_feature 0 1 2\nleaf 1 A 0 100000000 12\nleaf 2 L 0 100000000 12\n')
        unknown=run('M',100000,model)
        assert unknown['iterations']==baseline['iterations']
        assert unknown['discarded_iterations']==0
    print(f'{a.backend}: continue, cold L/A switch, global limit, short budget and missing-feature fallback passed')


if __name__=='__main__':
    main()
