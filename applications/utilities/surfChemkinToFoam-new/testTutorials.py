#!/usr/bin/env python3
"""Run converters on isolated output directories; never alter tutorial inputs."""
from pathlib import Path
import subprocess
import re
import json
import tempfile
import math
import os

repo=Path(os.environ.get('POROUS_PBR_REPO',str(Path(__file__).resolve().parents[3])))
root=repo/'tutorials/00_surfChemkinToFoam'
work=Path(tempfile.mkdtemp(prefix='surface-converter-regression-'))
results=[]
for tutorial in sorted(p for p in root.iterdir() if p.is_dir()):
    target=work/tutorial.name
    target.mkdir()
    status={}
    for suffix,exe in [('old','surfChemkinToFoam'),('new','surfChemkinToFoam-new')]:
        cmd=[exe,str(tutorial/'surf.inp'),str(tutorial/'surfTherm.dat'),str(tutorial/'surfaceThermoDict'),str(target/('reactions.'+suffix)),str(target/('thermo.'+suffix))]
        proc=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        (target/(suffix+'.log')).write_text(proc.stdout)
        status[suffix]=proc.returncode
        if proc.returncode: print(tutorial.name,suffix,proc.returncode,proc.stdout[-1000:])
    if not any(status.values()):
        def kinetics(path):
            text=path.read_text()
            values=[]
            for block in re.findall(r'un-named-reaction-\d+\s*\{(.*?)\n\s*\}',text,re.S):
                rate=re.search(r'type\s+(\w+);',block).group(1)
                equation=re.search(r'reaction\s+([^;]+);',block).group(1)
                # Chemical whitespace and coefficient form may be normalized.
                equation=re.sub(r'\s+','',equation)
                numbers=[float(re.search(r'\b'+key+r'\s+([^;]+);',block).group(1)) for key in ('A','beta','Ta')]
                cov=re.findall(r'\(([^\s()]+\(s\))\s+\(([^)]+)\)\)',block)
                values.append((rate,equation,numbers,[(s,[float(x) for x in c.split()]) for s,c in cov]))
            return values
        old=kinetics(target/'reactions.old');new=kinetics(target/'reactions.new')
        assert len(old)==len(new)>0,(tutorial.name,len(old),len(new))
        for a,b in zip(old,new):
            assert a[0:2]==b[0:2],(tutorial.name,a,b)
            assert all(math.isclose(x,y,rel_tol=1e-9,abs_tol=1e-12) for x,y in zip(a[2],b[2])),(tutorial.name,a,b)
            assert a[3]==b[3],(tutorial.name,'coverage mismatch')
        status['matching_reactions']=len(new)
        def properties(path):
            text=path.read_text()
            names=re.search(r'\bspecies\s*\d*\s*\((.*?)\);',text,re.S).group(1).split()
            result={}
            for name in names:
                match=re.search(r'(?m)^'+re.escape(name)+r'\s*\{',text)
                assert match,(path,name)
                start=match.end();depth=1;end=start
                while depth:
                    if text[end]=='{':depth+=1
                    if text[end]=='}':depth-=1
                    end+=1
                block=text[start:end-1]
                result[name]={key:float(re.search(r'\b'+key+r'\s+([^;]+);',block).group(1)) for key in ('molWeight','siteDensity','rho','Cp','Cv','Hf','mu','Pr')}
            return result
        a=properties(target/'thermo.old');b=properties(target/'thermo.new')
        assert set(a)==set(b),(tutorial.name,'species list mismatch')
        corrected=[]
        for name in a:
            for key in a[name]:
                if key=='molWeight' and name=='NI(s)':
                    assert math.isclose(b[name][key],58.69,rel_tol=1e-3)
                    corrected.append(name)
                else: assert math.isclose(a[name][key],b[name][key],rel_tol=1e-9,abs_tol=1e-12),(tutorial.name,name,key,a[name][key],b[name][key])
        status['matching_thermo_species']=len(b)
        if corrected:status['corrected_bare_metal_molecularWeights']=corrected
    results.append({'tutorial':tutorial.name,**status})
(work/'results.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
print('Regression outputs:',work)
assert all(item['new']==0 for item in results),'New converter did not convert all tutorials'
