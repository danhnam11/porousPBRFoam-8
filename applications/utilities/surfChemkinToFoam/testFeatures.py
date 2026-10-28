#!/usr/bin/env python3
"""Functional checks; all generated input/output files stay in a temporary directory."""
from pathlib import Path
import tempfile,subprocess,re,math,json,os

repo=Path(os.environ.get('POROUS_PBR_REPO',str(Path(__file__).resolve().parents[3])))
example=repo/'tutorials/00_surfChemkinToFoam/2011_Appari_Ni_NH3Decomposition'
work=Path(tempfile.mkdtemp(prefix='surface-converter-features-'))
thermo=(example/'surfTherm.dat').read_text()
defaults=(example/'surfaceThermoDict').read_text()
base='''SITE/Ni/ SDEN/2.6D-9/
Ni(s) H(s)
END
REACTIONS KJOULES/MOLE MOLES MWOFF
H2 + 2Ni(s) => 2H(s) 1.0D10 0.5 8.31446261815324
END
'''
passed=[]
def run(name,mechanism=base,therm=thermo,config=defaults,options=(),success=True,contains=None):
    directory=work/name;directory.mkdir()
    for file,text in [('surf.inp',mechanism),('therm.dat',therm),('defaults',config)]: (directory/file).write_text(text)
    command=['surfChemkinToFoam-new',*options,*[str(directory/f) for f in ('surf.inp','therm.dat','defaults','reactions.surf','thermo.solid')]]
    process=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (directory/'log').write_text(process.stdout)
    assert (process.returncode==0)==success,(name,process.stdout)
    if contains: assert contains in process.stdout,(name,process.stdout)
    if not success or '-checkOnly' in options: assert not (directory/'reactions.surf').exists() and not (directory/'thermo.solid').exists(),name
    passed.append(name)
    return directory
def rates(directory):
    t=(directory/'reactions.surf').read_text()
    return [float(re.search(r'\b'+k+r'\s+([^;]+);',t).group(1)) for k in ('A','beta','Ta')]
reference=run('kJ_D_exponent')
for name,unit,factor in [('joules','JOULES/MOLE',1000),('calories','CAL/MOLE',1000/4.184),('kilocalories','KCAL/MOLE',1/4.184),('kelvins','KELVINS',1000/8.31446261815324),('electronvolts','EVOLTS',1000/96485.33212331002)]:
    value=8.31446261815324*factor
    candidate=run(name,base.replace('KJOULES/MOLE',unit).replace('8.31446261815324',format(value,'.16g')))
    assert all(math.isclose(a,b,rel_tol=1e-9) for a,b in zip(rates(reference),rates(candidate))),name
# With gas order 1 and surface order 2, molecule prefactor divides by NA^2.
candidate=run('molecule_units',base.replace(' MOLES ',' MOLECULES ').replace('1.0D10',format(1e10/(6.02214076e23**2),'.16g')))
assert all(math.isclose(a,b,rel_tol=1e-9) for a,b in zip(rates(reference),rates(candidate)))
candidate=run('bare_site_mapping',base.replace('Ni(s)','*').replace('H(s)','H(*)'),thermo.replace('Ni(s)','*').replace('H(s)','H(*)'),options=['-siteSpecies','Ni(s)'])
assert rates(candidate)==rates(reference)
config=defaults+'''\nspeciesMap
{
    "Ni(s)" "Ni_site(s)";
}
surfaceThermoOverrides
{
    "Ni_site(s)" { molWeight 58.6934; rho 9000; }
    "H(s)" { Cp 123; Hf -456; }
}
'''
candidate=run('explicit_mapping_and_overrides',config=config)
result=(candidate/'thermo.solid').read_text()
assert 'Ni_site(s)' in result and 'rho 9000;' in result and 'Cp 123;' in result and 'Hf -456;' in result
run('dry_run',options=['-checkOnly'],contains='no output files written')
run('legacy_newFormat',options=['-newFormat'],contains='compatibility')
run('bad_number',base.replace('1.0D10','not_a_number'),success=False,contains='invalid number')
run('missing_thermo',therm=thermo.replace('Ni(s)','Other(s)'),success=False,contains='Missing thermodynamic record')
run('unknown_surface_species',base.replace('=> 2H(s)','=> 2X(s)'),success=False,contains='not declared in SITE')
run('bad_site_balance',base.replace('=> 2H(s)','=> H(s)'),success=False,contains='does not conserve')
run('unsupported_reversible',base.replace('=>','<=>'),success=False,contains='irreversible')
run('unsupported_FORD',base.replace('\nEND\n','\nEND\n',1).replace('8.31446261815324\n','8.31446261815324\nFORD/H2 0/\n'),success=False)
run('unsupported_units',base.replace('KJOULES/MOLE','BOGUS'),success=False,contains='unsupported REACTIONS option')
run('truncated_NASA',therm=thermo.splitlines()[0]+'\n300 1000 3000\nNi(s)\n',success=False,contains='truncated')
run('mapping_collision',config=defaults+'\nspeciesMap { "H(s)" "Ni(s)"; }\n',success=False,contains='mapping-colliding')
run('unknown_override',config=defaults+'\nsurfaceThermoOverrides { "X(s)" { Cp 10; } }\n',success=False,contains='not in the mapped SITE')
run('missing_reactions_END',mechanism=base.rsplit('END',1)[0],success=False,contains='missing END')
run('site_occupancy',mechanism=base.replace('Ni(s) H(s)','Ni(s)/2/ H(s)'),success=False,contains='multi-site occupancy')
run('negative_property',config=defaults+'\nsurfaceThermoOverrides { "H(s)" { rho -1; } }\n',success=False,contains='Invalid constant')
# Failed validation must leave pre-existing outputs untouched.
d=work/'failed_preserves_outputs';d.mkdir()
for f,s in [('surf.inp',base.replace('1.0D10','BAD')),('therm.dat',thermo),('defaults',defaults),('reactions.surf','old reactions'),('thermo.solid','old thermo')]: (d/f).write_text(s)
p=subprocess.run(['surfChemkinToFoam-new',*[str(d/f) for f in ('surf.inp','therm.dat','defaults','reactions.surf','thermo.solid')]],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
assert p.returncode and (d/'reactions.surf').read_text()=='old reactions' and (d/'thermo.solid').read_text()=='old thermo'
passed.append('failed_preserves_outputs')
# If installation of the second output is impossible, restore the first one.
d=work/'output_install_rollback';d.mkdir()
for f,s in [('surf.inp',base),('therm.dat',thermo),('defaults',defaults),('reactions.surf','previous reactions')]: (d/f).write_text(s)
(d/'thermo.solid').mkdir()
p=subprocess.run(['surfChemkinToFoam-new',*[str(d/f) for f in ('surf.inp','therm.dat','defaults','reactions.surf','thermo.solid')]],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
assert p.returncode and (d/'reactions.surf').read_text()=='previous reactions' and (d/'thermo.solid').is_dir()
assert not list(d.glob('*.tmp.*')) and not list(d.glob('*.backup.*'))
passed.append('output_install_rollback')
# Input/output alias protection, including symlinks.
d=work/'input_alias';d.mkdir()
for f,s in [('surf.inp',base),('therm.dat',thermo),('defaults',defaults)]: (d/f).write_text(s)
(d/'alias').symlink_to(d/'surf.inp')
p=subprocess.run(['surfChemkinToFoam-new',str(d/'surf.inp'),str(d/'therm.dat'),str(d/'defaults'),str(d/'alias'),str(d/'thermo.solid')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
assert p.returncode and 'aliases an input' in p.stdout and (d/'surf.inp').read_text()==base
passed.append('input_alias_protection')
# OpenFOAM itself must be able to read the generated dictionaries.
for f,entry in [('reactions.surf','reactions'),('thermo.solid','species')]:
    p=subprocess.run(['foamDictionary',str(reference/f),'-entry',entry,'-value'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    assert not p.returncode,p.stdout
passed.append('OpenFOAM_dictionary_read')
(work/'results.json').write_text(json.dumps(passed,indent=2))
print('Feature checks passed:',len(passed))
print('\n'.join(passed))
print('Feature outputs:',work)
