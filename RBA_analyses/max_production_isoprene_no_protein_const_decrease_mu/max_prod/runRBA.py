import sys
import os
import shutil
import csv


path_pycore = os.path.abspath('../pycore/')
if path_pycore not in sys.path:
    sys.path.append(path_pycore)

from simulate import get_GAMS_modelStat, RBA_result


from runRBA_options import *

#Decreasing growth rate and seeing resultants isoprene productions
# Measured growth rate is 0.359 h-
mus = [0.359, 0.357, 0.355, 0.35, 0.345, 0.34, 0.335, 0.33, 0.32, 0.28, 0.24, 0.20, 0.16, 0.12, 0.08, 0.05]

out_csv = 'mu_sweep_results.csv'

path_gams = '../GAMS/application/'
path_out = './'

#Copy run settings
abs_path_out = os.path.abspath(path_out)
shutil.copy(os.path.join(path_gams, 'runRBA_max_prod.gms'),
            os.path.join(abs_path_out, 'runRBA_max_prod.gms'))
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(abs_path_out, 'soplex.opt'))
shutil.copy(os.path.join(path_gams, 'runRBA_GAMS_settings.txt'),
            os.path.join(abs_path_out, 'runRBA_GAMS_settings.txt')) 

def solve_at_mu(mu_val):
    vglc = vglc0

    #removing old GAMS out 
    for stale in ['./runRBA.modelStat.txt', './runRBA.flux.txt']:
        if os.path.exists(stale):
            os.remove(stale)

    cmd = (f'module load gams; gams runRBA_max_prod.gms '
           f'--mu={mu_val} --vglc={vglc} --vprod="{vprod}" o=/dev/null')
    os.system(cmd)
    
    stat = get_GAMS_modelStat('./runRBA.modelStat.txt')


    iternum = 0
    while stat == 'need_rerun' and iternum < 100:
        iternum += 1
        vglc += 1e-3
        os.system(f'module load gams; gams runRBA_max_prod.gms '
                  f'--mu={mu_val} --vglc={vglc} --vprod="{vprod}" o=/dev/null')
        stat = get_GAMS_modelStat('./runRBA.modelStat.txt')

    if stat != 'optimal':
        return {'mu': mu_val, 'stat': stat, 'vglc': 0,
                'vprod': 0, 'molar_yield_pct': 0, 'mass_yield': 0}

    res = RBA_result(biom_id=biom_id)
    res.load_raw_flux(filepath='./runRBA.flux.txt')
    res.calculate_metabolic_flux()
    if vprod.split('-')[0] == 'RXNADD':
        res.metabolic_flux[vprod_coreid] = res.raw_flux[vprod]
    pflux = res.metabolic_flux[vprod_coreid]
    vglc_sim = -res.metabolic_flux['EX_glc_e']

    return {
        'mu': mu_val,
        'stat': stat,
        'vglc': vglc_sim,
        'vprod': pflux,
        'molar_yield_pct': 100.0 * pflux / vglc_sim,          
        'mass_yield': pflux * prod_mw / vglc_sim / 180.156,    
    }


rows = []
for m in mus:
    r = solve_at_mu(m)
    rows.append(r)

with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

for r in rows:
    print(f"{r['mu']:>6} {r['vglc']:>8.2f} {r['vprod']:>10.4f} {r['molar_yield_pct']:>9.3f}")
