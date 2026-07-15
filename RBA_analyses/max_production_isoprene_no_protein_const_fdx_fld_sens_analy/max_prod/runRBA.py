import sys
import os
import shutil

path_pycore = os.path.abspath('../pycore/')
if path_pycore not in sys.path:
    sys.path.append(path_pycore)


from simulate import get_GAMS_modelStat, RBA_result


from runRBA_options import *
import json
import pandas as pd


path_gams = '../GAMS/application/'
path_enz_mw = '../GAMS/inputs/enz_mw_g_per_mmol.txt'
path_pro_mw = '../GAMS/inputs/pro_mw_g_per_mmol.txt'
path_out = './'


#### Create directory and copy run settings
abs_path_out = os.path.abspath(path_out)

shutil.copy(os.path.join(path_gams, 'runRBA_max_prod.gms'),
            os.path.join(abs_path_out, 'runRBA_max_prod.gms'))
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(abs_path_out, 'soplex.opt'))
shutil.copy(os.path.join(path_gams, 'runRBA_GAMS_settings.txt'),
            os.path.join(abs_path_out, 'runRBA_GAMS_settings.txt')) 

# Initiate report
report = {k:None for k in ['stat', 'vglc', 'vprod', 'yield']}


# Remove old GAMS output files
for stale_file in ['./runRBA.modelStat.txt', './runRBA.flux.txt']:
    if os.path.exists(stale_file):
        os.remove(stale_file)

# Execute GAMS (Initial Run)
vglc = vglc0
gams_cmd = f'module load gams; gams runRBA_max_prod.gms --mu={mu} --vglc={vglc} --vprod="{vprod}"'
os.system(gams_cmd)

stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
if stat == 'infeasible':
    report['stat'] = stat
    report['vglc'] = 0
    report['vprod'] = 0
    report['yield'] = 0
        
elif stat == 'optimal':
    res = RBA_result(biom_id=biom_id)
    res.load_raw_flux(filepath='./runRBA.flux.txt')
    res.calculate_metabolic_flux()
    if vprod.split('-')[0] == 'RXNADD':
        res.metabolic_flux[vprod_coreid] = res.raw_flux[vprod]
    pflux = res.metabolic_flux[vprod_coreid]
    vglc_sim = - res.metabolic_flux['EX_glc_e']
        
    report['stat'] = stat
    report['vglc'] = vglc_sim
    report['vprod'] = pflux
    report['yield'] = pflux * prod_mw / vglc_sim / 180.156
        
elif stat == 'need_rerun':
    itermax = 100
    iternum = 0
    while stat == 'need_rerun' and iternum < itermax:
        iternum += 1
        vglc += 1e-3
        
   
        rerun_cmd = f'module load gams; gams runRBA_max_prod.gms --mu={mu} --vglc={vglc} --vprod="{vprod}" o=/dev/null'
        os.system(rerun_cmd)
        
        stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
            
        if stat == 'infeasible':
            report['stat'] = stat
            report['vglc'] = 0
            report['vprod'] = 0
            report['yield'] = 0
            
        
        elif stat == 'optimal':
            res = RBA_result(biom_id=biom_id)
            res.load_raw_flux(filepath='./runRBA.flux.txt')
            res.calculate_metabolic_flux()
            if vprod.split('-')[0] == 'RXNADD':
                res.metabolic_flux[vprod_coreid] = res.raw_flux[vprod]
            pflux = res.metabolic_flux[vprod_coreid]
            vglc_sim = - res.metabolic_flux['EX_glc_e']
            
            report['stat'] = stat
            report['vglc'] = vglc_sim
            report['vprod'] = pflux
            report['yield'] = pflux * prod_mw / vglc_sim / 180.156
        
        elif stat == 'need_rerun':
            report['stat'] = stat
        
        else:
            print('wtf')



# Write JSON results
if stat == 'optimal':
    # Load Enzyme MW
    enz_mw = dict()
    with open(path_enz_mw) as f:
        for line in f.read().splitlines():
            parts = line.split()  # Splits on ANY whitespace (tabs or spaces)
            if len(parts) >= 2:
                enz_mw[parts[0].strip()] = float(parts[1].strip())
        
    # Load Protein MW sa
    pro_mw = dict()
    with open(path_pro_mw) as f:
        for line in f.read().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                pro_mw[parts[0].strip()] = float(parts[1].strip())
        
    res.enzyme_mw = enz_mw
    res.protein_mw = pro_mw

    #calculate all capacities
    res.calculate_all()
    

    res.enzyme_mw = None
    res.protein_mw = None

    report['proteome_capacity'] = res.proteome_capacity_usage
    report['ribo_capacity'] = res.ribo_capacity_usage
    report['growth_rate'] = res.growth_rate

  
    # Write report text file
    text = []
    for k, v in report.items(): 
        text.append(f"{k}\t{v}")

    with open(os.path.join(abs_path_out, 'report.txt'), 'w') as f:
        f.write('\n'.join(text))