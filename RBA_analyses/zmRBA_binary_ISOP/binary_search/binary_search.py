import os
import shutil
import sys
import time 
import json
import pandas as pd
from binary_search_options import *



time0 = time.time()
base_dir = os.getcwd()

# Path definitions
path_pycore = os.path.abspath('../pycore/')
if path_pycore not in sys.path:
    sys.path.append(path_pycore)
from simulate import get_GAMS_modelStat, RBA_result

path_gams = os.path.abspath('../GAMS/application/')
path_enz_mw = os.path.abspath('../GAMS/inputs/enz_mw_g_per_mmol.txt')
path_pro_mw = os.path.abspath('../GAMS/inputs/pro_mw_g_per_mmol.txt')
path_out = os.path.abspath('./')

# Copy the GAMS settings to the base directory 
shutil.copy(os.path.join(path_gams, 'runRBA_binary_search.gms'), os.path.join(base_dir, 'runRBA_binary_search.gms'))
shutil.copy(os.path.join(path_gams, 'soplex.opt'), os.path.join(base_dir, 'soplex.opt'))
shutil.copy(os.path.join(path_gams, 'runRBA_GAMS_settings.txt'), os.path.join(base_dir, 'runRBA_GAMS_settings.txt')) 

# Reset variables
mu_min = mu_min0
mu_max = mu_max0
itercount = 0
vglc = 55.0

if mu_min0 < 0:
    print('Invalid negative value is set as the lower bound of mu. Mu is always positive. Revert to zero.')
    mu_min = 0

# Create directory for current vglc
folder_name = f"vglc_{vglc}"
abs_path_out = os.path.join(path_out, folder_name)
os.makedirs(abs_path_out, exist_ok=True)

# Define report file specifically mapped to the new folder
report_file = os.path.join(abs_path_out, 'run_log.txt')

# Initiate report
report = {k:None for k in ['stat', 'vglc']}
report['vglc'] = vglc

class bcolors:
    GREEN = '\033[92m' 
    RED = '\033[91m' 
    RESET = '\033[0m' 


# Test evaluation at zero
mu = 0.0  
gams_cmd = f'module load gams; gams runRBA_binary_search.gms --mu={mu} --vglc={vglc} o=/dev/null'
os.system(gams_cmd)

stat = get_GAMS_modelStat('./runRBA.modelStat.txt')

if stat == 'infeasible':
    text = f'Model is infeasible at mu = 0 for vglc = {vglc}, check model connectivity. Terminating.'
    print(f"{bcolors.RED}{text}{bcolors.RESET}")
    with open(report_file, 'w') as f:
        f.write(text + '\n')
    quit()
    
elif stat == 'optimal':
    with open(report_file, 'w') as f:
        f.write('Pass mu = 0 test, now proceed to binary search.\n')


# Start binary search
# Evaluate min feasibility
mu = mu_min
gams_cmd = f'module load gams; gams runRBA_binary_search.gms --mu={mu} --vglc={vglc} o=/dev/null'
os.system(gams_cmd)
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')

if stat == 'need_rerun':
    while stat == 'need_rerun':
        print(f"mu = {mu:.7f}, status = {stat}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")

        mu += mu_tol
        gams_cmd = f'module load gams; gams runRBA_binary_search.gms --mu={mu} --vglc={vglc} o=/dev/null'
        os.system(gams_cmd)
        stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
        
    if stat == 'optimal':
        mu_min = mu
        print(f"{bcolors.GREEN}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")
            
if stat == 'infeasible':
    mu_max = mu_min
    mu_min = 0
    print(f"{bcolors.RED}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
    with open(report_file, 'a') as f:
        f.write(f"mu = {mu:.7f}, status = {stat}\n")

# Evaluate max infeasibility
mu = mu_max
gams_cmd = f'module load gams; gams runRBA_binary_search.gms --mu={mu} --vglc={vglc} o=/dev/null'
os.system(gams_cmd)
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')

while stat == 'optimal':
    print(f"{bcolors.GREEN}mu = {mu:.7f}, status = {stat}{bcolors.RESET}\n")
    with open(report_file, 'a') as f:
        f.write(f"mu = {mu:.7f}, status = {stat}\n")
    mu_max = 1.5 * mu_max
    mu = mu_max
    gams_cmd = f'module load gams; gams runRBA_binary_search.gms --mu={mu} --vglc={vglc} o=/dev/null'
    os.system(gams_cmd)
    stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
# Update min-max
mu = float(mu_min + mu_max) / 2
final_res = RBA_result(biom_id=biom_id)

while mu_max - mu_min > mu_tol and itercount < maxiter:
    itercount += 1
        
    gams_cmd = f'module load gams; gams runRBA_binary_search.gms --mu={mu} --vglc={vglc} o=/dev/null'
    os.system(gams_cmd)
    stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
    if stat == 'need_rerun':
        while stat == 'need_rerun':
            print(f"mu = {mu:.7f}, status = {stat}")
            with open(report_file, 'a') as f:
                f.write(f"mu = {mu:.7f}, status = {stat}\n")
            mu += mu_tol
            gams_cmd = f'module load gams; gams runRBA_binary_search.gms --mu={mu} --vglc={vglc} o=/dev/null'
            os.system(gams_cmd)
            stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
    if stat == 'optimal':
        mu_min = mu
        final_res.load_raw_flux('./runRBA.flux.txt')
        print(f"{bcolors.GREEN}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")
        
    elif stat == 'infeasible':
        mu_max = mu
        print(f"{bcolors.RED}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")
            
    else:
        print("Error unknown!")
        quit()
        
    mu = float(mu_min + mu_max) / 2

# Final check at mu_min
mu = mu_min
gams_cmd = f'module load gams; gams runRBA_binary_search.gms --mu={mu} --vglc={vglc} o=/dev/null'
os.system(gams_cmd)
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
final_res.load_raw_flux('./runRBA.flux.txt')


# Write results
report['stat'] = stat
if stat == 'optimal':
    # Load Enzyme MW 
    enz_mw = dict()
    with open(path_enz_mw) as f:
        for line in f.read().splitlines():
            parts = line.split() 
            if len(parts) >= 2:
                enz_mw[parts[0].strip()] = float(parts[1].strip())
        
    # Load Protein MW 
    pro_mw = dict()
    with open(path_pro_mw) as f:
        for line in f.read().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                pro_mw[parts[0].strip()] = float(parts[1].strip())
        
    final_res.enzyme_mw = enz_mw
    final_res.protein_mw = pro_mw

    final_res.calculate_all()
    
    final_res.enzyme_mw = None
    final_res.protein_mw = None

    report['proteome_capacity'] = final_res.proteome_capacity_usage
    report['ribo_capacity'] = final_res.ribo_capacity_usage
    report['growth_rate'] = mu

    # Write report text file directly into the specific vglc folder
    text = []
    for k, v in report.items(): 
        text.append(f"{k}\t{v}")

    with open(os.path.join(abs_path_out, 'report.txt'), 'w') as f:
        f.write('\n'.join(text))

    # Save JSON directly into the specific vglc folder
    final_res.save_to_json(os.path.join(abs_path_out, 'RBA_result.json'))

dt = float(time.time() - time0) / 60
with open(report_file, 'a') as f:
    f.write(f"Run time = {dt:.2f} mins\n")
print(f"Run time = {dt:.2f} mins")

#copy the raw GAMS output texts into the specific folder
if os.path.exists('./runRBA.modelStat.txt'):
    shutil.copy('./runRBA.modelStat.txt', os.path.join(abs_path_out, 'runRBA.modelStat.txt'))
if os.path.exists('./runRBA.flux.txt'):
    shutil.copy('./runRBA.flux.txt', os.path.join(abs_path_out, 'runRBA.flux.txt'))