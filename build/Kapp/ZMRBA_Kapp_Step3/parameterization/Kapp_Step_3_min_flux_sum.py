path_gams = '../GAMS/'
path_out = './'


#### Create directory and copy run settings
import os,shutil


#### Load proteomics data and write protein translation fluxes
# Import from min_flux_violation outputs in GAMS

#### Append essential inactive reactions found in B2_min_flux_violation.py
fname = 'rxns_inactive.txt'
with open(fname) as f:
    rxns_inactive = f.read().split('\n')[1:-1]
rxns_inactive = [i[1:-1] for i in rxns_inactive]

fname = 'min_flux_violation.flux_essential_inactive_rxns.txt'
with open(fname) as f:
    rxns_essential_inactive = f.read().split('\n')
rxns_essential_inactive = [i.split('\t')[0] for i in rxns_essential_inactive if i != '']

rxns_inactive = [i for i in rxns_inactive if i not in rxns_essential_inactive]
rxns_inactive = ["'" + i + "'" for i in rxns_inactive]
rxns_inactive = ['/'] + rxns_inactive + ['/']
fname = os.path.join(path_out, 'rxns_inactive.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(rxns_inactive))

#### Simulation
import shutil
shutil.copy(os.path.join(path_gams, 'min_flux_sum.gms'),
            os.path.join(path_out, 'min_flux_sum.gms'));
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(path_out, 'soplex.opt'));
shutil.copy(os.path.join(path_gams, 'min_flux_sum_GAMS_settings.txt'),
            os.path.join(path_out, 'min_flux_sum_GAMS_settings.txt'));
 
 # Get the absolute path to the output directory
abs_path_out = os.path.abspath(path_out)

# Construct the commands dynamically
cmds = str('cd "' + abs_path_out + '"'), str('gams min_flux_sum.gms')
os.system('\n'.join(cmds))
           

#### Convert GAMS-scaled flux to actual flux
# All fluxes
fname = os.path.join(path_out, 'min_flux_sum.flux_gamsscaled.txt')
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
for i in fluxes:
    r,vtype,val = i.split('\t')
    fluxes_new.append('\t'.join([r, vtype, str(float(val) / 1e3)]))
fname = os.path.join(path_out, 'min_flux_sum.flux.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
