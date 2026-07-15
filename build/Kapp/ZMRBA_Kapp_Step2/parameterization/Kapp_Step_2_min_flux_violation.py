path_gams = '../GAMS/'
path_rxns_list = '../GAMS/RBA_rxns.txt'
path_out = './'
path_data = './Protein_ZM.xlsx'
path_enz_level = './enz_flux_calculation.txt'


pycore_path = '../pycore/'

#### Create directory 
import os,shutil


#### Load proteomics data and write protein translation fluxes
# Load data
import pandas as pd
df_data = pd.read_excel(path_data, sheet_name='RBA_proteins')
df_data.index = df_data['Gene'].to_list()
df_data = df_data[df_data['g/gDCW'] > 0]
# Excluding ribosome protein subunit (conflicting if fit to both enzymatic and ribosomal protein data)
df_data = df_data[(df_data.type == 'truedata_enz')]

data_val = []; data_idx = [];
for i in df_data.index:
        data_idx.append("'PROSYN-" + df_data.Gene[i] + "'")
        data_val.append("'PROSYN-" + df_data.Gene[i] + "' " + str(df_data.loc[i, 'vtrans (mmol/gDW/h)']))

data_val = ['/'] + data_val + ['/']
data_idx = ['/'] + data_idx + ['/']
with open(os.path.join(path_out, 'proteome_data.txt'), 'w') as f:
    f.write('\n'.join(data_val))
with open(os.path.join(path_out, 'proteome_data_set.txt'), 'w') as f:
    f.write('\n'.join(data_idx))

#### Determine active and inactive reactions informed by proteomics data and calculation in Kapp_Step1.py
# List out model reactions
with open(path_rxns_list) as f:
    idx = f.read().split('\n')[1:-1]
idx = [i[1:-1] for i in idx]
rxns_all = [i for i in idx if i.split('-')[0] == 'RXN']
rxns_enz_all = [i for i in idx if i.split('-')[0] == 'ENZLOAD']
rxns_enz_all = ['RXN-'+i[8:] for i in rxns_enz_all]
rxns_nonenz = [i for i in rxns_all if i not in rxns_enz_all]

# List active enzymatic reactions
with open(path_enz_level) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']

rxns_enz_active = []
for i in fluxes:
    enzid,v = i.split('\t')
    etype = enzid.split('-')[0]
    if etype == 'ENZLOAD':
        if float(v) > 0:
            rxns_enz_active.append('RXN-'+enzid[8:])

# List of inactive enzymatic reactions
rxns_inactive = [i for i in rxns_all if i not in rxns_enz_active + rxns_nonenz]
rxns_inactive = ["'" + i + "'" for i in rxns_inactive]
rxns_inactive = ['/'] + rxns_inactive + ['/']
fname = os.path.join(path_out, 'rxns_inactive.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(rxns_inactive))

#### Simulation
import shutil
shutil.copy(os.path.join(path_gams, 'min_flux_violation.gms'),
            os.path.join(path_out, 'min_flux_violation.gms'));
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(path_out, 'soplex.opt'));
shutil.copy(os.path.join(path_gams, 'min_flux_violation_GAMS_settings.txt'),
            os.path.join(path_out, 'min_flux_violation_GAMS_settings.txt'));

# Get the absolute path to the output directory
abs_path_out = os.path.abspath(path_out)

# Construct the commands dynamically
cmds = str('cd "' + abs_path_out + '"'), str('gams min_flux_violation.gms')
os.system('\n'.join(cmds))

#### Convert GAMS-scaled flux to actual flux
# All fluxes
fname = os.path.join(path_out, 'min_flux_violation.flux_gamsscaled.txt')
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
for i in fluxes:
    r,vtype,val = i.split('\t')
    fluxes_new.append('\t'.join([r, vtype, str(float(val) / 1e3)]))
fname = os.path.join(path_out, 'min_flux_violation.flux.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
    
# Only flux of reactions not associated with expressed enzymes
fname = os.path.join(path_out, 'min_flux_violation.flux_essential_inactive_rxns_gamsscaled.txt')
if os.path.exists(fname):
    with open(fname) as f:
        fluxes = f.read().split('\n')
    fluxes = [i for i in fluxes if i != '']
    fluxes_new = []
    for i in fluxes:
        r,vtype,val = i.split('\t')
        fluxes_new.append('\t'.join([r, vtype, str(float(val) / 1e3)]))
    fname = os.path.join(path_out, 'min_flux_violation.flux_essential_inactive_rxns.txt')
    with open(fname, 'w') as f:
        f.write('\n'.join(fluxes_new))


