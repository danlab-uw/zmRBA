path_gams = '../GAMS/'
path_out = './'
path_data = './Protein_ZM.xlsx'

pycore_path = '../pycore/'

#### Create directory and copy run settings
import os,shutil


#### Load data
import pandas as pd
df_all = pd.read_excel(path_data, sheet_name= 'RBA_proteins')
df_all.index = df_all['Gene'].to_list()

# Create a list of the ribosomal genes so we can exclude them entirely later (conflicting if fit to both enyzmatic and ribosomal protein data)
ribo_genes = df_all[df_all['type'] == 'truedata_ribo']['Gene'].to_list()

# filtering down to just the enzymatic proteins with valid data
df_data = df_all[df_all['g/gDCW'] > 0]
df_data = df_data[(df_data.type == 'truedata_enz')]


#### Process data
import sys
sys.path.append(pycore_path)
from utils import metabolites_dict_from_reaction_equation_RBA

with open(os.path.join(path_gams, 'pro_and_enz.txt')) as f:
    pro_list = f.read().split('\n')
pro_list = pro_list[1:-1]
pro_list = [i[1:-1] for i in pro_list]
pro_list = [i for i in pro_list if i.split('-')[0] == 'PRO']

data = []; pro_data = []; pro_nodata = []
for met in pro_list:
    _,sid = met.split('-', maxsplit=1)

    # If protein ribosomal skip it
    if sid in ribo_genes:
        continue
    #there is one ribosome included in model to make up full 50s that we had no concentration data for so making sure this doesn't get include in pro no data file
    if sid == 'ZMO1246':
        continue
    # checl for enzymatic proteins
    if sid in df_data.index:
        pro_data.append("'PROSYN-" + sid + "'")
        data.append("'PROSYN-" + sid + "' " + str(df_data.loc[sid, 'vtrans (mmol/gDW/h)']))
    else:
        # list of enzymatic proteins with no concentration data
        pro_nodata.append("'PROSYN-" + sid + "'")

data = ['/'] + data + ['/']
pro_data = ['/'] + pro_data + ['/']
pro_nodata = ['/'] + pro_nodata + ['/']

# Write out run files
with open(os.path.join(path_out, 'proteome_data.txt'), 'w') as f:
    f.write('\n'.join(data))
with open(os.path.join(path_out, 'rxns_pro_data.txt'), 'w') as f:
    f.write('\n'.join(pro_data))
with open(os.path.join(path_out, 'rxns_pro_nodata.txt'), 'w') as f:
    f.write('\n'.join(pro_nodata))

#### Simulation
shutil.copy(os.path.join(path_gams, 'new2.gms'),
            os.path.join(path_out, 'new2.gms'))
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(path_out, 'soplex.opt'))
shutil.copy(os.path.join(path_gams, 'enz_from_proteome_GAMS_settings.txt'),
            os.path.join(path_out, 'enz_from_proteome_GAMS_settings.txt' ))

from simulate import get_GAMS_modelStat
from collections import OrderedDict

with open(os.path.join(path_gams, 'rxns_enz.txt')) as f:
    enz_list = f.read().split('\n')
enz_list = enz_list[1:-1]
enz_list = [i[1:-1] for i in enz_list]

enzdict = OrderedDict()

for enz in enz_list:
    #os.system('gams new2.gms')
    #print (enz)
    abs_path_out = os.path.abspath(path_out)
    cmds = str('cd "' + abs_path_out + '"'), str('gams new2.gms --enzobj=') + enz
    print(cmds)
    os.system('\n'.join(cmds))
    fname = os.path.join(path_out,'enz_from_proteome.modelStat.txt')
    stat = get_GAMS_modelStat(fname)

    if stat == 'optimal':
        fname = os.path.join(path_out, 'enz_from_proteome.objval.txt')
        with open(fname) as f:
            v = float(f.read().split('\n')[0]) / 1e6
        enzdict[enz] = v
    else:
        enzdict[enz] = 0
        

enztext = [k+'\t'+str(v) for k,v in enzdict.items()]
with open(os.path.join(path_out, 'enz_flux_calculation.txt'), 'w') as f:
    f.write('\n'.join(enztext))
    
    #enzobj=enz
    #cmds = ['cd ' + 'gams new2.gms']
    #os.system('n/'.join(cmds))
            #'gams enz_from_proteome.gms --enzobj=' + enz]
    #print (cmds)
   
    
