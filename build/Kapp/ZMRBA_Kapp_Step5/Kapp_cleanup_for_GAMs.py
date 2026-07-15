import pandas as pd
import numpy as np

import sys
sys.path.append('../pycore/')
from utils import extract_details_from_rxnid

#Load info
df_kapp = '../ZMRBA_Kapp_Step4/Kapp comparison and compilation.xlsx'
kapp = pd.read_excel(df_kapp, sheet_name='Kapp for RBA')
kapp_df = kapp
kapp = kapp.set_index('rxnid')['kapp (1/s)']

df_stoich = '../../ZMRBA_outputs/RBA_Stoichiometry.xlsx'
df_stoich = pd.read_excel(df_stoich)
df_stoich.index = df_stoich.id.to_list()

df_enz = '../../ZMRBA_inputs/Enzyme_ZM.xlsx'
df_enz = pd.read_excel(df_enz)
df_enz.index = df_enz.id.to_list()

outfile = './Kapp_ZM.txt' 

# Get kapp median from just unique enzymes
kapp_df['rxnid'] = kapp_df['rxnid'].str.rsplit('-', n=1).str[-1]
kapp_df = kapp_df.drop_duplicates(subset='rxnid', keep='first').copy()
med = round(kapp_df['kapp (1/s)'].median(), 1)
print('Kapp median: ', med)

# Map kapp to ENZLOAD entries
idx = df_enz.index.to_list()
df_write = pd.DataFrame(index=idx, columns=['id', 'kapp (1/s)', 'source'])
df_write['id'] = df_write.index.to_list()

for i in df_enz.index:
    _,rxn,rdir,enz = extract_details_from_rxnid(i)
    if enz in ['SPONT', 'UNKNOWN']:
        df_write.loc[i, 'source'] = enz
        continue
    
    entry = i

    if entry in kapp.index:
        df_write.loc[i, 'kapp (1/s)'] = kapp[entry]
        df_write.loc[i, 'source'] = 'parameterization'
    else:
        df_write.loc[i, 'kapp (1/s)'] = med
        df_write.loc[i, 'source'] = 'kapp_median'


        
# Write to GAMS txt file
idx = df_stoich[df_stoich.coupling_type == 'rxn_enz'].index
kapp_list = []

for i in idx:
    lhs = "v('ENZLOAD-" + df_stoich.id[i][4:] + "') * " + "kapp('" + i + "')"
    kapp_list.append("'" + i + "' " + str(round(df_write.loc[i, 'kapp (1/s)'] * 3600, 6)))
    
kapp_list = ['/'] + kapp_list + ['/']
with open(outfile, 'w') as f:
    f.write('\n'.join(kapp_list))





