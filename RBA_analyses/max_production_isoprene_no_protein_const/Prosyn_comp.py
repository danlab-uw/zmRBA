import pandas as pd
import numpy as np


df_pro_path = './Protein_ZM.xlsx'
df_model_flux_path = './max_prod/runRBA.flux.txt' 


df_exp = pd.read_excel(df_pro_path, sheet_name='RBA_proteins')
df_exp = df_exp.rename(columns={'vtrans (mmol/gDW/h)': 'Exp_PROSYN'})


#flux data from model
fluxes = []
with open(df_model_flux_path, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            fluxes.append({'Reaction': parts[0], 'Flux': float(parts[2])})

df_flux = pd.DataFrame(fluxes)


#Extract Protein Synthesis
df_prosyn = df_flux[df_flux['Reaction'].str.startswith('PROSYN-')].copy()
df_prosyn['Gene'] = df_prosyn['Reaction'].str.replace('^PROSYN-', '', regex=True)
df_prosyn = df_prosyn.rename(columns={'Flux': 'Model_PROSYN'})

#Extract Protein Waste
df_prowaste = df_flux[df_flux['Reaction'].str.startswith('PROWASTE-')].copy()
df_prowaste['Gene'] = df_prowaste['Reaction'].str.replace('^PROWASTE-', '', regex=True)
df_prowaste = df_prowaste.rename(columns={'Flux': 'Model_PROWASTE'})


#merge experimental + model synthesis
analysis_df = pd.merge(df_exp[['Gene', 'Exp_PROSYN']], df_prosyn[['Gene', 'Model_PROSYN']], on='Gene', how='inner')

#model waster - left join because some genes might have 0 waste and no fluc
analysis_df = pd.merge(analysis_df, df_prowaste[['Gene', 'Model_PROWASTE']], on='Gene', how='left')
analysis_df['Model_PROWASTE'] = analysis_df['Model_PROWASTE'].fillna(0.0)



analysis_df['Active_Protein'] = analysis_df['Model_PROSYN'] - analysis_df['Model_PROWASTE']


#calculate utilization % 
analysis_df['Utilization_%'] = np.where(
    analysis_df['Model_PROSYN'] > 0, 
    (analysis_df['Active_Protein'] / analysis_df['Model_PROSYN']) * 100, 
    0.0
)

analysis_df['Model Use vs Exp'] = np.where(
    analysis_df['Model_PROSYN'] > 0, 
    analysis_df['Model_PROSYN']/analysis_df['Exp_PROSYN'],
    0.0
)
 

analysis_df = analysis_df.sort_values('Active_Protein', ascending=False)
analysis_df.to_excel('RBA_Protein_Syn_Analysis.xlsx', index=False)


