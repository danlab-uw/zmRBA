import gams.transfer as gt
import pandas as pd
import re


#SETUP & LOAD GDX DATA
Enzyme_cap_file = '../GAMS/inputs/RBA_enzCapacityConstraints_eqns_equality_version.txt'
gdx_file = 'model_results.gdx'

m = gt.Container()
m.read(gdx_file)


# EXTRACT FLUXES (Variable 'v')
df_fluxes = m.data['v'].records

# Filter to see only active fluxes (level > 0)
active_fluxes = df_fluxes[df_fluxes['level'] > 0]


# Descale flux levels scaled them in GAMS
df_fluxes['level_descaled'] = df_fluxes['level'] * 1e-3 


#EXTRACT METABOLITE MARGINALS (Equation 'Stoic')
df_stoic = m.data['Stoic'].records

#Sort by the most negative marginals
starving_metabolites = df_stoic.sort_values(by='marginal').head()
print(starving_metabolites[['uni', 'level', 'marginal']])


# REDOX CONSTRAINTS
redox_cap_list = []

# FLD RED
df_fld_red_cap = m.data['Eq_CarrierCap_fld_red'].records
bottlenecks_fld_red = df_fld_red_cap.sort_values(by='marginal').head().copy()
bottlenecks_fld_red['Carrier'] = 'fld_red'  
redox_cap_list.append(bottlenecks_fld_red)

# FLD OX
df_fld_ox_cap = m.data['Eq_CarrierCap_fld_ox'].records
bottlenecks_fld_ox = df_fld_ox_cap.sort_values(by='marginal').head().copy()
bottlenecks_fld_ox['Carrier'] = 'fld_ox'
redox_cap_list.append(bottlenecks_fld_ox)

# FDX RED
df_fdx_red_cap = m.data['Eq_CarrierCap_fdx_red'].records
bottlenecks_fdx_red = df_fdx_red_cap.sort_values(by='marginal').head().copy()
bottlenecks_fdx_red['Carrier'] = 'fdx_red'
redox_cap_list.append(bottlenecks_fdx_red)

# FDX OX
df_fdx_ox_cap = m.data['Eq_CarrierCap_fdx_ox'].records
bottlenecks_fdx_ox = df_fdx_ox_cap.sort_values(by='marginal').head().copy()
bottlenecks_fdx_ox['Carrier'] = 'fdx_ox'
redox_cap_list.append(bottlenecks_fdx_ox)


df_all_redox = pd.concat(redox_cap_list, ignore_index=True)


cols = ['Carrier'] + [c for c in df_all_redox.columns if c != 'Carrier']
df_all_redox = df_all_redox[cols]



# 5. EXTRACT ALL ENZYME CAPACITY MARGINALS
all_capacity_data = []

for symbol_name in m.data.keys():
    if symbol_name.startswith('EnzCap'):
        df_temp = m.data[symbol_name].records
        df_temp['Constraint_Name'] = symbol_name
        all_capacity_data.append(df_temp)

# Combine all the individual enzyme dataframes into one master table
if all_capacity_data:
    df_all_enzymes = pd.concat(all_capacity_data, ignore_index=True)
else:
    df_all_enzymes = pd.DataFrame()



#PARSE ENZYME NAMES FROM TEXT FILE
parsed_data = []

with open(Enzyme_cap_file, 'r') as f:
    for line in f:
        clean_line = line.strip()
        if not clean_line:
            continue
            
        # Extract the constraint name (EnzCap...) and the reaction ID
        eq_match = re.search(r'^([^.]+)\.\.', clean_line)
        rxn_match = re.search(r'\.\.(.*?)=e=', clean_line)

        if eq_match and rxn_match:
            eq_name = eq_match.group(1).strip()
            rxn_id = rxn_match.group(1).strip()
            
            parsed_data.append({
                'Constraint_Name': eq_name,
                'Reaction_Enzyme_ID': rxn_id
            })

enzymecap_info = pd.DataFrame(parsed_data)



#MERGE MARGINALS WITH ENZYME NAMES
if not df_all_enzymes.empty and not enzymecap_info.empty:
    
    # Left join to attach the text file info to the GAMS data
    final_enzyme_data = pd.merge(df_all_enzymes, enzymecap_info, on='Constraint_Name', how='left')

    # Reorder columns for Excel
    columns_to_keep = ['Constraint_Name', 'Reaction_Enzyme_ID', 'level', 'marginal']
    other_cols = [col for col in final_enzyme_data.columns if col not in columns_to_keep]
    final_enzyme_data = final_enzyme_data[columns_to_keep + other_cols]

    final_enzyme_data = final_enzyme_data.sort_values(by='Constraint_Name', ascending=True)




#EXPORT TO EXCEL
with pd.ExcelWriter('RBA_Marginal_Results.xlsx') as writer:
    df_fluxes.to_excel(writer, sheet_name='Fluxes', index=False)
    df_stoic.to_excel(writer, sheet_name='Metabolite_Marginals', index=False)
    df_all_redox.to_excel(writer, sheet_name='Redox_Caps', index=False)
    final_enzyme_data.to_excel(writer, sheet_name='Enzyme_Caps', index=False)


