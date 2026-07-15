import pandas as pd
import numpy as np
import sys
sys.path.append('../pycore/')
from simulate import RBA_result
from utils import extract_details_from_rxnid

# Load enzyme info
df_enz = pd.read_excel('../input/Enzyme_ZM.xlsx')


biom_id = 'BIOSYN-BIODIL'

res_metab = RBA_result(biom_id=biom_id)
res_metab.load_raw_flux('../input/min_flux_sum.flux.RNF.txt')
#res_metab.load_raw_flux('../input/min_flux_sum.flux.FNOR.FLDR.txt')
res_metab.calculate_metabolic_flux()

res_esyn = RBA_result(biom_id=biom_id, twocol_format=True)
res_esyn.load_raw_flux('../input/enz_flux_calculation.txt')



mu = res_metab.growth_rate
print('Growth rate:', mu)


#### Map rxn to enz
rxndict = {k:[] for k,v in res_metab.metabolic_flux.items() if abs(v) > 0}

enzdict = dict()
for k,v in res_esyn.raw_flux.items():
    if k.split('-')[0] == 'ENZSYN' and v > 0:
        _,enz = k.split('-', maxsplit=1)
        enzdict[enz] = []
        
for i in df_enz.index:
    rxn = df_enz.rxn_src[i]
    enz = df_enz.enz[i]
    if rxn in rxndict.keys() and enz in ['SPONT', 'UNKNOWN']:
        rxndict[rxn].append('zeroCost')
        
    if rxn in rxndict.keys() and enz in enzdict.keys():
        rxndict[rxn].append(enz)
        enzdict[enz].append(rxn)
        
rxndict = {k:set(v) for k,v in rxndict.items()}
rxndict = {k:v for k,v in rxndict.items() if v != {'zeroCost'}}
enzdict = {k:set(v) for k,v in enzdict.items()}


# Find enzymes whose together carry a total load of reactions
x = {k:v for k,v in rxndict.items() if len(v) > 1.5}
enz_share_rxn_load = set().union(*[v for v in x.values()])

# Find enzymes whose individually carry loads of multiple reactions
x = {k:v for k,v in enzdict.items() if len(v) > 1.5}
enz_multiload = set(x.keys())

# Set 1: Enzyme-reaction one-to-one load mapping
set1 = set(enzdict.keys()) - enz_share_rxn_load - enz_multiload
set1 = set([i for i in set1 if len(enzdict[i]) > 0.5])

# Set 2: Enzyme-reaction one-to-many
set2 = enz_multiload - enz_share_rxn_load

# Set 3: Enzyme-reaction many-to-one
set3 = enz_share_rxn_load - enz_multiload

kapp = dict()

### Set 1: Enzyme-reaction one-to-one load mapping
### Note that when calculating the Kapp mulitplying by growth rate (basically dividing enzyme synthetis flux by mu) to get enzyme concentration in the denominator
for enz in set1:
    rxn = [i for i in enzdict[enz]][0]
    rval = res_metab.metabolic_flux[rxn]
    if rval > 0:
        rdir = 'FWD'
    elif rval < 0:
        rdir = 'REV'
    else:
        print('rval == 0, check enzyme ' + enz + ' and reaction ' + rxn)
    
    rid = 'RXN-' + rxn + '_' + rdir + '-' + enz
    
    enzval = res_esyn.raw_flux['ENZSYN-' + enz]
    kapp[rid] = mu * abs(rval) / enzval / 3600


### Set 2: Enzyme-reaction one-to-many
for enz in set2:
    rids = []; rvalsum = 0;
    for rxn in enzdict[enz]:
        rval = res_metab.metabolic_flux[rxn]
        if rval > 0:
            rdir = 'FWD'
        elif rval < 0:
            rdir = 'REV'
        else:
            print('rval == 0, check enzyme ' + enz + ' and reaction ' + rxn)

        rids.append('RXN-' + rxn + '_' + rdir + '-' + enz)
        rvalsum += abs(rval)
        
    enzval = res_esyn.raw_flux['ENZSYN-' + enz]
    
    for rid in rids:
        kapp[rid] = mu * rvalsum / enzval / 3600
        
### Set 3: Enzyme-reaction many-to-one
rxns = set().union(*[enzdict[enz] for enz in set3])
for rxn in rxns:
    rval = res_metab.metabolic_flux[rxn]
    if rval > 0:
        rdir = 'FWD'
    elif rval < 0:
        rdir = 'REV'
    else:
        # Avoid error in loop if flux is zero
        print('rval == 0, check reaction ' + rxn)
        continue # Skip this reaction
        
    rids = []; enzvals = [];
    for enz in rxndict[rxn]:
        # Only include enzymes that are actually synthesized
        if 'ENZSYN-' + enz in res_esyn.raw_flux:
            rids.append('RXN-' + rxn + '_' + rdir + '-' + enz)
            enzvals.append(res_esyn.raw_flux['ENZSYN-' + enz])

    # Sum the synthesis fluxes of all contributing isoenzymes
    enzval = sum(enzvals)

    # Avoid division by zero if no active enzyme was found
    if enzval > 1e-12: 
        for rid in rids:
            kapp[rid] = mu * abs(rval) / enzval / 3600
    else:
        print(f"Warning: Zero total enzyme synthesis for active reaction {rxn}. Skipping kapp calculation.")


        

# Setting the calcualted kpp for RNF , ispG/H equal for both fld and fdx versions (in pfba model only chooses one)
iso_pairs_to_sync = [
    ('RXN-RNFfld_FWD-ZMOenz44',      'RXN-RNFfdx_FWD-ZMOenz44'),
    ('RXN-DMPPS_fld_FWD-ZMO0875',    'RXN-DMPPS_fdx_FWD-ZMO0875'),
    ('RXN-IPDPS_fld_FWD-ZMO0875',    'RXN-IPDPS_fdx_FWD-ZMO0875'),
    ('RXN-MECDPDH5_fld_FWD-ZMO0180', 'RXN-MECDPDH5_fdx_FWD-ZMO0180')
]

# Loop through each pair and sync the Kapp value from the active one to the inactive one
for r_fld, r_fdx in iso_pairs_to_sync:
    if r_fld in kapp and r_fdx not in kapp:
        # fld is active, fdx is inactive
        kapp[r_fdx] = kapp[r_fld]
    elif r_fdx in kapp and r_fld not in kapp:
        # fdx is active, fld is inactive
        kapp[r_fld] = kapp[r_fdx]



#Calculate the median kapp from the valid, active reactions
if kapp.values():
    median_kapp = np.median(list(kapp.values()))
    print(f"Calculated Median kapp: {median_kapp}")
else:
    median_kapp = 0
    print("No kapp values were calculated. Cannot impute.")


       

texts = ['rxnid\tkapp (1/s)']
for k,v in kapp.items():
    texts.append(k + '\t' + str(v))

with open('./kapps_in_vivo_RNF.txt', 'w') as f:
    f.write('\n'.join(texts))

# with open('./kapps_in_vivo_FNOR_FLDR.txt', 'w') as f:
#     f.write('\n'.join(texts))

        



