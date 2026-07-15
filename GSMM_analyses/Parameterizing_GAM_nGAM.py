import cobra
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np


model = cobra.io.read_sbml_model("i_ZM4_489.xml")


# removing old ATP values from the model
biomass_rxn = model.reactions.get_by_id("BIOMASS_ZM")
atp_maint_rxn = model.reactions.get_by_id("ATPM")
glc_exchange = model.reactions.get_by_id("EX_glc_e")
G6PDH2r = model.reactions.get_by_id('G6PDH2r')
G6PDH2r.bounds = (0,1000)
G6PDH2xr = model.reactions.get_by_id('G6PDH2xr')
G6PDH2xr.bounds = (0,1000)

rnf_fld = model.reactions.get_by_id('RNFfld')
rnf_fld.bounds = (0,1000)
rnf_fdx = model.reactions.get_by_id('RNFfdx')
rnf_fdx.bounds = (0,1000)

# Unlock Maintenance
atp_maint_rxn.bounds = (0, 1000) 


# Zero out old GAM coefficients in Biomass
met_atp = model.metabolites.get_by_id("atp_c")
met_adp = model.metabolites.get_by_id("adp_c")
met_pi  = model.metabolites.get_by_id("pi_c")
met_h2o = model.metabolites.get_by_id("h2o_c")
met_h   = model.metabolites.get_by_id("h_c")

for met in [met_atp, met_adp, met_pi, met_h2o, met_h]:
    if met in biomass_rxn.metabolites:
        old_coeff = biomass_rxn.metabolites[met]
        biomass_rxn.add_metabolites({met: -old_coeff}, combine=True)



# calculating GAM and NGAM using batch ZMM growth and uptake data 
# data from 2011 Widiastuti et al. - see 'Parameterizing GAM and NGAM in updated GSM.xlsx'
experimental_data = [
    (0.1574, 22.9960),
    (0.1869, 27.0433),
    (0.2254, 32.40),
    (0.2870, 40.39),
    (0.3880, 52.02)
    #(0.36, 41.9)              #this is a single data point from separate paper Jacobson et al. 2019
]

results_mu = []
results_qATP = []

print(f"\n{'Mu (1/h)':<10} | {'q_Glc_Exp':<10} | {'q_Glc_Used':<10} | {'Calc q_ATP'}")
print("-" * 55)

for mu_val, q_glc_exp in experimental_data:
    with model:
        # Constrain growth to experimental measurements
        biomass_rxn.bounds = (mu_val, mu_val)
        
        # Constrain glucose uptake based on experimental measurements
        #Way wrote it allows for slack if needed (i.e., allow model to take up to 5% more or less glucose if needed)
        # Experimental error is often 5-10%
        # ONLY allow for slack if doesn't solve for time
        upper_limit = -q_glc_exp * 1.0     #e.g., * 0.95
        lower_limit = -q_glc_exp * 1.0     # e.g., * 1.05
        
        glc_exchange.bounds = (lower_limit, upper_limit)
        
        # Solve model for optimizing ATP
        model.objective = atp_maint_rxn
        solution = model.optimize()
        
        if solution.status == 'optimal':
            q_atp = solution.objective_value
            # Get actual glucose used - if had to apply slack
            actual_glc = abs(solution.fluxes[glc_exchange.id])
            
            print(f"{mu_val:<10} | {q_glc_exp:<10} | {actual_glc:.3f}      | {q_atp:.3f}")
            results_mu.append(mu_val)
            results_qATP.append(q_atp)
        else:
            # If slack used isn't enough, just report failure
            print(f"{mu_val:<10} | {q_glc_exp:<10} | {'---':<10} | Infeasible (>10% dev)")


# Regression analysis and plot to determine GAM and NGAM (from Pirt 1965 and Pirt 1982)
if len(results_mu) > 2:
    slope, intercept, r_value, p_value, std_err = stats.linregress(results_mu, results_qATP)
    
    print("\n" + "="*30)
    print(f"FINAL PARAMETERS (N={len(results_mu)})")
    print("="*30)
    print(f"R-squared:            {r_value**2:.4f}")
    print(f"New GAM (Slope):      {slope:.4f} mmol ATP/gDW")
    print(f"New nGAM (Intercept): {intercept:.4f} mmol ATP/gDW/h")
    
    # Simple ASCII Plot check
    print("\nCheck: Is the nGAM positive? " + ("YES" if intercept > 0 else "NO (Warning!)"))

else:
    print("\nStill not enough points. Try increasing the slack .")