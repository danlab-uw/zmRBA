#Trying cobrapy pfba and seeing what fluxes give
import cobra
from cobra.flux_analysis import pfba 


# Load the model
model = cobra.io.read_sbml_model('i_ZM4_489_isoprene.xml')

# Set solver
model.solver = 'gurobi'

# Set constraints
model.reactions.EX_glc_e.bounds = (-55.1,-46.9)
model.reactions.BIOMASS_ZM.bounds = (0.359, 0.359)
#model.reactions.EX_etoh_e.bounds = (73.3, 83.1)

#enforcing anaerobic conditions
model.reactions.EX_o2_e.bounds = (0, 0)
#model.reactions.LDH_D.bounds = (0,0)
rnf_fld = model.reactions.get_by_id('RNFfld')
rnf_fld.bounds = (0,1000)
rnf_fdx = model.reactions.get_by_id('RNFfdx')
rnf_fdx.bounds = (0,1000)

# f = model.reactions.get_by_id('FNOR')
# fd = model.reactions.get_by_id('FLDR2')
# f.bounds = (0,0)
# fd.bounds = (0,0)

# Run pFBA
print("Running pFBA...")
solution = pfba(model)

# Check the solution status
print(f"\nSolution Status: {solution.status}")

if solution.status == 'optimal':
    print(f"Objective Value (Min Sum of Fluxes): {solution.objective_value:.4f}")
    
    print("\n--- Key Fluxes ---")
    print(f"Biomass: {solution.fluxes['BIOMASS_ZM']:.4f}")
    print(f"Glucose: {solution.fluxes['EX_glc_e']:.4f}")
    print(f"Ethanol: {solution.fluxes['EX_etoh_e']:.4f}")


    print("\n--- All Fluxes ---")
    
    for rxn_id, flux_value in solution.fluxes.items():
        #filrer out reactions with zero flux
        if abs(flux_value) > 1e-6: 
             print(f"{rxn_id}: {flux_value:.4f}")
    
else:
    print("pFBA failed to find a solution.")


