import cobra


# Load the model
model = cobra.io.read_sbml_model('i_ZM4_489_isoprene.xml')

# Set solver
model.solver = 'glpk'

# Set constraints
model.reactions.EX_glc_e.bounds = (-55.6,-41.9)
model.reactions.BIOMASS_ZM.bounds = (0.359, 0.359)
#model.reactions.EX_etoh_e.bounds = (73.3,83.1)

rnf_fld = model.reactions.get_by_id('RNFfld')
rnf_fld.bounds = (0,1000)
rnf_fdx = model.reactions.get_by_id('RNFfdx')
rnf_fdx.bounds = (0,1000)

# f = model.reactions.get_by_id('FNOR')
# fd = model.reactions.get_by_id('FLDR2')
# f.bounds = (0,0)
# fd.bounds = (0,0)
#enforcing anaerobic conditions
model.reactions.EX_o2_e.bounds = (0, 0)
model.reactions.LDH_D.bounds = (0,0)


model.objective = "EX_ISOP_e"
solution = model.optimize()
print(solution)


if solution.status == 'optimal':
    print(f"Objective Value: {solution.objective_value:.4f}")
    
    print("\n--- Key Fluxes ---")
    print(f"Biomass: {solution.fluxes['BIOMASS_ZM']:.4f}")
    print(f"Glucose: {solution.fluxes['EX_glc_e']:.4f}")
    print(f"Ethanol: {solution.fluxes['EX_etoh_e']:.4f}")


    print("\n--- All Fluxes ---")
    
    # Option 2: Loop and print them neatly
    for rxn_id, flux_value in solution.fluxes.items():
        #Filtering out fluxes of 0
        if abs(flux_value) > 1e-6: 
             print(f"{rxn_id}: {flux_value:.4f}")

