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



# #Escher map
# from escher import Builder
# import json
# import os
# flux_data = solution.fluxes.to_dict() 


# map_filename = 'i_ZM4_489.json' 
    
# # Check if a local map file exists
# if os.path.exists(map_filename):
#     print(f"Loading local map: {map_filename}")
#     with open(map_filename, 'r') as f:
#         local_map_json = json.load(f)
            
#     builder = Builder(
#         p_json=local_map_json,  # Load local map structure
#         reaction_data=flux_data,  # Overlay fluxes
#         model=model               # Pass actual model object for metadata
#         )
# else:
#     print(f"\n[!] Warning: Local map '{map_filename}' not found.")
#     print("Attempting to load standard map (this may look empty if IDs don't match)...")
        

# # Customize the appearance
# builder.reaction_scale = [
#     {'type': 'min', 'color': '#cccccc', 'size': 10},  # Grey for low flux
#     {'type': 'mean', 'color': '#0000ff', 'size': 20}, # Blue for medium
#      {'type': 'max', 'color': '#ff0000', 'size': 30}   # Red for high flux
#     ]


# # Save to HTML
# output_file = 'pfba_map.html'
# builder.save_html(output_file)
# print(f"Map saved to {output_file}")