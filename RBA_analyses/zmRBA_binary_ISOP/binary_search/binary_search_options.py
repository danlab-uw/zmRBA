path_gams = '../../GAMS/'
path_pycore = '../../pycore/'
path_enz_mw = '../../application/input/GAMS_model_application/enz_mw_g_per_mmol_ZM.txt'
path_pro_mw = '../../application/input/GAMS_model_application/pro_mw_g_per_mmol_ZM.txt'

report_file = './binary_search_report.txt' # Text file recording binary search process
mu_tol = 1e-5; # Tolerance of upper and lower bound gap to tolerance search
maxiter = 100; # Maximum number of iteration
mu_min0 = 0.0; mu_max0 = 0.05; # User-set initial upper and lower bounds of mu
biom_id = 'BIOSYN-BIODIL' # Set the ID of the biomass reaction
