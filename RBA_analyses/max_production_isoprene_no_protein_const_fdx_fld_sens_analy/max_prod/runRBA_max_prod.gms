************************* Run RBA model ********************
************************* Run RBA model ********************
*       Author: Hoang Dinh
************************************************************

$INLINECOM /*  */
$include "./runRBA_GAMS_settings.txt"
* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal nscale 1e3
$setGlobal nscaleback 1e-3


* For second run when the protein capacity is corrected, overwrite this default value of 1
* Value of 1 means the protein capacity is used freely, which leads to cheap proteins that
* catalyzed flux cycle to be produced.
* Flux cycle consists of the forward and reverse directions of a reversible reaction
* catalyzed by a protein that is cheaper than the dummy protein. Why is it cheaper?
* It is because the composition of that protein has higher fractions of cheap
* amino acids compared to the dummy protein which has the composition of the
* experimentally observations.
* This happens to system where protein availability
* is not the limiting factor (e.g., S. cerevisiae where nutrient or rRNA availability is
* the limiting factor).
*$setGlobal corrected_protein_capacity_percentage 1


* Ribosome efficiency
$setGlobal kribonuc 9.2*3600

* Enforce part of proteome allocate to non-modeled protein
$setGlobal nonmodeled_proteome_allocation 0.41


options
    LP = cplex /*Solver selection*/
    limrow = 0 /*number of equations listed, 0 is suppresed*/
    limcol = 0 /*number of variables listed, 0 is suppresed*/
    iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
    decimals = 8 /*decimal places for display statement*/
    reslim = 1000000 /*wall-clock time limit for solver in seconds*/
    sysout = on /*solver status file report option*/
    solprint = on /*solution printing option*/
        
Sets
i
$include "%species_path%"
j
$include "%rxns_path%"
prosyn(j)
$include "%prosyn_path%"
prowaste(j)
$include "%prowaste_path%"
uptake(j)
$include "%uptake_path%"
media(j) /*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"
;

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
kapp(j)
$include "%kapp_path%"
;

Variables
z, v(j)
;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e3 * %nscale%;


* Simulation top-level settings
* Enable or disable wasteful protein production, disabled by default (to solve faster)
* Note in solving: Enable protein waste flux might cause error for solver
* This is because enabling protein waste introduces several thousands more free variable to the system
* Thus, protein waste should only be implemented with actual data to constrain the free variable
v.fx(j)$prowaste(j) = 0;

* --- LOAD EMPIRICAL PROTEOMICS DATA ---
* This forces the proteins to your exact measured amounts and frees their waste variables
*$include %proconst_path%
v.up('PROWASTE-ZMO0860') = 1e3 * %nscale%;

v.up('PROWASTE-ZMO0456') = 1e3 * %nscale%;

v.up('PROWASTE-ZMO1818') = 1e3 * %nscale%;

v.up('PROWASTE-ZMO0220') = 1e3 * %nscale%;

v.up('PROWASTE-ZMO1851') = 1e3 * %nscale%;


* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1e3 * %nscale%;


$include %phenotype_path%
* Substrate and oxygenation. Set in phenotype.txt
*v.up('RXN-EX_glc__D_e_REV-SPONT') = 1e3 * %nscale%;
*v.fx('RXN-EX_glc__D_e_FWD-SPONT') = 0;

*v.up('RXN-EX_o2_e_REV-SPONT') = 1e3 * %nscale%;
*v.up('RXN-EX_o2_e_FWD-SPONT') = 0;

* Set your NGAM in phenotype.txt since NGAM value depends on growth-condition
* Set your biomass composition in phenotype.txt since biomass composition depends on growth-condition

Scalars
    kapp_fld / 17280.0 /  /* Global turnover/collision limit for flavodoxin */
    kapp_fdx / 17280.0 /  /* Global turnover/collision limit for ferredoxin */
;

Sets
    rxns_fld_red(j) "Reactions that oxidize reduced flavodoxin"
    / 'RXN-DMPPS_fld_FWD-ZMO0875',
      'RXN-FLDR2_REV-ZMO1753',
      'RXN-IPDPS_fld_FWD-ZMO0875',
      'RXN-MECDPDH5_fld_FWD-ZMO0180',
      'RXN-RNFfld_REV-ZMOenz44',
      'RXN-RNTR1c2_fld_FWD-ZMO1025',
      'RXN-RNTR2c2_fld_FWD-ZMO1025',
      'RXN-RNTR3c3_fld_FWD-ZMO1025',
      'RXN-RNTR4c2_fld_FWD-ZMO1025',
      'RXN-HPN1_fld_FWD-ZMO0874',
      'RXN-NIT1b_fld_FWD-ZMOenz45'
    /
    
    rxns_fld_ox(j) "Reactions that reduce oxidized flavodoxin (e.g., RNF, FNOR)"
    / 'RXN-FLDR2_FWD-ZMO1753',
      'RXN-RNFfld_FWD-ZMOenz44',
      'RXN-NIT1b_fld_REV-ZMOenz45'/

    rxns_fdx_red(j) "Reactions that oxidize reduced ferredoxin"
    / 'RXN-FNOR_FWD-ZMO1753',
      'RXN-DMPPS_fdx_FWD-ZMO0875',
      'RXN-IPDPS_fdx_FWD-ZMO0875',
      'RXN-MECDPDH5_fdx_FWD-ZMO0180',
      'RXN-RNFfdx_REV-ZMOenz44',
      'RXN-NIT1b_fdx_FWD-ZMOenz45',
      'RXN-RNTR1c2_fdx_FWD-ZMO1025',
      'RXN-RNTR2c2_fdx_FWD-ZMO1025',
      'RXN-RNTR3c2_fdx_FWD-ZMO1025',
      'RXN-RNTR4c2_fdx_FWD-ZMO1025',
      'RXN-HPN1_fdx_FWD-ZMO0874'/

    rxns_fdx_ox(j) "Reactions that reduce oxidized ferredoxin (e.g., RNF, FNOR)"
    / 'RXN-FNOR_REV-ZMO1753',
    'RXN-RNFfdx_FWD-ZMOenz44',
    'RXN-NIT1b_fdx_REV-ZMOenz45'/
;
*
*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic, RiboCapacityNuc, NonModelProtAllo,
*Obj, Stoic, RiboCapacityNuc, NonModelProtAllo, ModelProtAlloCorrection
Eq_CarrierCap_fld_red, Eq_CarrierCap_fdx_red,
Eq_CarrierCap_fld_ox, Eq_CarrierCap_fdx_ox    
$include %enz_cap_declares_path%
;

Obj..           z =e= -v('%vprod%');
Stoic(i)..      sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityNuc..   v('RIBOSYN-ribonuc') * %kribonuc% =g= %mu% * sum(j$prosyn(j), NAA(j) * v(j));
NonModelProtAllo..  v('BIOSYN-PROTMODELED') =e= (1 - %nonmodeled_proteome_allocation%) * (v('BIOSYN-PROTMODELED') + v('BIOSYN-PROTDUMMY'));

*NonModelProtAllo..  v('BIOSYN-PROTMODELED') =e= (1 - %nonmodeled_proteome_allocation%) * (v('BIOSYN-PROTMODELED') + v('BIOSYN-PROTDUMMY'));
*ModelProtAlloCorrection..   v('BIOSYN-PROTDUMMY2') =g= (1 - %corrected_protein_capacity_percentage%) * (1 - %nonmodeled_proteome_allocation%) * (v('BIOSYN-PROTMODELED') + v('BIOSYN-PROTDUMMY'));

* 1. Bottlenecks for reactions using the REDUCED pool
Eq_CarrierCap_fld_red(j)$rxns_fld_red(j)..
   %mu% * v(j) =L= (kapp_fld) * v('DILUTION-fld_red') ;

Eq_CarrierCap_fdx_red(j)$rxns_fdx_red(j)..
    %mu% * v(j) =L= (kapp_fdx) * v('DILUTION-fdx_red') ;

* 2. Bottlenecks for reactions using the OXIDIZED pool (RNF, FNOR, etc.)
Eq_CarrierCap_fld_ox(j)$rxns_fld_ox(j)..
   %mu% * v(j) =L= (kapp_fld) * v('DILUTION-fld_ox') ;

Eq_CarrierCap_fdx_ox(j)$rxns_fdx_ox(j)..
   %mu% * v(j) =L= (kapp_fdx) * v('DILUTION-fdx_ox') ;

$include %enz_cap_eqns_path%

*** BUILD OPTIMIZATION MODEL ***
Model rba
/Obj, Stoic, RiboCapacityNuc, NonModelProtAllo,
Eq_CarrierCap_fld_red, Eq_CarrierCap_fdx_red, 
 Eq_CarrierCap_fld_ox, Eq_CarrierCap_fdx_ox
*/Obj, Stoic, RiboCapacityNuc, NonModelProtAllo
$include %enz_cap_declares_path%

/;
rba.optfile = 1;

*** SOLVE (baseline, kapp_fld and kapp_fdx at their default medium value) ***
Solve rba using lp minimizing z;

****SETUP PARSIMONIOUS RBA (pRBA) ***
** Declare a new variable for the total sum of all fluxes
*Variables
*    v_sum
*;
*
** Declare a new equation to calculate this sum
*Equations
*    pRBA_Obj
*;
*
** Because variables are all >= 0 (FWD and REV are split), 
**  simply sum them up to get the total network flux
*pRBA_Obj..  v_sum =e= sum(j, v(j));
*
** Create a new model that includes everything from base model and the new equation
*Model rba_prba / rba, pRBA_Obj /;
*rba_prba.optfile = 1;

* allow a tiny 0.1% relaxation (multiplying by 0.999) to prevent numerical 
* infeasibility issues in the solver during the second pass.
*v.lo('%vprod%') = v.l('%vprod%') * 0.999;


* capture the single baseline solve above, before the sweep runs.
execute_unload "model_results.gdx";

file ff /runRBA.modelStat.txt/;
put ff;
put rba.modelStat/;
putclose ff;

file ff2 /runRBA.flux.txt/;
put ff2;
loop(j,
    if ( (v.l(j) gt 0),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j) * %nscaleback%):0:15/;
    );
);
putclose ff2;

*** ================= Sensitivity analysis of Kapps for fld and fdx ==================================
* quantify how much the assumed medium kapp value matters to isoprene production

Set steps / step1*step17 / ;

Parameters
    mult(steps)          'multiplier applied to the medium kapp'
        / step1 0.10, step2 0.25, step3 0.5, step4 1, step5 1.5, step6 2.0, step7 5.0, step8 10.0, step9 50.0, step10 100.0, step11 250.0, step12 500.0,
        step13 750.0, step14 1000.0, step15 5000.0, step16 7500.0, step17 10000.0 /
    obj_history(steps)   'isoprene production at each multiplier'
    stat_history(steps)  'model status (1 = optimal, anything else = suspect)'
;

Scalars base_fld, base_fdx ;
base_fld = kapp_fld ;
base_fdx = kapp_fdx ;

loop(steps,
    kapp_fld = base_fld * mult(steps) ;
    kapp_fdx = base_fdx * mult(steps) ;
    solve rba using lp minimizing z ;
    obj_history(steps)  = -z.l ;
    stat_history(steps) = rba.modelStat ;
);


* Restore the medium kapp so the model is left in its baseline state
kapp_fld = base_fld ;
kapp_fdx = base_fdx ;

*** ================= SWEEP OUTPUT FILE ==================================
file ff3 /runRBA.kappSensitivity.txt/;
put ff3;
put 'multiplier', system.tab, 'production', system.tab, 'modelStat'/;
loop(steps,
    put mult(steps):0:6, system.tab,
        obj_history(steps):0:15, system.tab,
        stat_history(steps):0:0/;
);
putclose ff3;
