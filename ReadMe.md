## Notes for Ramachandran & Gadgil 2023 PBPK human model:
- Oral Dose Folder:
  - There is a driver function for each of the figures to recreate from the paper: PBPK_human_Fig3.m, PBPK_human_Fig4.m, PBPK_human_Fig5.m, PBPK_human_Fig6.m
  - Equations are humanEqns.m
  - loadPhysiology.m
  - loadPartitionCoefficients.m

The code from the authors is on Box (Code from Gadgil PBPK…)

- Important to note that if we want to use this for other drugs besides the 4 in the paper, you’d have to have data to compare to make sure that this model can properly capture that specific drug dynamic. The code we were given also contains the equations to calculate the partition coefficients. 

- Currently the folder for lung dose is just a copy of the same scripts from the oral dose with NO modifications yet. To do a lung dose, compartment 18 would need to be removed. The drug would be added to compartment 3 (index for lung in our version of the equations). Things to consider – would the ka (absorption rate) be the same in the lung as in the gut? Would the gut (eqn 14) and the gut lumen (eq 17) need to be modified for a lung dose compared to an oral dose? 

- Worth noting that the tissue partition coefficients are calculated with rat values (fraction of neutral lipids, neutral phospholipids, extracellular water, intracellular water, concAP, Ratio albumin to plasma and phIW).
### Differences between supplement vs in Matlab code:
- The liver equation in the supplement doesn’t match the form the authors have in their scripts. But I believe their scripts are correct given what the general form is from the Lyons paper. Specifically, the numerator for the last term should be the same as the first 3 terms in the equations: QLa*CA +Qsp*CVsp + (Ggu-Lgu)*Cgu. Currently it is: QLa*CA +Qsp*Csp + Ggu*Cgu.

- Figure 4 pleura concentrations for PYZ dose is 2000 mg in Matlab scripts, but the paper caption says 1500mg. When trying 1500 mg couldn’t match the graph, but when trying 2000 mg that matched the graph in Figure 4 for PYZ. 

- Notice that the time steps vary for drugs. For Fig 3: RIF is 0.01, while EMB and INH are 0.1, and no step for PYZ. For Fig 5: RIF is 0.01, INH is 0.1 and no step for PYZ. 

### Differences in parameters in supplement vs in Matlab code:
- Using the ones from their code is what matches the graphs. I did a check for concentrations between our model and theirs to make sure for each figure. 

| RIF	| Supplement | 	Code | 
| --- | ---------- | ----- |
| K.Others | 1.0047 | 1.047 (this is the correct number calculated as median of the other K values) | 
|ka	| 1.07	| 1.08
|fR	| 0.07	| 0.1830
|CL	| 7.79	| 7.86

| INH Fast	| Supplement | 	Code | 
| --- | ---------- | ----- |
|ka	 | 2.86 | 	2.89|
|CL |	24.56	| 24.34|

| INH Slow	| Supplement | 	Code | 
| --- | ---------- | ----- |
|CL	|9.16|	9.17|

| EMB	| Supplement | 	Code | 
| --- | ---------- | ----- |
|CL	|49.99|	49.93|

| PYZ	| Supplement | 	Code | 
| --- | ---------- | ----- |
| CL	| 4.10	| 4.14| 
| ka	| 1.36	| 1.39|

## Notes for Reali et al 2023 minimal PBPK mouse model:
- Main driver function is PBPK_Mouse.m
- Equations are mouseEqns.m
- loadPhysiologyMouse.m loads in the blood flow rates and tissue volumes I could find. See below for more details. 
- After loading the physiology the script has 4 options for loading partition coefficients for tissues:
  1.	loadPartitionCoefficients.m - For the missing partition coefficients, I used a mix of some directly from Ramachandran and the Kplu from the Reali supplement Table 3, along with some random numbers that I can’t figure out where I came up with them. This also had two methods for the Other compartment. This one gives the correct shape, although the values seem off. 
  2.	loadPartitionCoefficients_Lyons.m - For the missing partition coefficients, I used those in Lyons et al 2013. These are for a rifampin not rifampicin though and maybe that’s why the graph is incorrect shape (exponential growth of drug in lung, etc). 
  3.	loadPartitionCoefficients_Ramachandran.m - This just uses the partition coefficients from Ramachandran as is. This gives the incorrect shape (exponential growth of drug in lung, etc).
  4.	This option calculates the partition coefficients with the scripts from Ramachandran. It’s identical to #3 currently. Ramachandran calculated the partition coefficients using rat data for their human PBPK. Ruark et al 2014 has the fraction of neutral lipids and fraction of neutral phospholipids for various species, including rat and mice, but I’m not sure how to use the other tables for the other parameters needed for the partition coefficient calculations. This gives the incorrect shape (exponential growth of drug in lung, etc). This option uses: RIF_Properties.m, RatTissueComposition.m, calculateTissuePartition.m, and BloodProperties.m.

- Parameters given (for RIF): 
  - fup = 0.029   - Fraction unbound in plasma (Supp Table 1)
  - CLr = 19 - Fraction of renal clearance fraction (Supp Table 1)
  - F = 1 - % Drug bioavailability (Supp Table 1)
  - ka = 0.3713 - Rate of absorption (1/h) (Supp Table 3) It says the 'best fitted parameter'
  - CL = 0.037 - Total body clearance (L/h) (Supp Table 3) It says the 'best fitted parameter'
  - BP = 0.9 - Blood:plasma ratio

- Parameters need to be identified with a potential number and a reference:
  - All tissue volumes
  - All tissue blood flows
  - All tissue partition coefficients
- I tried to identify the tissue volumes and blood flows from Lee 2020 Supplement and Ruark 2014 Supplement as well as Brown 1997 (https://journals.sagepub.com/doi/epdf/10.1177/074823379701300401)  but unsure how they are doing the 'other' compartment. 

The lumping strategy is a little confusing. “At each step, parallel connected tissues not involved in ADME processes—namely, adipose, bone, brain, heart, muscle, gonads, and skin—were incrementally lumped and the reduced model was validated via a visual predictive comparison with the original model outputs. The lumping criterium was data-driven and applied to the above TB-unessential tissues and organs regardless of their perfusion rates, body proximity, or functionality (Ryu et al., 2022). The connections involving the liver, spleen and pancreas were pooled into the splenic compartment to simplify the chain of compartments. Model simulations supported the definition of a lumped compartment, named “other”, that summarizes all eight tissues, allowing for a 77% reduction in the number of physiological variables in the ODE system (Stader et al., 2019).”
•	Seven compartments are noted in red text. But they mention eight tissues lumped? What is the eighth compartment? The pancreas? 

## Lung Sim ODEs from the Kirschner lab from Wessler et al 2020
- There are 2 folders: Original and Adapted.
- I couldn't just run their scripts so I had to modify it and those scripts are in the Adapted folder
  - AdaptedFuncRunODEs.m is the driver function that reads the parameters from ODEParams.csv and then creates an LHS matrix using LHS_Call.m from the Kirschner labs PRCC paper (Marino et al 2008) then executes the ODEs in ODEquations.m. UnpackODEVars.m and UnpackODEParams.m are called from ODEquations.m.
  - The last part of the script plots 2 figures. The first plots all the variables in a 4x4. The second plots a 1x3 to recreate Figure 4 in Wessler et al 2020. 
