## Notes for Ramachandran & Gadgil 2023 PBPK human model:
- Oral Dose Folder:
  - There is a driver function for each of the figures to recreate from the paper: PBPK_human_Fig3.m, PBPK_human_Fig4.m, PBPK_human_Fig5.m, PBPK_human_Fig6.m
  - Equations are humanEqns.m
  - loadPhysiology.m
  - loadPartitionCoefficients.m

The code from the authors is on Box (Code from Gadgil PBPK…)

- The code in the lung dose folder is just copied from the code for an oral dose. Compartment 18 (absorption) needs to be removed, and the equation for the gut needs to be altered. The drug needs to be added to compartment 3 (lungs), with a corresponding ka (or the equivalent if absorption is not simply first-order).

- Tissue partition coefficients are calculated with rat values (fraction of neutral lipids, neutral phospholipids, extracellular water, intracellular water, concAP, Ratio albumin to plasma and phIW).

### Differences between supplement vs in Matlab code:
- The liver equation in the supplement doesn’t match the form the authors have in their scripts. But I believe their scripts are correct given what the general form is from the Lyons paper. Specifically, the numerator for the last term should be the same as the first 3 terms in the equations: QLa*CA +Qsp*CVsp + (Ggu-Lgu)*Cgu. Currently it is: QLa*CA +Qsp*Csp + Ggu*Cgu.

- Notice that the time steps vary for drugs. For Fig 3: RIF is 0.01. For Fig 5: RIF is 0.01.

### Differences in parameters in supplement vs in Matlab code:
- Using the ones from their code is what matches the graphs. I did a check for concentrations between our model and theirs to make sure for each figure. 

| RIF	| Supplement | 	Code | 
| --- | ---------- | ----- |
| K.Others | 1.0047 | 1.047 (this is the correct number calculated as median of the other K values) | 
|ka	| 1.07	| 1.08
|fR	| 0.07	| 0.1830
|CL	| 7.79	| 7.86