# Inhaled Tuberculosis Treatment Population Pharmacokinetic/Pharmacodynamic Model

The model presented here predicts the performance of inhaled rifampin (RIF) against the standard oral dose for the treatment of tuberculosis (TB). Toxicity in the liver and kidneys is also taken into consideration. Scripts and functions for the other first-line TB medications (isoniazid, pyrazinamide, and ethambutol) are currently in development. Drug-specific scripts and functions are contained within their respective folders, and general methods functions are in Methods/.

The model for each drug is separated into two scripts - regimen_analysis_[DRUG].m and run_popPK_[DRUG]. The former simulates the performance of several possible inhaled regimens of the drug using mean parameter values, calculates the PK metrics AUC and C<sub>max</sub>, and compares them to the standard oral dose for that drug. Once an optimal regimen is chosen, the latter script may be used to simulate population-level pharmacokinetics and pharmacodynamics (PK/PD), returning probability of target attainment (PTA) curves and cumulative fraction of response (CFR) predictions.

Dependencies:
- MATLAB parallel computing toolbox
- MATLAB statistics and machine learning toolbox

Model structure/ODEs were adapted from Ramachandran & Gadgil, 2023, and Himstedt et al., 2022. All model parameters have been altered. Inputs have been scaled down for quick example runs.

Ramachandran & Gadgil DOI: 10.1002/psp4.13008
Himstedt et al. DOI: 10.1093/jac/dkac240

Contact nstrawha@purdue.edu with questions.

## Drug Scripts

### regimen_analysis_[DRUG].m (WIP for all drugs)

Uses a combinatorial analysis approach to generate contour plots and data tables of percent differences between AUC and C<sub>max</sub> for several possible inhaled dose regimens and the standard oral dose. Used to select for an optimal dose for simulation in the run_popPK_[DRUG].m script.

Important inputs:
- n_days_[DRUG]: Days to simulate a dose for. The last day of dosing will be used to calculate PK metrics.
- relevant_compts_[DRUG]: Compartments to track PK metrics for.
- oral_dose_[DRUG]: Standard oral dose amount (mg/dose) to compare all inhaled drug regimens to.
- oral_dose_freq_[DRUG]: Standard oral dose frequency (doses/day) to compare all inhaled drug regimens to.
- lung_dose_min/max_[DRUG]: The minimum/maximum inhaled dose amount to simulate (mg/dose).
- lung_dose_inc_[DRUG]: The increment (mg) to increase the amount of drug given in a certain dose from the minimum to the maximum.
- lung_dose_freq_min/max_[DRUG]: The minimum/maximum inhaled dose frequency to simulate (doses/day). Increment is 1 dose/day.

### run_popPK_[DRUG].m

Conducts a population-level PK/PD simulation of both an oral and an inhaled dose of the drug for TB treatment and toxicity comparison. Generates figures and tables of PK metric comparisons, PTA curves using a drug-TB MIC distribution from the literature, and CFR predictions.

Important inputs:
- n_days_[DRUG]: Days to simulate a dose for. The last day of dosing will be used to calculate PK metrics.
- relevant_compts_[DRUG]: Compartments to track PK metrics for.
- oral/lung_dose_[DRUG]: Oral/inhaled dose amount (mg/dose) to simulate.
- oral/lung_dose_freq_[DRUG]: Oral/inhaled dose frequency (doses/day) to simulate.

## Methods Functions

### [`getParamPDs.m`](methods/getParamPDs.m)

Sets up and returns probability distributions for most model parameters for later sampling. Assumes all blood flow parameters and compartment volume parameters have identical CVs, respectively, which are taken as input. Also sets up and returns empty tables for patient parameter storage.

### [`loadPhysParams.m`](methods/loadPhysParams.m)

Uses the probability distributions returned by [`getParamPDs.m`](methods/getParamPDs.m) to sample and return a set of random parameters, representing an individual patient, that are later used to solve the model's systems of ODEs.

### [DRUG]Oral/LungODEs.m

Used when calling ode15s to set up the system of ODEs for an oral/lung dose. 

### [`solveODEs.m`](methods/solveODEs.m)

Used to link between the main drug scripts and the functions containing the sytsem of ODEs to be solved for the oral and inhaled dose. Handles formatting of ode15s input and output for each system, taking into account the number of days the dose is to be simulated for, the dose amount, and the dose frequency.

### [`plotTimeCourses.m`](methods/plotTimeCourses.m)

Uses the output of [`solveODEs.m`](methods/solveODEs.m) to plot concentration-time courses for all patients individually (for each administration method), as well as percentile plots. 

### [`trackPKMetrics.m`](methods/trackPKMetrics.m)

Uses the output of [`solveODEs.m`](methods/solveODEs.m) to calculate and compare AUC and C<sub>max</sub> for each dose administration method.

### [`plotPTAs.m`](methods/plotPTAs.m)

Uses the output of [`trackPKMetrics.m`](methods/trackPKMetrics.m), as well as a drug-TB MIC distribution and target values from the literature, to generate PTA curves and CFR predictions for each administration method.