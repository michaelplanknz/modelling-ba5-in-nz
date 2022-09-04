# modelling-ba5-in-nz

MODELLING THE IMPACT OF THE OMICRON BA.5 SUBVARIANT IN NEW ZEALAND



Matlab code
===========

The top-level folder contains the Matlab .m files to run the model and reproduce the results in the paper

The key files are:
- main.m - runs the ABC model fitting algorithm, generates model results for a specified set of scenarios, and saves them in /results/
- plotGraphs.m - reads the saved results for a specified scenario, plots graphs and saves them as .png files in /plots/




Data
====

Data files read in by the model are stored in /data/

These are:
- epidataPartial_DD-MMM-YYYY.csv - Ministry of Health data as at DD-MMM-YYYY showing New Zealand's total national number of reported daily cases, new daily hospital admissions, daily deaths, daily infections per 100k in a routinely tested cohort of border works, proportion of new reported cases in over 60s, daily hopsital discharges, proportion of new reported cases that are potential reinfecitons (have already reported a positive test at least 28 days previously). 
- nzcontmatrix.xlsx - 16x16 matrix containing the relative contact rates between 16 five-year age bands taken based on results of Prem et al. and modified by Steyn et al. for the New Zealand population
- popsize_national.xlsx - New Zealand population size in five-year age bands
- reshaped_b2_projecitons_final_2022-07-13.csv - projections of the future uptake of fourth vaccine doses in age bands over 50 years old supplied by the Ministry of Health
- vaccine_data_national.csv - cumulative number of nth doses given to people in 16 five-year age bands by date




