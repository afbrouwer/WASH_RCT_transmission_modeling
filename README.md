# WASH_RCT_transmission_modeling
This repository contains the data and code for the papers using the SISE-RCT model. The SISE-RCT model is a compartmental infectious disease transission model for enteric pathogens that uses a susceptible-infectious-susceptible (SIS) framework and environmental (E) transmission. It is used to predict the effectiveness of an intervention in a water, sanitation, and hygiene (WASH) randomized controlled trial (RCT) as a function of the intervention coverage, compliance, the fraction of people not in the study who have WASH conditions comparable to the intervention, the intervention efficacy of reducting transmission or shedding, the transmission potential measured by the basic reproduction number, and the fraction of transmission that is along pathways that can be intervened on. 

A Shiny app for exploring how intervention effectiveness is affected by these WASH factors is available here: https://umich-biostatistics.shinyapps.io/sise_rct/


## Directories:
### Support for SISE-RCT Shiny app
The R code supporting the app is given in "SISE-RCT_for_shiny_app.R"

### WASH Benefits Bangladesh data
Contains the data and codebook used in this project. The original data may be accessed at https://osf.io/tprw2/

Data - "washmodeldata_bangladesh.csv"

Data codebook - "WASH Bangladesh Dataset Codebook.docx"

### Original scenario
This folder includes the files supporting "Leveraging infectious disease models to interpret 
randomized controlled trials: controlling enteric pathogen transmission through water, sanitation, and hygiene interventions"  (https://doi.org/10.1371/journal.pcbi.1010748).

Code - "Code.R"

This code can reproduce the figures in the paper. Several intermediate data objects are saved to facilitate code modularity.

"sample_and_NLL_cluster_coverage.RDS"

"resample.RDS"

"prevalences.RDS"

### Counterfactuals
This folder contains the additional code needed to reproduce the results and figures from "Improving the effectiveness of water, sanitation, and hygiene interventions: a simulation approach to generalizing the outcomes of intervention trials" (https://doi.org/10.1101/2022.11.15.22282349). It includes the following files:

Simulate counterfactuals - "Code_counterfactuals_Bangladesh"

Analyze and plot counterfactual results - "Plot_counterfactuals_Bangladesh"

Counterfactual results - All files starting with "Bprevalences"
