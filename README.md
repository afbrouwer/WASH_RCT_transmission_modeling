# WASH_RCT_transmission_modeling

This repository contains the data and code for the papers using the susceptible, infections, susceptible transmission model with an environmental pathway applied to randomized controlled trials (SISE-RCT model).

A public Shiny app for the code is available here: https://umich-biostatistics.shinyapps.io/sise_rct/

The R code supporting the app is given in "SISE-RCT_for_shiny_app.R"


The main directory includes files supporting "Leveraging infectious disease models to interpret 
randomized controlled trials: controlling enteric pathogen transmission through water, sanitation, and hygiene interventions."

Data - "washmodeldata_bangladesh.csv"

Data codebook - "WASH Bangladesh Dataset Codebook.docx"

Code - "Code.R"

This code can reproduce the figures in the paper. Several intermediate data objects are saved to facilitate code modularity.

"sample_and_NLL_cluster_coverage.RDS"

"resample.RDS"

"prevalences.RDS"

The subdirectory Counterfactuals contains the additional code needed to reproduce the results and figrures from "Improving the effectiveness of water, sanitation, and hygiene interventions: a simulation approach to generalizing the outcomes of intervention trials." It includes the following files:

Simulate counterfactuals - "Code_counterfactuals_Bangladesh"

Analyze and plot counterfactual results - "Plot_counterfactuals_Bangladesh"

Counterfactual results - All files starting with "Bprevalences"
