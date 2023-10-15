This folder should consist of one .RData file and two .R files. 

In this folder:

1. asfr.RData: the age-specific fertility rate (ASFR) data set.

2. Functions.R: includes the main functions to impletment the proposed methods (NNBB, NNSS, NNBR & NNSR), along with two existing models (FoS and FAM), and perform hyperparameter tuning for all approaches. This file needs to be sourced.

3. RealApplication.R: contains the R codes for analyzing the age-specific fertility rate data set using the proposed models (NNBB, NNSS, NNBR & NNSR) and the existing approaches (FoS and FAM). For each model, we run 20 replicates and then evalute the model's predictive performance for further comparison. 
