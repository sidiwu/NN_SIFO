# NN_SIFO

This repo contains the data set and R code for the real application described in the manuscript "Neural Networks for Scalar Input and Functional Output". The manuscript is available at https://link.springer.com/article/10.1007/s11222-023-10287-3.

You should find 3 files in this repo, including:

- **asfr.RData**: the age-specific fertility rate (ASFR) data set.

- **Functions.R**: includes the main functions to impletment the proposed methods (NNBB, NNSS, NNBR & NNSR), along with two existing models (FoS and FAM), and perform hyperparameter tuning for all approaches.
  
- **RealApplication.R**: contains the R code for analyzing the age-specific fertility rate data set using the proposed models (NNBB, NNSS, NNBR & NNSR) and the existing approaches (FoS and FAM). For each model, we run 20 replicates and then evalute the model's predictive performance for further comparison. 
