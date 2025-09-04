Code for **Fast phase prediction of charged copolymer blends by white-box machine learning surrogates**

This repository contains all the code and data to reproduce all numerical results from the paper:

`Phase_Functions.R`:  All auxiliary functions for RPA simulation, uncertainty-aware prediction, etc.
`Datasets_Jun2.R`:  Code for generating all simulated training/testing sets

`FormFactorPreds_July11.R`: Code for generating direct predictions of form factor entries (Figures 2, S2, S3)
`MethodComparison_July10.R`: Code for comparing white-box surrogate model with conventional black-box models (Figure 3)
`Timing_July11.R`: Code for comparing runtime between white-box predictions and RPA simulation (Figure 4)
`CaseStudies_July11.R`: Code for generating true/predicted phase diagrams of selected case studies (Figure 5)
`PhaseBoundaries_July14.R`: Code for generating phase boundary plots (Figure 6)
`Uncertainty_July15.R`: Code for generating plots relevant to model uncertainty (Figure 7)

The `All_Datasets` folder contains all relevant training/testing datasets.

This software is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 2, April 2013.
