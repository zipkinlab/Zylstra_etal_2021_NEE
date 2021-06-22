# ReducedAnnualCycleModel_1994-2003

## R code for 1994-2003 models 
1. [ReducedModel_1994-2003.R](ReducedAnnualCycleModel_1994-2003/ReducedModel_1994-2003.R): R code used to model population dynamics of monarch butterflies in eastern North America between 1994 and 2003. This is the model used for inferences that includes all seasonal covariates.
2. [HierarchicalPartitioning_1994-2003.R](ReducedAnnualCycleModel_1994-2003/HierarchicalPartitioning_1994-2003.R): R code used to run all models needed to evaluate the relative importance of covariates in the summer submodel. Includes code to run hierarchical partitioning analysis.

## Stan model files for 1994-2003
1. [ReducedModel_1994-2003.stan](ReducedAnnualCycleModel_1994-2003/ReducedModel_1994-2003.stan): Stan model file, used in ReducedModel_1994-2003.R.
2. [ReducedModel_HierPart.stan](ReducedAnnualCycleModel_1994-2003/ReducedModel_HierPart.stan): Stan model file, used in HierarchicalPartitioning_1994-2003.R.
