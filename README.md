#  Monarch Population Dynamics, 1994-2018

### Erin R. Zylstra, et al.
_______________________________________________________________________________________________________________________________________

## R code for 2004-2018 models 
1. [FullModel_2004-2018.R](FullAnnualCycleModel_2004-2018/FullModel_2004-2018.R): R code used to model population dynamics of monarch butterflies (*Danaus plexippus*) in eastern North America between 2004 and 2018. This is the model used for inferences that includes all seasonal covariates.
2. [HierarchicalPartitioning_Summer_2004-2018.R](FullAnnualCycleModel_2004-2018/HierarchicalPartitioning_Summer_2004-2018.R): R code used to run all models needed to evaluate the relative importance of covariates in the summer submodel. Includes code to run hierarchical partitioning analysis.
3. [HierarchicalPartitioning_EarlyWinter_2004-2018.R](FullAnnualCycleModel_2004-2018/HierarchicalPartitioning_EarlyWinter_2004-2018.R): R code used to run all models needed to evaluate the relative importance of covariates in the winter submodel. Includes code to run hierarchical partitioning analysis.

## Stan model files for 2004-2018
1. [FullModel.stan](FullAnnualCycleModel_2004-2018/FullModel.stan): Stan model file, used in FullModel_2004-2018.R.
2. [FullModel_HierPart_Summer.stan](FullAnnualCycleModel_2004-2018/FullModel_HierPart_Summer.stan): Stan model file, used in HierarchicalPartitioning_Summer_2004-2018.R.
3. [FullModel_HierPart_Winter_0FE_NoSummer.stan](FullAnnualCycleModel_2004-2018/FullModel_HierPart_Winter_0FE_NoSummer.stan): Stan model file, used in HierarchicalPartitioning_EarlyWinter_2004-2018.R.
4. [FullModel_HierPart_Winter_0FE_InclSummer.stan](FullAnnualCycleModel_2004-2018/FullModel_HierPart_Winter_0FE_InclSummer.stan): Stan model file, used in HierarchicalPartitioning_EarlyWinter_2004-2018.R.
5. [FullModel_HierPart_Winter_1FE_NoSummer.stan](FullAnnualCycleModel_2004-2018/FullModel_HierPart_Winter_1FE_NoSummer.stan): Stan model file, used in HierarchicalPartitioning_EarlyWinter_2004-2018.R.
6. [FullModel_HierPart_Winter_1FE_InclSummer.stan](FullAnnualCycleModel_2004-2018/FullModel_HierPart_Winter_1FE_InclSummer.stan): Stan model file, used in HierarchicalPartitioning_EarlyWinter_2004-2018.R.
7. [FullModel_HierPart_Winter_2FE_NoSummer.stan](FullAnnualCycleModel_2004-2018/FullModel_HierPart_Winter_2FE_NoSummer.stan): Stan model file, used in HierarchicalPartitioning_EarlyWinter_2004-2018.R.
8. [FullModel_HierPart_Winter_2FE_InclSummer.stan](FullAnnualCycleModel_2004-2018/FullModel_HierPart_Winter_2FE_InclSummer.stan): Stan model file, used in HierarchicalPartitioning_EarlyWinter_2004-2018.R.

## R code for 1994-2003 models 
1. [ReducedModel_1994-2003.R](ReducedAnnualCycleModel_1994-2003/ReducedModel_1994-2003.R): R code used to model population dynamics of monarch butterflies in eastern North America between 1994 and 2003. This is the model used for inferences that includes all seasonal covariates.
2. [HierarchicalPartitioning_1994-2003.R](ReducedAnnualCycleModel_1994-2003/HierarchicalPartitioning_1994-2003.R): R code used to run all models needed to evaluate the relative importance of covariates in the summer submodel. Includes code to run hierarchical partitioning analysis.

## Stan model files for 1994-2003
1. [ReducedModel.stan](ReducedAnnualCycleModel_1994-2003/ReducedModel.stan): Stan model file, used in ReducedModel_1994-2003.R.
2. [ReducedModel_HierPart.stan](ReducedAnnualCycleModel_1994-2003/ReducedModel.stan): Stan model file, used in HierarchicalPartitioning_1994-2003.R.

## Data
Will be provided when paper is published.
