# Data

## Data files for 2004-2018 (full annual-cycle model)
All monarch data from the overwintering grounds and covariate data are publicly available.  Monarch data from the summer breeding grounds are proprietary and are therefore not publicly available (though we provide descriptions of those files here).  
1. [Monarchs_winter.csv](Data/Monarchs_winter.csv): Data on overwintering monarch aggregations, 2004-2018. 
    - supercolony.ind: a unique index identifying each supercolony (1:13)
    - supercolony: name of each supercolony
    - reserve: indicates whether a supercolony is inside (1) or outside (0) of the butterfly reserve
    - yr: year 
    - area: area occupied (ha) by monarch butterflies
    - forest.dense.p: percent of surrounding area with dense forest cover
2. [Covariates_County.txt](Data/Covariates_County.txt): Covariate values associated with each county on the summer breeding range, in the U.S. and Canada.
    - county.ind: a unique index identifying each county (1:545)
    - state.county: combined state-county FIPS code
    - state: FIPS code for each state
    - county: FIPS code for each county in a state
    - state.name: 2-letter code for each state
    - county.name: name of each county
    - area.land.sqmi: land area of each county (sq. mi)
    - Xcentroid: longitude of county centroid
    - Ycentroid: latitude of county centroid
    - avgGDD: average of annual GDD values accumulated between 3 May and 15 Aug (weeks 10-24), 2004-2018
    - avgPCP: average of annual, cumulative precipitation (mm) between Apr and Aug, 2004-2018
    - perc.open: percent of each county that is unforested
    - perc.crop: percent of each county associated with agricultural crops
    - surveyed: indicates whether one or more monarch surveys were conducted in the county between 2004 and 2018 or not (1 or 0, respectively)
3. [Covariates_Year.csv](Data/Covariates_Year.csv): Annual covariate values, 2004-2018.
    - yr: year
    - feb: total area occupied (ha) by monarchs on the overwintering grounds in late February
    - spGDD.east: annual GDD values accumulated between 22 Mar and 2 May (weeks 4-9) in eastern Texas.
    - spPCP.east: cumulative precipitation (mm) between Feb and Apr in eastern Texas.
    - NDVIR1: mean NDVI for the first (northern) half of autumn migration.
4. [Covariates_CountyYear.csv](Data/Covariates_CountyYear.csv): Annual covariate values associated with each county on the summer breeding range, in the U.S. and Canada, 2004-2018.
    - county.ind: a unique index identifying each county (1:545)
    - state.county: combined state-county FIPS code
    - yr: year
    - glyphosate: estimated proportion of corn and soy crops in each county sprayed with glyphosate
    - diffPCP: difference between annual precipitation (cumulative, Apr-Aug) and average precipitation (average of annual values, 2004-2018) in each county
5. [Covariates_CountyWeek.csv](Data/Covariates_CountyWeek.csv): Weekly covariate values associated with each county on the summer breeding range, in the U.S. and Canada, 2004-2018.
    - county.ind: a unique index identifying each county (1:545)
    - state.county: combined state-county FIPS code
    - yr: year
    - wk: week (16-24)
    - diffGDD: difference between GDD and average GDD (2004-2018) for that week and county
6. Monarchs_summer.csv: Data from monarch surveys on the summer breeding grounds, 2004-2018 (not publicly available). 
    - program: name of monitoring program
    - state.county: combined state-county FIPS code
    - lat: latitude of survey location
    - long: longitude of survey location
    - yr: year
    - wk: week (16-24)
    - monarch: total number of adult monarchs observed during survey
    - duration: number of person/party hours spent surveying
    - site.ind: a unique index identifying each survey location
    - county.ind: a unique index identifying each county
    - perc.open: percent of area immediately surrounding each survey location that is unforested 

## Data files for reduced annual-cycle models (1994-2003 or 2004-2018)
1. [Monarchs_winter_1994-2018.csv](Data/Monarchs_winter_1994-2018.csv): Data on overwintering monarch aggregations, 1994-2018. 
    - yr: year 
    - dec: total area occupied (ha) by monarchs on the overwintering grounds in late December
    - dense: percent of surrounding area with dense forest cover, averaged among supercolonies
2. [Covariates_County_1994-2018.txt](Data/Covariates_County_1994-2018.txt): Covariate values associated with each county on the summer breeding range in the U.S.
    - county.ind: a unique index identifying each county (1:502)
    - state.county: combined state-county FIPS code
    - state: FIPS code for each state
    - county: FIPS code for each county in a state
    - state.name: 2-letter code for each state
    - county.name: name of each county
    - area.land.sqmi: land area of each county (sq. mi)
    - Xcentroid: longitude of county centroid
    - Ycentroid: latitude of county centroid
    - avgGDD: average of annual GDD values accumulated between 3 May and 15 Aug (weeks 10-24), 1994-2018
    - avgGDD.9403: average of annual GDD values accumulated between 3 May and 15 Aug (weeks 10-24), 1994-2003
    - avgGDD.0418: average of annual GDD values accumulated between 3 May and 15 Aug (weeks 10-24), 2004-2018
    - avgPCP: average of annual, cumulative precipitation (mm) between Apr and Aug, 1994-2018
    - avgPCP.9403: average of annual, cumulative precipitation (mm) between Apr and Aug, 1994-2003
    - avgPCP.0418: average of annual, cumulative precipitation (mm) between Apr and Aug, 2004-2018
    - perc.open: percent of each county that is unforested
    - perc.crop: percent of each county associated with agricultural crops
    - surveyed: indicates whether one or more monarch surveys were conducted in the county between 2004 and 2018 or not (1 or 0, respectively)
3. [Covariates_Year_1994-2018.csv](Data/Covariates_Year_1994-2018.csv): Annual covariate values, 1994-2018. File structure is the same as that of Covariates_Year.csv except that February estimates of area occupied on the overwintering grounds (feb) and NDVI values in autumn (NDVIR1) are not included because they are unavailable during this time period.  
4. [Covariates_CountyYear_1994-2018.csv](Data/Covariates_CountyYear_1994-2018.csv): Annual covariate values associated with each county on the summer breeding range, in the U.S., 1994-2018. 
    - county.ind: a unique index identifying each county (1:545)
    - state.county: combined state-county FIPS code
    - yr: year
    - glyphosate: estimated proportion of corn and soy crops in each county sprayed with glyphosate
    - diffPCP: difference between annual precipitation (cumulative, Apr-Aug) and average precipitation (average of annual values, 1994-2018) in each county
    - diffPCP.9403: difference between annual precipitation (cumulative, Apr-Aug) and average precipitation (average of annual values, 1994-2003) in each county
    - diffPCP.0418: difference between annual precipitation (cumulative, Apr-Aug) and average precipitation (average of annual values, 2004-2018) in each county
5. [Covariates_CountyWeek_1994-2018.csv](Data/Covariates_CountyWeek_1994-2018.csv): Weekly covariate values associated with each county on the summer breeding range, in the U.S., 1994-2018.
    - county.ind: a unique index identifying each county (1:545)
    - state.county: combined state-county FIPS code
    - yr: year
    - wk: week (16-24)
    - diffGDD: difference between GDD and average GDD (1994-2018) for that week and county
    - diffGDD.9403: difference between GDD and average GDD (1994-2003) for that week and county
    - diffGDD.0418: difference between GDD and average GDD (2004-2018) for that week and county
6. Monarchs_summer_1994-2018.csv: Data from monarch surveys on the summer breeding grounds, 1994-2018 (not publicly available). File structure is the same as that of Monarchs_summer.csv. 
