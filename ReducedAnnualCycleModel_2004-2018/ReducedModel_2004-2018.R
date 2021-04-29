###################################################################################################################################
#Reduced annual-cycle model for monarch butterflies in eastern North America
  #Years included: 2004-2018
  #Model that includes winter, spring, summer covariates

#Not using this model for inferences.  
#Running this to verify that estimates are consistent with those from the full annual-cycle model 
#using data from the same period. 
###################################################################################################################################

#------------------------------------------------------------#
# Set working directory and load packages
#------------------------------------------------------------#

setwd()
# rm(list=ls()) #clear workspace if needed

library(plyr)
library(reshape2)
library(rstan)
rstan_options(auto_write = TRUE)

#------------------------------------------------------------#
# Read in data (assuming these files are in a "Data" subfolder)
#------------------------------------------------------------#

#Summer monarch data
  summer <- read.csv('Data/Monarchs_summer_1994-2018.csv',header=TRUE,stringsAsFactors=FALSE)

#Winter monarch data
  winter <- read.csv('Data/Monarchs_winter_1994-2018.csv',header=TRUE,stringsAsFactors=FALSE)

#Covariates: year
  cov.y <- read.csv('Data/Covariates_Year_1994-2018.csv',header=TRUE)

#Covariates: county
  cov.c <- read.table('Data/Covariates_County_1994-2018.txt',sep='\t',header=TRUE,quote='\"',
                      colClasses=c(rep('numeric',2),rep('character',4),rep('numeric',12)))

#Covariates: county * year
  cov.cy <- read.csv('Data/Covariates_CountyYear_1994-2018.csv',header=TRUE)

#Covariates: county * year * week
  cov.cw <- read.csv('Data/Covariates_CountyWeek_1994-2018.csv',header=TRUE)

#------------------------------------------------------------#
# Format summer survey data
#------------------------------------------------------------#    

#Subset data by year
  yr.min <- 2004
  yr.max <- 2018
  
  summer <- summer[summer$yr %in% yr.min:yr.max,]
  winter <- winter[winter$yr %in% yr.min:yr.max,]
  cov.y <- cov.y[cov.y$yr %in% yr.min:yr.max,]
  cov.cy <- cov.cy[cov.cy$yr %in% yr.min:yr.max,]
  cov.cw <- cov.cw[cov.cw$yr %in% yr.min:yr.max,]

#Specify number of weeks, years, counties, sites, etc.
  uyears <- sort(unique(summer$yr))
  n_years <- length(uyears)
  uweeks <- sort(unique(summer$wk))
  n_weeks <- length(uweeks)
  ucounties <- cov.c$county.ind
  n_counties <- length(ucounties)

#Need to reset site index (since it was built for 1994-2018 data)
  siteind.match <- data.frame(usiteID=sort(unique(summer$usiteID)),site.ind=1:length(unique(summer$usiteID)))
  usites <- siteind.match$site.ind
  n_sites <- length(usites)  
  summer$site.ind <- siteind.match$site.ind[match(summer$usiteID,siteind.match$usiteID)]
  
#Add indicators for state monitoring programs (using NABA as a reference level)
  summer$ia.ind <- ifelse(summer$program=='Iowa',1,0)
  summer$il.ind <- ifelse(summer$program=='Illinois',1,0)
  summer$mi.ind <- ifelse(summer$program=='Michigan',1,0)
  summer$oh.ind <- ifelse(summer$program=='Ohio',1,0)

#Scale effort by the mean
  effortmean <- mean(summer$duration)  #2.17 party hours
  summer$effort.sc <- summer$duration/effortmean  

#Standardize estimates of percent non-forested in immediate survey area (open)
  openS.m <- mean(summer$perc.open)
  openS.sd <- sd(summer$perc.open)
  summer$openS.st <- (summer$perc.open - openS.m)/openS.sd

#------------------------------------------------------------#
# Standardizing annual covariates for summer model
#------------------------------------------------------------#   

#GDD in spring, eastern Texas (spGDD)
  spGDD.m <- mean(cov.y$spGDD.east)
  spGDD.sd <- sd(cov.y$spGDD.east)
  cov.y$spGDD.st <- (cov.y$spGDD.east - spGDD.m)/spGDD.sd
  #Quadratic
  cov.y$spGDD.st2 <- cov.y$spGDD.st*cov.y$spGDD.st
  
#Precipitation in spring, eastern Texas (spPCP)
  spPCP.m <- mean(cov.y$spPCP.east)
  spPCP.sd <- sd(cov.y$spPCP.east)
  cov.y$spPCP.st <- (cov.y$spPCP.east - spPCP.m)/spPCP.sd
  #Quadratic
  cov.y$spPCP.st2 <- cov.y$spPCP.st*cov.y$spPCP.st

#------------------------------------------------------------#
# Standardize county covariates
#------------------------------------------------------------#     

#First, make sure dataframe is sorted by county index
  cov.c <- cov.c[with(cov.c,order(county.ind)),]
  
#Average GDD for weeks 10-24 across all years of the study (avgGDD)
  avgGDD.m <- mean(cov.c$avgGDD.0418)
  avgGDD.sd <- sd(cov.c$avgGDD.0418)
  cov.c$avgGDD.st <- (cov.c$avgGDD.0418 - avgGDD.m)/avgGDD.sd
  
#Average precipitation for AMJJA across all years of the study (avgPCP)
  avgPCP.m <- mean(cov.c$avgPCP.0418)
  avgPCP.sd <- sd(cov.c$avgPCP.0418)
  cov.c$avgPCP.st <- (cov.c$avgPCP.0418 - avgPCP.m)/avgPCP.sd

#Percent crop (crop)
  cropC.m <- mean(cov.c$perc.crop)
  cropC.sd <- sd(cov.c$perc.crop)
  cov.c$cropC.st <- (cov.c$perc.crop - cropC.m)/cropC.sd 

#------------------------------------------------------------#
# Standardize county*year covariates
#------------------------------------------------------------#

#First, make sure dataframe is sorted by county and year index
  cov.cy <- cov.cy[with(cov.cy,order(yr,county.ind)),]  
  
#Glyphosate use (proportion of corn/soy crops sprayed) (gly)
  gly.m <- mean(cov.cy$glyphosate)
  gly.sd <- sd(cov.cy$glyphosate)
  cov.cy$gly.st <- (cov.cy$glyphosate-gly.m)/gly.sd
  
#Summer precipitation (annual deviations from 2004-2018 means for AMJJA) (diffPCP)
  diffPCP.m <- mean(cov.cy$diffPCP.0418)
  diffPCP.sd <- sd(cov.cy$diffPCP.0418)
  cov.cy$diffPCP.st <- (cov.cy$diffPCP.0418 - diffPCP.m)/diffPCP.sd

#------------------------------------------------------------#
# Standardize county*week covariates
#------------------------------------------------------------#

#First, make sure dataframe is sorted by county, week, and year index
  cov.cw <- cov.cw[with(cov.cw,order(yr,wk,county.ind)),]   
  
#Difference between GDD and average GDD for that week and county across all years of the study (diffGDD)
  diffGDD.m <- mean(cov.cw$diffGDD.0418)
  diffGDD.sd <- sd(cov.cw$diffGDD.0418)
  cov.cw$diffGDD.st <- (cov.cw$diffGDD.0418-diffGDD.m)/diffGDD.sd  
  
#------------------------------------------------------------#
# Combine ecological covariates for summer model in long form (nrows = counties*years*weeks)
#------------------------------------------------------------# 
  cyw <- expand.grid(wk=uweeks,county.ind=ucounties,yr=uyears)
  cyw$yr.ind <- cyw$yr-uyears[1]+1
    #Standarize week
    wk.m <- mean(cyw$wk)
    wk.sd <- sd(cyw$wk)
    cyw$wk.st <- (cyw$wk-wk.m)/wk.sd
  cyw <- join(cyw,cov.y[,c('yr','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2')],by='yr',type='left')
  cyw <- join(cyw,cov.c[,c('county.ind','avgGDD.st','avgPCP.st','cropC.st')],by='county.ind',type='left')
  cyw <- join(cyw,cov.cy[,c('county.ind','yr','diffPCP.st','gly.st')],by=c('county.ind','yr'),type='left')
  cyw <- join(cyw,cov.cw[,c('county.ind','yr','wk','diffGDD.st')],by=c('county.ind','yr','wk'),type='left')

  cyw$diffGDD.st2 <- cyw$diffGDD.st*cyw$diffGDD.st
  cyw$diffavgGDD <- cyw$diffGDD.st*cyw$avgGDD.st
  cyw$diffPCP.st2 <- cyw$diffPCP.st*cyw$diffPCP.st
  cyw$diffavgPCP <- cyw$diffPCP.st*cyw$avgPCP.st
  cyw$glycrop <- cyw$gly.st*cyw$cropC.st
  
  cyw <- cyw[,c('county.ind','yr','yr.ind','wk','wk.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2',
                'avgGDD.st','diffGDD.st','diffGDD.st2','diffavgGDD',
                'avgPCP.st','diffPCP.st','diffPCP.st2','diffavgPCP',
                'gly.st','cropC.st','glycrop')]
  names(cyw)[6:20] <- c('spGDD','spGDD2','spPCP','spPCP2',
                        'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                        'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                        'gly','crop','glycrop')
  
  #Need to sort in a particular way to create the model-based index in STAN
    #Put first 5 weeks (16-20) on top
    #Then sort week-fastest, year-slowest: first 4 rows would be county=1, yr=1994, wk=21-24, next 4 rows = county=2, yr=1994, wk=21-24
    cyw15 <- cyw[cyw$wk %in% 16:20,]
    cyw69 <- cyw[cyw$wk %in% 21:24,]
    # head(cyw69); tail(cyw69)
    cyw <- rbind(cyw15,cyw69)
  
  #Create an index for unique combinations of county, year, and week
    cyw$cyw.ind <- 1:nrow(cyw)
  
  #Attach indices to the survey data
    summer <- join(summer,cyw[,c('county.ind','wk','yr','cyw.ind')],by=c('county.ind','wk','yr'),type='left')

#------------------------------------------------------------#
# Standardizing annual covariates for winter model
#------------------------------------------------------------#   

#Percent dense forest cover (forest)
  dense.m <- mean(winter$dense)
  dense.sd <- sd(winter$dense)
  winter$dense.st <- (winter$dense - dense.m)/dense.sd    

#------------------------------------------------------------#
# For each county, calculate weights based on area non-forested land
#------------------------------------------------------------#

#Area non-forested land
  area.open <- as.vector(cov.c$area.land.sqmi*cov.c$perc.open/100)
  weights.open <- area.open/sum(area.open)

#------------------------------------------------------------#
# Calculating observed values for PPC
#------------------------------------------------------------#  

#Summer counts
  #Counts divided by scaled effort (note: includes 0 counts)
  sccounts <- summer$monarch/summer$effort.sc
  sccount.mn <- mean(sccounts)
  sccount.sd <- sd(sccounts)
  
#Winter areas
  area.mn <- mean(winter$dec)
  area.sd <- sd(winter$dec)

#------------------------------------------------------------#
# Package up data
#------------------------------------------------------------#
  
#Survey-level covariates in matrix
  X_survey <- as.matrix(summer[,c('ia.ind','il.ind','mi.ind','oh.ind','openS.st')]) 
  
#County-year-week covariates 
  X_county <- as.matrix(cyw[,c('spGDD','spGDD2','spPCP','spPCP2',
                               'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                               'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                               'gly','crop','glycrop')]) 
  
#Covariates in winter model  
  X_winter <- as.vector(winter$dense.st)
  
#Bundle data
  standata <- list(n_years=n_years,
                   n_counties=n_counties,
                   n_cyw=nrow(cyw),
                   n_sites=n_sites,
                   n_surveys=nrow(summer), 
                   year_id=cyw$yr.ind,
                   county_id=cyw$county.ind,
                   site_id=summer$site.ind,
                   cyw_id=summer$cyw.ind,
                   n_cov_alpha=ncol(X_county),
                   n_cov_beta=ncol(X_survey),
                   X_county=X_county,
                   week_st=cyw$wk.st,
                   X_survey=X_survey,
                   effort=summer$effort.sc,
                   X_winter=X_winter,
                   y_count=summer$monarch,
                   area=winter$dec,
                   weights=weights.open,
                   ind1=seq(1,n_counties*n_years*4,by=4),
                   ind2=seq(1,n_counties*n_years,by=n_counties),
                   start_peak=which(cyw$wk %in% 21:24)[1],
                   n_peak=sum(cyw$wk %in% 21:24),
                   n_cy=n_counties*n_years)
                   
#------------------------------------------------------------#
# MCMC parameters, initial values, parameters to monitor
#------------------------------------------------------------#  

# MCMC parameters
    ni <- 12500     # No. iterations (including warmup)
    nb <- 10000     # No. burn-in iterations to discard (ie, warmup)
    nt <- 1         # Thin rate
    nc <- 3         # No. chains
  
#Initial values
  set.seed(123)
  inits <- lapply(1:nc, function(i)
           list(alpha0=runif(1,1,4),
                alphaFE=runif(ncol(X_county),-0.5,0.5),
                alphaRE_week=runif(1,0,1),
                alphaRE_week2=runif(1,-1,0),
                betaFE=runif(ncol(X_survey),-1,1),
                gamma0=runif(1,0,2),
                gamma_sum=runif(1,0,0.5),
                gammaFE=runif(1,-0.5,0.5),
                sd_county=runif(1,0,1),
                sd_week=runif(1,0,1),
                sd_week2=runif(1,0,1),
                sd_site=runif(1,0,1),
                r_count=runif(1,0,2),
                shape=runif(1,0,2)))

#Parameters to monitor
  params <- c('alpha0','alphaFE','alphaRE_week','alphaRE_week2','betaFE',
              'gamma0','gamma_sum','gammaFE','sd_county','sd_week','sd_week2',
              'sd_site','r_count','shape','pred_orig','pred_sum',
              'y_pred_prop0','sccount_pred_mn','sccount_pred_max','sccount_pred_sd',
              'area_pred_mn','area_pred_max','area_pred_sd',
              'sum_log_lik','win_log_lik')

#------------------------------------------------------------#
# Call STAN from R
#------------------------------------------------------------#
  
  out <- stan('ReducedModel_2004-2018.stan',
              control=list(adapt_delta=0.99),
              data=standata, init=inits, pars=params,
              chains=nc, iter=ni, warmup=nb, thin=nt,
              seed=1, cores=3, open_progress=FALSE)

  print(out,probs=c(0.025,0.5,0.975),digits=3)
