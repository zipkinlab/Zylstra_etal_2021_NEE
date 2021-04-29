###########################################################################################################################################
#Full-annual-cycle model for monarch butterflies in eastern North America
  #Years included: 2004-2018

#Running models with different combinations of covariates in the winter submodel
#Using hierarchical partitioning methods to calculate the percent of explained variation in the area occupied
  #by monarchs in early winter that is attributable to each covariate.
###########################################################################################################################################

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
  summer <- read.csv('Data/Monarchs_summer.csv',header=TRUE,stringsAsFactors=FALSE)

#Winter monarch data
  winter <- read.csv('Data/Monarchs_winter.csv',header=TRUE,stringsAsFactors=FALSE)

#Covariates: year
  cov.y <- read.csv('Data/Covariates_Year.csv',header=TRUE)

#Covariates: county
  cov.c <- read.table('Data/Covariates_County.txt',sep='\t',header=TRUE,quote='\"',
                      colClasses=c(rep('numeric',2),rep('character',4),rep('numeric',8)))

#Covariates: county * year
  cov.cy <- read.csv('Data/Covariates_CountyYear.csv',header=TRUE)

#Covariates: county * year * week
  cov.cw <- read.csv('Data/Covariates_CountyWeek.csv',header=TRUE)

#------------------------------------------------------------#
# Format summer survey data
#------------------------------------------------------------#    

#Specify number of weeks, years, counties, sites, etc.
  uyears <- sort(unique(summer$yr))
  n_years <- length(uyears)
  uweeks <- sort(unique(summer$wk))
  n_weeks <- length(uweeks)
  ucounties <- cov.c$county.ind
  n_counties <- length(ucounties)
  usites <- sort(unique(summer$site.ind))
  n_sites <- length(usites)

#Add indicators for state monitoring programs (using NABA as a reference level)
  summer$ia.ind <- ifelse(summer$program=='Iowa',1,0)
  summer$il.ind <- ifelse(summer$program=='Illinois',1,0)
  summer$mi.ind <- ifelse(summer$program=='Michigan',1,0)
  summer$oh.ind <- ifelse(summer$program=='Ohio',1,0)

#Scale effort by the mean
  effortmean <- mean(summer$duration)  #2.85 party hours
  summer$effort.sc <- summer$duration/effortmean  

#Standardize estimates of percent non-forested in immediate survey area (open)
  openS.m <- mean(summer$perc.open)
  openS.sd <- sd(summer$perc.open)
  summer$openS.st <- (summer$perc.open - openS.m)/openS.sd

#------------------------------------------------------------#
# Standardize annual covariates for summer model
#------------------------------------------------------------#   

#Size of monarch population in late February (Feb)
  feb.m <- mean(cov.y$feb)
  feb.sd <- sd(cov.y$feb)
  cov.y$feb.st <- (cov.y$feb - feb.m)/feb.sd
  
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
  avgGDD.m <- mean(cov.c$avgGDD)
  avgGDD.sd <- sd(cov.c$avgGDD)
  cov.c$avgGDD.st <- (cov.c$avgGDD - avgGDD.m)/avgGDD.sd
  
#Average precipitation for AMJJA across all years of the study (avgPCP)
  avgPCP.m <- mean(cov.c$avgPCP)
  avgPCP.sd <- sd(cov.c$avgPCP)
  cov.c$avgPCP.st <- (cov.c$avgPCP - avgPCP.m)/avgPCP.sd

#Percent crop (crop)
  cropC.m <- mean(cov.c$perc.crop)
  cropC.sd <- sd(cov.c$perc.crop)
  cov.c$cropC.st <- (cov.c$perc.crop - cropC.m)/cropC.sd 

#------------------------------------------------------------#
# Standardize county*year covariates
#------------------------------------------------------------#
  
#First, make sure dataframe is sorted by county and year index
  cov.cy <- cov.cy[with(cov.cy,order(yr,county.ind)),]  
  
#Proportion of corn/soy crops sprayed with glyphosate (gly)
  gly.m <- mean(cov.cy$glyphosate)
  gly.sd <- sd(cov.cy$glyphosate)
  cov.cy$gly.st <- (cov.cy$glyphosate-gly.m)/gly.sd
  
#Summer precipitation (annual deviations from 2004-2018 means for AMJJA) (diffPCP)
  diffPCP.m <- mean(cov.cy$diffPCP)
  diffPCP.sd <- sd(cov.cy$diffPCP)
  cov.cy$diffPCP.st <- (cov.cy$diffPCP - diffPCP.m)/diffPCP.sd

#------------------------------------------------------------#
# Standardize county*week covariates
#------------------------------------------------------------#
  
#First, make sure dataframe is sorted by county, week, and year index
  cov.cw <- cov.cw[with(cov.cw,order(yr,wk,county.ind)),]   
  
#Difference between GDD and average GDD for that week and county across all years of the study (diffGDD)
  diffGDD.m <- mean(cov.cw$diffGDD)
  diffGDD.sd <- sd(cov.cw$diffGDD)
  cov.cw$diffGDD.st <- (cov.cw$diffGDD-diffGDD.m)/diffGDD.sd  
  
#------------------------------------------------------------#
# Combine ecological covariates for summer model in long form (nrows = counties*years*weeks)
#------------------------------------------------------------# 

  cyw <- expand.grid(wk=uweeks,county.ind=ucounties,yr=uyears)
  cyw$yr.ind <- cyw$yr-2003
    #Standarize week
    wk.m <- mean(cyw$wk)
    wk.sd <- sd(cyw$wk)
    cyw$wk.st <- (cyw$wk-wk.m)/wk.sd
  cyw <- join(cyw,cov.y[,c('yr','feb.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2')],by='yr',type='left')
  cyw <- join(cyw,cov.c[,c('county.ind','avgGDD.st','avgPCP.st','cropC.st')],by='county.ind',type='left')
  cyw <- join(cyw,cov.cy[,c('county.ind','yr','diffPCP.st','gly.st')],by=c('county.ind','yr'),type='left')
  cyw <- join(cyw,cov.cw[,c('county.ind','yr','wk','diffGDD.st')],by=c('county.ind','yr','wk'),type='left')
  
  cyw$diffavgGDD <- cyw$diffGDD.st*cyw$avgGDD.st
  cyw$diffGDD.st2 <- cyw$diffGDD.st*cyw$diffGDD.st
  cyw$diffavgPCP <- cyw$diffPCP.st*cyw$avgPCP.st
  cyw$diffPCP.st2 <- cyw$diffPCP.st*cyw$diffPCP.st
  cyw$glycrop <- cyw$gly.st*cyw$cropC.st
  
  cyw <- cyw[,c('county.ind','yr','yr.ind','wk','wk.st','feb.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2',
                'avgGDD.st','diffGDD.st','diffGDD.st2','diffavgGDD',
                'avgPCP.st','diffPCP.st','diffPCP.st2','diffavgPCP',
                'gly.st','cropC.st','glycrop')]
  names(cyw)[6:21] <- c('feb','spGDD','spGDD2','spPCP','spPCP2',
                        'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                        'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                        'gly','crop','glycrop')
  
  #Need to sort in a particular way to create the model-based index in STAN
    #Put first 5 weeks (16-20) on top
    #Then sort week-fastest, year-slowest: first 4 rows would be county=1, yr=2004, wk=21-24, next 4 rows = county=2, yr=2004, wk=21-24
    cyw15 <- cyw[cyw$wk %in% 16:20,]
    cyw69 <- cyw[cyw$wk %in% 21:24,]
    # head(cyw69); tail(cyw69)
    cyw <- rbind(cyw15,cyw69)
  
  #Create an index for unique combinations of county, year, and week
    cyw$cyw.ind <- 1:nrow(cyw)
  
  #Attach indices to the survey data
    summer <- join(summer,cyw[,c('county.ind','wk','yr','cyw.ind')],by=c('county.ind','wk','yr'),type='left')

#------------------------------------------------------------#
# Format winter data
#------------------------------------------------------------#     

#Sort dataframe
  winter <- winter[with(winter,order(supercolony.ind,yr)),]
  
#Specify number of supercolonies
  ucolonies <- sort(unique(winter$supercolony.ind))
  n_colonies <- length(ucolonies)

#Create an index for monarch presence
  winter$present <- ifelse(winter$area>0,1,0)
  index_present <- which(winter$present==1)

#------------------------------------------------------------#
# Standardize annual covariates for winter model
#------------------------------------------------------------#   

#NDVI in northern half of autumn migration corridor (nectar)
  winter$nectar1 <- cov.y$NDVIR1[match(winter$yr,cov.y$yr)]
  nectar1.m <- mean(winter$nectar1)
  nectar1.sd <- sd(winter$nectar1)
  winter$nectar1.st <- (winter$nectar1 - nectar1.m)/nectar1.sd    

#------------------------------------------------------------#
# Standardize supercolony*year covariate
#------------------------------------------------------------#   

#Percent dense forest cover (forest)
  dense.m <- mean(winter$forest.dense.p)
  dense.sd <- sd(winter$forest.dense.p)
  winter$dense.st <- (winter$forest.dense.p-dense.m)/dense.sd
  
#------------------------------------------------------------#
# For each county, calculate weights based on area non-forested land
#------------------------------------------------------------#

#Area non-forested land
  area.open <- as.vector(cov.c$area.land.sqmi*cov.c$perc.open/100)
  weights.open <- area.open/sum(area.open)

#------------------------------------------------------------#
# Calculating observed values for posterior predictive checks
#------------------------------------------------------------#  

#Summer counts
  #Counts divided by scaled effort (note: includes 0 counts)
  sccounts <- summer$monarch/summer$effort.sc
  sccount.mn <- mean(sccounts)
  sccount.sd <- sd(sccounts)
  
#Winter areas
  #Non-zero areas
  area.n0 <- winter$area[winter$area>0]
  area.mn <- mean(area.n0)
  area.sd <- sd(area.n0)
  
#-------------------------------------------------------------#
# Package up data, initial values, parameters, and run in STAN
#-------------------------------------------------------------#

#Survey-level covariates in matrix
  X_survey <- as.matrix(summer[,c('ia.ind','il.ind','mi.ind','oh.ind','openS.st')]) 
  
#County-year-week covariates 
  X_county <- as.matrix(cyw[,c('feb','spGDD','spGDD2','spPCP','spPCP2',
                               'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                               'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                               'gly','crop','glycrop')]) 
  
#MCMC parameter
  ni <- 3000      # No. iterations (including warmup)
  nb <- 500       # No. burn-in iterations to discard (ie, warmup)
  nt <- 1         # Thin rate
  nc <- 3         # No. chains  
  
#-Null model----------------------------------------------------------------#
  
  #Covariates in winter model  
    #X_winter <- as.matrix(winter[,c('nectar1.st','dense.st')])
    #X_winter <- as.vector(winter[,c('nectar1.st')])
    #X_winter <- as.vector(winter[,c('dense.st')])
    
  #Bundle data
    standata <- list(n_years=n_years,
                     n_counties=n_counties,
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     n_colonies=n_colonies,
                     n_obs=nrow(winter),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     colony_id=winter$supercolony.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     #n_cov_gamma=ncol(X_winter),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     #X_winter=X_winter,
                     y_count=summer$monarch,
                     area=winter$area,
                     index_present=index_present)
                     
  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  pocc_mn=runif(n_colonies,0,1),
                  gamma0=runif(n_colonies,-3,1),
                  #gamma_sum=runif(1,-1,1),
                  #gammaFE=runif(ncol(X_winter),-1,1),
                  shape=runif(1,0,2),
                  phi=runif(1,0,5)))
  
  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'pocc_mn','gamma0','shape',
                'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
                'sum_log_lik','win_log_lik')
  
  out.null <- stan('FullModel_HierPart_Winter_0FE_NoSummer.stan',
                   control=list(adapt_delta=0.8), data=standata, init=inits, pars=params,
                   chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)

#-Summer model----------------------------------------------------------------#
  
  #Covariates in winter model  
    #X_winter <- as.matrix(winter[,c('nectar1.st','dense.st')])
    #X_winter <- as.vector(winter[,c('nectar1.st')])
    #X_winter <- as.vector(winter[,c('dense.st')])
    
  #Bundle data
    standata <- list(n_years=n_years,
                     n_counties=n_counties,
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     n_colonies=n_colonies,
                     n_obs=nrow(winter),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     colony_id=winter$supercolony.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     #n_cov_gamma=ncol(X_winter),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     #X_winter=X_winter,
                     y_count=summer$monarch,
                     area=winter$area,
                     weights=weights.open,
                     ind1=seq(1,n_counties*n_years*4,by=4),
                     ind2=seq(1,n_counties*n_years,by=n_counties),
                     start_peak=which(cyw$wk %in% 21:24)[1],
                     n_peak=sum(cyw$wk %in% 21:24),
                     n_cy=n_counties*n_years,
                     index_present=index_present)
                     
  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  pocc_mn=runif(n_colonies,0,1),
                  gamma0=runif(n_colonies,-3,-1),
                  gamma_sum=runif(1,-0.5,0.5),
                  #gammaFE=runif(ncol(X_winter),-0.5,0.5),
                  shape=runif(1,0,2),
                  phi=runif(1,0,5)))
  
  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'pocc_mn','gamma0','gamma_sum','shape','pred_orig','pred_sum',
                'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
                'sum_log_lik','win_log_lik')
  
  out.S <- stan('FullModel_HierPart_Winter_0FE_InclSummer.stan',
                control=list(adapt_delta=0.8), data=standata, init=inits, pars=params,
                chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)

#-Nectar model----------------------------------------------------------------#
  
  #Covariates in winter model  
    #X_winter <- as.matrix(winter[,c('nectar1.st','dense.st')])
    X_winter <- as.vector(winter[,c('nectar1.st')])
    #X_winter <- as.vector(winter[,c('dense.st')])
    
  #Bundle data
    standata <- list(n_years=n_years,
                     n_counties=n_counties,
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     n_colonies=n_colonies,
                     n_obs=nrow(winter),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     colony_id=winter$supercolony.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     #n_cov_gamma=ncol(X_winter),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     X_winter=X_winter,
                     y_count=summer$monarch,
                     area=winter$area,
                     index_present=index_present)
                     
  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  pocc_mn=runif(n_colonies,0,1),
                  gamma0=runif(n_colonies,-3,-1),
                  #gamma_sum=runif(1,-0.5,0.5),
                  gammaFE=runif(1,-0.5,0.5),
                  shape=runif(1,0,2),
                  phi=runif(1,0,5)))
  
  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'pocc_mn','gamma0','gammaFE','shape',
                'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
                'sum_log_lik','win_log_lik')
  
  out.N <- stan('FullModel_HierPart_Winter_1FE_NoSummer.stan',
                control=list(adapt_delta=0.8), data=standata, init=inits, pars=params,
                chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)
  
#-Forest cover model----------------------------------------------------------------#
  
  #Covariates in winter model  
    #X_winter <- as.matrix(winter[,c('nectar1.st','dense.st')])
    #X_winter <- as.vector(winter[,c('nectar1.st')])
    X_winter <- as.vector(winter[,c('dense.st')])
    
  #Bundle data
    standata <- list(n_years=n_years,
                     n_counties=n_counties,
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     n_colonies=n_colonies,
                     n_obs=nrow(winter),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     colony_id=winter$supercolony.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     #n_cov_gamma=ncol(X_winter),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     X_winter=X_winter,
                     y_count=summer$monarch,
                     area=winter$area,
                     index_present=index_present)
                     
  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  pocc_mn=runif(n_colonies,0,1),
                  gamma0=runif(n_colonies,-3,-1),
                  #gamma_sum=runif(1,-0.5,0.5),
                  gammaFE=runif(1,-0.5,0.5),
                  shape=runif(1,0,2),
                  phi=runif(1,0,5)))
  
  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'pocc_mn','gamma0','gammaFE','shape',
                'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
                'sum_log_lik','win_log_lik')
  
  out.C <- stan('FullModel_HierPart_Winter_1FE_NoSummer.stan',
                control=list(adapt_delta=0.8), data=standata, init=inits, pars=params,
                chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)
  
#-Summer, nectar model----------------------------------------------------------------#
  
  #Covariates in winter model  
    #X_winter <- as.matrix(winter[,c('nectar1.st','dense.st')])
    X_winter <- as.vector(winter[,c('nectar1.st')])
    #X_winter <- as.vector(winter[,c('dense.st')])
    
  #Bundle data
    standata <- list(n_years=n_years,
                     n_counties=n_counties,
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     n_colonies=n_colonies,
                     n_obs=nrow(winter),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     colony_id=winter$supercolony.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     #n_cov_gamma=ncol(X_winter),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     X_winter=X_winter,
                     y_count=summer$monarch,
                     area=winter$area,
                     weights=weights.open,
                     ind1=seq(1,n_counties*n_years*4,by=4),
                     ind2=seq(1,n_counties*n_years,by=n_counties),
                     start_peak=which(cyw$wk %in% 21:24)[1],
                     n_peak=sum(cyw$wk %in% 21:24),
                     n_cy=n_counties*n_years,
                     index_present=index_present)
                     
  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  pocc_mn=runif(n_colonies,0,1),
                  gamma0=runif(n_colonies,-3,-1),
                  gamma_sum=runif(1,-0.5,0.5),
                  gammaFE=runif(1,-0.5,0.5),
                  shape=runif(1,0,2),
                  phi=runif(1,0,5)))
  
  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'pocc_mn','gamma0','gamma_sum','gammaFE','shape','pred_orig','pred_sum',
                'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
                'sum_log_lik','win_log_lik')
  
  out.SN <- stan('FullModel_HierPart_Winter_1FE_InclSummer.stan',
                 control=list(adapt_delta=0.8), data=standata, init=inits, pars=params,
                 chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)
  
#-Summer, forest cover model----------------------------------------------------------------#
  
  #Covariates in winter model  
    #X_winter <- as.matrix(winter[,c('nectar1.st','dense.st')])
    #X_winter <- as.vector(winter[,c('nectar1.st')])
    X_winter <- as.vector(winter[,c('dense.st')])
    
  #Bundle data
    standata <- list(n_years=n_years,
                     n_counties=n_counties,
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     n_colonies=n_colonies,
                     n_obs=nrow(winter),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     colony_id=winter$supercolony.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     #n_cov_gamma=ncol(X_winter),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     X_winter=X_winter,
                     y_count=summer$monarch,
                     area=winter$area,
                     weights=weights.open,
                     ind1=seq(1,n_counties*n_years*4,by=4),
                     ind2=seq(1,n_counties*n_years,by=n_counties),
                     start_peak=which(cyw$wk %in% 21:24)[1],
                     n_peak=sum(cyw$wk %in% 21:24),
                     n_cy=n_counties*n_years,
                     index_present=index_present)
                     
  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  pocc_mn=runif(n_colonies,0,1),
                  gamma0=runif(n_colonies,-3,-1),
                  gamma_sum=runif(1,-0.5,0.5),
                  gammaFE=runif(1,-0.5,0.5),
                  shape=runif(1,0,2),
                  phi=runif(1,0,5)))
  
  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'pocc_mn','gamma0','gamma_sum','gammaFE','shape','pred_orig','pred_sum',
                'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
                'sum_log_lik','win_log_lik')
  
  out.SC <- stan('FullModel_HierPart_Winter_1FE_InclSummer.stan',
                 control=list(adapt_delta=0.8), data=standata, init=inits, pars=params,
                 chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)
  
#-Nectar, forest cover model----------------------------------------------------------------#
  
  #Covariates in winter model  
    X_winter <- as.matrix(winter[,c('nectar1.st','dense.st')])
    #X_winter <- as.vector(winter[,c('nectar1.st')])
    #X_winter <- as.vector(winter[,c('dense.st')])
    
  #Bundle data
    standata <- list(n_years=n_years,
                     n_counties=n_counties,
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     n_colonies=n_colonies,
                     n_obs=nrow(winter),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     colony_id=winter$supercolony.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     n_cov_gamma=ncol(X_winter),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     X_winter=X_winter,
                     y_count=summer$monarch,
                     area=winter$area,
                     index_present=index_present)
                     
  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  pocc_mn=runif(n_colonies,0,1),
                  gamma0=runif(n_colonies,-3,-1),
                  #gamma_sum=runif(1,-0.5,0.5),
                  gammaFE=runif(ncol(X_winter),-0.5,0.5),
                  shape=runif(1,0,2),
                  phi=runif(1,0,5)))
  
  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'pocc_mn','gamma0','gammaFE','shape',
                'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
                'sum_log_lik','win_log_lik')
  
  out.NC <- stan('FullModel_HierPart_Winter_2FE_NoSummer.stan',
                 control=list(adapt_delta=0.8), data=standata, init=inits, pars=params,
                 chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)
  
#-Summer, nectar, forest cover model----------------------------------------------------------------#
  
  #Covariates in winter model  
    X_winter <- as.matrix(winter[,c('nectar1.st','dense.st')])
    #X_winter <- as.vector(winter[,c('nectar1.st')])
    #X_winter <- as.vector(winter[,c('dense.st')])
    
  #Bundle data
    standata <- list(n_years=n_years, 
                     n_counties=n_counties,
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     n_colonies=n_colonies,
                     n_obs=nrow(winter),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     colony_id=winter$supercolony.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     n_cov_gamma=ncol(X_winter),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     X_winter=X_winter,
                     y_count=summer$monarch,
                     area=winter$area,
                     weights=weights.open,
                     ind1=seq(1,n_counties*n_years*4,by=4),
                     ind2=seq(1,n_counties*n_years,by=n_counties),
                     start_peak=which(cyw$wk %in% 21:24)[1],
                     n_peak=sum(cyw$wk %in% 21:24),
                     n_cy=n_counties*n_years,
                     index_present=index_present)
                     
  #Initial values
    set.seed(1234)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  pocc_mn=runif(n_colonies,0,1),
                  gamma0=runif(n_colonies,-3,-1),
                  gamma_sum=runif(1,-0.5,0.5),
                  gammaFE=runif(ncol(X_winter),-0.5,0.5),
                  shape=runif(1,0,2),
                  phi=runif(1,0,5)))
  
  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'pocc_mn','gamma0','gamma_sum','gammaFE','shape','pred_orig','pred_sum',
                'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
                'sum_log_lik','win_log_lik')
  
  out.SNC <- stan('FullModel_HierPart_Winter_2FE_InclSummer.stan',
                  control=list(adapt_delta=0.8), data=standata, init=inits, pars=params,
                  chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)
  
#------------------------------------------------------------#
# Hierarchical partitioning
#------------------------------------------------------------#
  
  outlist.w <- list(out.null=out.null,
                    out.S=out.S,
                    out.N=out.N,
                    out.C=out.C,
                    out.SN=out.SN,
                    out.SC=out.SC,
                    out.NC=out.NC,
                    out.SNC=out.SNC)
  
  df.w <- data.frame(model=names(outlist.w),stringsAsFactors=FALSE)
  df.w$S <- ifelse(grepl('S',df.w$model),1,0)
  df.w$N <- ifelse(grepl('N',df.w$model),1,0)
  df.w$C <- ifelse(grepl('C',df.w$model),1,0)
  df.w$logL <- sapply(outlist.w,function(x) round(summary(x)$summary['win_log_lik','50%'],2)) 
  df.w$maxR <- sapply(outlist.w,function(x) round(max(summary(x)$summary[,'Rhat']),4))
  df.w[order(-df.w$logL),]

#Level of hierarchy (0-3 = 0-3 covariate groups, respectively)
  df.w$hier <- apply(df.w[,c('S','N','C')],1,sum)
  
#Binary string that indicates which covariates in each model
  df.w$modelc <- paste0(df.w$S,df.w$N,df.w$C)
  
#Create dataframe with all nested model pairs
  allpairs.w <- expand.grid(model1=df.w$modelc,model0=df.w$modelc)
  allpairs.w <- allpairs.w[allpairs.w$model1!=allpairs.w$model0,]
  allpairs.w$hier1 <- df.w$hier[match(allpairs.w$model1,df.w$modelc)]
  allpairs.w$hier0 <- df.w$hier[match(allpairs.w$model0,df.w$modelc)]
  allpairs.w <- allpairs.w[allpairs.w$hier1-allpairs.w$hier0==1,]
  allpairs.w$S.1 <- as.numeric(substr(allpairs.w$model1,1,1))
  allpairs.w$N.1 <- as.numeric(substr(allpairs.w$model1,2,2))
  allpairs.w$C.1 <- as.numeric(substr(allpairs.w$model1,3,3))
  allpairs.w$S.0 <- as.numeric(substr(allpairs.w$model0,1,1))
  allpairs.w$N.0 <- as.numeric(substr(allpairs.w$model0,2,2))
  allpairs.w$C.0 <- as.numeric(substr(allpairs.w$model0,3,3))
  allpairs.w$S <- allpairs.w$S.1 + allpairs.w$S.0
  allpairs.w$N <- allpairs.w$N.1 + allpairs.w$N.0
  allpairs.w$C <- allpairs.w$C.1 + allpairs.w$C.0
  allpairs.w$nBoth <- apply(allpairs.w[,c('S','N','C')],1,function(x) sum(x==2))
  allpairs.w <- allpairs.w[allpairs.w$hier0==allpairs.w$nBoth,]

#Calculate difference in logL for each pair
  allpairs.w$logL1 <- df.w$logL[match(allpairs.w$model1,df.w$modelc)]
  allpairs.w$logL0 <- df.w$logL[match(allpairs.w$model0,df.w$modelc)]
  allpairs.w$diff <- allpairs.w$logL1 - allpairs.w$logL0
  
#Identify covariate group that's different in each pair
  allpairs.w$param.S <- ifelse(allpairs.w$S.1==1 & allpairs.w$S.0==0,1,0)
  allpairs.w$param.N <- ifelse(allpairs.w$N.1==1 & allpairs.w$N.0==0,1,0)
  allpairs.w$param.C <- ifelse(allpairs.w$C.1==1 & allpairs.w$C.0==0,1,0)  
  
#Average logL differences for each covariate group and level of hierarchy
  hpl.w <- data.frame(expand.grid(param=c('S','N','C'),hier=unique(allpairs.w$hier1),stringsAsFactors=FALSE))
  for(i in 1:nrow(hpl.w)){
    hpl.w$avg.level[i] <- mean(allpairs.w$diff[allpairs.w$hier1==hpl.w$hier[i] & allpairs.w[,paste0('param.',hpl.w$param[i])]==1]) 
  }
  
#Mean of those averages for each covariate group  
  hp.w <- data.frame(param=unique(hpl.w$param))
  for(i in 1:nrow(hp.w)){
    hp.w$IC[i] <- mean(hpl.w$avg.level[hpl.w$param==hp.w$param[i]])
    hp.w$IC.corr[i] <- max(hp.w$IC[i],0)
  }

#Relativize values (divide by total)
  hp.w$IC.perc <- round(hp.w$IC.corr/sum(hp.w$IC.corr)*100,2)
  hp.w
