###################################################################################################################################
#Reduced annual-cycle model for monarch butterflies in eastern North America
  #Years included: 1994-2003
  #All models for hierarchical partitioning of explained variation in summer counts
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
  yr.min <- 1994
  yr.max <- 2003
  
  summer <- summer[summer$yr %in% yr.min:yr.max,]
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

#Need to reset site.ind (since it was built for 1994-2018 data)
  siteind.match <- data.frame(usiteID=sort(unique(summer$usiteID)),site.ind=1:length(unique(summer$usiteID)))
  usites <- siteind.match$site.ind
  n_sites <- length(usites)  
  summer$site.ind <- siteind.match$site.ind[match(summer$usiteID,siteind.match$usiteID)]
  
#Add indicators for state monitoring programs (using NABA as a reference level)
  #summer$ia.ind <- ifelse(summer$program=='Iowa',1,0)
  summer$il.ind <- ifelse(summer$program=='Illinois',1,0)
  #summer$mi.ind <- ifelse(summer$program=='Michigan',1,0)
  summer$oh.ind <- ifelse(summer$program=='Ohio',1,0)

#Scale effort by the mean
  effortmean <- mean(summer$duration)  #2.55 party hours
  summer$effort.sc <- summer$duration/effortmean  

#Standardize estimates of percent non-forested in immediate survey area (open)
  openS.m <- mean(summer$perc.open)
  openS.sd <- sd(summer$perc.open)
  summer$openS.st <- (summer$perc.open - openS.m)/openS.sd

#------------------------------------------------------------#
# Standardize annual covariates for summer model
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
  avgGDD.m <- mean(cov.c$avgGDD.9403)
  avgGDD.sd <- sd(cov.c$avgGDD.9403)
  cov.c$avgGDD.st <- (cov.c$avgGDD.9403 - avgGDD.m)/avgGDD.sd
  
#Average precipitation for AMJJA across all years of the study (avgPCP)
  avgPCP.m <- mean(cov.c$avgPCP.9403)
  avgPCP.sd <- sd(cov.c$avgPCP.9403)
  cov.c$avgPCP.st <- (cov.c$avgPCP.9403 - avgPCP.m)/avgPCP.sd

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
  
#Summer precipitation (deviation from 1994-2003 means, AMJJA) (diffPCP)
  diffPCP.m <- mean(cov.cy$diffPCP.9403)
  diffPCP.sd <- sd(cov.cy$diffPCP.9403)
  cov.cy$diffPCP.st <- (cov.cy$diffPCP.9403-diffPCP.m)/diffPCP.sd

#------------------------------------------------------------#
# Standardize county*week covariates
#------------------------------------------------------------#
  
#First, make sure dataframe is sorted by county, week, and year index
  cov.cw <- cov.cw[with(cov.cw,order(yr,wk,county.ind)),]   
  
#Difference between GDD and average GDD for that week and county across all years of the study (diffGDD)
  diffGDD.m <- mean(cov.cw$diffGDD.9403)
  diffGDD.sd <- sd(cov.cw$diffGDD.9403)
  cov.cw$diffGDD.st <- (cov.cw$diffGDD.9403-diffGDD.m)/diffGDD.sd  
  
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
# Calculating observed values for posterior predictive checks
#------------------------------------------------------------#  

#Summer counts
  #Counts divided by scaled effort (note: includes 0 counts)
  sccounts <- summer$monarch/summer$effort.sc
  sccount.mn <- mean(sccounts)
  sccount.sd <- sd(sccounts)

#-------------------------------------------------------------#
# Package up data, initial values, parameters, and run in STAN
#-------------------------------------------------------------#
  
#Survey-level covariates in matrix
  X_survey <- as.matrix(summer[,c('il.ind','oh.ind','openS.st')]) 

# MCMC parameters
  ni <- 3000      # No. iterations (including warmup)
  nb <- 500       # No. burn-in iterations to discard (ie, warmup)
  nt <- 1         # Thin rate
  nc <- 3         # No. chains  
  

#-Null model----------------------------------------------------------------#  
    
  #County-year-week covariates 
    X_county <- as.matrix(cyw[,c('avgGDD','avgPCP','crop'
                                 #,'spGDD','spGDD2','spPCP','spPCP2'
                                 #,'diffGDD','diffGDD2','diffavgGDD','diffPCP','diffPCP2','diffavgPCP'
                                 #,'gly','glycrop'
                                 )]) 
  #Bundle data
    standata <- list(n_years=n_years,                 
                     n_counties=n_counties,           
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     y_count=summer$monarch)

  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  phi=runif(1,0,5)))

  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'sccount_pred_mn','sccount_pred_sd','sum_log_lik')
    
  #Run model
    out.null <- stan('ReducedModel_HierPart.stan',
                     control=list(adapt_delta=0.8),data=standata, init=inits, pars=params,
                     chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)

#-Spring weather model----------------------------------------------------------------#  
    
  #County-year-week covariates 
    X_county <- as.matrix(cyw[,c('avgGDD','avgPCP','crop'
                                 ,'spGDD','spGDD2','spPCP','spPCP2'
                                 #,'diffGDD','diffGDD2','diffavgGDD','diffPCP','diffPCP2','diffavgPCP'
                                 #,'gly','glycrop'
                                 )]) 
  #Bundle data
    standata <- list(n_years=n_years,                 
                     n_counties=n_counties,           
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     y_count=summer$monarch)

  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  phi=runif(1,0,5)))

  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'sccount_pred_mn','sccount_pred_sd','sum_log_lik')
    
  #Run model
    out.spC <- stan('ReducedModel_HierPart.stan',
                    control=list(adapt_delta=0.8),data=standata, init=inits, pars=params,
                    chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)

#-Summer weather model----------------------------------------------------------------#  
    
  #County-year-week covariates 
    X_county <- as.matrix(cyw[,c('avgGDD','avgPCP','crop'
                                 #,'spGDD','spGDD2','spPCP','spPCP2'
                                 ,'diffGDD','diffGDD2','diffavgGDD','diffPCP','diffPCP2','diffavgPCP'
                                 #,'gly','glycrop'
                                 )]) 
  #Bundle data
    standata <- list(n_years=n_years,                 
                     n_counties=n_counties,           
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     y_count=summer$monarch)

  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  phi=runif(1,0,5)))

  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'sccount_pred_mn','sccount_pred_sd','sum_log_lik')
    
  #Run model
    out.suC <- stan('ReducedModel_HierPart.stan',
                    control=list(adapt_delta=0.8),data=standata, init=inits, pars=params,
                    chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)

#-Summer land-use model----------------------------------------------------------------#  
    
  #County-year-week covariates 
    X_county <- as.matrix(cyw[,c('avgGDD','avgPCP','crop'
                                 #,'spGDD','spGDD2','spPCP','spPCP2'
                                 #,'diffGDD','diffGDD2','diffavgGDD','diffPCP','diffPCP2','diffavgPCP'
                                 ,'gly','glycrop'
                                 )]) 
  #Bundle data
    standata <- list(n_years=n_years,                 
                     n_counties=n_counties,           
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     y_count=summer$monarch)

  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  phi=runif(1,0,5)))

  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'sccount_pred_mn','sccount_pred_sd','sum_log_lik')
    
  #Run model
    out.suL <- stan('ReducedModel_HierPart.stan',
                    control=list(adapt_delta=0.8),data=standata, init=inits, pars=params,
                    chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)
    
 
#-Spring weather, summer weather model----------------------------------------------------------------#  
    
  #County-year-week covariates 
    X_county <- as.matrix(cyw[,c('avgGDD','avgPCP','crop'
                                 ,'spGDD','spGDD2','spPCP','spPCP2'
                                 ,'diffGDD','diffGDD2','diffavgGDD','diffPCP','diffPCP2','diffavgPCP'
                                 #,'gly','glycrop'
                                 )]) 
  #Bundle data
    standata <- list(n_years=n_years,                 
                     n_counties=n_counties,           
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     y_count=summer$monarch)

  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  phi=runif(1,0,5)))

  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'sccount_pred_mn','sccount_pred_sd','sum_log_lik')
    
  #Run model
    out.spCsuC <- stan('ReducedModel_HierPart.stan',
                       control=list(adapt_delta=0.8),data=standata, init=inits, pars=params,
                       chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE) 
    
#-Spring weather, summer land-use model----------------------------------------------------------------#  
    
  #County-year-week covariates 
    X_county <- as.matrix(cyw[,c('avgGDD','avgPCP','crop'
                                 ,'spGDD','spGDD2','spPCP','spPCP2'
                                 #,'diffGDD','diffGDD2','diffavgGDD','diffPCP','diffPCP2','diffavgPCP'
                                 ,'gly','glycrop'
                                 )]) 
  #Bundle data
    standata <- list(n_years=n_years,                 
                     n_counties=n_counties,           
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     y_count=summer$monarch)

  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  phi=runif(1,0,5)))

  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'sccount_pred_mn','sccount_pred_sd','sum_log_lik')
    
  #Run model
    out.spCsuL <- stan('ReducedModel_HierPart.stan',
                       control=list(adapt_delta=0.8),data=standata, init=inits, pars=params,
                       chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)

#-Summer weather, summer land-use model----------------------------------------------------------------#  
    
  #County-year-week covariates 
    X_county <- as.matrix(cyw[,c('avgGDD','avgPCP','crop'
                                 #,'spGDD','spGDD2','spPCP','spPCP2'
                                 ,'diffGDD','diffGDD2','diffavgGDD','diffPCP','diffPCP2','diffavgPCP'
                                 ,'gly','glycrop'
                                 )]) 
  #Bundle data
    standata <- list(n_years=n_years,                 
                     n_counties=n_counties,           
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     y_count=summer$monarch)

  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  phi=runif(1,0,5)))

  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'sccount_pred_mn','sccount_pred_sd','sum_log_lik')
    
  #Run model
    out.suCsuL <- stan('ReducedModel_HierPart.stan',
                       control=list(adapt_delta=0.8),data=standata, init=inits, pars=params,
                       chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)
    
#-Spring weather, summer weather, summer land-use model----------------------------------------------------------------#  
    
  #County-year-week covariates 
    X_county <- as.matrix(cyw[,c('avgGDD','avgPCP','crop'
                                 ,'spGDD','spGDD2','spPCP','spPCP2'
                                 ,'diffGDD','diffGDD2','diffavgGDD','diffPCP','diffPCP2','diffavgPCP'
                                 ,'gly','glycrop'
                                 )]) 
  #Bundle data
    standata <- list(n_years=n_years,                 
                     n_counties=n_counties,           
                     n_cyw=nrow(cyw),
                     n_surveys=nrow(summer),
                     year_id=cyw$yr.ind,
                     cyw_id=summer$cyw.ind,
                     n_cov_alpha=ncol(X_county),
                     n_cov_beta=ncol(X_survey),
                     X_county=X_county,
                     week_st=cyw$wk.st,
                     X_survey=X_survey,
                     effort=summer$effort.sc,
                     y_count=summer$monarch)

  #Initial values
    set.seed(123)
    inits <- lapply(1:nc, function(i)
             list(alpha0=runif(1,1,4),
                  alphaFE=runif(ncol(X_county),-1,1),
                  alpha_week=runif(n_years,-1,1),
                  alpha_week2=runif(n_years,-1,1),
                  betaFE=runif(ncol(X_survey),-1,1),
                  phi=runif(1,0,5)))

  #Parameters to monitor
    params <- c('alpha0','alphaFE','alpha_week','alpha_week2','betaFE','phi',
                'sccount_pred_mn','sccount_pred_sd','sum_log_lik')
    
  #Run model
    out.spCsuCsuL <- stan('ReducedModel_HierPart.stan',
                          control=list(adapt_delta=0.8),data=standata, init=inits, pars=params,
                          chains=nc, iter=ni, warmup=nb, thin=nt, seed=1, cores=3, open_progress=FALSE)    

#------------------------------------------------------------#
# Hierarchical partitioning
#------------------------------------------------------------#

  outlist.s <- list(out.null=out.null,
                    out.spC=out.spC,
                    out.suC=out.suC,
                    out.suL=out.suL,
                    out.spCsuC=out.spCsuC,
                    out.spCsuL=out.spCsuL,
                    out.suCsuL=out.suCsuL,
                    out.spCsuCsuL=out.spCsuCsuL)
  
  df.s <- data.frame(model=names(outlist.s),stringsAsFactors=FALSE)
  df.s$spC <- ifelse(grepl('spC',df.s$model),1,0)
  df.s$suC <- ifelse(grepl('suC',df.s$model),1,0)
  df.s$suL <- ifelse(grepl('suL',df.s$model),1,0)
  df.s$logL <- sapply(outlist.s,function(x) round(summary(x)$summary['sum_log_lik','50%'])) 
  df.s$maxR <- sapply(outlist.s,function(x) round(max(summary(x)$summary[,'Rhat']),3))
  df.s[order(-df.s$logL),]

#Level of hierarchy (0-3 = 0-3 covariate groups, respectively)
  df.s$hier <- apply(df.s[,c('spC','suC','suL')],1,sum)
  
#Binary string that indicates which covariates in each model
  df.s$modelc <- paste0(df.s$spC,df.s$suC,df.s$suL)
  
#Create dataframe with all nested model pairs
  allpairs.s <- expand.grid(model1=df.s$modelc,model0=df.s$modelc)
  allpairs.s <- allpairs.s[allpairs.s$model1!=allpairs.s$model0,]
  allpairs.s$hier1 <- df.s$hier[match(allpairs.s$model1,df.s$modelc)]
  allpairs.s$hier0 <- df.s$hier[match(allpairs.s$model0,df.s$modelc)]
  allpairs.s <- allpairs.s[allpairs.s$hier1-allpairs.s$hier0==1,]
  allpairs.s$spC.1 <- as.numeric(substr(allpairs.s$model1,1,1))
  allpairs.s$suC.1 <- as.numeric(substr(allpairs.s$model1,2,2))
  allpairs.s$suL.1 <- as.numeric(substr(allpairs.s$model1,3,3))
  allpairs.s$spC.0 <- as.numeric(substr(allpairs.s$model0,1,1))
  allpairs.s$suC.0 <- as.numeric(substr(allpairs.s$model0,2,2))
  allpairs.s$suL.0 <- as.numeric(substr(allpairs.s$model0,3,3))
  allpairs.s$spC <- allpairs.s$spC.1 + allpairs.s$spC.0
  allpairs.s$suC <- allpairs.s$suC.1 + allpairs.s$suC.0
  allpairs.s$suL <- allpairs.s$suL.1 + allpairs.s$suL.0
  allpairs.s$nBoth <- apply(allpairs.s[,c('spC','suC','suL')],1,function(x) sum(x==2))
  allpairs.s <- allpairs.s[allpairs.s$hier0==allpairs.s$nBoth,]

#Calculate difference in logL for each pair
  allpairs.s$logL1 <- df.s$logL[match(allpairs.s$model1,df.s$modelc)]
  allpairs.s$logL0 <- df.s$logL[match(allpairs.s$model0,df.s$modelc)]
  allpairs.s$diff <- allpairs.s$logL1 - allpairs.s$logL0
  
#Identify covariate group that's different in each pair
  allpairs.s$param.spC <- ifelse(allpairs.s$spC.1==1 & allpairs.s$spC.0==0,1,0)
  allpairs.s$param.suC <- ifelse(allpairs.s$suC.1==1 & allpairs.s$suC.0==0,1,0) 
  allpairs.s$param.suL <- ifelse(allpairs.s$suL.1==1 & allpairs.s$suL.0==0,1,0) 
  
#Average logL differences for each covariate group and level of hierarchy
  hpl.s <- data.frame(expand.grid(param=c('spC','suC','suL'),hier=unique(allpairs.s$hier1),stringsAsFactors=FALSE))
  for(i in 1:nrow(hpl.s)){
      hpl.s$avg.level[i] <- mean(allpairs.s$diff[allpairs.s$hier1==hpl.s$hier[i] & allpairs.s[,paste0('param.',hpl.s$param[i])]==1]) 
  }
  
#Mean of those averages for each covariate group  
  hp.s <- data.frame(param=unique(hpl.s$param))
  for(i in 1:nrow(hp.s)){
      hp.s$IC[i] <- mean(hpl.s$avg.level[hpl.s$param==hp.s$param[i]])
  }

#Relativize values (divide by total)
  hp.s$IC.perc <- round(hp.s$IC/sum(hp.s$IC)*100,1)
  hp.s
 
