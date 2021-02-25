###################################################################################################################################
#Full annual-cycle model for monarch butterflies in eastern North America
  #Years included: 2004-2018
  #Model for inference that includes all winter, spring, summer, autumn covariates
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
# Read in data
#------------------------------------------------------------#

#Summer monarch data
  summer <- read.csv('~/Data/Monarchs_summer.csv',header=TRUE,stringsAsFactors=FALSE)

#Winter monarch data
  winter <- read.csv('~/Data/Monarchs_winter.csv',header=TRUE,stringsAsFactors=FALSE)

#Covariates: year
  cov.y <- read.csv('~/]Data/Covariates_Year.csv',header=TRUE)

#Covariates: county
  cov.c <- read.table('~/=Data/Covariates_County.txt',sep='\t',header=TRUE,quote='\"',
                      colClasses=c(rep('numeric',2),rep('character',4),rep('numeric',8)))

#Covariates: county * year
  cov.cy <- read.csv('~/Data/Covariates_CountyYear.csv',header=TRUE)

#Covariates: county * year * week
  cov.cw <- read.csv('~/Data/Covariates_CountyWeek.csv',header=TRUE)

#------------------------------------------------------------#
# Format summer survey data
#------------------------------------------------------------#    

#Specify number of weeks, years, counties, sites, etc.
  uyears <- sort(unique(summer$yr))
  n_years <- length(uyears)
  uweeks <- sort(unique(summer$wk))
  n_weeks <- length(uweeks)
  ucounties <- cov.c$county.ind
  n_counties <- length(ucounties) #total number of counties in study area (regardless of whether surveyed or not)
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
  cov.cy$gly.st <- (cov.cy$glyphosate - gly.m)/gly.sd

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
  cyw$diff2avgGDD <- cyw$diffGDD.st2*cyw$avgGDD.st
  cyw$diffavgPCP <- cyw$diffPCP.st*cyw$avgPCP.st
  cyw$diffPCP.st2 <- cyw$diffPCP.st*cyw$diffPCP.st
  cyw$diff2avgPCP <- cyw$diffPCP.st2*cyw$avgPCP.st
  cyw$glycrop <- cyw$gly.st*cyw$cropC.st
  
  cyw <- cyw[,c('county.ind','yr','yr.ind','wk','wk.st','feb.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2',
                'avgGDD.st','diffGDD.st','diffGDD.st2','diffavgGDD','diff2avgGDD',
                'avgPCP.st','diffPCP.st','diffPCP.st2','diffavgPCP','diff2avgPCP',
                'gly.st','cropC.st','glycrop')]
  names(cyw)[6:23] <- c('feb','spGDD','spGDD2','spPCP','spPCP2',
                        'avgGDD','diffGDD','diffGDD2','diffavgGDD','diff2avgGDD',
                        'avgPCP','diffPCP','diffPCP2','diffavgPCP','diff2avgPCP',
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
# Standardize annual covariate for winter model
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
# Calculate observed values for posterior predictive checks
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

#------------------------------------------------------------#
# Package up data
#------------------------------------------------------------#

#Survey-level covariates in matrix
  X_survey <- as.matrix(summer[,c('ia.ind','il.ind','mi.ind','oh.ind','openS.st')]) 
  
#County-year-week covariates 
  X_county <- as.matrix(cyw[,c('feb','spGDD','spGDD2','spPCP','spPCP2',
                               'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                               'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                               'gly','crop','glycrop')]) 

#Covariates in winter model  
  X_winter <- as.matrix(winter[,c('reserve','dense.st','nectar1.st')])
  
#Bundle data
  standata <- list(n_years=n_years,
                   n_counties=n_counties,
                   n_cyw=nrow(cyw),
                   n_sites=n_sites, 
                   n_surveys=nrow(summer),
                   n_colonies=n_colonies,
                   n_obs=nrow(winter),
                   year_id=cyw$yr.ind,
                   county_id=cyw$county.ind,
                   site_id=summer$site.ind,
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
                   
#------------------------------------------------------------#
# MCMC parameters, initial values, parameters to monitor
#------------------------------------------------------------#  

# MCMC parameters
  
    ni <- 3000       # No. iterations (including warmup)
    nb <- 500        # No. burn-in iterations to discard (ie, warmup)
    nt <- 1          # Thin rate
    nc <- 3          # No. chains
  
#Initial values
  set.seed(1234)
  inits <- lapply(1:nc, function(i)
           list(alpha0=runif(1,2,3),
                alphaFE=runif(ncol(X_county),-0.5,0.5),
                alphaRE_week=runif(1,-0.5,0.5),
                alphaRE_week2=runif(1,-0.5,0.5),
                betaFE=runif(ncol(X_survey),-0.5,0.5),
                pocc_mn=runif(1,0.5,1),
                gamma0=runif(1,-3,-1),
                gamma_sum=runif(1,0,0.5),
                gammaFE=runif(ncol(X_winter),-0.2,0.6),
                sd_county=runif(1,0,1),
                sd_week=runif(1,0,1),
                sd_week2=runif(1,0,1),
                sd_site=runif(1,0,1),
                sd_colony_l=runif(1,0,1),
                sd_colony_g=runif(1,0,1),
                r_count=runif(1,0,2),
                shape=runif(1,0,2)))

#Parameters to monitor
  params <- c('alpha0','alphaFE','alphaRE_week','alphaRE_week2','betaFE','pocc_mn',
              'gamma0','gamma_sum','gammaFE','sd_county','sd_week','sd_week2',
              'sd_site','sd_colony_l','sd_colony_g','r_count','shape','pred_orig','pred_sum',
              'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
              'sum_log_lik','win_log_lik','mu_win')

#------------------------------------------------------------#
# Call STAN from R
#------------------------------------------------------------#

  out <- stan('FullModel.stan',
              control=list(adapt_delta=0.8),
              data=standata, init=inits, pars=params,
              chains=nc, iter=ni, warmup=nb, thin=nt,
              seed=1,cores=3,open_progress=FALSE)
  
  print(out,probs=c(0.025,0.5,0.975),digits=2,include=TRUE)
  posterior <- as.matrix(out)  
  
#------------------------------------------------------------#
# Posterior predictive checks
#------------------------------------------------------------#
  
  #Mean of summer counts (count/effort.sc)
    sccount.mn.post <- posterior[,'sccount_pred_mn']
    # hist(sccount.mn.post,breaks=25,col='gray80')
    # abline(v=sccount.mn,col='dodgerblue3',lwd=2)
    sccount.mn.p <- sum(sccount.mn.post>=sccount.mn)/length(sccount.mn.post)
    
  #SD of summer counts (count/effort.sc)
    sccount.sd.post <- posterior[,'sccount_pred_sd']
    # hist(sccount.sd.post,breaks=25,col='gray80')
    # abline(v=sccount.sd,col='dodgerblue3',lwd=2)
    sccount.sd.p <- sum(sccount.sd.post>=sccount.sd)/length(sccount.sd.post)
    
  #Mean of area occupied in winter, conditional on presence
    area.mn.post <- posterior[,'area_pred_mn']
    # hist(area.mn.post,breaks=25,col='gray80')
    # abline(v=area.mn,col='dodgerblue3',lwd=2)
    area.mn.p <- sum(area.mn.post>=area.mn)/length(area.mn.post)

  #SD of area occupied in winter, conditional on presence
    area.sd.post <- posterior[,'area_pred_sd']
    # hist(area.sd.post,breaks=25,col='gray80')
    # abline(v=area.sd,col='dodgerblue3',lwd=2)
    area.sd.p <- sum(area.sd.post>=area.sd)/length(area.sd.post)
