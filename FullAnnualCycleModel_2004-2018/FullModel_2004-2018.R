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
library(lme4)
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

#Specify number of weeks, years, counties, sites
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

#Scale survey effort by the mean
  effortmean <- mean(summer$duration)  #2.85 party hours
  summer$effort.sc <- summer$duration/effortmean  

#Standardize estimates of percent non-forested area in immediate survey region (open)
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
# Combine ecological covariates for summer model in long form (nrows = n_counties * n_years * n_weeks)
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
  
  #Need to sort in a particular way to create the model-based peak summer index in STAN
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
  
  #Can include 'countyyr' in parameter list to output expected counts of monarchs in each county each year 
  #(on a typical NABA survey with average effort), but note that this adds 545*15 parameters.

#------------------------------------------------------------#
# Call STAN from R
#------------------------------------------------------------#

  out <- stan('FullModel.stan',
              control=list(adapt_delta=0.8),
              data=standata, init=inits, pars=params,
              chains=nc, iter=ni, warmup=nb, thin=nt,
              seed=1,cores=3,open_progress=FALSE)
  
  print(out,probs=c(0.025,0.5,0.975),digits=3)
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
    
#------------------------------------------------------------#
# Marginal effects of spring weather (Fig. 3c)
#------------------------------------------------------------#

  library(plot.matrix)
    
  spring <- cov.y[,c('yr','spGDD.east','spPCP.east')]
  names(spring)[2:3] <- c('gdd','pcp')
	
	#Posterior samples from relevant parameters 
		amat.sp <- posterior[,c('alpha0','alphaRE_week','alphaRE_week2','alphaFE[2]','alphaFE[3]','alphaFE[4]','alphaFE[5]')]  
 
  #Standardized week values
    wk.r <- sort(unique(cyw$wk.st))

  #Marginal effect of spring climate (on counts in week 21, all other covariates at mean)
    spGDD.r.orig <- range(spring$gdd[spring$yr %in% 2004:2018])  #332-469
    spGDD.orig <- seq(330,470,by=5) 
    # length(spGDD.orig)  #29 columns in heat map
    spGDD.z <- (spGDD.orig - mean(spring$gdd[spring$yr %in% 2004:2018]))/sd(spring$gdd[spring$yr %in% 2004:2018])
    
    spPCP.r.orig <- range(spring$pcp[spring$yr %in% 2004:2018])  #67.9-326.17
    spPCP.orig <- seq(65,330,by=5)
    # length(spPCP.orig)   #54 rows in heat map
    spPCP.z <- (spPCP.orig - mean(spring$pcp[spring$yr %in% 2004:2018]))/sd(spring$pcp[spring$yr %in% 2004:2018])

    #Matrix of GDD*PCP combinations
    spclim <- expand.grid(gdd=spGDD.z,pcp=spPCP.z)
    X.sp <- data.frame(int=1,wk=wk.r[6],wk2=wk.r[6]^2,spGDD=spclim$gdd,spGDD2=spclim$gdd^2,
                       spPCP=spclim$pcp,spPCP2=spclim$pcp^2)
    predl.sp <- as.matrix(X.sp) %*% t(amat.sp) 
    pred.sp.avgeff <- exp(predl.sp)
    pred.sp <- pred.sp.avgeff/effortmean
    #Calculate mean/median, CRIs
    cent.sp <- apply(pred.sp,1,median)
    cri.sp <- apply(pred.sp,1,quantile,probs=c(0.025,0.975))
    #Attach to covariate values:
    sp <- data.frame(gdd.z=spclim$gdd,pcp.z=spclim$pcp,ct.cent=cent.sp,ct.lwr=cri.sp[1,],ct.upr=cri.sp[2,])
    sp$gdd <- sp$gdd.z*sd(spring$gdd[spring$yr %in% 2004:2018]) + mean(spring$gdd[spring$yr %in% 2004:2018])
    sp$pcp <- sp$pcp.z*sd(spring$pcp[spring$yr %in% 2004:2018]) + mean(spring$pcp[spring$yr %in% 2004:2018])
    
    #Create matrix
    sp <- sp[with(sp,order(-pcp,gdd)),]
    ct.mat <- matrix(sp$ct.cent,nrow=length(unique(sp$pcp)),ncol=length(unique(sp$gdd)),byrow=TRUE)    
    
    #Green color ramp
    my_palette <- colorRampPalette(c("#F1F7F1","#508250"))    

    #Converting x/y on plot to GDD/PCP values:
    plot.gdd <- data.frame(gdd=sort(unique(sp$gdd)),x=1:length(unique(sp$gdd)))
    plot.pcp <- data.frame(pcp=sort(unique(sp$pcp)),y=1:length(unique(sp$pcp)))
    
    plot.gddi <- data.frame(gdd=seq(min(plot.gdd$gdd),max(plot.gdd$gdd)))
    plot.gddi$x <- plot.gdd$x[match(plot.gddi$gdd,plot.gdd$gdd)]
    plot.gddi$x <- seq(1,max(plot.gddi$x,na.rm=T),by=0.2)
    
    plot.pcpi <- data.frame(pcp=seq(min(plot.pcp$pcp),max(plot.pcp$pcp)))
    plot.pcpi$y <- plot.pcp$y[match(plot.pcpi$pcp,plot.pcp$pcp)]
    plot.pcpi$y <- seq(1,max(plot.pcpi$y,na.rm=T),by=0.2)
    
    spring2 <- spring[spring$yr %in% 2004:2018,]
    spring2$pcp <- round(spring2$pcp)
    spring2$x <- plot.gddi$x[match(spring2$gdd,plot.gddi$gdd)]
    spring2$y <- plot.pcpi$y[match(spring2$pcp,plot.pcpi$pcp)]   
    
    #Creating colorscale for legend    
    plot.col <- data.frame(ct.cent=seq(min(sp$ct.cent),max(sp$ct.cent),length=101))
    ct.cent100 <- rep(NA,100)
    for(i in 1:100){ct.cent100[i] <- (plot.col$ct.cent[i]+plot.col$ct.cent[i+1])/2}
    plot.col100 <- data.frame(n=1:100,ct.cent=ct.cent100,col=my_palette(100))
    #Out of 100 cells with color gradient, find ones with ct.cent closest to 1,2,3,4,5
    refcols <- plot.col100$col[c(7,26,45,64,83)]
    refdist <- c(7,26,45,64,83)/100*10  #10 is height of legend    
    
    #Plot marginal effects of spring weather (Fig. 3c)
    par(mfrow=c(1,1),mar=c(3,4.5,0.5,3.5),xpd=NA)
    plot(ct.mat,col=my_palette,breaks=seq(min(sp$ct.cent),max(sp$ct.cent),length=100),main='',border=NA,
         axis.row=NULL,axis.col=NULL,xlab='',ylab='',key=NULL)
      axis(1,at=c(0,29.5),tcl=0,labels=F,pos=0)
      axis(1,at=which(plot.gdd$gdd %in% c(360,400,440)),tcl=-0.25,labels=c(360,400,440),pos=0,mgp=c(1.5,0.4,0))
      axis(2,at=c(0,54.5),tcl=0,labels=F,pos=0)
      axis(2,at=which(plot.pcp$pcp %in% seq(100,300,by=50)),tcl=-0.25,labels=seq(100,300,by=50),pos=0,las=1,mgp=c(1.5,0.5,0))
      mtext(expression(paste("GDD (",degree,"C)")),side=1,line=1.8)
      mtext('Precipitation (mm)',side=2,line=3.0)
      arrows(y0=0.5,y1=54.5,x0=plot.gddi$x[which(plot.gddi$gdd==round(mean(spring2$gdd)))],
             x1=plot.gddi$x[which(plot.gddi$gdd==round(mean(spring2$gdd)))],length=0,col='gray90',lty=2)  
      arrows(x0=0.5,x1=29.5,y0=plot.pcpi$y[which(plot.pcpi$pcp==round(mean(spring2$pcp)))],
             y1=plot.pcpi$y[which(plot.pcpi$pcp==round(mean(spring2$pcp)))],length=0,col='gray90',lty=2)  
      text(x=spring2$x,y=spring2$y,as.character(spring2$yr),adj=c(0.5,0.5))
      segments(x0=30.1,x1=30.1,y0=seq(25,34.9,by=0.1),y1=seq(25.1,35,by=0.1),col=my_palette(100),lwd=15)
      segments(x0=30.5,x1=30.7,y0=25+refdist,y1=25+refdist,col='black')
      text(x=30.8,y=25+refdist,adj=c(0,0.5),labels=seq(1,5))
      text(x=29.8,y=36,'Count',adj=c(0,0))
    
#------------------------------------------------------------#
# Marginal effects of summer temperatures (GDD Fig. 4c)
#------------------------------------------------------------#
	#Posterior samples from relevant parameters 
		amat.gdd <- posterior[,c('alpha0','alphaRE_week','alphaRE_week2','alphaFE[6]','alphaFE[7]','alphaFE[8]','alphaFE[9]')]
      
  #Predict effects of diffGDD for a cool/average/warm county
    avgGDD.q <- as.vector(quantile(cyw$avgGDD,probs=c(0.1,0.5,0.9)))
    (avgGDD.q.orig <- avgGDD.q*sd(cov.c$avgGDD) + mean(cov.c$avgGDD))  #701, 909, 1064
    
    #Create sequence of standardized diffGDD values (for all weeks)
    diffGDD.r <- round(range(cyw$diffGDD),1)
    diffGDD <- rep(seq(diffGDD.r[1],diffGDD.r[2],length=100),3)
    
    #Reference avgGDD values on standardized scale
    avgGDD <- rep(avgGDD.q,each=100)

    X.sugdd <- data.frame(int=1,wk=wk.r[6],wk2=wk.r[6]^2,avgGDD=avgGDD,diffGDD=diffGDD,diffGDD2=diffGDD*diffGDD,diffavgGDD=diffGDD*avgGDD)
    predl.sugdd <- as.matrix(X.sugdd) %*% t(amat.gdd)
    pred.sugdd.avgeff <- exp(predl.sugdd)
    pred.sugdd <- pred.sugdd.avgeff/effortmean
    #Calculate medians and credible intervals
    cent.sugdd1 <- apply(pred.sugdd[1:100,],1,median)
    cri.sugdd1 <- apply(pred.sugdd[1:100,],1,quantile,probs=c(0.025,0.975))
    cent.sugdd2 <- apply(pred.sugdd[101:200,],1,median)
    cri.sugdd2 <- apply(pred.sugdd[101:200,],1,quantile,probs=c(0.025,0.975))
    cent.sugdd3 <- apply(pred.sugdd[201:300,],1,median)
    cri.sugdd3 <- apply(pred.sugdd[201:300,],1,quantile,probs=c(0.025,0.975))

    plot.sugdd <- X.sugdd$diffGDD[1:100]*sd(cov.cw$diffGDD) + mean(cov.cw$diffGDD)
      #Note: want to restrict predictions to range of diffGDD values observed at cool/average/warm locations in week 21
      #Using range of GDD differences in week 21 for sites with avgGDD within 5 degree days of cool/average/warm quantiles
      cov.cw <- join(cov.cw,cov.c[,c('state.county','avgGDD')],by='state.county',type='left')
      cold.diff <- range(cov.cw$diffGDD[cov.cw$wk==21 & abs(cov.cw$avgGDD-avgGDD.q.orig[1])<5])
      avg.diff <- range(cov.cw$diffGDD[cov.cw$wk==21 & abs(cov.cw$avgGDD-avgGDD.q.orig[2])<5])
      hot.diff <- range(cov.cw$diffGDD[cov.cw$wk==21 & abs(cov.cw$avgGDD-avgGDD.q.orig[3])<5])  
      cold.ind <- which(plot.sugdd>=cold.diff[1] & plot.sugdd<=cold.diff[2])
      avg.ind <- which(plot.sugdd>=avg.diff[1] & plot.sugdd<=avg.diff[2])
      hot.ind <- which(plot.sugdd>=hot.diff[1] & plot.sugdd<=hot.diff[2])

    #3-color summer pallet 
    mycol3 <- col2rgb(c('salmon1','salmon3','salmon4'))
    col3.1 <- rgb(mycol3[1,1],mycol3[2,1],mycol3[3,1],alpha=255,max=255)
    col3.1p <- rgb(mycol3[1,1],mycol3[2,1],mycol3[3,1],alpha=0.4*255,max=255)
    col3.2 <- rgb(mycol3[1,2],mycol3[2,2],mycol3[3,2],alpha=255,max=255)
    col3.2p <- rgb(mycol3[1,2],mycol3[2,2],mycol3[3,2],alpha=0.4*255,max=255)
    col3.3 <- rgb(mycol3[1,3],mycol3[2,3],mycol3[3,3],alpha=255,max=255)
    col3.3p <- rgb(mycol3[1,3],mycol3[2,3],mycol3[3,3],alpha=0.4*255,max=255)

    #Plot marginal effects of summer GDD (Fig. 4c)
    plot(cent.sugdd1[cold.ind]~plot.sugdd[cold.ind],type='l',lty=1,bty='l',lwd=1.2,col=col3.1,
         xlim=c(-140,200),ylim=c(-0.8,20),axes=F,xlab='',ylab='',xaxs='i',yaxs='i')
      axis(1,at=c(par('usr')[1],200),tcl=0,labels=F)
      axis(1,at=seq(-120,200,by=80),tcl=-0.25,labels=seq(-120,200,by=80))
      axis(2,at=par('usr')[3:4],tcl=0,labels=F)
      axis(2,at=seq(0,20,by=5),tcl=-0.25,labels=seq(0,20,by=5),las=1)  
      polygon(c(plot.sugdd[cold.ind],rev(plot.sugdd[cold.ind])),c(cri.sugdd1[1,cold.ind],rev(cri.sugdd1[2,cold.ind])),col=col3.1p,border=NA)
      lines(cent.sugdd2[avg.ind]~plot.sugdd[avg.ind],type='l',lty=2,bty='l',lwd=1.2,col=col3.2)
      polygon(c(plot.sugdd[avg.ind],rev(plot.sugdd[avg.ind])),c(cri.sugdd2[1,avg.ind],rev(cri.sugdd2[2,avg.ind])),col=col3.2p,border=NA)
      lines(cent.sugdd3[hot.ind]~plot.sugdd[hot.ind],type='l',lty=3,bty='l',lwd=1.2,col=col3.3)
      polygon(c(plot.sugdd[hot.ind],rev(plot.sugdd[hot.ind])),c(cri.sugdd3[1,hot.ind],rev(cri.sugdd3[2,hot.ind])),col=col3.3p,border=NA)
      legend(x=-140,y=20,legend=c('Cool','Average','Warm'),lty=c(1,2,3),col=c(col3.1,col3.2,col3.3),bty='n',cex=0.9)  
      mtext(expression(paste("GDD (",degree,"C), deviation")),side=1,line=2,cex=0.9)
      mtext('Monarchs count (adults/hr)',side=2,las=0,line=1.8,cex=0.9)

#------------------------------------------------------------#
# Marginal effects of summer precipitation (Fig. 4d)
#------------------------------------------------------------#
	#Posterior samples from relevant parameters 
		amat.pcp <- posterior[,c('alpha0','alphaRE_week','alphaRE_week2','alphaFE[10]','alphaFE[11]','alphaFE[12]','alphaFE[13]')]

  #Predict effects of diffPCP for a dry/average/wet county
    avgPCP.q <- as.vector(quantile(cyw$avgPCP,probs=c(0.1,0.5,0.9)))
    (avgPCP.q.orig <- avgPCP.q*sd(cov.c$avgPCP) + mean(cov.c$avgPCP))  #427, 532, 626
      #The high value reflects the 0.9 quantile for 2004-2018 values, which is at the extreme of avgPCP values for 1994-2003,
      #and the 1994-2003 values are on the x-axis for the adjacent panel with precipitation trends. So, adjusting the highest value a bit
    avgPCP.q.orig[3] <- 600
    avgPCP.q[3] <- (avgPCP.q.orig[3] - mean(cov.c$avgPCP))/sd(cov.c$avgPCP)
    
    #Create sequence of standardized diffPCP values
    diffPCP.r <- round(range(cyw$diffPCP),1)
    diffPCP <- rep(seq(diffPCP.r[1],diffPCP.r[2],length=100),3)
    
    #Reference avgPCP values on standardized scale
    avgPCP <- rep(avgPCP.q,each=100)

    X.supcp <- data.frame(int=1,wk=wk.r[6],wk2=wk.r[6]^2,avgPCP=avgPCP,diffPCP=diffPCP,diffPCP2=diffPCP*diffPCP,diffavgPCP=diffPCP*avgPCP)
    predl.supcp <- as.matrix(X.supcp) %*% t(amat.pcp)
    pred.supcp.avgeff <- exp(predl.supcp)
    pred.supcp <- pred.supcp.avgeff/effortmean
    #Calculate medians, CRIs
    cent.supcp1 <- apply(pred.supcp[1:100,],1,median)
    cri.supcp1 <- apply(pred.supcp[1:100,],1,quantile,probs=c(0.025,0.975))
    cent.supcp2 <- apply(pred.supcp[101:200,],1,median)
    cri.supcp2 <- apply(pred.supcp[101:200,],1,quantile,probs=c(0.025,0.975))
    cent.supcp3 <- apply(pred.supcp[201:300,],1,median)
    cri.supcp3 <- apply(pred.supcp[201:300,],1,quantile,probs=c(0.025,0.975))    

    plot.supcp <- X.supcp$diffPCP[1:100]*sd(cov.cy$diffPCP) + mean(cov.cy$diffPCP)
      #Note: want to restrict predictions to range of diffPCP values observed at dry/avg/wet locations
      #Using range of PCP differences for sites with avgPCP within 10 mm of dry/avg/wet quantiles
      cov.cy <- join(cov.cy,cov.c[,c('state.county','avgPCP')],by='state.county',type='left')
      dry.diff <- range(cov.cy$diffPCP[abs(cov.cy$avgPCP-avgPCP.q.orig[1])<10])
      avgp.diff <- range(cov.cy$diffPCP[abs(cov.cy$avgPCP-avgPCP.q.orig[2])<10])
      wet.diff <- range(cov.cy$diffPCP[abs(cov.cy$avgPCP-avgPCP.q.orig[3])<10])  
      dry.ind <- which(plot.supcp>=dry.diff[1] & plot.supcp<=dry.diff[2])
      avgp.ind <- which(plot.supcp>=avgp.diff[1] & plot.supcp<=avgp.diff[2])
      wet.ind <- which(plot.supcp>=wet.diff[1] & plot.supcp<=wet.diff[2])

    #Plot marginal effects of summer PCP (Fig. 4d)
    plot(cent.supcp1[dry.ind]~plot.supcp[dry.ind],type='l',lty=1,bty='l',lwd=1.2,col=col3.1,
         xlim=c(-360,525),ylim=c(-0.8,20),axes=F,xlab='',ylab='',xaxs='i',yaxs='i')
      axis(1,at=c(par('usr')[1],520),tcl=0,labels=F)
      axis(1,at=seq(-300,500,by=200),tcl=-0.25,labels=seq(-300,500,by=200))
      axis(2,at=par('usr')[3:4],tcl=0,labels=F)
      axis(2,at=seq(0,20,by=5),tcl=-0.25,labels=seq(0,20,by=5),las=1)    
      polygon(c(plot.supcp[dry.ind],rev(plot.supcp[dry.ind])),c(cri.supcp1[1,dry.ind],rev(cri.supcp1[2,dry.ind])),col=col3.1p,border=NA)
      lines(cent.supcp2[avgp.ind]~plot.supcp[avgp.ind],type='l',lty=2,bty='l',lwd=1.2,col=col3.2)
      polygon(c(plot.supcp[avgp.ind],rev(plot.supcp[avgp.ind])),c(cri.supcp2[1,avgp.ind],rev(cri.supcp2[2,avgp.ind])),col=col3.2p,border=NA)
      lines(cent.supcp3[wet.ind]~plot.supcp[wet.ind],type='l',lty=3,bty='l',lwd=1.2,col=col3.3)
      polygon(c(plot.supcp[wet.ind],rev(plot.supcp[wet.ind])),c(cri.supcp3[1,wet.ind],rev(cri.supcp3[2,wet.ind])),col=col3.3p,border=NA)
      legend(x=-340,y=20,legend=c('Dry','Average','Wet'),lty=c(1,2,3),col=c(col3.1,col3.2,col3.3),bty='n',cex=0.9)  
      mtext("Precipitation (mm), difference",side=1,line=1.8,cex=0.9) 
      mtext('Monarchs count (adults/hr)',side=2,las=0,line=1.8,cex=0.9)

#------------------------------------------------------------#
# Correlations between seasonal population indices (Fig. 5)
#------------------------------------------------------------#
  
  #Late winter [February] (note that 2004 is imputed based on linear regression: Feb~Dec for 2005-2018)
  popdf <- cov.y[,c('yr','feb')]
  
  #Peak summer (index values from 2004-2018 FACM model)
  summ <- posterior[,grep('pred_orig',colnames(posterior))]
  #pred_orig is the expected count for a NABA survey with average effort (=2.85 hours)
  #Calculating expected count per person hour
  popdf$sum.md <- apply(summ,2,median)/2.85
  popdf$sum.lcl <- apply(summ,2,quantile,probs=0.025)/2.85
  popdf$sum.ucl <- apply(summ,2,quantile,probs=0.975)/2.85
  
  #Early winter [December]
  popdf$dec <- ddply(winter,.(yr),summarize,dec=sum(area))[,2]

  #Regression: Summer ~ February (excluding 2004 since Feb value was imputed)
    library(jagsUI) 
    library(runjags)
    jagsdata <- list(nobs=nrow(popdf[popdf$yr!=2004,]),
                     x=popdf$feb[popdf$yr!=2004],
                     y=popdf$sum.md[popdf$yr!=2004])

    sink("LinearRegression.txt")
	  cat("
		  model{
			  for(i in 1:nobs){
				  y[i] ~ dnorm(mu[i],tau)
		  		mu[i] <- b0 + b1*x[i]
		  	 }

			  b0 ~ dnorm(0,0.0001)
			  b1 ~ dnorm(0,0.0001)
			  sigma ~ dunif(0,100)
			  tau <- 1/(sigma*sigma)
		  } #model
	  ",fill=TRUE)
	  sink()

  	ni <- 2000; na <- 1000; nb <- 8000; nc <- 3; ni.tot <- ni + nb
  	params <- c('b0','b1','sigma')
    inits <- lapply(1:nc, function(i)
             list(b0=runif(1,-2,2),
                  b1=runif(1,-2,2),
                  sigma=runif(1,1,10)))	
    fit.sumfeb <- jags(data=jagsdata, inits=inits, parameters.to.save=params,
                       model.file='LinearRegression.txt',
                       n.chains=nc, n.adapt=na, n.iter=ni.tot, n.burnin=nb,
                       parallel=T, n.cores=3, DIC=F)
    print(fit.sumfeb)
		fit.sumfeb.s <- fit.sumfeb$samples
		comb.sumfeb <- combine.mcmc(fit.sumfeb.s)
    Xsumfeb <- matrix(c(rep(1,nrow(popdf[popdf$yr!=2004,])),popdf$feb[popdf$yr!=2004]),byrow=FALSE,ncol=2)     
    
    #Function to calculate Bayesian R2 values based on Gelman et al. 2019 (Am Stat)
    #For linear regressions with one predictor variable (X)
    #Xmat = matrix with first column of 1's and second column = X
    #comb = matrix with posterior samples from JAGS fit: nrow = nsamples, cols = b0, b1, sigma
    bayes_R2 <- function(Xmat,comb){
		  y_pred <- Xmat %*% t(comb[,1:2])
      var_fit <- apply(y_pred,2,var)
      var_res <- comb[,3]^2
      bayesR2 <- var_fit / (var_fit + var_res)		  
		}    
    
    rsq_bayes <- bayes_R2(Xsumfeb,comb.sumfeb)
    median(rsq_bayes)  #0.31
    predfeb <- seq(min(popdf$feb[popdf$yr!=2004]),max(popdf$feb[popdf$yr!=2004]),length=100)
    sumfeb.preds <- matrix(c(rep(1,100),predfeb),ncol=2,byrow=FALSE) %*% t(comb.sumfeb[,1:2])
    cent.sumfeb <- apply(sumfeb.preds,1,median)
    cri.sumfeb <- apply(sumfeb.preds,1,quantile,probs=c(0.025,0.975))
  
  #Regression: December ~ Summer -- using median of summer index
  	jagsdata <- list(nobs=nrow(popdf),
                     x=popdf$sum.md,
                     y=popdf$dec)    
		fit.decsum <- jags(data=jagsdata, inits=inits, parameters.to.save=params,
                       model.file='LinearRegression.txt',
                       n.chains=nc, n.adapt=na, n.iter=ni.tot, n.burnin=nb,
                       parallel=T, n.cores=3, DIC=F)
		print(fit.decsum)    
    fit.decsum.s <- fit.decsum$samples
		comb.decsum <- combine.mcmc(fit.decsum.s)
    Xdecsum <- matrix(c(rep(1,nrow(popdf)),popdf$sum.md),byrow=FALSE,ncol=2)     
    rsq_bayes <- bayes_R2(Xdecsum,comb.decsum)
    median(rsq_bayes)  #0.67
    predsum <- seq(min(popdf$sum.lcl),max(popdf$sum.ucl),length=100)
    decsum.preds <- matrix(c(rep(1,100),predsum),ncol=2,byrow=FALSE) %*% t(comb.decsum[,1:2])
    cent.decsum <- apply(decsum.preds,1,median)
    cri.decsum <- apply(decsum.preds,1,quantile,probs=c(0.025,0.975))
 
  #Regression: February ~ December
    feb1 <- popdf$feb[popdf$yr %in% 2005:2018]
    dec1 <- popdf$dec[popdf$yr %in% 2004:2017]    
  	jagsdata <- list(nobs=length(dec1),x=dec1,y=feb1)
		fit.febdec <- jags(data=jagsdata, inits=inits, parameters.to.save=params,
                       model.file='LinearRegression.txt',
                       n.chains=nc, n.adapt=na, n.iter=ni.tot, n.burnin=nb,
                       parallel=T, n.cores=3, DIC=F)
		print(fit.febdec)
		fit.febdec.s <- fit.febdec$samples
		comb.febdec <- combine.mcmc(fit.febdec.s)
    Xfebdec <- matrix(c(rep(1,length(dec1)),dec1),byrow=FALSE,ncol=2)     
    rsq_bayes <- bayes_R2(Xfebdec,comb.febdec)
    median(rsq_bayes)  #0.43
    preddec <- seq(min(dec1),max(dec1),length=100)
    febdec.preds <- matrix(c(rep(1,100),preddec),ncol=2,byrow=FALSE) %*% t(comb.febdec[,1:2])
    cent.febdec <- apply(febdec.preds,1,median)
    cri.febdec <- apply(febdec.preds,1,quantile,probs=c(0.025,0.975))

  #3-panel figure, with each index plotted against index from previous season (Fig. 5)
  cex.mt <- 0.95
  cext <- 1.1

  par(mfrow=c(3,1),mar=c(3.6,3.5,0,1),oma=c(0,0,1,0),cex=0.8)
  plot(dec~sum.md,data=popdf,type='n',xlim=c(0.7,11.5),ylim=c(-0.7,12),axes=F,xlab='',ylab='',xaxs='i',yaxs='i')
    arrows(x0=popdf$sum.lcl,x1=popdf$sum.ucl,y0=popdf$dec,y1=popdf$dec,length=0,col='gray50')
    points(dec~sum.md,data=popdf,pch=19,col='grey40')
    lines(cent.decsum~predsum,col='grey40')
    polygon(c(predsum,rev(predsum)),c(cri.decsum[1,],rev(cri.decsum[2,])),col=rgb(0,0,0,0.1),border=NA)
    axis(1,at=par('usr')[1:2],tcl=0,labels=F)
    axis(1,at=seq(2,10,2),tcl=-0.25,labels=seq(2,10,2),cex=0.9,mgp=c(0,0.6,0))
    axis(2,at=par('usr')[3:4],tcl=0,labels=F)
    axis(2,at=seq(0,12,3),tcl=-0.25,labels=seq(0,12,3),las=1,mgp=c(0,0.7,0))
    mtext('Early-winter population (ha)',side=2,line=2.0,cex=cex.mt,col='steelblue2')
    mtext('Summer population (mon/hr)',side=1,line=1.9,cex=cex.mt,col='salmon2')
    text(x=par('usr')[2]*0.98,y=-0.7+0.14*12.7,adj=c(1,0),'y = 0.38 + 0.76x',cex=cext)
    text(x=par('usr')[2]*0.98,y=-0.7+0.04*12.7,adj=c(1,0),expression(paste(italic(R^2),' = 0.67')),cex=cext)
    text(x=par('usr')[1]+0.03*(par('usr')[2]-par('usr')[1]),
         y=0.95*(par('usr')[4]-par('usr')[3])+par('usr')[3],'a',cex=1.4,font=2)
  plot(feb1~dec1,pch=19,col='grey40',xlim=c(0.5,6.9),ylim=c(-0.5,8),axes=F,xlab='',ylab='',xaxs='i',yaxs='i')
    lines(cent.febdec~preddec,col='grey40')
    polygon(c(preddec,rev(preddec)),c(cri.febdec[1,],rev(cri.febdec[2,])),col=rgb(0,0,0,0.1),border=NA)
    axis(1,at=par('usr')[1:2],tcl=0,labels=F)
    axis(1,at=seq(2,6,2),tcl=-0.25,labels=seq(2,6,2),cex=0.9,mgp=c(0,0.6,0))
    axis(2,at=par('usr')[3:4],tcl=0,labels=F)
    axis(2,at=seq(0,8,2),tcl=-0.25,labels=seq(0,8,2),las=1,mgp=c(0,0.7,0))
    mtext('Late-winter population (ha)',side=2,line=2.0,cex=cex.mt,col='steelblue4')
    mtext('Early-winter population (ha)',side=1,line=1.9,cex=cex.mt,col='steelblue2')
    text(x=par('usr')[2]*0.98,y=-0.5+0.14*8.5,adj=c(1,0),'y = 0.33 + 0.58x',cex=cext)
    text(x=par('usr')[2]*0.98,y=-0.5+0.04*8.5,adj=c(1,0),expression(paste(italic(R^2),' = 0.43')),cex=cext)
    text(x=par('usr')[1]+0.03*(par('usr')[2]-par('usr')[1]),
         y=0.95*(par('usr')[4]-par('usr')[3])+par('usr')[3],'b',cex=1.4,font=2)
  plot(sum.md~feb,data=popdf[popdf$yr %in% 2005:2018,],type='n',xlim=c(0.6,6.5),ylim=c(-0.7,12),axes=F,xlab='',ylab='',xaxs='i',yaxs='i')
    arrows(x0=popdf$feb,x1=popdf$feb,y0=popdf$sum.lcl,y1=popdf$sum.ucl,length=0,col='gray50')
    points(sum.md~feb,data=popdf[popdf$yr %in% 2005:2018,],pch=19,col='gray40')
    lines(cent.sumfeb~predfeb,col='grey40')
    polygon(c(predfeb,rev(predfeb)),c(cri.sumfeb[1,],rev(cri.sumfeb[2,])),col=rgb(0,0,0,0.1),border=NA)
    axis(1,at=par('usr')[1:2],tcl=0,labels=F)
    axis(1,at=seq(2,6,2),tcl=-0.25,labels=seq(2,6,2),cex=0.9,mgp=c(0,0.6,0))
    axis(2,at=par('usr')[3:4],tcl=0,labels=F)
    axis(2,at=seq(0,12,3),tcl=-0.25,labels=seq(0,12,3),las=1,mgp=c(0,0.7,0))
    mtext('Summer population (mon/hr)',side=2,line=2.0,cex=cex.mt,col='salmon2')
    mtext('Late-winter population (ha)',side=1,line=1.9,cex.mt,col='steelblue4')
    text(par('usr')[2]*0.98,-0.7+0.14*12.7,adj=c(1,0),'y = 2.41 + 0.82x',cex=cext)
    text(par('usr')[2]*0.98,-0.7+0.04*12.7,adj=c(1,0),expression(paste(italic(R^2),' = 0.31')),cex=cext)
    text(x=par('usr')[1]+0.03*(par('usr')[2]-par('usr')[1]),
         y=0.95*(par('usr')[4]-par('usr')[3])+par('usr')[3],'c',cex=1.4,font=2)

#------------------------------------------------------------#
# Residuals from winter submodel (Ext Data Fig. 4)
#------------------------------------------------------------#  

  #Focusing only on those supercolony-year combinations when monarchs were present
    winterp <- winter[index_present,]
    mu_win <- posterior[,grep('mu_win',colnames(posterior))]
    mu.md <- apply(mu_win[,index_present],2,median)
    winterp$resid.md <- winterp$area - mu.md
  
  #Post-hoc regression (residual ~ yr)
    winterp$yr0 <- winterp$yr-2004
  	jagsdata <- list(nobs=nrow(winterp),
                     x=winterp$yr0,
                     y=winterp$resid.md)  
    ni <- 2000; na <- 1000; nb <- 8000; nc <- 3; ni.tot <- ni + nb
  	params <- c('b0','b1','sigma')
    inits <- lapply(1:nc, function(i)
             list(b0=runif(1,-2,2),
                  b1=runif(1,-2,2),
                  sigma=runif(1,1,10)))	
    
    fit <- jags(data=jagsdata, inits=inits, parameters.to.save=params,
                model.file='LinearRegression.txt',
                n.chains=nc, n.adapt=na, n.iter=ni.tot, n.burnin=nb,
                parallel=T, n.cores=3, DIC=F)
  	print(fit)
  	fits <- fit$samples
  	comb <- combine.mcmc(fits)
    predyr <- seq(min(winterp$yr0),max(winterp$yr0),length=100)
  	pred <- matrix(c(rep(1,100),predyr),ncol=2,byrow=FALSE) %*% t(comb[,1:2])
    cent <- apply(pred,1,median)
    cri <- apply(pred,1,quantile,probs=c(0.025,0.975))
    plotyr <- predyr + min(winterp$yr)
    
    #Plot residuals over time, with regression line (Extended Data Fig. 4)
    par(mfrow=c(1,1),mar=c(3,3,0.5,0.7),cex=0.8)
    plot(cent~plotyr,type='n',xlab='',ylab='',ylim=c(-0.8,1.5),axes=F,yaxs='i',bty='l')
      axis(1,at=par('usr')[1:2],tcl=0,labels=F)
      axis(1,at=seq(2005,2017,by=4),tcl=-0.25,labels=seq(2005,2017,by=4),cex=0.9,mgp=c(0,0.6,0))
      axis(2,at=par('usr')[3:4],tcl=0,labels=F)
      axis(2,at=seq(-0.5,1,by=0.5),tcl=-0.25,labels=seq(-0.5,1,by=0.5),las=1,mgp=c(0,0.7,0))
      polygon(c(plotyr,rev(plotyr)),c(cri[1,],rev(cri[2,])),col='gray80',border=NA)
      lines(cent~plotyr,col='gray50')
      points(resid.md~yr,data=winterp,pch=19,col='steelblue2')
      mtext('Residual',side=2,line=2.0,cex=0.9)
      mtext('Year',side=1,line=1.9,cex=0.9)
      
#---------------------------------------------------------------------------------------------------------------#
# Spatiotemporal variation in glyphosate use; relationships between glyphosate use and trends in monarch counts 
#---------------------------------------------------------------------------------------------------------------#  

  #Note: need to output "countyyr" parameter to evaluate county-level trends in expected monarch counts
  
  #Spatiotemporal variation in glyphosate use   
    gly.cy <- cov.cy[,c('yr','state.county','glyphosate')]
    gly.cy$yr <- as.factor(gly.cy$yr)
    gly.cy$state.county <- as.factor(gly.cy$state.county)
    #treating county and year as random effects (only one obs per yr:county combination)
    gly1 <- lmer(glyphosate ~ 1 + (1|yr) + (1|state.county), data=gly.cy)      
    summary(gly1) 
    #Percent of total variance associated with county (i.e., spatial variation)
    re.dat <- as.data.frame(VarCorr(gly1)) 
    re.dat[re.dat$grp=='state.county','vcov']/sum(re.dat[,'vcov'])*100  #74.3%
    
  #Calculate mean glyphosate use between 2004-2018 for each county
    gly <- ddply(cov.cy,.(state.county),summarize,gly=mean(glyphosate))
    
  #Extract countyyr parameters from fitted model
    cy.ind <- grep('countyyr',colnames(posterior))
    cy <- cyw[cyw$wk==21,c('county.ind','yr')]
    #Extract median of posterior distributions
    cy$pred.md <- apply(posterior[,cy.ind],2,median)
  
  #Attach county IDs, lat/long, and mean glyphosate use
    cy <- join(cy,cov.c[,c('county.ind','state.county','state.name','Xcentroid','Ycentroid')],
               by='county.ind',type='left')
    cy <- join(cy,gly,by='state.county',type='left')
    
  #Evaluating whether mean monarch counts vary with glyphosate application rates
    cy$gly.z <- (cy$gly-mean(cy$gly))/sd(cy$gly)
    cy$yr0 <- cy$yr - min(cy$yr)
    #Analysis that includes all counties
    summary(mean.gly <- lmer(pred.md ~ gly.z + (1 | state.county), data=cy))
    #Analysis that excludes counties where glyphosate use = 0
    summary(mean.gly1 <- lmer(pred.md ~ gly.z + (1 | state.county), data=cy[cy$gly>0,])) 
    
  #Evaluating whether temporal trends in monarch counts vary with glyphosate application rate 
  #(excluding counties with no glyphosate use)
    summary(trend.gly1 <- lmer(pred.md ~ gly.z + yr0 + gly.z*yr0 + (1 | state.county), data=cy[cy$gly>0,])) 
    #Negative trend in counts (year coefficient = -0.24)
    #But no evidence that trends vary among counties depending on glyphosate use (interaction coefficient = -0.04, SE = 0.05)
    