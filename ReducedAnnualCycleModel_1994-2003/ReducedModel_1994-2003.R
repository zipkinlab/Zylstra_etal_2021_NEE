###################################################################################################################################
#Reduced annual-cycle model for monarch butterflies in eastern North America
  #Years included: 1994-2003
  #Model for inference that includes all winter, spring, summer covariates
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
  yr.min <- 1994
  yr.max <- 2003

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
# Standardize annual covariates for winter model
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
# Calculate observed values for PPC
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
  X_survey <- as.matrix(summer[,c('il.ind','oh.ind','openS.st')]) 
  
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
    ni <- 3000      # No. iterations (including warmup)
    nb <- 500       # No. burn-in iterations to discard (ie, warmup)
    nt <- 1         # Thin rate
    nc <- 3         # No. chains
  
#Initial values
  set.seed(123)
  inits <- lapply(1:nc, function(i)
           list(alpha0=runif(1,1,4),
                alphaFE=runif(ncol(X_county),-1,1),
                alphaRE_week=runif(1,-1,1),
                alphaRE_week2=runif(1,-1,1),
                betaFE=runif(ncol(X_survey),-1,1),
                gamma0=runif(1,1,3),
                gamma_sum=runif(1,-0.5,0.5),
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
              'sccount_pred_mn','sccount_pred_sd','area_pred_mn','area_pred_sd',
              'sum_log_lik','win_log_lik','mu_win')

#------------------------------------------------------------#
# Call STAN from R
#------------------------------------------------------------#
  
  out <- stan('ReducedModel_1994-2003.stan',
              control=list(adapt_delta=0.8),
              data=standata, init=inits, pars=params,
              chains=nc, iter=ni, warmup=nb, thin=nt,
              seed=1, cores=3, open_progress=FALSE)
  
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
    
  #Mean of area occupied in winter
    area.mn.post <- posterior[,'area_pred_mn']
    # hist(area.mn.post,breaks=25,col='gray80')
    # abline(v=area.mn,col='dodgerblue3',lwd=2)
    area.mn.p <- sum(area.mn.post>=area.mn)/length(area.mn.post)

  #SD of area occupied in winter
    area.sd.post <- posterior[,'area_pred_sd']
    # hist(area.sd.post,breaks=25,col='gray80')
    # abline(v=area.sd,col='dodgerblue3',lwd=2)
    area.sd.p <- sum(area.sd.post>=area.sd)/length(area.sd.post)
    
#--------------------------------------------------------------------#
# Marginal effects of summer temperatures (GDD; Extended Data Fig. 5a)
#--------------------------------------------------------------------#
  
  #Posterior samples from relevant parameters  
    amat.gdd <- posterior[,c('alpha0','alphaRE_week','alphaRE_week2','alphaFE[5]','alphaFE[6]','alphaFE[7]','alphaFE[8]')]
   
  #Standardized week values
    wk.r <- sort(unique(cyw$wk.st))

  #Predict effects of diffGDD for a cool/average/warm county
    avgGDD.q <- as.vector(quantile(cyw$avgGDD,probs=c(0.1,0.5,0.9)))
    (avgGDD.q.orig <- avgGDD.q*sd(cov.c$avgGDD.9403) + mean(cov.c$avgGDD.9403))  #711, 898, 1033
    
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
    
    plot.sugdd <- X.sugdd$diffGDD[1:100]*sd(cov.cw$diffGDD.9403) + mean(cov.cw$diffGDD.9403)
      #Note: want to restrict predictions to range of diffGDD values observed at cool/average/warm locations in week 21
      #Using range of GDD differences in week 21 for sites with avgGDD within 5 degree days of cool/average/warm quantiles
      cov.cw <- join(cov.cw,cov.c[,c('state.county','avgGDD.9403')],by='state.county',type='left')
      cold.diff <- range(cov.cw$diffGDD.9403[cov.cw$wk==21 & abs(cov.cw$avgGDD.9403-avgGDD.q.orig[1])<5])
      avg.diff <- range(cov.cw$diffGDD.9403[cov.cw$wk==21 & abs(cov.cw$avgGDD.9403-avgGDD.q.orig[2])<5])
      hot.diff <- range(cov.cw$diffGDD.9403[cov.cw$wk==21 & abs(cov.cw$avgGDD.9403-avgGDD.q.orig[3])<5])  
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
    
    #Plot marginal effects of summer GDD (Extended Data Fig. 5a)
    plot(cent.sugdd1[cold.ind]~plot.sugdd[cold.ind],type='l',lty=1,bty='l',lwd=1.2,col=col3.1,
         xlim=c(-120,160),ylim=c(-0.8,20),axes=F,xlab='',ylab='',xaxs='i',yaxs='i')
      axis(1,at=c(par('usr')[1],160),tcl=0,labels=F)
      axis(1,at=seq(-100,100,by=100),tcl=-0.25,labels=seq(-100,100,by=100))
      axis(2,at=par('usr')[3:4],tcl=0,labels=F)
      axis(2,at=seq(0,20,by=5),tcl=-0.25,labels=seq(0,20,by=5),las=1)  
      polygon(c(plot.sugdd[cold.ind],rev(plot.sugdd[cold.ind])),c(cri.sugdd1[1,cold.ind],rev(cri.sugdd1[2,cold.ind])),col=col3.1p,border=NA)
      lines(cent.sugdd2[avg.ind]~plot.sugdd[avg.ind],type='l',lty=2,bty='l',lwd=1.2,col=col3.2)
      polygon(c(plot.sugdd[avg.ind],rev(plot.sugdd[avg.ind])),c(cri.sugdd2[1,avg.ind],rev(cri.sugdd2[2,avg.ind])),col=col3.2p,border=NA)
      lines(cent.sugdd3[hot.ind]~plot.sugdd[hot.ind],type='l',lty=3,bty='l',lwd=1.2,col=col3.3)
      polygon(c(plot.sugdd[hot.ind],rev(plot.sugdd[hot.ind])),c(cri.sugdd3[1,hot.ind],rev(cri.sugdd3[2,hot.ind])),col=col3.3p,border=NA)
      legend(x=-80,y=20,legend=c('Cool','Average','Warm'),lty=c(1,2,3),col=c(col3.1,col3.2,col3.3),bty='n',cex=0.9)  
      mtext(expression(paste("GDD (",degree,"C), deviation")),side=1,line=2,cex=0.9)
      mtext('Monarchs count (adults/hr)',side=2,las=0,line=1.8,cex=0.9)    
      
#--------------------------------------------------------------------#
# Marginal effects of summer precipitation (Extended Data Fig. 5b)
#--------------------------------------------------------------------#

  #Posterior samples from relevant parameters  
    amat.pcp <- posterior[,c('alpha0','alphaRE_week','alphaRE_week2','alphaFE[9]','alphaFE[10]','alphaFE[11]','alphaFE[12]')]

  #Predict effects of diffPCP for a dry/average/wet county
    avgPCP.q <- as.vector(quantile(cyw$avgPCP,probs=c(0.1,0.5,0.9)))
    (avgPCP.q.orig <- avgPCP.q*sd(cov.c$avgPCP.9403) + mean(cov.c$avgPCP.9403))  #422, 525, 578
      
    #Create sequence of standardized diffPCP values
    diffPCP.r <- round(range(cyw$diffPCP),1)
    diffPCP <- rep(seq(diffPCP.r[1],diffPCP.r[2],length=100),3)
    
    #Reference avgPCP values on standardized scale
    avgPCP <- rep(avgPCP.q,each=100)
      
    X.supcp <- data.frame(int=1,wk=wk.r[6],wk2=wk.r[6]^2,avgPCP=avgPCP,diffPCP=diffPCP,diffPCP2=diffPCP*diffPCP,diffavgPCP=diffPCP*avgPCP)
    predl.supcp <- as.matrix(X.supcp) %*% t(amat.pcp)
    pred.supcp.avgeff <- exp(predl.supcp)
    pred.supcp <- pred.supcp.avgeff/effortmean
    #Calculate medians and credible intervals
    cent.supcp1 <- apply(pred.supcp[1:100,],1,median)
    cri.supcp1 <- apply(pred.supcp[1:100,],1,quantile,c(0.025,0.975))
    cent.supcp2 <- apply(pred.supcp[101:200,],1,median)
    cri.supcp2 <- apply(pred.supcp[101:200,],1,quantile,c(0.025,0.975))
    cent.supcp3 <- apply(pred.supcp[201:300,],1,median)
    cri.supcp3 <- apply(pred.supcp[201:300,],1,quantile,c(0.025,0.975))    

    plot.supcp <- X.supcp$diffPCP[1:100]*sd(cov.cy$diffPCP.9403) + mean(cov.cy$diffPCP.9403)
      #Note: want to restrict predictions to range of diffPCP values observed at dry/avg/wet locations
      #Using range of PCP differences for sites with avgPCP within 10 mm of dry/avg/wet quantiles
      cov.cy <- join(cov.cy,cov.c[,c('state.county','avgPCP.9403')],by='state.county',type='left')
      dry.diff <- range(cov.cy$diffPCP.9403[abs(cov.cy$avgPCP.9403-avgPCP.q.orig[1])<10])
      avgp.diff <- range(cov.cy$diffPCP.9403[abs(cov.cy$avgPCP.9403-avgPCP.q.orig[2])<10])
      wet.diff <- range(cov.cy$diffPCP.9403[abs(cov.cy$avgPCP.9403-avgPCP.q.orig[3])<10])  
      dry.ind <- which(plot.supcp>=dry.diff[1] & plot.supcp<=dry.diff[2])
      avgp.ind <- which(plot.supcp>=avgp.diff[1] & plot.supcp<=avgp.diff[2])
      wet.ind <- which(plot.supcp>=wet.diff[1] & plot.supcp<=wet.diff[2])
    
    #Plot marginal effects of summer precipitation (Extended Data Fig. 5b)
    plot(cent.supcp1[dry.ind]~plot.supcp[dry.ind],type='l',lty=1,bty='l',lwd=1.2,col=col3.1,
         xlim=c(-250,325),ylim=c(-0.8,20),axes=F,xlab='',ylab='',xaxs='i',yaxs='i')
      axis(1,at=c(par('usr')[1],325),tcl=0,labels=F)
      axis(1,at=seq(-200,200,by=200),tcl=-0.25,labels=seq(-200,200,by=200))
      axis(2,at=par('usr')[3:4],tcl=0,labels=F)
      axis(2,at=seq(0,20,by=5),tcl=-0.25,labels=seq(0,20,by=5),las=1)    
      polygon(c(plot.supcp[dry.ind],rev(plot.supcp[dry.ind])),c(cri.supcp1[1,dry.ind],rev(cri.supcp1[2,dry.ind])),col=col3.1p,border=NA)
      lines(cent.supcp2[avgp.ind]~plot.supcp[avgp.ind],type='l',lty=2,bty='l',lwd=1.2,col=col3.2)
      polygon(c(plot.supcp[avgp.ind],rev(plot.supcp[avgp.ind])),c(cri.supcp2[1,avgp.ind],rev(cri.supcp2[2,avgp.ind])),col=col3.2p,border=NA)
      lines(cent.supcp3[wet.ind]~plot.supcp[wet.ind],type='l',lty=3,bty='l',lwd=1.2,col=col3.3)
      polygon(c(plot.supcp[wet.ind],rev(plot.supcp[wet.ind])),c(cri.supcp3[1,wet.ind],rev(cri.supcp3[2,wet.ind])),col=col3.3p,border=NA)
      legend(x=-150,y=20,legend=c('Dry','Average','Wet'),lty=c(1,2,3),col=c(col3.1,col3.2,col3.3),bty='n',cex=0.9)  
      mtext("Precipitation (mm), difference",side=1,line=1.8,cex=0.9) 
      mtext('Monarchs count (adults/hr)',side=2,las=0,line=1.8,cex=0.9)         
