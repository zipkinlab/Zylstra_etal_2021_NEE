data {

	// dimensions of data
	int<lower=1> n_counties;    // number of counties included in summer breeding region
	int<lower=1> n_years;       // number of years
	int<lower=1> n_cyw;         // n_counties * n_years * n_weeks
	int<lower=1> n_surveys;     // number of summer surveys
	int<lower=1> n_colonies;    // number of supercolonies
	int<lower=1> n_obs;         // number of supercolony-year combinations

	// indices
	int<lower=1,upper=n_years>    year_id[n_cyw];            // year id (1:15)
	int<lower=1,upper=n_cyw>      cyw_id[n_surveys];         // id for colony-year-week combination
	int<lower=1,upper=n_colonies> colony_id[n_obs];          // supercolony id (1:13)
	int<lower=1>                  index_present[140];        // index indicating when monarchs were present in supercolonies

	// number of regression coefficients
	int<lower=1> n_cov_alpha;   // number of fixed effects in county-level model (summer)
	int<lower=1> n_cov_beta;    // number of fixed effects in survey level model (summer)
	int<lower=1> n_cov_gamma;   // number of fixed effects in winter (gamma) model MINUS 1 for summer index
	
	// covariate data
	matrix[n_cyw,n_cov_alpha] X_county;    // matrix of covariates for summer county-level model (all standardized)
	vector[n_cyw] week_st;                 // week, in county-level model (standardized)
	matrix[n_surveys,n_cov_beta] X_survey; // matrix of covariates for summer survey-level model (all standardized)
	vector<lower=0>[n_surveys] effort;     // party hours spent on each survey, scaled by the mean
	matrix[n_obs,n_cov_gamma] X_winter;    // matrix of covariates for winter model (all standardized)

	// response data
	int<lower=0> y_count[n_surveys];   // monarch count on each survey
	vector<lower=0>[n_obs] area;       // area occupied by monarchs in each supercolony, each year

}

parameters {
	
	// regression coefficients
	vector[n_cov_alpha] alphaFE;   // betas associated with fixed effects in county-level model
	vector[n_years] alpha_week;    // mean effect of week
	vector[n_years] alpha_week2;   // mean effect of week^2
	vector[n_cov_beta] betaFE;     // betas associated with fixed effects in survey-level model
	vector[n_cov_gamma] gammaFE;   // betas associated with fixed effects in winter model (not including coef for summer index)

	// intercepts
	real alpha0;                                  // intercept in county-level model (mean expected count on NABA survey, on log scale)
	vector<lower=0,upper=1>[n_colonies] pocc_mn;  // intercepts in winter, logistic model (mean probability monarchs present at each colony)
	vector[n_colonies] gamma0;                    // intercepts in winter (different for each colony), gamma model (mean area occupied when monarchs present, on log scale)
	
	// overdispersion parameter for negative binomial
	real <lower=0> phi;  
	
	// shape parameter in Gamma distribution in winter model
	real<lower=0> shape;   	

}

transformed parameters {

	// expected count (on an average NABA survey with average effort) in each county, year, week
	vector<lower=0>[n_cyw] mu_county; 

	// expected count on each survey (as a function of survey type, effort, local land use)
	vector<lower=0>[n_surveys] lambda;
	
	// probability that monarchs are present in supercolony each year 
	vector<lower=0,upper=1>[n_obs] pocc; 
	
	// mean area occupied in December, conditional on monarch presence
	vector<lower=0>[n_obs] mu_win;
	
	// rate parameter for Gamma model (winter)
	vector<lower=0>[n_obs] rate;

	// summer: county-level model
	for(i in 1:n_cyw)
		mu_county[i] = exp(alpha0 + alpha_week[year_id[i]] * week_st[i]
		               + alpha_week2[year_id[i]] * week_st[i] * week_st[i]
		               + X_county[i] * alphaFE); 
									 
	// summer: survey-level model
	for(i in 1:n_surveys)
		lambda[i] = exp(log(mu_county[cyw_id[i]]) + log(effort[i]) 
		            + X_survey[i] * betaFE);
	
	// winter: logistic part of hurdle model
	for(i in 1:n_obs)
		pocc[i] = pocc_mn[colony_id[i]];
		
	for(i in 1:n_obs)
		mu_win[i] = exp(gamma0[colony_id[i]] + X_winter[i] * gammaFE);		

	rate = shape ./ mu_win;
}

model {
	
  // priors on intercepts and regression coefficients (implicit U(0,1) prior on pocc_mn) 
  target += normal_lpdf(alphaFE       | 0, sqrt(1000));
  target += normal_lpdf(betaFE        | 0, sqrt(1000));
  target += normal_lpdf(gammaFE       | 0, sqrt(1000));
	target += normal_lpdf(alpha0        | 0, sqrt(1000));	
  target += normal_lpdf(alpha_week    | 0, sqrt(1000));
  target += normal_lpdf(alpha_week2   | 0, sqrt(1000)); 
	target += normal_lpdf(gamma0        | 0, sqrt(1000));
	
	// prior for overdisperion parameter in neg-binom
	target += uniform_lpdf(phi | 0, 20);
	
	// prior for shape parameter in Gamma model (winter)
	target += gamma_lpdf(shape | 2, 2);
	
	// modeling counts in summer (Negative binomial distribution)
	target += neg_binomial_2_lpmf(y_count | lambda, phi);	
	
	// modeling area occupied in winter (Gamma hurdle model)
	for(i in 1:n_obs)
	{
		if (area[i] == 0)
			target += log1m(pocc[i]);
		else
			target += log(pocc[i]) + gamma_lpdf(area[i] | shape, rate[i]);
	}
}

generated quantities {

	vector[n_surveys] log_lik_summer;
	real sum_log_lik;
	vector[n_obs] log_lik_winter;
	real win_log_lik;
	vector<lower=0>[n_surveys] y_pred;
	vector<lower=0>[n_surveys] sccount_pred;
	real<lower=0> sccount_pred_mn;
	real<lower=0> sccount_pred_sd;
	vector<lower=0>[n_obs] area_pred;
	real<lower=0> area_pred_mn;
	real<lower=0> area_pred_sd;

	// calculating log-likelihood for model comparisons
	for(i in 1:n_surveys)
	{
		log_lik_summer[i] = neg_binomial_2_lpmf(y_count[i] | lambda[i], phi);
	}
	for(i in 1:n_obs)
	{
		if (area[i] == 0)
			log_lik_winter[i] = log1m(pocc[i]);
		else
			log_lik_winter[i] = log(pocc[i]) + gamma_lpdf(area[i] | shape, rate[i]);
	}
	sum_log_lik = sum(log_lik_summer);
	win_log_lik = sum(log_lik_winter);
		
	// fit stats: generating new counts (summer model)
	for(i in 1:n_surveys)
	{
		if(lambda[i]>150){
			y_pred[i] = neg_binomial_2_rng(150, phi);
	  } else {
			y_pred[i] = neg_binomial_2_rng(lambda[i], phi);
		}
	}
	sccount_pred = y_pred ./ effort;
	sccount_pred_mn = mean(sccount_pred);
	sccount_pred_sd = sd(sccount_pred);
	
	// fit stats: generating new areas occupied (winter model)
	// only comparing areas for sites where monarchs were observed 2004-2018 
	for(i in 1:n_obs)
	{
		area_pred[i] = gamma_rng(shape, rate[i]);
	}
	area_pred_mn = mean(area_pred[index_present]);
	area_pred_sd = sd(area_pred[index_present]);
	
}


