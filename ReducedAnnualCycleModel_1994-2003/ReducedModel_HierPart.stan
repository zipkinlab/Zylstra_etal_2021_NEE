data {

	// dimensions of data
	int<lower=1> n_counties;    // number of counties included in summer breeding region
	int<lower=1> n_years;       // number of years
	int<lower=1> n_cyw;         // n_counties * n_years * n_weeks
	int<lower=1> n_surveys;     // number of summer surveys

	// indices
	int<lower=1,upper=n_years>    year_id[n_cyw];            // year id (1:15)
	int<lower=1,upper=n_cyw>      cyw_id[n_surveys];         // id for colony-year-week combination

	// number of regression coefficients
	int<lower=1> n_cov_alpha;   // number of fixed effects in county-level model (summer)
	int<lower=1> n_cov_beta;    // number of fixed effects in survey level model (summer)
	
	// covariate data
	matrix[n_cyw,n_cov_alpha] X_county;    // matrix of covariates for summer county-level model (all standardized)
	vector[n_cyw] week_st;                 // week, in county-level model (standardized)
	matrix[n_surveys,n_cov_beta] X_survey; // matrix of covariates for summer survey-level model (all standardized)
	vector<lower=0>[n_surveys] effort;     // party hours spent on each survey, scaled by the mean

	// response data
	int<lower=0> y_count[n_surveys];   // monarch count on each survey

}

parameters {
	
	// regression coefficients
	vector[n_cov_alpha] alphaFE;   // betas associated with fixed effects in county-level model
	vector[n_years] alpha_week;    // mean effect of week
	vector[n_years] alpha_week2;   // mean effect of week^2
	vector[n_cov_beta] betaFE;     // betas associated with fixed effects in survey-level model

	// intercepts
	real alpha0;                    // intercept in county-level model (mean expected count on NABA survey, on log scale)
	
	// overdispersion parameter for negative binomial
	real <lower=0> phi;  
	
}

transformed parameters {

	// expected count (on an average NABA survey with average effort) in each county, year, week
	vector<lower=0>[n_cyw] mu_county; 

	// expected count on each survey (as a function of survey type, effort, local land use)
	vector<lower=0>[n_surveys] lambda;
	
	// Summer: County-level model
	for(i in 1:n_cyw)
		mu_county[i] = exp(alpha0 + alpha_week[year_id[i]] * week_st[i]
		               + alpha_week2[year_id[i]] * week_st[i] * week_st[i]
		               + X_county[i] * alphaFE); 
									 
	// Summer: Survey-level model
	for(i in 1:n_surveys)
		lambda[i] = exp(log(mu_county[cyw_id[i]]) + log(effort[i]) 
		            + X_survey[i] * betaFE);
	
}

model {

	// priors on variance components
	
	// priors on intercepts and regression coefficients (implicit U(0,1) prior on pocc_mn) 
	target += normal_lpdf(alphaFE       | 0, sqrt(1000));
	target += normal_lpdf(betaFE        | 0, sqrt(1000));
	target += normal_lpdf(alpha0        | 0, sqrt(1000));	
	target += normal_lpdf(alpha_week    | 0, sqrt(1000));
	target += normal_lpdf(alpha_week2   | 0, sqrt(1000)); 

	// prior for overdisperion parameter in neg-binom
	target += uniform_lpdf(phi | 0, 20);
	
	// Modeling counts in summer (Negative binomial distribution)
	target += neg_binomial_2_lpmf(y_count | lambda, phi);	
	
}

generated quantities {

	vector[n_surveys] log_lik_summer;
	real sum_log_lik;
	vector<lower=0>[n_surveys] y_pred;
	vector<lower=0,upper=1>[n_surveys] y_pred_0;
	real<lower=0,upper=1> y_pred_prop0; 
	vector<lower=0>[n_surveys] sccount_pred;
	real<lower=0> sccount_pred_mn;
	real<lower=0> sccount_pred_max;
	real<lower=0> sccount_pred_sd;

	// calculating log-likelihood for model comparisons
	for(i in 1:n_surveys)
	{
		log_lik_summer[i] = neg_binomial_2_lpmf(y_count[i] | lambda[i], phi);
	}
	sum_log_lik = sum(log_lik_summer);
		
	// fit stats: generating new counts (summer model)  // adding a limit for lambda to avoid super high predicted values
	for(i in 1:n_surveys)
	{
		if(lambda[i]>150){
			y_pred[i] = neg_binomial_2_rng(150, phi);
	  } else {
			y_pred[i] = neg_binomial_2_rng(lambda[i], phi);
		}
		y_pred_0[i] = !(y_pred[i]);
	}
	y_pred_prop0 = sum(y_pred_0) / n_surveys;
	sccount_pred = y_pred ./ effort;
	sccount_pred_mn = mean(sccount_pred);
	sccount_pred_max = max(sccount_pred);
	sccount_pred_sd = sd(sccount_pred);
	
}


