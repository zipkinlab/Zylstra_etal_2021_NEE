data {

	// dimensions of data
	int<lower=1> n_counties;    // number of counties included in summer breeding region
	int<lower=1> n_years;       // number of years
	int<lower=1> n_cyw;         // n_counties * n_years * n_weeks
	int<lower=1> n_surveys;     // number of summer surveys
	int<lower=1> n_peak;        // number of county-year-week combinations during peak
	int<lower=1> n_cy;          // n_counties * n_years	
	int<lower=1> n_colonies;    // number of metacolonies
	int<lower=1> n_obs;         // number of metacolony-year combinations

	// indices
	int<lower=1,upper=n_years>    year_id[n_cyw];            // year id (1:15)
	int<lower=1,upper=n_cyw>      cyw_id[n_surveys];         // id for colony-year-week combination
	int<lower=1,upper=n_colonies> colony_id[n_obs];          // metacolony id (1:13)
	int<lower=1>                  ind1[n_counties*n_years];  // used to calculate summer pop index
	int<lower=1>                  ind2[n_years];             // used to calculate summer pop index
	int<lower=1>                  start_peak;                // first row of X_county/week_st that has data for summer peak (wks 6-9)
	int<lower=1>                  index_present[140];        // index indicating when monarchs were present in metacolonies

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
	vector<lower=0>[n_obs] area;       // area occupied by monarchs in each metacolony, each year
	
	// county weights (for calculating the summer index)
	vector<lower=0,upper=1>[n_counties] weights;  // based on non-forested area

}

parameters {
	
	// regression coefficients
	vector[n_cov_alpha] alphaFE;   // betas associated with fixed effects in county-level model
	vector[n_years] alpha_week;    // mean effect of week
	vector[n_years] alpha_week2;   // mean effect of week^2
	vector[n_cov_beta] betaFE;     // betas associated with fixed effects in survey-level model
	real gamma_sum;                // effect of summer population size

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
	
	// probability that monarchs are present in metacolony each year 
	vector<lower=0,upper=1>[n_obs] pocc; 
	
	// mean area occupied in December, conditional on monarch presence
	vector<lower=0>[n_obs] mu_win;
	
	// rate parameter for Gamma model (winter)
	vector<lower=0>[n_obs] rate;
	
	// objects associated with index of peak summer population size
	vector<lower=0>[n_peak] wk69;
	vector<lower=0>[n_counties*n_years] countyyr;
	vector<lower=0>[n_years] pred_orig;
	vector[n_years] pred_sum;
	vector[n_years*n_colonies] pred_sumC;

	// Summer: County-level model
	for(i in 1:n_cyw)
		mu_county[i] = exp(alpha0 + alpha_week[year_id[i]] * week_st[i]
		               + alpha_week2[year_id[i]] * week_st[i] * week_st[i]
		               + X_county[i] * alphaFE); 
									 
	// Getting index of "peak summer abundance"
	wk69 = segment(mu_county,start_peak,n_peak); // extracting weeks 6-9 (start_peak is first row with data from weeks 6-9)
	
	for(i in 1:n_cy)  // get means for each county-year combination
		countyyr[i] = mean(segment(wk69,ind1[i],4));

	for(i in 1:n_years)  // get weighted mean of counties in each year
		pred_orig[i] = sum(weights .* segment(countyyr,ind2[i],n_counties));

	pred_sum = (pred_orig - 15.7) ./ 8.5;

	// Summer: Survey-level model
	for(i in 1:n_surveys)
		lambda[i] = exp(log(mu_county[cyw_id[i]]) + log(effort[i]) 
		            + X_survey[i] * betaFE);
	
	// Winter: Logistic part of hurdle model
	for(i in 1:n_obs)
		pocc[i] = pocc_mn[colony_id[i]];
		
	// Winter: Gamma part of hurdle model
	pred_sumC = append_row(pred_sum,(append_row(pred_sum,(append_row(pred_sum,(append_row(pred_sum,
	            append_row(pred_sum,(append_row(pred_sum,(append_row(pred_sum,(append_row(pred_sum,
				      append_row(pred_sum,(append_row(pred_sum,(append_row(pred_sum,(append_row(pred_sum,pred_sum)))))))))))))))))))));

	for(i in 1:n_obs)
		mu_win[i] = exp(gamma0[colony_id[i]] + gamma_sum * pred_sumC[i]);		

	rate = shape ./ mu_win;
}

model {
	
  // priors on intercepts and regression coefficients (implicit U(0,1) prior on pocc_mn) 
  target += normal_lpdf(alphaFE       | 0, sqrt(1000));
  target += normal_lpdf(betaFE        | 0, sqrt(1000));
	target += normal_lpdf(alpha0        | 0, sqrt(1000));	
  target += normal_lpdf(alpha_week    | 0, sqrt(1000));
  target += normal_lpdf(alpha_week2   | 0, sqrt(1000)); 
	target += normal_lpdf(gamma0        | 0, sqrt(1000));
	target += normal_lpdf(gamma_sum     | 0, sqrt(1000));
	
	// prior for overdisperion parameter in neg-binom
	target += uniform_lpdf(phi | 0, 20);
	
	// prior for shape parameter in Gamma model (winter)
	target += gamma_lpdf(shape | 2, 2);
	
	// Modeling counts in summer (Negative binomial distribution)
	target += neg_binomial_2_lpmf(y_count | lambda, phi);	
	
	// Modeling area occupied in winter (Gamma hurdle model)
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
	vector<lower=0,upper=1>[n_surveys] y_pred_0;
	real<lower=0,upper=1> y_pred_prop0; 
	vector<lower=0>[n_surveys] sccount_pred;
	real<lower=0> sccount_pred_mn;
	real<lower=0> sccount_pred_max;
	real<lower=0> sccount_pred_sd;
	vector<lower=0,upper=1>[n_obs] p_unocc;
	vector<lower=0,upper=1>[n_obs] area_pred_0;
	real<lower=0,upper=1> area_pred_prop0;
	vector<lower=0>[n_obs] area_pred;
	real<lower=0> area_pred_mn;
	real<lower=0> area_pred_max;
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
		y_pred_0[i] = !(y_pred[i]);
	}
	y_pred_prop0 = sum(y_pred_0) / n_surveys;
	sccount_pred = y_pred ./ effort;
	sccount_pred_mn = mean(sccount_pred);
	sccount_pred_max = max(sccount_pred);
	sccount_pred_sd = sd(sccount_pred);
	
	// fit stats: generating new areas occupied (winter model)
	// only comparing areas for sites where monarchs were observed 2004-2018 (can't have dynamic vector length)
	p_unocc = 1 - pocc;
	for(i in 1:n_obs)
	{
		area_pred_0[i] = bernoulli_rng(p_unocc[i]);
		area_pred[i] = gamma_rng(shape, rate[i]);
	}
	area_pred_prop0 = sum(area_pred_0) / n_obs;		
	area_pred_mn = mean(area_pred[index_present]);
	area_pred_max = max(area_pred[index_present]);
	area_pred_sd = sd(area_pred[index_present]);
	
}


