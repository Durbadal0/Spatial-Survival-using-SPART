stan_code <- "data {
  int<lower=0> N;        // Number of observations
  int<lower=0> n;       // Number of counties
  real<lower=0> time[N]; // Time to event
  int<lower=1> p; // number of covariates
  int<lower=0, upper=1> death[N]; // Censoring indicator (1 for censored, 0 for uncensored)
   real x[N,p];     // Covariate matrix
   int<lower=0> L;
  real<lower=0> a[L];         // Interval breakpoints
  int<lower=1,upper=n> county[N];
  matrix[n, n] D;                // Data matrix
  matrix[n, n] A;
  real v1;
  real v2;
  
}

parameters {
   
  real b[L,p];              // Hazard coefficients
   vector[n] W;
  real<lower=0> sigma;           // Standard deviation of the multivariate normal prior
  real<lower=1/v1, upper=1/v2> rho;  // Correlation parameter
}



model {
  
  
   for(l in 1:L){
   for(k in 1:p){
  b[l,k] ~ normal(0, 1);
   }
 }
 
 W ~ multi_normal(rep_vector(0, n), sigma^2 * (D - rho * A));
  sigma ~ gamma(0.001, 0.001);
  // You can adjust these limits based on your knowledge
  rho ~ uniform(1/v1, 1/v2);
  
  for (i in 1:N) {
    int interval = 0;

    for (l in 1:L) {
      if (time[i] >= a[l] && time[i] < a[l + 1]) {
        interval = l;
     
      }
    }

   real hazard = 0;
   real interval_term = 0;
    if(interval >1){
    for (l in 1:interval-1) {
     
      
      for (k in 1:p) {
        interval_term += x[i, k] * b[l, k];
      }
      hazard += exp(interval_term)*(a[l+1]-a[l]);
    }


    } 
    real interval_term2 = 0;
    for (k in 1:p) {
        interval_term2 += x[i, k] * b[interval, k];
      }
     hazard+= exp(W[county[i]])* exp(interval_term2)*(time[i]-a[interval]);    
         
    if (death[i] == 0) {
      target +=  -hazard; // Using the CDF of normal distribution
    } else {
      target += W[county[i]]+interval_term2 - hazard;
    }
  }


}

  


generated quantities {
  // You can add any additional quantities you want to compute here.
}
"






data=Surv_data[which(Surv_data$Race==1 & Surv_data$Stage==2),]
stan_data <- list(
  N = nrow(data),  # Number of observations
  n=67,
  p = 6,
  time = data$as.numeric.date_diff.,  # Time to event
  death = data$death,  # Censoring indicator
  x = data[,c(3,5,6,7,8,10)],
  L=101,# Assuming your covariate data is in columns 4 to 103
  a = seq(0,6000,60),  # Vector of interval breakpoints
  county=as.integer(data$county),
  D = D,    # Replace with the D matrix
  A = Adj.M,  
  v1=v1,
  v2=v2# Replace with the A matrix
)

fit <- stan(model_code = stan_code, data = stan_data,
            chains = 1,
            iter = 2000,
            warmup = 1000,
            thin = 1)




# Extract posterior samples
posterior_samples <- rstan::extract(fit)

# View the summary of posterior samples
summary(posterior_samples)




trace_plot <- stan_trace(fit)  # Replace with your posterior samples
plot(trace_plot)



rhats <- stan_rhat(fit)  # Replace with your posterior samples


geweke_test <- rstan::geweke(stan_model)  # Replace with the parameter of interest


# Assuming 'posterior_samples' is a matrix or list of MCMC samples
mcmc_samples <- as.mcmc.list(posterior_samples)

acf(posterior_samples$b)


dim(posterior_samples$p)

smu_hat_W=colMeans(posterior_samples$p)

posterior_samples_SMU_A=extract(stan_model_A)

posterior_samples_SMU_W=extract(stan_model_W)

save(posterior_samples_SMU_A,file="posterior_samples_SMU_A.RData")
save(posterior_samples_SMU_W,file="posterior_samples_SMU_W.RData")
