data {
  int<lower=0> n;
  int<lower=0> p;
  vector[n] y;
  matrix[n,p] X;
  array[n] vector[2] coords;
  
  vector<lower=0>[p] scale_theta;
  real<lower=0> scale_sigma;
  real<lower=0> scale_tau;
  
  real<lower = 0> a; // Shape parameters in the prior for ell
  real<lower = 0> b; // Scale parameters in the prior for ell
  
}

transformed data{
  
}

parameters {
  vector[p] theta_std;
  real<lower=0> ell;
  real<lower=0> sigma_std;
  real<lower=0> tau_std;
}

transformed parameters{
  vector[p] theta = scale_theta .* theta_std;
  real sigma = scale_sigma * sigma_std;
  real tau = scale_sigma * tau_std;
}

model {
  theta_std ~ std_normal();
  ell ~ inv_gamma(a,b);
  sigma_std ~ std_normal();
  tau_std ~ std_normal();
  vector[n] mu = X*theta;
  matrix[n,n] Sigma = gp_matern32_cov(coords, sigma, ell);
  matrix[n,n] L = cholesky_decompose(add_diag(Sigma, square(tau)));
  y ~ multi_normal_cholesky(mu, L);
}


