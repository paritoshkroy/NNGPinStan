data {
  int<lower=0> n;
  int<lower=0> p;
  vector[n] y;
  matrix[n,p] X;
  array[n] vector[2] coords;
  
  vector<lower=0>[p] scale_theta;
  real<lower=0> scale_sigma;
  real<lower=0> scale_tau;
  
  real<lower=0> a;
  real<lower=0> b;
}

transformed data{
  
}

parameters {
  vector[p] theta_std;
  real<lower=0> ell;
  real<lower=0> sigma_std;
  real<lower=0> tau_std;
  vector[n] noise;
}

transformed parameters{
  vector[p] theta = scale_theta .* theta_std;
  real sigma = scale_sigma * sigma_std;
  real tau = scale_sigma * tau_std;
  //vector[n] z = cholesky_decompose(add_diag(gp_matern32_cov(coords, sigma, ell), 1e-8)) * noise;
  vector[n] z = cholesky_decompose(add_diag(gp_exponential_cov(coords, sigma, ell), 1e-8)) * noise;
  
  //matrix[n,n] Sigma = gp_exponential_cov(coords, sigma, phi);
  //matrix[n,n] V = add_diag(V, 1e-8);
  //matrix[n,n] L = cholesky_decompose(V);
  //vector[n] z = L * noise;
}

model {
  theta_std ~ std_normal();
  ell ~ inv_gamma(a,b);
  sigma_std ~ std_normal();
  tau_std ~ std_normal();
  noise ~ std_normal();
  vector[n] mu = X*theta;
  
  y ~ normal(mu + z, tau);
}

generated quantities{
  
}


