functions {

real responseNNGP_matern32_lpdf(vector y, vector mu, real sigma, real tau, real lscale, matrix site2neiDist, matrix neiDistMat, array[,] int neiID, int n, int m) {
    
    vector[n] V;
    vector[n] resid = y - mu;
    vector[n] U = resid;
    real tausq = square(tau);
    real sigmasq = square(sigma);
    
    real variance_ratio_plus_1 = tausq * inv(sigmasq) + 1; // variance ratio plus 1
    int h;
    for (i in 2:n) {
      int dim = (i < (m + 1))? (i - 1) : m;
      matrix[dim, dim] neiCorMat;
      matrix[dim, dim] neiCorChol;
      vector[dim] site2neiCor;
      vector[dim] v;
      row_vector[dim] v2;
      
      if(dim == 1){
        neiCorMat[1, 1] = variance_ratio_plus_1;
        } else {
          h = 0;
          for (j in 1:(dim - 1)){
            for (k in (j + 1):dim){
              h = h + 1;
              neiCorMat[j, k] = (1 + sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale)) * exp(-sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale));
              neiCorMat[k, j] = neiCorMat[j, k];
              }
            }
            for(j in 1:dim){
              neiCorMat[j, j] = variance_ratio_plus_1;
            }
        }

        neiCorChol = cholesky_decompose(add_diag(neiCorMat, 1e-7));
        site2neiCor = to_vector((1 + sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)) .* exp(-sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)));
        v = mdivide_left_tri_low(neiCorChol, site2neiCor);
        V[i] = variance_ratio_plus_1 - dot_self(v);
        v2 = mdivide_right_tri_low(v', neiCorChol);
        U[i] = U[i] - v2 * resid[neiID[(i - 1), 1:dim]];
        }
        V[1] = variance_ratio_plus_1;
        return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) + sum(log(V)) + n * log(sigmasq));
      }
 
}


data{
  int n;
  int m;
  int p;
  vector[n] y;
  matrix[n,p] X;
  
  vector<lower=0>[p] scale_theta;
  real<lower = 0> scale_tau;
  real<lower = 0> scale_sigma;
  
  real<lower = 0> a; // Shape parameters in the prior for ell
  real<lower = 0> b; // Scale parameters in the prior for ell
  
  array[n-1, m] int neiID;
  matrix[n-1, m] site2neiDist;
  matrix[n-1, (m*(m-1))%/%2] neiDistMat;
  array[n-1] int nNei;

}

transformed data {
}

parameters {
  vector[p] theta_std;
  real<lower=0> tau_std;
  real<lower=0> ell;
  real<lower=0> sigma_std;
}

transformed parameters{
  vector[p] theta = scale_theta .* theta_std;
  real tau = scale_tau * tau_std;
  real sigma = scale_sigma * sigma_std;
}


model {
  theta_std ~ std_normal();
  tau_std ~ std_normal();
  sigma_std ~ std_normal();
  ell ~ inv_gamma(a,b);
  
  vector[n] mu = X*theta;
  target += responseNNGP_matern32_lpdf(y | mu, sigma, tau, ell, site2neiDist, neiDistMat, neiID, n, m);
  
}

generated quantities{
  
}





