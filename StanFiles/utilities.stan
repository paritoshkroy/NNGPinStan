functions{

    /*Exposed Exponential Covariance Matrix*/
    matrix exposed_exponential_cov(array[] vector coords, real sigma, real length_scale) {
    return gp_exponential_cov(coords, sigma, length_scale);
    }
    
    /*Exposed Exponential Covariance Matrix for Prediction Location to Observed Data Location*/
    matrix exposed_pred2obs_exponential_cov(array[] vector coords1, array[] vector coords2, real sigma, real length_scale) {
    return gp_exponential_cov(coords1, coords2, sigma, length_scale);
    }
    
    
    /*Exposed Matern3/2 Covariance Matrix*/
    matrix exposed_matern32_cov(array[] vector coords, real sigma, real length_scale) {
    return gp_matern32_cov(coords, sigma, length_scale);
    }
    
    /*Exposed Matern3/2 Covariance Matrix for Prediction Location to Observed Data Location*/
    matrix exposed_pred2obs_matern32_cov(array[] vector coords1, array[] vector coords2, real sigma, real length_scale) {
    return gp_matern32_cov(coords1, coords2, sigma, length_scale);
    }
    
   /*Stan function for backsolve*/
   vector mdivide_left_tri_upper(matrix U, vector b) {
     int n = rows(U);
     vector[n] x = b; 
     real cs;
     array[n] int ids = sort_desc(linspaced_int_array(n, 1, n));
     for (i in ids){
       x[i] = x[i]/U[i,i];
       cs = 0;
       for (j in (i+1):n){
         cs = cs + x[j]*U[i,j];
         }
       x[i] = x[i] - cs/U[i,i];
       }
    return x;
 }
 
  /*Generate from conditional multivariate normal given 
    partition component of a multivariate normal*/
  vector multi_conditional_normal_rng(vector y2, vector mu1, vector mu2, matrix V11, matrix V12, matrix V22){
   
   int n = size(y2);
   int q = size(mu1);
   vector[q] yc;
   vector[q] mc;
   matrix[q,q] vc;
   matrix[n, n] ch;
   
   ch = cholesky_decompose(V22);
   mc = mu1 + V12*mdivide_left_tri_upper(ch',mdivide_left_tri_low(ch, y2-mu2));
   vc = V11 - crossprod(mdivide_left_tri_low(ch, V12'));
   
   yc = multi_normal_rng(mc, vc);
   
   return yc;
 }

  /*Generate from conditional multivariate normal given 
    partition component of a multivariate normal*/
  real uni_conditional_normal_rng(vector y2, real mu1, vector mu2, real V11, vector V12, matrix V22){
   
   int n = size(y2);
   real yc;
   real mc;
   real vc;
   matrix[n, n] ch;
   
   ch = cholesky_decompose(V22);
   mc = mu1 + dot_product(V12,mdivide_left_tri_upper(ch',mdivide_left_tri_low(ch, y2-mu2)));
   vc = V11 - dot_self(mdivide_left_tri_low(ch, V12));
   
   yc = normal_rng(mc, sqrt(vc));
   
   return yc;
 }

  /*Recovery after fitting marginal linear geostatistical model*/
  vector revoery_latent_gp_exponential_cov_rng(array[] vector coords, vector y, vector mu, real sigma, real phi, real tau){
   
   int n = size(y);
   vector[n] z;
   
   vector[n] estd = (y-mu)*inv_square(tau);
   matrix[n,n] L = cholesky_decompose(add_diag(inverse(gp_exponential_cov(coords, sigma, phi)), inv_square(tau)));
   vector[n] mn = mdivide_right_tri_low(mdivide_left_tri_low(L, estd)',L)';
   z = multi_normal_cholesky_rng(mn, L);
   return z;
 }


array[] vector predict_responseNNGP_rng(
                vector y, matrix X, matrix pred_X, 
                array[] vector coords, 
                array[] vector pred2obs_nei_dist,
                array[,] int pred2obs_nei_id, 
                array[] vector theta,
                vector sigma, vector lscale, vector tau, 
                int nsize, 
                int psize, 
                int L, 
                int print_interval){
    
    array[L] vector[psize] out;
    int nprint = L %/% print_interval;
    int m = dims(pred2obs_nei_dist)[2];
    int p = dims(X)[2];
    
    for(l in 1:L) {
      
      if(l%nprint == 0) print("Starts for prediction location : ", l);
      
      for(i in 1:psize) {
        
        
        // Correlation between a prediction location and m nearest neighbors
        vector[m] ds = pred2obs_nei_dist[i] * inv(lscale[l]) * sqrt(3.0);
        vector[m] pred2obs_nei_cor = (1 + ds) .* exp(-ds);
        
        // Cholesky decomposition of correlation between m nearest neighbors
        matrix[m,m] obs_nei_chol = cholesky_decompose(add_diag(gp_matern32_cov(coords[pred2obs_nei_id[i,1:m]], 1, lscale[l]), rep_vector(square(tau[l])*inv_square(sigma[l]), m)));
        
        // Conditional mean for the predictive distribution
        real conditional_mean =  pred_X[i,1:p]*theta[l] + mdivide_left_tri_low(obs_nei_chol, pred2obs_nei_cor)' * mdivide_left_tri_low(obs_nei_chol, y[pred2obs_nei_id[i,1:m]] - X[pred2obs_nei_id[i,1:m],1:p] * theta[l]);
        
        // Conditional variance for the predictive distribution
        real conditional_variance = square(sigma[l]) * (1 + square(tau[l]) * inv_square(sigma[l]) - dot_self(mdivide_left_tri_low(obs_nei_chol, pred2obs_nei_cor)));
        
        // Generate from predictive distribution
        out[l][i] = normal_rng(conditional_mean, sqrt(conditional_variance));
      }
    }
    return out;
  }


}



data {
  
}

parameters {
}

model {

}

