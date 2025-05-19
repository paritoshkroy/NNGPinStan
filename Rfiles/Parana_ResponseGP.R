rm(list=ls())
graphics.off()
library(sf)
library(tidyverse)

rainfall <- read_csv(file = "Data/ParanaRainfall.csv")
border <- read_csv(file = "Data/ParanaBorder.csv")
pred_data <- read_csv(file = "Data/ParanaPredLocs.csv")

parana_poly <- st_sfc(st_polygon(list(cbind(border$east, border$north))))
ggplot(parana_poly) + geom_sf(fill = NA) + theme_bw()

distvec <- as.vector(dist(rainfall[,c("east","north")]))
hist(distvec)
a <- 5
b <- median(distvec)
hist(1/rgamma(n = 1000, shape = a, rate = b))

input <- list(n = length(rainfall$rainfall), 
              p = 3, 
              y = rainfall$rainfall, 
              X = unname(cbind(1, scale(rainfall$east)[,1], scale(rainfall$north)[,1])),
              coords = unname(cbind(rainfall$east, rainfall$north)),
              scale_theta = c(100,1,1),
              scale_sigma = sd(rainfall$rainfall)/2.58,
              scale_tau = sd(rainfall$rainfall)/2.58,
              a = a,
              b = b)
str(input)

library(cmdstanr)
stan_file <- "StanFiles/ResponseGP.stan"
mod <- cmdstan_model(stan_file, compile = FALSE)
mod$check_syntax(pedantic = TRUE)
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$print()
fit <- mod$sample(data = input, 
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          init = 0.25)
fit$cmdstan_diagnose()
sampler_diag <- fit$sampler_diagnostics(format = "df")
str(sampler_diag)

elapsed_time <- fit$time()
elapsed_time
elapsed_time$total/60

pars <- c(paste0("theta[",1:3,"]"), "sigma", "ell", "tau")
quantile2.5 <- function(x) quantile(x, prob = 0.025)
quantile50 <- function(x) quantile(x, prob = 0.50)
quantile97.5 <- function(x) quantile(x, prob = 0.975)
fixed_summary <- fit$summary(NULL, c("mean","sd","quantile50","quantile2.5","quantile97.5","rhat","ess_bulk","ess_tail"))
fixed_summary <- fixed_summary %>% filter(variable %in% pars)
fixed_summary %>% print(digits = 3)

draws_df <- fit$draws(format = "df")
str(draws_df)

library(bayesplot)
color_scheme_set("brewer-Spectral")
mcmc_trace(draws_df,  pars = pars, facet_args = list(ncol = 3)) + 
  facet_text(size = 15)

###############################################################
# Extract Posterior Samples
###############################################################
size_post_samples <- nrow(draws_df); size_post_samples
post_theta <- as_tibble(draws_df) %>% dplyr::select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
post_sigma <- draws_df$sigma; str(post_sigma)
post_ell <- draws_df$ell; str(post_ell)
post_tau <- draws_df$tau; str(post_tau)

########################################################################
# Spatial Prediction/Interpolation/Kriging
# To make predictions at unobserved locations, you must have the same covariates (used in model building) available for those locations
#######################################################################
rstan::expose_stan_functions(stanmodel = "StanFiles/utilities.stan")

## Coordinates of prediction locations
pred_data
pred_coords <- unname(as.matrix(pred_data[,c("east","north")]))
pred_coords

## Matrix of covariate at prediction locations
## Standardize the covariate **"east"** using the `scale` function, ensuring that the same **mean (center)** and **standard deviation (scale)** from the observed data (used for model fitting) are applied.
pred_x1 <- scale(x = pred_data$east, 
                center = mean(rainfall$east), 
                scale = sd(rainfall$east))

pred_x2 <- scale(x = pred_data$north, 
                center = mean(rainfall$north), 
                scale = sd(rainfall$north))

pred_X <- cbind(1,pred_x1,pred_x2)
pred_X

## The outcome vector is represented as **c(y1, y2)**, where **y1** corresponds to the outcomes at prediction locations and **y2** corresponds to the outcomes at observed locations. Similarly, the covariance matrix is partitioned following the same order.

## The following function generates samples from the conditional distribution of **y1** given **y2**, which represents the **predictive distribution** of **y1** given **y2**. The mean vectors **mu1** and **mu2** depend on the **mean regression coefficients** and the **covariates** at the prediction and observed data locations, respectively.
args(multi_conditional_normal_rng)

l <- 1
#######################################
## Joint Prediction
#######################################
draws_pred_y_list <- lapply(1:size_post_samples, function(l){
  
  ## data used to model fitting
  y2 <- input$y
  ## mean vector at prediction locations
  mu1 <- drop(pred_X %*% post_theta[l,])
  ## mean vector at observation locations
  mu2 <- drop(input$X %*% post_theta[l,])
  
  ## Covariance among the outcomes at prediction locations: Cov(y1)
  V11 <- exposed_exponential_cov(
    coords = lapply(1:nrow(pred_coords), function(i) pred_coords[i,]),
    sigma = post_sigma[l], length_scale = post_ell[l]) + diag(x = post_tau[l]^2, nrow = 4, ncol = 4)
  
  ## Covariance among the outcomes at observation locations: Cov(y2)
  V22 <- exposed_exponential_cov(
    coords = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    sigma = post_sigma[l], length_scale = post_ell[l]) + diag(x = post_tau[l]^2, nrow = input$n, ncol = input$n)
  
  ## Covariance between the outcomes at prediction locations and observation locations: Cov(y1,y2)
  V12 <- exposed_pred2obs_exponential_cov(
    coords1 = lapply(1:nrow(pred_coords), function(i) pred_coords[i,]),
    coords2 = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    sigma = post_sigma[l], length_scale = post_ell[l])
  
  predy <- multi_conditional_normal_rng(y2 = y2, mu1 = mu1, mu2 = mu2, V11 = V11, V12 = V12,  V22 = V22)
  
  return(predy)
})

draws_predy <- do.call(rbind, draws_pred_y_list)
str(draws_predy)
draws_predy_df <- as_tibble(draws_predy, .name_repair = ~paste0(1:4))
draws_predy_df                 

draws_predy_df %>%
  gather(pred_loc_id, pred) %>%
  ggplot(aes(x = pred)) + 
  geom_density() + 
  facet_wrap(~pred_loc_id, scales = "free_y", labeller = label_bquote(Prediction~Location~.(pred_loc_id))) +
  theme_bw() +
  theme(strip.background = element_blank())

#######################################
## Univariate Prediction
## When the number of prediction locations is very large, computing predictions for all locations at once can be computationally challenging. In such cases, predicting one location at a time is often preferable. This approach is known as **univariate prediction**.
#######################################

l <- 1
draws_uni_predy_list <- lapply(1:size_post_samples, function(l){
  
  V11 <- post_sigma[l]^2 + post_tau[l]^2
  
  V22 <- exposed_exponential_cov(
    coords = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    sigma = post_sigma[l], length_scale = post_ell[l]) + diag(x = post_tau[l]^2, nrow = input$n, ncol = input$n)
  
  predy <- unlist(lapply(1:4, function(k){
    
    V12 <- exposed_pred2obs_exponential_cov(
      coords1 = list(pred_coords[k,]),
      coords2 = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
      sigma = post_sigma[l], length_scale = post_ell[l])
    
    out <- uni_conditional_normal_rng(
      y2 = input$y, 
      mu1 = drop(pred_X[k,] %*% post_theta[l,]), 
      mu2 = drop(input$X %*% post_theta[l,]), 
      V11 = V11, V12 = V12,  V22 = V22)
    
    return(out)
    }))
  
  return(predy)
})

draws_uni_predy <- do.call(rbind, draws_uni_predy_list)
str(draws_uni_predy)
draws_uni_predy_df <- as_tibble(draws_uni_predy, .name_repair = ~paste0(1:4))
draws_uni_predy_df                 

draws_uni_predy_df %>%
  gather(pred_loc_id, pred) %>%
  ggplot(aes(x = pred)) + 
  geom_density() + 
  facet_wrap(~pred_loc_id, scales = "free_y", labeller = label_bquote(Prediction~Location~.(pred_loc_id))) +
  theme_bw() +
  theme(strip.background = element_blank())

####################################################################
## !!!!!!!!!!!!!!!!!!!!!!!!! Caution !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## To obtain smooth surface you need more grid points
## which is computationally time demanding
## Therefore univariate prediction may be a faster approach.
## Also you need associated covariates for all grid locations
## Prediction method is the same as before
####################################################################

## For example, consider a grid of 10 x 10 resolutions
## If we consider 100 x 100, it will take longer time to run

prid_grid <- st_make_grid(x = parana_poly, n = c(10,10), what = "centers")
ggplot(prid_grid) + geom_sf(size = 0.25) + 
  geom_sf(data = parana_poly, fill = NA) + theme_bw()
prid_grid_within <- st_intersection(parana_poly, prid_grid)
ggplot(prid_grid_within) + geom_sf(size = 0.25) + 
  geom_sf(data = parana_poly, fill = NA) + theme_bw()

pred_coords <- unname(st_coordinates(prid_grid_within))
head(pred_coords)
nrow(pred_coords)

pred_X1 <- scale(x = pred_coords[,1], 
                center = mean(rainfall$east), 
                scale = sd(rainfall$east))
pred_X2 <- scale(x = pred_coords[,2],
                center = mean(rainfall$north), 
                scale = sd(rainfall$north))
pred_X <- cbind(1,pred_X1,pred_X2)
head(pred_X)
nrow(pred_X)

l <- 1
## The following perform a joint prediction, however, univariate prediction can be done similarly as before.
draws_pred_y_list <- lapply(1:size_post_samples, function(l){
    
    V11 <- exposed_exponential_cov(
      coords = lapply(1:nrow(pred_coords), function(i) pred_coords[i,]),
      sigma = post_sigma[l], length_scale = post_ell[l]) + diag(x = post_tau[l]^2, nrow = nrow(pred_coords), ncol = nrow(pred_coords))
    
    V22 <- exposed_exponential_cov(
      coords = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
      sigma = post_sigma[l], length_scale = post_ell[l]) + diag(x = post_tau[l]^2, nrow = input$n, ncol = input$n)
    
    V12 <- exposed_pred2obs_exponential_cov(
      coords1 = lapply(1:nrow(pred_coords), function(i) pred_coords[i,]),
      coords2 = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
      sigma = post_sigma[l], length_scale = post_ell[l])
    
    predy <- multi_conditional_normal_rng(y2 = input$y, mu1 = drop(pred_X %*% post_theta[l,]), mu2 = drop(input$X %*% post_theta[l,]), V11 = V11, V12 = V12,  V22 = V22)
    
    return(predy)
  })
  
draws_pred_y <- do.call(rbind, draws_pred_y_list)
str(draws_pred_y)
draws_pred_y_df <- as_tibble(draws_pred_y, .name_repair = ~paste0(1:nrow(pred_coords)))
draws_predy_df                 
  
draws_pred_y_df %>%
  gather(pred_loc_id, pred) %>%
  mutate(pred_loc_id = as.numeric(pred_loc_id)) %>%
  filter(pred_loc_id %in% sample(unique(pred_loc_id), size = 5, replace = FALSE)) %>%
  ggplot(aes(x = pred)) + 
  geom_density() + 
  facet_wrap(~pred_loc_id, scales = "free_y", 
             labeller = label_bquote(Prediction~Location~.(pred_loc_id))) +
  theme_bw() +
  theme(strip.background = element_blank())


pred_y_post_summary <- draws_pred_y_df %>% 
  gather(pred_loc_id, pred) %>% 
  mutate(pred_loc_id = as.numeric(pred_loc_id)) %>% 
  group_by(pred_loc_id) %>%
  summarise(post.mean = mean(pred),
            post.sd = sd(pred),
            post.q2.5 = quantile2.5(pred),
            post.q50 = quantile50(pred),
            post.q97.5 = quantile97.5(pred)) %>%
  mutate(east = pred_coords[,1]) %>%
  mutate(north = pred_coords[,2])


pred_y_post_summary_sf <- st_as_sf(pred_y_post_summary, coords = c("east", "north"))
ggplot() + 
  geom_sf(data = parana_poly, fill = NA) + 
  geom_sf(data = pred_y_post_summary_sf, aes(col = post.mean)) +
  viridis::scale_color_viridis() +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(legend.title = element_blank())

#######################################
## Recovery of z
#######################################
l <- 1
args(revoery_latent_gp_exponential_cov_rng)
draw_z_list <- lapply(1:size_post_samples, function(l){
  
  recovered_z <- revoery_latent_gp_exponential_cov_rng(
    coords = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    y = input$y,
    mu = drop(input$X %*% post_theta[l,]),
    sigma = post_sigma[l],
    phi = post_ell[l],
    tau = post_tau[l])
  
  return(recovered_z)
})

draws_z <- do.call(rbind, draw_z_list)
str(draws_z)
draws_z_df <- as_tibble(draws_z, .name_repair = ~paste0(1:input$n))
draws_z_df                 

z_post_summary <- draws_z_df %>% 
  gather(locID, z) %>% 
  mutate(locID = as.numeric(locID)) %>% 
  group_by(locID) %>%
  summarise(post.mean = mean(z),
            post.sd = sd(z),
            post.q2.5 = quantile2.5(z),
            post.q50 = quantile50(z),
            post.q97.5 = quantile97.5(z)) %>%
  mutate(east = rainfall$east) %>%
  mutate(north = rainfall$north)


z_post_summary_sf <- st_as_sf(z_post_summary, coords = c("east", "north"))
ggplot() + 
  geom_sf(data = parana_poly, fill = NA) + 
  geom_sf(data = z_post_summary_sf, aes(size = post.mean), shape = 1) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(legend.title = element_blank()) + 
  theme_bw() +
  theme(legend.title = element_blank())
