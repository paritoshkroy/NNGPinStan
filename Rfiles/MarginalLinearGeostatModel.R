rm(list=ls())
library(sf)
library(tidyverse)

rainfall <- read_csv(file = "data/ParanaRainfall.csv")
border <- read_csv(file = "data/ParanaBorder.csv")
predLocs <- read_csv(file = "data/ParanaPredLocs.csv")

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
              X = unname(cbind(1, scale(rainfall$east), scale(rainfall$north))),
              coords = unname(cbind(rainfall$east, rainfall$north)),
              scale_beta = c(100,1,1),
              scale_sigma = sd(rainfall$rainfall)/2.58,
              scale_tau = sd(rainfall$rainfall)/2.58,
              a = a,
              b = b)
str(input)

library(cmdstanr)
stan_file <- "MarginalLinGeostatModel.stan"
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

pars <- c(paste0("beta[",1:3,"]"),"sigma","phi","tau")
quantile2.5 <- function(x) quantile(x, prob = 0.025)
quantile50 <- function(x) quantile(x, prob = 0.50)
quantile97.5 <- function(x) quantile(x, prob = 0.975)
fixed_summary <- fit$summary(NULL, c("mean","sd","quantile50","quantile2.5","quantile97.5","rhat","ess_bulk","ess_tail"))
fixed_summary <- fixed_summary %>% filter(variable %in% pars)
fixed_summary %>% print(digits = 3)
# A tibble: 6 × 9
#  variable     mean      sd    `50%` `2.5%` `97.5%`  rhat ess_bulk ess_tail
#1 phi      294.     123.    269.     138.    611.    1.00    2194.    2551.
#2 beta[1]  232.      30.3   235.     164.    286.    1.00    2911.    1874.
#3 beta[2]   -0.0999   0.994  -0.0844  -2.03    1.83  1.00    4670.    3111.
#4 beta[3]   -0.331    1.01   -0.330   -2.31    1.63  1.00    3902.    3276.
#5 sigma     46.3      7.17   45.5     34.8    63.1   1.00    2156.    2106.
#6 tau       18.3      1.94   18.2     14.5    22.0   1.00    3641.    2236.

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
post_beta <- as_tibble(draws_df) %>% dplyr::select(starts_with("beta[")) %>% as.matrix() %>% unname(); str(post_beta)
post_sigma <- draws_df$sigma; str(post_sigma)
post_phi <- draws_df$phi; str(post_phi)
post_tau <- draws_df$tau; str(post_tau)

########################################################################
# Spatial Prediction/Interpolation/Kriging
# To make predictions at unobserved locations, you must have the same covariates (used in model building) available for those locations
#######################################################################
rstan::expose_stan_functions(stanmodel = "utilities.stan")

## Coordinates of prediction locations
predLocs
predCoords <- unname(as.matrix(predLocs[,c("east","north")]))
predCoords

## Matrix of covariate at prediction locations
## Standardize the covariate **"east"** using the `scale` function, ensuring that the same **mean (center)** and **standard deviation (scale)** from the observed data (used for model fitting) are applied.
predX1 <- scale(x = predLocs$east, 
                center = mean(rainfall$east), 
                scale = sd(rainfall$east))

predX2 <- scale(x = predLocs$north, 
                center = mean(rainfall$north), 
                scale = sd(rainfall$north))

predX <- cbind(1,predX1,predX2)
predX

## The outcome vector is represented as **c(y1, y2)**, where **y1** corresponds to the outcomes at prediction locations and **y2** corresponds to the outcomes at observed locations. Similarly, the covariance matrix is partitioned following the same order.

## The following function generates samples from the conditional distribution of **y1** given **y2**, which represents the **predictive distribution** of **y1** given **y2**. The mean vectors **mu1** and **mu2** depend on the **mean regression coefficients** and the **covariates** at the prediction and observed data locations, respectively.
args(multi_conditional_normal_rng)

l <- 1
#######################################
## Joint Prediction
#######################################
draws_predy_list <- lapply(1:size_post_samples, function(l){
  
  ## data used to model fitting
  y2 <- input$y
  ## mean vector at prediction locations
  mu1 <- drop(predX %*% post_beta[l,])
  ## mean vector at observation locations
  mu2 <- drop(input$X %*% post_beta[l,])
  
  ## Covariance among the outcomes at prediction locations: Cov(y1)
  V11 <- exposed_exponential_cov(
    coords = lapply(1:nrow(predCoords), function(i) predCoords[i,]),
    sigma = post_sigma[l], length_scale = post_phi[l]) + diag(x = post_tau[l]^2, nrow = 4, ncol = 4)
  
  ## Covariance among the outcomes at observation locations: Cov(y2)
  V22 <- exposed_exponential_cov(
    coords = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    sigma = post_sigma[l], length_scale = post_phi[l]) + diag(x = post_tau[l]^2, nrow = input$n, ncol = input$n)
  
  ## Covariance between the outcomes at prediction locations and observation locations: Cov(y1,y2)
  V12 <- exposed_pred2obs_exponential_cov(
    coords1 = lapply(1:nrow(predCoords), function(i) predCoords[i,]),
    coords2 = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    sigma = post_sigma[l], length_scale = post_phi[l])
  
  predy <- multi_conditional_normal_rng(y2 = y2, mu1 = mu1, mu2 = mu2, V11 = V11, V12 = V12,  V22 = V22)
  
  return(predy)
})

draws_predy <- do.call(rbind, draws_predy_list)
str(draws_predy)
draws_predy_df <- as_tibble(draws_predy, .name_repair = ~paste0(1:4))
draws_predy_df                 

draws_predy_df %>%
  gather(predLocID, predy) %>%
  ggplot(aes(x = predy)) + 
  geom_density() + 
  facet_wrap(~predLocID, scales = "free_y", labeller = label_bquote(Prediction~Location~.(predLocID))) +
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
    sigma = post_sigma[l], length_scale = post_phi[l]) + diag(x = post_tau[l]^2, nrow = input$n, ncol = input$n)
  
  predy <- unlist(lapply(1:4, function(k){
    
    V12 <- exposed_pred2obs_exponential_cov(
      coords1 = list(predCoords[k,]),
      coords2 = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
      sigma = post_sigma[l], length_scale = post_phi[l])
    
    out <- uni_conditional_normal_rng(
      y2 = input$y, 
      mu1 = drop(predX[k,] %*% post_beta[l,]), 
      mu2 = drop(input$X %*% post_beta[l,]), 
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
  gather(predLocID, predy) %>%
  ggplot(aes(x = predy)) + 
  geom_density() + 
  facet_wrap(~predLocID, scales = "free_y", labeller = label_bquote(Prediction~Location~.(predLocID))) +
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

predCoords <- unname(st_coordinates(prid_grid_within))
head(predCoords)
nrow(predCoords)

predX1 <- scale(x = predCoords[,1], 
                center = mean(rainfall$east), 
                scale = sd(rainfall$east))
predX2 <- scale(x = predCoords[,2],
                center = mean(rainfall$north), 
                scale = sd(rainfall$north))
predX <- cbind(1,predX1,predX2)
head(predX)
nrow(predX)

l <- 1
## The following perform a joint prediction, however, univariate prediction can be done similarly as before.
draws_predy_list <- lapply(1:size_post_samples, function(l){
    
    V11 <- exposed_exponential_cov(
      coords = lapply(1:nrow(predCoords), function(i) predCoords[i,]),
      sigma = post_sigma[l], length_scale = post_phi[l]) + diag(x = post_tau[l]^2, nrow = nrow(predCoords), ncol = nrow(predCoords))
    
    V22 <- exposed_exponential_cov(
      coords = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
      sigma = post_sigma[l], length_scale = post_phi[l]) + diag(x = post_tau[l]^2, nrow = input$n, ncol = input$n)
    
    V12 <- exposed_pred2obs_exponential_cov(
      coords1 = lapply(1:nrow(predCoords), function(i) predCoords[i,]),
      coords2 = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
      sigma = post_sigma[l], length_scale = post_phi[l])
    
    predy <- multi_conditional_normal_rng(y2 = input$y, mu1 = drop(predX %*% post_beta[l,]), mu2 = drop(input$X %*% post_beta[l,]), V11 = V11, V12 = V12,  V22 = V22)
    
    return(predy)
  })
  
draws_predy <- do.call(rbind, draws_predy_list)
str(draws_predy)
draws_predy_df <- as_tibble(draws_predy, .name_repair = ~paste0(1:nrow(predCoords)))
draws_predy_df                 
  
draws_predy_df %>%
  gather(predLocID, predy) %>%
  filter(predLocID %in% sample(unique(predLocID), size = 5, replace = FALSE)) %>%
  ggplot(aes(x = predy)) + 
  geom_density() + 
  facet_wrap(~predLocID, scales = "free_y", 
             labeller = label_bquote(Prediction~Location~.(predLocID))) +
  theme_bw() +
  theme(strip.background = element_blank())


predy_post_summary <- draws_predy_df %>% 
  gather(locID, y) %>% 
  mutate(locID = as.numeric(locID)) %>% 
  group_by(locID) %>%
  summarise(post.mean = mean(y),
            post.sd = sd(y),
            post.q2.5 = quantile2.5(y),
            post.q50 = quantile50(y),
            post.q97.5 = quantile97.5(y)) %>%
  mutate(east = predCoords[,1]) %>%
  mutate(north = predCoords[,2])


predy_post_summary_sf <- st_as_sf(predy_post_summary, coords = c("east", "north"))
ggplot() + 
  geom_sf(data = parana_poly, fill = NA) + 
  geom_sf(data = predy_post_summary_sf, aes(col = post.mean)) +
  colorspace::scale_color_continuous_sequential(palette = "Viridis") +
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
    mu = drop(input$X %*% post_beta[l,]),
    sigma = post_sigma[l],
    phi = post_phi[l],
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
  theme(legend.title = element_blank())
theme_bw()
