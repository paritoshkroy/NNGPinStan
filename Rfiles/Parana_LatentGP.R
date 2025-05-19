rm(list=ls())
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
              X = unname(cbind(1, scale(rainfall$east), scale(rainfall$north))),
              coords = unname(cbind(rainfall$east, rainfall$north)),
              scale_theta = c(100,1,1),
              scale_sigma = sd(rainfall$rainfall)/2.58,
              scale_tau = sd(rainfall$rainfall)/2.58,
              a = a,
              b = b)
str(input)

library(cmdstanr)
stan_file <- "StanFiles/LatentGP.stan"
mod <- cmdstan_model(stan_file, compile = FALSE)
mod$check_syntax(pedantic = TRUE)
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$print()
fit <- mod$sample(data = input, 
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 200,
                          iter_sampling = 200,
                          init = 0.25)
fit$cmdstan_diagnose()
sampler_diag <- fit$sampler_diagnostics(format = "df")
str(sampler_diag)

elapsed_time <- fit$time()
elapsed_time
elapsed_time$total/60

pars <- c(paste0("theta[",1:3,"]"),"sigma","ell","tau","z[1]","z[50]")
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
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
post_sigma <- draws_df$sigma; str(post_sigma)
post_ell <- draws_df$ell; str(post_ell)
post_tau <- draws_df$tau; str(post_tau)
post_z <- as_tibble(draws_df) %>% select(starts_with("z[")) %>% as.matrix() %>% unname(); str(post_z)

###############################################################
# Spatial Prediction/Interpolation/Kriging
###############################################################
rstan::expose_stan_functions(stanmodel = "StanFiles/utilities.stan")

## Joint Prediction
pred_data
pred_coords <- unname(as.matrix(pred_data[,c("east","north")]))
pred_coords

pred_x1 <- scale(x = pred_data$east, 
                center = mean(rainfall$east), 
                scale = sd(rainfall$east))

pred_x2 <- scale(x = pred_data$north, 
                center = mean(rainfall$north), 
                scale = sd(rainfall$north))

pred_X <- cbind(1,pred_x1,pred_x2)
pred_X

args(multi_conditional_normal_rng)

l <- 1

## Prediction of latent value
draws_predz_list <- lapply(1:size_post_samples, function(l){
  
  V11 <- exposed_exponential_cov(
    coords = lapply(1:nrow(pred_coords), function(i) pred_coords[i,]),
    sigma = post_sigma[l], length_scale = post_ell[l])
  
  V22 <- exposed_exponential_cov(
    coords = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    sigma = post_sigma[l], length_scale = post_ell[l])
  
  V12 <- exposed_pred2obs_exponential_cov(
    coords1 = lapply(1:nrow(pred_coords), function(i) pred_coords[i,]),
    coords2 = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    sigma = post_sigma[l], length_scale = post_ell[l])
  
  predz <- multi_conditional_normal_rng(y2 = post_z[l,], mu1 = rep(0,4), mu2 = rep(0,input$n), V11 = V11, V12 = V12,  V22 = V22)
  
  return(predz)
})

draws_predz <- do.call(rbind, draws_predz_list)
str(draws_predz)
draws_predz_df <- as_tibble(draws_predz, .name_repair = ~paste0(1:4))
draws_predz_df

draws_pred_y_list <- lapply(1:size_post_samples, function(l){
  
  mean_pred <- drop(pred_X %*% post_theta[l,]) + draws_predz_list[[l]]
  pred <- rnorm(n = 4, mean =  mean_pred, sd = post_tau[l])
  
  return(pred)
})

draws_pred_y <- do.call(rbind, draws_pred_y_list)
str(draws_pred_y)
draws_pred_y_df <- as_tibble(draws_pred_y, .name_repair = ~paste0(1:4))
draws_pred_y_df                 

draws_pred_y_df %>%
  gather(pred_loc_id, pred) %>%
  ggplot(aes(x = pred)) + 
  geom_density() + 
  facet_wrap(~pred_loc_id, scales = "free_y", labeller = label_bquote(Prediction~Location~.(pred_loc_id))) +
  theme_bw() +
  theme(strip.background = element_blank())



