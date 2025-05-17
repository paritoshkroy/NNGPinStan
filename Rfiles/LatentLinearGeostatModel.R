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
stan_file <- "LatentLinGeostatModel.stan"
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

pars <- c(paste0("beta[",1:3,"]"),"sigma","phi","tau","z[1]","z[50]")
quantile2.5 <- function(x) quantile(x, prob = 0.025)
quantile50 <- function(x) quantile(x, prob = 0.50)
quantile97.5 <- function(x) quantile(x, prob = 0.975)
fixed_summary <- fit$summary(NULL, c("mean","sd","quantile50","quantile2.5","quantile97.5","rhat","ess_bulk","ess_tail"))
fixed_summary <- fixed_summary %>% filter(variable %in% pars)
fixed_summary %>% print(digits = 3)
# A tibble: 6 × 9
# variable     mean     sd   `50%` `2.5%` `97.5%`  rhat ess_bulk ess_tail
#1 phi      295.     119.   274.    138.    594.    1.00     999.    1264.
#2 beta[1]  233.      29.0  235.    169.    284.    1.00    1301.    1416.
#3 beta[2]   -0.0845   1.00  -0.114  -1.96    1.89  1.00    4403.    1490.
#4 beta[3]   -0.311    1.03  -0.303  -2.36    1.61  1.00    4067.    1509.
#5 sigma     46.3      6.97  45.6    35.1    62.3   1.00    1894.    1438.
#6 tau       18.4      1.86  18.4    14.8    22.1   1.01     516.     840.
#7 z[1]      79.0     31.8   76.4    21.4   145.    1.00    1452.    1334.
#8 z[50]     66.8     31.6   65.2     9.40  134.    1.00    1433.    1422.
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
post_beta <- as_tibble(draws_df) %>% select(starts_with("beta[")) %>% as.matrix() %>% unname(); str(post_beta)
post_sigma <- draws_df$sigma; str(post_sigma)
post_phi <- draws_df$phi; str(post_phi)
post_tau <- draws_df$tau; str(post_tau)
post_z <- as_tibble(draws_df) %>% select(starts_with("z[")) %>% as.matrix() %>% unname(); str(post_z)

###############################################################
# Spatial Prediction/Interpolation/Kriging
###############################################################
source("expose_cmdstanr_functions.R")
utils <- expose_cmdstanr_functions(model_path = "utilities.stan")

## Joint Prediction
predLocs
predCoords <- unname(as.matrix(predLocs[,c("east","north")]))
predCoords

predX1 <- scale(x = predLocs$east, 
                center = mean(rainfall$east), 
                scale = sd(rainfall$east))

predX2 <- scale(x = predLocs$north, 
                center = mean(rainfall$north), 
                scale = sd(rainfall$north))

predX <- cbind(1,predX1,predX2)
predX

args(utils$multi_conditional_normal_rng)

l <- 1

## Prediction of latent value
draws_predz_list <- lapply(1:size_post_samples, function(l){
  
  V11 <- utils$exposed_exponential_cov(
    coords = lapply(1:nrow(predCoords), function(i) predCoords[i,]),
    sigma = post_sigma[l], length_scale = post_phi[l])
  
  V22 <- utils$exposed_exponential_cov(
    coords = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    sigma = post_sigma[l], length_scale = post_phi[l])
  
  V12 <- utils$exposed_pred2obs_exponential_cov(
    coords1 = lapply(1:nrow(predCoords), function(i) predCoords[i,]),
    coords2 = lapply(1:nrow(input$coords), function(i) input$coords[i,]),
    sigma = post_sigma[l], length_scale = post_phi[l])
  
  predz <- utils$multi_conditional_normal_rng(y2 = post_z[l,], mu1 = rep(0,4), mu2 = rep(0,input$n), V11 = V11, V12 = V12,  V22 = V22)
  
  return(predz)
})

draws_predz <- do.call(rbind, draws_predz_list)
str(draws_predz)
draws_predz_df <- as_tibble(draws_predz, .name_repair = ~paste0(1:4))
draws_predz_df

draws_predy_list <- lapply(1:size_post_samples, function(l){
  
  mean_predy <- drop(predX %*% post_beta[l,]) + draws_predz_list[[l]]
  predy <- rnorm(n = 4, mean =  mean_predy, sd = post_tau[l])
  
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






