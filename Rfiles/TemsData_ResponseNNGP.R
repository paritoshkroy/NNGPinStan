rm(list=ls())
graphics.off()
load("~/git/NNGPinStan/Data/AllSatelliteTemps.RData")
ls()
## It is a dataset of size 150000
str(all.sat.temps)
table(is.na(all.sat.temps$MaskTemp[!is.na(all.sat.temps$MaskTemp)]))
table(is.na(all.sat.temps$TrueTemp))

## Take a random sample of size 5000
temp_data <- all.sat.temps %>% 
  filter(!is.na(MaskTemp)) %>% # drop the missing cases
  sample_n(size = 5000) ## Change the size as you want

## Draw a random sample of size 500 among the selected 5000
train_data <- temp_data %>% sample_n(size = 500) ## Change the size as you want
str(train_data)
nsize <- nrow(train_data)
p <- 3
y <- train_data$MaskTemp
X <- unname(cbind(1, scale(train_data$Lon), scale(train_data$Lat)))
coords <- unname(cbind(train_data$Lon, train_data$Lat))

#########################################################
# NNGP preparation
#########################################################
m <- 10
library(BayesNSGP)
ord <- order(coords[,1])
ord_coords <- coords[ord,]
obs_nei_id <- determineNeighbors(coords = ord_coords, k = m)[-1,]
str(obs_nei_id)
nei_size <- rowSums(obs_nei_id!=-1)

site2nei_dist <- matrix(-1, nrow = nsize-1, ncol = m); str(site2nei_dist)
for(i in 2:nsize){
  site2nei_dist[i-1,1:nei_size[i-1]] <- rdist(rbind(ord_coords[i,], ord_coords[obs_nei_id[i-1,1:nei_size[i-1]],]))[-1,1]
}
head(site2nei_dist)
nei_dist_mat <- matrix(-1, nrow = nsize-1, ncol = m*(m-1)/2)
str(nei_dist_mat)

for(i in 2:nsize){
  dim <- nei_size[i-1]
  nei_dist_mat[i-1,1:(0.5*dim*(dim-1))] <- as.vector(dist(ord_coords[obs_nei_id[i-1,1:dim],]))
}
str(nei_dist_mat)

##########################################################
## Import some utility functions
##########################################################
source("Rfiles/Utilities.R")
##########################################################
## Order the input
##########################################################
ord_y <- y[ord]; str(ord_y)
ord_X <- X[ord,]; str(ord_X)
input <- list(n = nsize,
              m = m,
              p = 3,
              y = ord_y,
              X = ord_X,
              theta_multiplier = c(100, 1, 1),
              tau_multiplier = sd(y)/qhalfnorm(p = 0.99,sigma = 1),
              sigma_multiplier = sd(y)/qhalfnorm(p = 0.99,sigma = 1),
              a = 2,
              b = 0.1,
              neiID = obs_nei_id,
              site2neiDist = site2nei_dist,
              neiDistMat = nei_dist_mat,
              nNei = nei_size)
str(input)

stan_file <- "StanFiles/ResponseNNGP.stan"
mod <- cmdstan_model(stan_file, compile = FALSE)
mod$check_syntax(pedantic = TRUE)
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$print()

fit <- mod$sample(data = input, 
                  chains = 4,
                  parallel_chains = 4,
                  iter_warmup = 500,
                  iter_sampling = 500)
fit$time()
fit$time()$total/60
fit$cmdstan_diagnose()
fit$summary()

sampler_diag <- fit$sampler_diagnostics(format = "df")
str(sampler_diag)

q2.5 <- function(x) quantile(x, prob = 0.025)
q25 <- function(x) quantile(x, prob = 0.25)
q50 <- function(x) quantile(x, prob = 0.50)
q75 <- function(x) quantile(x, prob = 0.75)
q97.5 <- function(x) quantile(x, prob = 0.975)

fit_summary <- fit$summary(variables = c("theta","sigma","ell","tau"), c("mean","sd","q2.5","q25","q50","q75","q97.5","rhat","ess_bulk","ess_tail"))
fit_summary

## Posterior draws
draws_df <- fit$draws(format = "df")
draws_df

library(bayesplot)
color_scheme_set("brewer-Spectral")
mcmc_trace(draws_df,  pars = c("theta[1]", "theta[2]", "theta[3]", "sigma", "ell","tau"), facet_args = list(ncol = 3)) + facet_text(size = 15)

###############################################################
# Extract Posterior Samples
###############################################################
size_post_samples <- nrow(draws_df); size_post_samples
post_theta <- as_tibble(draws_df) %>% dplyr::select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
post_sigma <- draws_df$sigma; str(post_sigma)
post_ell <- draws_df$ell; str(post_ell)
post_tau <- draws_df$tau; str(post_tau)

########################################################################
# NNGP Spatial Prediction/Interpolation/Kriging
# By definition NNGP spatial prediction is univariate prediction
# Meaning that prediction is performed for each location separately
#######################################################################

## Preparing the prediction settings from cross-validation data
pred_data <- anti_join(temp_data, train_data, by = join_by(Lon, Lat, MaskTemp, TrueTemp))
str(pred_data)

## Prediction location
pred_coords <- unname(as.matrix(pred_data[,c("Lon","Lat")]))
str(pred_coords)
psize <- nrow(pred_coords)
psize

## Covariate matrix at prediction locations
pred_x1 <- scale(x = pred_data$Lon, center = mean(train_data$Lon), scale = sd(train_data$Lon))[,1]
pred_x2 <- scale(x = pred_data$Lat, center = mean(train_data$Lat), scale = sd(train_data$Lat))[,1]
pred_X <- unname(cbind(1, pred_x1, pred_x2))
str(pred_X)

## Nearest neighbors for prediction location
nei_info_pred <- FNN::get.knnx(ord_coords, pred_coords, k = m)
str(nei_info_pred)
pred2obs_nei_id <- nei_info_pred$nn.index
str(pred2obs_nei_id)

pred2obs_nei_dist <- nei_info_pred$nn.dist
str(pred2obs_nei_dist)

## Exposed Stan function requires to make prediction
rstan::expose_stan_functions(stanmodel = "StanFiles/utilities.stan")
args(predict_responseNNGP_rng)

l <- 1
draws_pred_y_list <- predict_responseNNGP_rng(
  y = input$y, 
  X = input$X, 
  pred_X = pred_X, 
  coords = lapply(1:nrow(ord_coords), function(i) ord_coords[i,]), 
  pred2obs_nei_dist = lapply(1:nrow(pred2obs_nei_dist), function(i) pred2obs_nei_dist[i,]),
  pred2obs_nei_id = lapply(1:nrow(pred2obs_nei_id), function(i) pred2obs_nei_id[i,]),
  theta = lapply(1:nrow(post_theta), function(i) post_theta[i,]), 
  sigma = post_sigma, 
  lscale = post_ell, 
  tau = post_tau, 
  nsize = nsize, 
  psize = psize, 
  L = size_post_samples, 
  print_interval = 10)

draws_pred_y <- do.call(rbind, draws_pred_y_list)
str(draws_pred_y)
draws_pred_y_dt <- as_tibble(draws_pred_y, .name_repair = ~paste0(1:psize)) %>%
  gather(pred_loc_id, pred) %>%
  mutate(pred_loc_id = as.numeric(pred_loc_id))
draws_pred_y_dt 


## Compute posterior summaries
pred_summary_dt <- draws_pred_y_dt %>%
  group_by(pred_loc_id) %>%
  summarise(post_mean = mean(pred),
            post_sd = sd(pred),
            post_q2.5 = quantile(pred, prob = 0.025),
            post_q50 = quantile(pred, prob = 0.50),
            post_q97.5 = quantile(pred, prob = 0.975)) %>%
  ungroup()

## Join the summaries of the predictive distributions and cross-validation data set
pred_summary_dt %>% nrow()
pred_summary_dt <- inner_join(pred_summary_dt, pred_data %>% mutate(pred_loc_id = row_number()), by = join_by(pred_loc_id))
pred_summary_dt %>% nrow()

## Compare the true values and predicted summaries
pred_summary_dt %>%
  drop_na() %>%
  ggplot(aes(x = TrueTemp, y = post_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = post_q2.5, ymax = post_q97.5)) +
  tune::coord_obs_pred() +
  theme_bw() +
  xlab("Observed") +
  ylab("Predicted mean")
