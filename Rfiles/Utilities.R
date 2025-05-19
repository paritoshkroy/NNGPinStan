getIGamma <- function(par, lRange, uRange, prob) {
  
  lProb <- (1 - prob)/2
  uProb <- 1 - lProb
  
  invlRange <- 1/lRange
  invuRange <- 1/uRange
  
  c(qgamma(p = lProb, shape = par[1], rate = par[2], lower.tail = TRUE) - invuRange, qgamma(p = uProb, shape = par[1], rate = par[2], lower.tail = TRUE) - invlRange)
}


dinvgamma <- function (x, shape, scale, log = FALSE) {
  ldensity <- dgamma(1/x, shape, rate = scale, log = TRUE) - 2 * log(x)
  if (log){
    out <- ldensity
  } else{
    out <- exp(ldensity)
  }
  return(out)
}

qinvgamma <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE){
  1/qgamma(1 - p, shape, rate = scale, lower.tail = lower.tail, log.p = log.p)
}

rinvgamma <- function(n, shape, scale){
  1/rgamma(n, shape, rate = scale)
}

pinvgamma <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE){
  pgamma(1/q, shape, rate = scale, lower.tail = !lower.tail, log.p = log.p)
}


## Half normal
dhalfnorm <- function(x, sigma, log = FALSE){
  if(log){
    d <- ifelse(x < 0, -Inf, log(2) + dnorm(x, mean = 0, sd = sigma, log = TRUE))
  }else{
    d <- ifelse(x < 0, 0, 2 * dnorm(x, mean = 0, sd = sigma))
  }
  return(d)
}

phalfnorm <- function(q, sigma, lower.tail = TRUE, log.p = FALSE){
  p <- ifelse(q < 0, 0, 2 * pnorm(q, mean = 0, sd = sigma) - 1)
  if(lower.tail == FALSE){
    p <- 1-p
  }
  if(log.p){
    p <- log(p)
  }
  return(p)
}

qhalfnorm <- function(p, sigma, lower.tail = TRUE, log.p = FALSE){
  if(log.p) p <- exp(p)
  if(lower.tail == FALSE) p <- 1-p
  q <- ifelse(p<0, NaN, qnorm((p+1)/2, mean = 0, sd = sigma))
  return(q)
}

rhalfnorm <- function(n, sigma){
  return(abs(rnorm(n, mean = 0, sd = sigma)))
}
