library(tidyverse)
library(lme4)

### functions
test_lin <- function(y,G1,G2,weight = c(1/3,1/3,1/3)){
  a <- weight
  data <- data.frame(y = y, G1 = G1, G2 = G2)
  p <- suppressMessages(data %>% group_by(G1,G2) %>% summarise(p = mean(y),n = n()))
  compute_v <- function(p,n){
    ((1/dnorm(qnorm(p)))^2)*(p*(1 - p))/n
  }
  compute_z <- function(p){
    qnorm(p)
  }
  v00 <- compute_v(p$p[p$G1 == 0 & p$G2 == 0],p$n[p$G1 == 0 & p$G2 == 0])
  v01 <- compute_v(p$p[p$G1 == 0 & p$G2 == 1],p$n[p$G1 == 0 & p$G2 == 1])
  v02 <- compute_v(p$p[p$G1 == 0 & p$G2 == 2],p$n[p$G1 == 0 & p$G2 == 2])
  v10 <- compute_v(p$p[p$G1 == 1 & p$G2 == 0],p$n[p$G1 == 1 & p$G2 == 0])
  v11 <- compute_v(p$p[p$G1 == 1 & p$G2 == 1],p$n[p$G1 == 1 & p$G2 == 1])
  v12 <- compute_v(p$p[p$G1 == 1 & p$G2 == 2],p$n[p$G1 == 1 & p$G2 == 2])
  v20 <- compute_v(p$p[p$G1 == 2 & p$G2 == 0],p$n[p$G1 == 2 & p$G2 == 0])
  v21 <- compute_v(p$p[p$G1 == 2 & p$G2 == 1],p$n[p$G1 == 2 & p$G2 == 1])
  v22 <- compute_v(p$p[p$G1 == 2 & p$G2 == 2],p$n[p$G1 == 2 & p$G2 == 2])
  z00 <- compute_z(p$p[p$G1 == 0 & p$G2 == 0])
  z01 <- compute_z(p$p[p$G1 == 0 & p$G2 == 1])
  z02 <- compute_z(p$p[p$G1 == 0 & p$G2 == 2])
  z10 <- compute_z(p$p[p$G1 == 1 & p$G2 == 0])
  z11 <- compute_z(p$p[p$G1 == 1 & p$G2 == 1])
  z12 <- compute_z(p$p[p$G1 == 1 & p$G2 == 2])
  z20 <- compute_z(p$p[p$G1 == 2 & p$G2 == 0])
  z21 <- compute_z(p$p[p$G1 == 2 & p$G2 == 1])
  z22 <- compute_z(p$p[p$G1 == 2 & p$G2 == 2])
  S1 <- a[1] * (z10 - z00) + a[2] * (z11 - z01) + a[3] * (z12 - z02)
  varS1 <- (a[1]^2)*(v10 + v00) + (a[2]^2)*(v11 + v01) + (a[3]^2)*(v12 + v02)
  S2 <- a[1] * (z20 - z10) + a[2] * (z21 - z11) + a[3] * (z22 - z12)
  varS2 <- (a[1]^2)*(v20 + v10) + (a[2]^2)*(v21 + v11) + (a[3]^2)*(v22 + v12)
  C <- -((a[1]^2)*v10 + (a[2]^2)*v11 + (a[3]^2)*v12)
  TestStats <- (((S1 - S2)^2)/(varS1 + varS2 - 2*C))
  1 - pchisq(TestStats,df = 1)
}
simulate_test <- function(n,p1,p2,beta = c(-1.5,0.8,0.7,1), inter = F, true_weight = F, num_trial = 1000, meanE = 1, sdE = 6){
  n <- n
  p1 <- p1
  q1 <- 1 - p1
  p2 <- p2
  q2 <- 1 - p2
  beta0 <- beta[1]
  beta1 <- beta[2]
  beta2 <- beta[3]
  beta3 <- beta[4]
  G1 = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
  G2 = apply(X = rmultinom(n,1,prob = c(p2^2,2*p2*q2,q2^2)) > 0, FUN = "which",MARGIN = 2) - 1
  E <- rnorm(n, mean = meanE, sd = sdE)
  p_val <- c()
  for (i in 1:num_trial) {
    G1 = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
    G2 = apply(X = rmultinom(n,1,prob = c(p2^2,2*p2*q2,q2^2)) > 0, FUN = "which",MARGIN = 2) - 1
    if (inter == T) latent_y <- beta0 + beta1*G1 + beta2*G2 + beta3*G1*E + rnorm(n = n)
    else latent_y <- beta0 + beta1*G1 + beta2*G2 + rnorm(n = n) 
    y <- ifelse(latent_y > 0,1,0)
    if (true_weight == T) p_val[i] <- test_lin(y,G1,G2, weight = c(p2^2,2*p2*q2,q2^2))
    else p_val[i] <- test_lin(y,G1,G2)
  }
  p_val
}

### simulations
set.seed(1234)
n = 3000
p1 <- 0.2
q1 <- 1 - p1

p2 <- 0.5
q2 <- 1 - p2

G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
Z = apply(X = rmultinom(n,1,prob = c(p2^2,2*p2*q2,q2^2)) > 0, FUN = "which",MARGIN = 2) - 1
E <- rnorm(n, mean = 1, sd = 6)

beta0 <- -1.2
beta1 <- 0.8
beta2 <- 0.3
beta3 <- 0.6

latent_y <- beta0 + beta1*G + beta2*Z + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G = G, Z = Z)
p <- data %>% group_by(G,Z) %>% summarise(p = mean(y))
test_lin(y,G,Z)

latent_y <- beta0 + beta1*G + beta2*Z + beta3*G*E + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G = G, Z = Z)
p <- data %>% group_by(G,Z) %>% summarise(p = mean(y))
test_lin(y,G,Z)
test_lin(y,Z,G)


#### Large interaction effect
p_values <- simulate_test(n = 3000, p1 = 0.2, p2 = 0.5, inter = FALSE, true_weight = T,num_trial = 1000,beta = c(-1.2,0.8,0.3,0.6))
hist(p_values,freq = F,breaks = 30)
p_values1 <- simulate_test(n = 3000, p1 = 0.2, p2 = 0.5, inter = T, true_weight = T,num_trial = 1000,beta = c(-1.2,0.8,0.3,0.6))
hist(p_values1,freq = F,breaks = 30)
p_values2 <- simulate_test_forG2(n = 3000, p1 = 0.2, p2 = 0.5, inter = T, true_weight = T,num_trial = 1000,beta = c(-1.2,0.8,0.3,0.6))
hist(p_values2,freq = F,breaks = 30)
#### Very small interaction effect
p_values <- simulate_test(n = 3000, p1 = 0.2, p2 = 0.5, inter = FALSE, true_weight = T,num_trial = 5000,beta = c(-1.2,0.8,0.3,0.003))
hist(p_values,freq = F,breaks = 30)
p_values1 <- simulate_test(n = 3000, p1 = 0.2, p2 = 0.5, inter = T, true_weight = T,num_trial = 5000,beta = c(-1.2,0.8,0.3,0.08))
hist(p_values1,freq = F,breaks = 30)
p_values2 <- simulate_test_forG2(n = 3000, p1 = 0.2, p2 = 0.5, inter = T, true_weight = T,num_trial = 5000,beta = c(-1.2,0.8,0.3,0.003))
hist(p_values2,freq = F,breaks = 30)
