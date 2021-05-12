library(tidyverse)
library(lme4)

##############
n = 30000

set.seed(123)

p1 <- 0.7
q1 <- 0.3

p2 <- 0.7
q2 <- 0.3


##### Generate random genotype for G1 and G2, and a normal environmental factor that is unknown:
G1 = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2))>0, FUN = "which",MARGIN = 2) - 1
G2 = apply(X = rmultinom(n,1,prob = c(p2^2,2*p2*q2,q2^2))>0, FUN = "which",MARGIN = 2) - 1
E <- rnorm(n, mean = 1, sd = 6)
# E <- rt(n = n, df = 1)
# E <- rexp(n, rate = 2) 

### Case 1: If the true model is nicely additive without interaction (Assuming probit model, inverse normal CDF as link function)
beta0 <- -1.5
beta1 <- 0.8
beta2 <- 0.7
beta3 <- 1
latent_y <- beta0 + beta1*G1 + beta2*G2 + rnorm(n = n)
y <- ifelse(latent_y>0,1,0)
data <- data.frame(y = y, G1 = G1, G2 = G2)

test_lin <- function(y,G1,G2){
  data <- data.frame(y = y, G1 = G1, G2 = G2)
  p <- suppressMessages(data %>% group_by(G1,G2) %>% summarise(p = mean(y),n=n()))
  compute_v <- function(p,n){
    ((1/dnorm(qnorm(p)))^2)*(p*(1-p))/n
  }
  v00 <- compute_v(p$p[p$G1 == 0 & p$G2 == 0],p$n[p$G1 == 0 & p$G2 == 0])
  v01 <- compute_v(p$p[p$G1 == 1 & p$G2 == 0],p$n[p$G1 == 1 & p$G2 == 0])
  v11 <- compute_v(p$p[p$G1 == 1 & p$G2 == 1],p$n[p$G1 == 1 & p$G2 == 1])
  v12 <- compute_v(p$p[p$G1 == 2 & p$G2 == 1],p$n[p$G1 == 2 & p$G2 == 1])
  s11 <- qnorm(p$p[p$G1 == 2 & p$G2 == 1]) - qnorm(p$p[p$G1 == 1 & p$G2 == 1])
  s00 <- qnorm(p$p[p$G1 == 1 & p$G2 == 0]) - qnorm(p$p[p$G1 == 0 & p$G2 == 0])
  t <- ((s11-s00)^2)/(v00+v01+v11+v12)
  1 - pchisq(t,1)
}

test_lin(y,G1,G2)
test_lin(y,G2,G1)



beta0 <- -1.5
beta1 <- 0.8
beta2 <- 0.7
beta3 <- 1
latent_y <- beta0 + beta1*G1 + beta2*G2 + beta3 * G1*E + rnorm(n = n)
y <- ifelse(latent_y>0,1,0)
test_lin(y,G1,G2)
test_lin(y,G2,G1)







