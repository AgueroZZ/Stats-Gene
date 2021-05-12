library(tidyverse)
library(lme4)
library(doParallel)
library(foreach)
library(parallel)

set.seed(123)
n = 6000
p1 <- 0.5
q1 <- 1 - p1

p2 <- 0.5
q2 <- 1 - p2

G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
Z = apply(X = rmultinom(n,1,prob = c(p2^2,2*p2*q2,q2^2)) > 0, FUN = "which",MARGIN = 2) - 1
E <- rnorm(n, mean = 0, sd = 3)

G1 <- ifelse(G == 1, 1,0)
G2 <- ifelse(G == 2, 1,0)
G0 <- ifelse(G1 == G2, 1,0)


beta0 <- -1.2
beta1 <- 0.8
beta2 <- 0.3
beta3 <- 0.5
beta4 <- -0.5
betaz <- 0.5
betaE <- 0.5

latent_y <- beta0 + beta1*G1 + beta2*G2 + betaE*E + beta3*G1*E + beta4*G2*E  + betaz*Z + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G),Z=Z)
data$OLRE1 <- 1:nrow(data)
data$OLRE1[data$G != 1] <- 0
data$OLRE2 <- 1:nrow(data)
data$OLRE2[data$G != 2] <- 0

model21 <- glmer(y~ G + Z + (-1+G1|OLRE1) + (-1+G2|OLRE2) ,family = binomial(link = "probit"), data = data,  nAGQ = 1L)
summary(model21)
model22 <- glm(y~ G + Z ,family = binomial(link = "probit"), data = data)
summary(model22)
as.numeric((1 - pchisq(2*(stats::logLik(model21) - stats::logLik(model22)), df = 2))/4) + as.numeric((1 - pchisq(2*(stats::logLik(model21) - stats::logLik(model22)), df = 1))/2)




##### An approximation: using one random effect only


data2 <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G), Z=Z, Gsum = G1 + G2)
data2$OLRE1 <- 1:nrow(data2)
data2$OLRE1[data2$Gsum != 1] <- 0

model11 <- glmer(y~ G + Z + (-1+Gsum|OLRE1), family = binomial(link = "probit"), data = data2,  nAGQ = 8L)
model12 <- glm(y~ G + Z ,family = binomial(link = "probit"), data = data2)
as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)












#### Without Z:
set.seed(123)
n = 3000
p1 <- 0.5
q1 <- 1 - p1

G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
E <- rnorm(n, mean = 0, sd = 3)

G1 <- ifelse(G == 1, 1,0)
G2 <- ifelse(G == 2, 1,0)
G0 <- ifelse(G1 == G2, 1,0)


beta0 <- -1.2
beta1 <- 0.8
beta2 <- 0.3
beta3 <- 0.5
beta4 <- 0.5
betaE <- 0

latent_y <- beta0 + beta1*G1 + beta2*G2 + betaE*E + beta3*G1*E + beta4*G2*E  + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G))
data$OLRE1 <- 1:nrow(data)
data$OLRE1[data$G != 1] <- 0
data$OLRE2 <- 1:nrow(data)
data$OLRE2[data$G != 2] <- 0

model21 <- glmer(y~ G + (-1+G1|OLRE1) + (-1+G2|OLRE2) ,family = binomial(link = "probit"), data = data,  nAGQ = 1L)
summary(model21)
model22 <- glm(y~ G ,family = binomial(link = "probit"), data = data)
summary(model22)
as.numeric((1 - pchisq(2*(stats::logLik(model21) - stats::logLik(model22)), df = 2))/4) + as.numeric((1 - pchisq(2*(stats::logLik(model21) - stats::logLik(model22)), df = 1))/2)






##### An approximation: using one random effect only
data2 <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G), Gsum = G1 + G2)
data2$OLRE1 <- 1:nrow(data2)
data2$OLRE1[data2$Gsum != 1] <- 0

model11 <- glmer(y~ G + (-1+Gsum|OLRE1), family = binomial(link = "probit"), data = data2,  nAGQ = 20L)
model12 <- glm(y~ G, family = binomial(link = "probit"), data = data2)
as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)























NonAdditive_Simulator_With_Z <- function(beta0 = -1.2, betaG = c(0.8,0.3), betaE = 0.6, betaGE = c(0.6,0.4), muE = 0, sdE = 3, p1 = 0.5, p2 = 0.5, size = 3000, num_trails = 1000, parallel = F, nAGQ = 1L, betaZ = 1, SingleApproxi = F){
  n <- size
  q1 <- 1 - p1
  q2 <- 1 - p2
  beta1 <- betaG[1]
  beta2 <- betaG[2]
  betaGE1 <- betaGE[1]
  betaGE2 <- betaGE[2]
  
  if(parallel == F){
    p_val <- c()
    for (i in 1:num_trails) {
      G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
      
      G1 <- ifelse(G == 1, 1,0)
      G2 <- ifelse(G == 2, 1,0)
      G0 <- ifelse(G1 == G2, 1,0)
      Z <- apply(X = rmultinom(n,1,prob = c(p2^2,2*p2*q2,q2^2)) > 0, FUN = "which",MARGIN = 2) - 1
      E <- rnorm(n,muE,sdE)
      
      ylat <- beta0 + beta1*G1 + beta2*G2 +  betaE*E + betaGE1*E*G1 + betaGE2*E*G2 + betaZ * Z + rnorm(n)
      y <- ifelse(ylat > 0,1,0)
      
      
      data <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G),Z=Z)
      data$OLRE1 <- 1:nrow(data)
      data$OLRE1[data$G != 1] <- 0
      data$OLRE2 <- 1:nrow(data)
      data$OLRE2[data$G != 2] <- 0
      if(SingleApproxi == F){
      model11 <- suppressMessages(glmer(y~ G + Z + (-1+G1|OLRE1) + (-1+G2|OLRE2) ,family = binomial(link = "probit"), data = data,  nAGQ = nAGQ))
      model12 <- glm(y~ G + Z ,family = binomial(link = "probit"), data = data)
      p_val[i] <- as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 2))/4) + as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
      else{
        data$Gsum <- G1 + G2
        data$OLRE <- 1:nrow(data)
        data$OLRE[data$Gsum == 0] <- 0
        model11 <- GLMMadaptive::mixed_model(fixed = y ~ G + Z, random = ~ -1 + Gsum|OLRE, data = data, family = binomial(link = 'probit'), control = list(nAGQ = nAGQ))
        model12 <- glm(y~ G + Z ,family = binomial(link = "probit"), data = data)
        p_val[i] <- as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)

      }
    }
    p_val
  }
  
  #setup parallel backhand to use many processors
  if(parallel == T){
    cores= parallel::detectCores()
    cl <- parallel::makeCluster(cores[1]-2) #not to overload your computer
    doParallel::registerDoParallel(cl)
    p_val <- foreach(i=1:num_trails, .combine='c') %dopar% {
      G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
      
      G1 <- ifelse(G == 1, 1,0)
      G2 <- ifelse(G == 2, 1,0)
      G0 <- ifelse(G1 == G2, 1,0)
      
      E <- rnorm(n,muE,sdE)
      Z = apply(X = rmultinom(n,1,prob = c(p2^2,2*p2*q2,q2^2)) > 0, FUN = "which",MARGIN = 2) - 1
      ylat <- beta0 + beta1*G1 + beta2*G2 +  betaE*E + betaGE1*E*G1 + betaGE2*E*G2 + betaZ * Z + rnorm(n)
      y <- ifelse(ylat > 0,1,0)
      
      
      data <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G),Z=Z)
      data$OLRE1 <- 1:nrow(data)
      data$OLRE1[data$G != 1] <- 0
      data$OLRE2 <- 1:nrow(data)
      data$OLRE2[data$G != 2] <- 0
      
      if(SingleApproxi == F){
        model11 <- lme4::glmer(y~ G + Z + (-1+G1|OLRE1) + (-1+G2|OLRE2) ,family = binomial(link = "probit"), data = data,  nAGQ = nAGQ)
        model12 <- glm(y~ G + Z ,family = binomial(link = "probit"), data = data)
        as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 2))/4) + as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
      else{
        data$Gsum <- G1 + G2
        data$OLRE <- 1:nrow(data)
        data$OLRE[data$Gsum == 0] <- 0
        model11 <- GLMMadaptive::mixed_model(fixed = y ~ G + Z, random = ~ -1 + Gsum|OLRE, data = data, family = binomial(link = 'probit'), control = list(nAGQ = nAGQ))
        model12 <- glm(y~ G + Z ,family = binomial(link = "probit"), data = data)
        as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
    }
    #stop cluster
    parallel::stopCluster(cl)
    p_val
  }
  p_val
}






NonAdditive_Simulator_Without_Z <- function(beta0 = -1.2, betaG = c(0.8,0.3), betaE = 0.6, betaGE = c(0.6,0.4), muE = 0, sdE = 3, p1 = 0.5, size = 3000, num_trails = 1000, parallel = F, nAGQ = 1L, SingleApproxi = F){
  n <- size
  q1 <- 1 - p1
  beta1 <- betaG[1]
  beta2 <- betaG[2]
  betaGE1 <- betaGE[1]
  betaGE2 <- betaGE[2]
  
  if(parallel == F){
    p_val <- c()
    for (i in 1:num_trails) {
      G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
      
      G1 <- ifelse(G == 1, 1,0)
      G2 <- ifelse(G == 2, 1,0)
      G0 <- ifelse(G1 == G2, 1,0)
      
      E <- rnorm(n,muE,sdE)
      ylat <- beta0 + beta1*G1 + beta2*G2 +  betaE*E + betaGE1*E*G1 + betaGE2*E*G2  + rnorm(n)
      y <- ifelse(ylat > 0,1,0)
      
      
      data <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G))
      data$OLRE1 <- 1:nrow(data)
      data$OLRE1[data$G != 1] <- 0
      data$OLRE2 <- 1:nrow(data)
      data$OLRE2[data$G != 2] <- 0
      
      if(SingleApproxi == F){
        model11 <- suppressMessages(glmer(y~ G + (-1+G1|OLRE1) + (-1+G2|OLRE2) ,family = binomial(link = "probit"), data = data,  nAGQ = nAGQ))
        model12 <- glm(y~ G, family = binomial(link = "probit"), data = data)
        p_val[i] <- as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 2))/4) + as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
      else{
        data$Gsum <- G1 + G2
        data$OLRE <- 1:nrow(data)
        data$OLRE[data$Gsum == 0] <- 0
        model11 <- GLMMadaptive::mixed_model(fixed = y ~ G, random = ~ -1 + Gsum|OLRE, data = data, family = binomial(link = 'probit'), control = list(nAGQ = nAGQ))
        model12 <- glm(y~ G, family = binomial(link = "probit"), data = data)
        p_val[i] <- as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
    }
    p_val
  }
  
  #setup parallel backhand to use many processors
  if(parallel == T){
    cores= parallel::detectCores()
    cl <- parallel::makeCluster(cores[1]-2) #not to overload your computer
    doParallel::registerDoParallel(cl)
    p_val <- foreach(i=1:num_trails, .combine='c') %dopar% {
      G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
      
      G1 <- ifelse(G == 1, 1,0)
      G2 <- ifelse(G == 2, 1,0)
      G0 <- ifelse(G1 == G2, 1,0)
      
      E <- rnorm(n,muE,sdE)
      ylat <- beta0 + beta1*G1 + beta2*G2 +  betaE*E + betaGE1*E*G1 + betaGE2*E*G2 + rnorm(n)
      y <- ifelse(ylat > 0,1,0)
      
      
      data <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G))
      data$OLRE1 <- 1:nrow(data)
      data$OLRE1[data$G != 1] <- 0
      data$OLRE2 <- 1:nrow(data)
      data$OLRE2[data$G != 2] <- 0
      
      if(SingleApproxi == F){
        model11 <- lme4::glmer(y~ G + (-1+G1|OLRE1) + (-1+G2|OLRE2) ,family = binomial(link = "probit"), data = data,  nAGQ = nAGQ)
        model12 <- glm(y~ G , family = binomial(link = "probit"), data = data)
        as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 2))/4) + as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
      else{
        data$Gsum <- G1 + G2
        data$OLRE <- 1:nrow(data)
        data$OLRE[data$Gsum == 0] <- 0
        model11 <- GLMMadaptive::mixed_model(fixed = y ~ G, random = ~ -1 + Gsum|OLRE, data = data, family = binomial(link = 'probit'), control = list(nAGQ = nAGQ))
        model12 <- glm(y~ G ,family = binomial(link = "probit"), data = data)
        as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
    }
    #stop cluster
    parallel::stopCluster(cl)
    p_val
  }
  p_val
}





################ in most cases, better to use single approximation, but theoretically not valid if E has an effect (quite robust in simulation though)

#### With Z:
p_val11 <- NonAdditive_Simulator_With_Z(betaGE = c(0,0), size = 6000, num_trails = 100, parallel = T, betaE = 0.1, betaZ = 0)
p_val21 <- NonAdditive_Simulator_With_Z(betaGE = c(0.5,-0.5), size = 6000, num_trails = 100, parallel = T, betaE = 0.1, betaZ = 0)

p_val31 <- NonAdditive_Simulator_With_Z(betaGE = c(0,0), size = 3000, num_trails = 100, parallel = T, betaE = 0.5, betaZ = 0, nAGQ = 10L, SingleApproxi = TRUE)
p_val41 <- NonAdditive_Simulator_With_Z(betaGE = c(0.7,0.1), size = 3000, num_trails = 100, parallel = T, betaE = 0.5, betaZ = 0, nAGQ = 10L, SingleApproxi = TRUE)



#### Without Z: It seems like without Z variable, this method will not work at all...

p_val12 <- NonAdditive_Simulator_Without_Z(betaGE = c(0,0), size = 3000, num_trails = 100, parallel = T, betaE = 0.1)
p_val22 <- NonAdditive_Simulator_Without_Z(betaGE = c(0.1,-0.2), size = 3000, num_trails = 100, parallel = T, betaE = 0.1)

p_val32 <- NonAdditive_Simulator_Without_Z(betaGE = c(0,0), size = 3000, num_trails = 100, parallel = T, betaE = 0.5, nAGQ = 10L, SingleApproxi = TRUE)
p_val42 <- NonAdditive_Simulator_Without_Z(betaGE = c(0.1,-0.2), size = 3000, num_trails = 100, parallel = T, betaE = 0.5, nAGQ = 10L, SingleApproxi = TRUE)










