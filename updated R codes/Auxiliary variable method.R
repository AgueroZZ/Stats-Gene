library(tidyverse)



set.seed(123)
n = 6000
p1 <- 0.5
q1 <- 1 - p1



G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
Z = rnorm(n,0,3)
E <- rnorm(n, mean = 0, sd = 3)

G1 <- ifelse(G == 1, 1,0)
G2 <- ifelse(G == 2, 1,0)
G0 <- ifelse(G1 == G2, 1,0)


### No interaction:

beta0 <- -1.2
betaG1 <- 0.8
betaG2 <- 0.3
betaG1E <- 0
betaG2E <- 0
betaZ <- 0.5
betaE <- 0.5

latent_y <- beta0 + betaG1*G1 + betaG2*G2 + betaE*E + betaG1E*G1*E + betaG2E*G2*E  + betaZ*Z + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G),Z=Z)


mod1 <- glm(y~G+G*Z+Z, family = binomial(link = "probit"), data=data)
mod2 <- glm(y~G+Z, family = binomial(link = "probit"), data=data)
1 - pchisq(2*(logLik(mod1) - logLik(mod2)), df = 2)



### With interaction:
beta0 <- -1.2
betaG1 <- 0.8
betaG2 <- 0.3
betaG1E <- 0.1
betaG2E <- -0.1
betaZ <- 0.5
betaE <- 0.5

latent_y <- beta0 + betaG1*G1 + betaG2*G2 + betaE*E + betaG1E*G1*E + betaG2E*G2*E  + betaZ*Z + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G),Z=Z)


mod1 <- glm(y~G+G*Z+Z, family = binomial(link = "probit"), data=data)
mod2 <- glm(y~G+Z, family = binomial(link = "probit"), data=data)
1 - pchisq(2*(logLik(mod1) - logLik(mod2)), df = 2)







#### Simulator:
Aux_simulator <- function(beta0 = -1.2, betaZ = 0.6, betaG = c(0.8, 0.3), betaE = 0.6, betaGE = c(0.6,0.4), muZ = 0, muE = 0, sdE = 1, p1 = 0.5, sdZ = 1, size = 3000, num_trails = 1000){
  betaG1 <- betaG[1]
  betaG2 <- betaG[2]
  betaG1E <- betaGE[1]
  betaG2E <- betaGE[2]
  n = size
  q1 <- 1 - p1
  p_val <- c()
  for (i in 1:num_trails) {
    G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
    Z = rnorm(n,muZ,sdZ)
    E <- rnorm(n, mean = muE, sd = sdE)
    
    G1 <- ifelse(G == 1, 1,0)
    G2 <- ifelse(G == 2, 1,0)
    G0 <- ifelse(G1 == G2, 1,0)
    
    
    latent_y <- beta0 + betaG1*G1 + betaG2*G2 + betaE*E + betaG1E*G1*E + betaG2E*G2*E  + betaZ*Z + rnorm(n = n)
    y <- ifelse(latent_y > 0,1,0)
    data <- data.frame(y = y, G1 = G1, G2 = G2, G = as.factor(G),Z=Z)
    
    
    mod1 <- glm(y~G+G*Z+Z, family = binomial(link = "probit"), data=data)
    mod2 <- glm(y~G+Z, family = binomial(link = "probit"), data=data)
    p_val[i] <- 1 - pchisq(2*(logLik(mod1) - logLik(mod2)), df = 2)
  }
  p_val
}






###: When the auxiliary variable has a weak effect: almost uniform p-value
p1 <- Aux_simulator(betaZ = 0.1, betaGE = c(0.3,0.2),sdZ = 1)
hist(p1, breaks = 30)
sum(p1 <= 0.05)/length(p1)


###: When the auxiliary variable has a stronger effect: higher power
p2 <- Aux_simulator(betaZ = 0.4, betaGE = c(0.3,0.2),sdZ = 1)
hist(p2, breaks = 30)
sum(p2 <= 0.05)/length(p2)


###: When the auxiliary variable has a stronger effect and higher standard deviation: higher power
p3 <- Aux_simulator(betaZ = 0.4, betaGE = c(0.3,0.2),sdZ = 3)
hist(p3, breaks = 30)
sum(p3 <= 0.05)/length(p3)


###: When the auxiliary variable has a stronger effect and higher standard deviation, and the interaction effects differ strongly between groups: much higher power
p4 <- Aux_simulator(betaZ = 0.4, betaGE = c(0.3,-0.2),sdZ = 3)
hist(p4, breaks = 30)
sum(p4 <= 0.05)/length(p3)

