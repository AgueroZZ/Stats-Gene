library(tidyverse)
library(lme4)

set.seed(123)
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



###### Method 2: Random Slope:
set.seed(12345)
latent_y <- beta0 + beta1*G + beta2*Z + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G = G, Z = Z)
data$OLRE <- 1:nrow(data)
model11 <- glmer(y~ G + Z + (-1+G|OLRE) ,family = binomial(link = "probit"), data = data,  nAGQ = 25L)
summary(model11)
model12 <- glm(y~ G + Z ,family = binomial(link = "probit"), data = data)
summary(model12)
as.numeric((1 - pchisq(2*(logLik(model11) - logLik(model12)), df = 1))/2)




set.seed(12345)
latent_y <- beta0 + beta1*G + beta2*Z + beta3*G*E + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G = G, Z = Z)
data$OLRE <- 1:nrow(data)
model21 <- glmer(y~ G + Z + (-1+G|OLRE) ,family = binomial(link = "probit"), data = data, nAGQ = 25L)
summary(model21)
model22 <- glm(y~ G + Z ,family = binomial(link = "probit"), data = data)
summary(model22)
as.numeric((1 - pchisq(2*(logLik(model21) - logLik(model22)), df = 1))/2)

