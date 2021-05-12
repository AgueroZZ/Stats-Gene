library(tidyverse)

###### For a simulation of size n:

n = 300000

set.seed(123)



##### Generate random genotype for G1 and G2, and a normal environmental factor that is unknown:
G1 = apply(X = rmultinom(n,1,prob = c(0.2,0.5,0.3))>0, FUN = "which",MARGIN = 2)
G2 = apply(X = rmultinom(n,1,prob = c(0.3,0.5,0.2))>0, FUN = "which",MARGIN = 2)
E <- rnorm(n, mean = 3, sd = 1)

### Case 1: If the true model is nicely additive without interaction (Assuming probit model, inverse normal CDF as link function)
beta0 <- -5.5
beta1 <- 1.5
beta2 <- 1
beta3 <- 1

latent_y <- beta0 + beta1*G1 + beta2*G2 + rnorm(n = n)
y <- ifelse(latent_y>0,1,0)
data <- data.frame(y = y, G1 = G1, G2 = G2)

### Estimated probability of cases in different level of G1:
p <- data %>% group_by(G1) %>% summarise(p = mean(y))
knitr::kable(p)
### Put the estimated probability of being case into Inverse Normal CDF:
qnorm(p$p)
diff(qnorm(p$p)) #### Approximately linear increase with rate 1.22

### Estimated probability of cases in different level of G2:
p <- data %>% group_by(G2) %>% summarise(p = mean(y))
knitr::kable(p)
### Put the estimated probability of being case into Inverse Normal CDF:
qnorm(p$p)
diff(qnorm(p$p)) #### Approximately linear increase with rate 0.66
#### We don't really care how accurate are those rates compared to true beta, we care about whether they are linear. 


### Fit a probit model for the data:
model0 <- glm(y~., family = binomial(link = "probit"), data = data)
summary(model0)
### The estimated coefficients are very accurate as we expect. 



### Case 2: If the true model contains a interaction term with the unknown variable E:
set.seed(123)
latent_y <- beta0 + beta1*G1 + beta2*G2 + beta3*G1*E + rnorm(n = n)
y <- ifelse(latent_y>0,1,0)
data <- data.frame(y = y, G1 = G1, G2 = G2)


### Estimated probability of cases in different level of G1:
p <- data %>% group_by(G1) %>% summarise(p = mean(y))
knitr::kable(p)
### Put the estimated probability of being case into Inverse Normal CDF:
qnorm(p$p)
diff(qnorm(p$p)) #### Definitely not linear, one is 1.7431254 another is 0.7921436


### Estimated probability of cases in different level of G2:
p <- data %>% group_by(G2) %>% summarise(p = mean(y))
knitr::kable(p)
### Put the estimated probability of being case into Inverse Normal CDF:
qnorm(p$p)
diff(qnorm(p$p)) #### Not super linear, since 0.41 is a bit different from 0.51


### A finer grouping for the difference of proportion to be estimated more accurately:
### A weighted sum for G2's difference:
p <- data %>% group_by(G1,G2) %>% summarise(p = mean(y))
a1 <- diff(qnorm(p$p[1:3]))
a2 <- diff(qnorm(p$p[4:6]))
a3 <- diff(qnorm(p$p[7:9]))
((table(G1)[1]) * a1 + (table(G1)[2]) * a2 + (table(G1)[3]) * a3)/n
#### Seems to be linear as 0.465 close to 0.481, G2 is fine

p <- data %>% group_by(G2,G1) %>% summarise(p = mean(y))
a1 <- diff(qnorm(p$p[1:3]))
a2 <- diff(qnorm(p$p[4:6]))
a3 <- diff(qnorm(p$p[7:9]))
### A weighted sum:
((table(G2)[1]) * a1 + (table(G2)[2]) * a2 + (table(G2)[3]) * a3)/n
#### Definitely not linear as 1.7777863 much larger than 0.7721195, G1 has interaction

### Fit a misspecified model anyway:
model0 <- glm(y~., family = binomial(link = "probit"), data = data)
summary(model0)
### All the coefficients are biased.







###### Consider a random slope type of method to capture interaction:
library(lme4)
data$OLRE <- 1:nrow(data)
model1 <- glmer(y~ G1 + (0+G1|OLRE) ,family = binomial(link = "probit"), data = data)
summary(model1)
#### The variance of the random slope seems to be very large, suggesting the presence of interaction


### If the true model does not have interaction
latent_y <- beta0 + beta1*G1 + beta2*G2 + rnorm(n = n)
y <- ifelse(latent_y>0,1,0)
data <- data.frame(y = y, G1 = G1, G2 = G2)
data$OLRE <- 1:nrow(data)
model2 <- glmer(y~ G1 + (0+G1|OLRE) ,family = binomial(link = "probit"), data = data)
summary(model2)
#### The variance of the random slope seems to be very small, suggesting no interaction is present.



























