library(dplyr)
library(nlme)
library(aod)

set.seed(123)
n = 3000
p1 <- 0.2
q1 <- 1 - p1

G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
E <- rnorm(n, mean = 0, sd = 2)

beta0 <- -1.2
beta1 <- 0.8
beta2 <- 0.3
beta3 <- 0.8

latent_y <- beta0 + beta1*G + beta2*E + beta3*G*E + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G = G)


##### Stage 1:

model_1 <- lm(y~factor(G),data = data)
p <- data %>% group_by(G) %>% summarise(p = mean(y))

data$d <- qnorm(model_1$fitted.values)

#### Stage 2:
v <- model_1$fitted.values*(1-model_1$fitted.values)/(dnorm(qnorm(model_1$fitted.values))^2)
# for (i in 1:length(v)) {
# v[i] <- v[i]/sum(data$G == data$G[i])
# }

model_2 <- lm(d~factor(G),data = data, x = TRUE)
summary(model_2)
Sig_inv <- diag(1/v,ncol = length(v))
X <- as.matrix(model_2$x)
covmatrix <- solve((t(X) %*% Sig_inv %*% X))[-1,-1]

print(wald.test(covmatrix, b = model_2$coefficients[-1], H0 = matrix(0,nrow = 1,ncol = 1), L = matrix(c(2,-1),nrow = 1)),,digits = 4)

print(Test_one_G(data,"y","G"),digits = 4)




###############
model_try <- glm(y~factor(G), data = data, family = binomial(link = "probit"))







