set.seed(100)
N <- 1000
G <- sample(c(0,1,2),size = N, replace = T, prob = c(0.09,0.42,0.49))
Z <- rnorm(N,sd = 3)

## Eur:
beta0 <- -0.5
betaZ <- 0.8
betaG <- 0.3

ylat_Eur <- beta0 + betaG*G + betaZ*Z + rlogis(N)
y_Eur <- ifelse(ylat_Eur >=0, 1, 0)
table(y_Eur)
### Theoretical Power
mod_Eur <- glm(y_Eur~Z + G, family = binomial(link = "logit"))
#### Get the design matrix:
X <- cbind(rep(1,N),mod_Eur$model[,-1])
### Compute the weight matrix W:
beta <- c(beta0,betaZ,betaG)
#beta <- as.numeric(mod_Eur$coefficients)
w <- c()
for (i in 1:N) {
  si <- as.numeric(as.numeric(X[i,]) %*% beta)
  w[i] <- (dlogis(si)^2)/(plogis(si)*(1-plogis(si)))
}
I <- as.matrix(t(X)) %*% diag(w,nrow = N,ncol = N) %*% as.matrix(X)
#### Invert to get the true covariance matrix 
V <- solve(I)
### Compute the power function:
delta <- sqrt(1/V[3,3])*(0-beta[3])
alpha <- 0.05
Power_EU <- 1- pnorm(delta - qnorm(alpha/2)) + pnorm(delta + qnorm(alpha/2))















## Asia:
beta0 <- -0.5
betaZ <- 0.1
betaG <- 0.3

ylat_As <- beta0 + betaG*G + betaZ*Z + rlogis(N)
y_As <- ifelse(ylat_As >=0, 1, 0)
table(y_As)
### Theoretical Power
mod_As <- glm(y_As~ Z + G, family = binomial(link = "logit"))
#### Get the design matrix:
X <- cbind(rep(1,N),mod_As$model[,-1])
### Compute the weight matrix W:
beta <- c(beta0,betaZ,betaG)
#beta <- as.numeric(mod_As$coefficients)
w <- c()
for (i in 1:N) {
  si <- as.numeric(as.numeric(X[i,]) %*% beta)
  w[i] <- (dlogis(si)^2)/(plogis(si)*(1-plogis(si)))
}
I <- as.matrix(t(X)) %*% diag(w,nrow = N,ncol = N) %*% as.matrix(X)
#### Invert to get the true covariance matrix 
V <- solve(I)
### Compute the power function:
delta <- sqrt(1/V[3,3])*(0-beta[3])
alpha <- 0.05
Power_AS <- 1- pnorm(delta - qnorm(alpha/2)) + pnorm(delta + qnorm(alpha/2))




Power_EU
Power_AS







#### Illustrate that this power is correct:
set.seed(100)
p1 <- c()
for (i in 1:1000) {
  ylat <- beta0 + betaG*G + betaE*E + betaGE*G*E + rnorm(N)
  y <- ifelse(ylat >=0, 1, 0)
  mod <- glm(y~G+E+I(G*E), family = binomial(link = "probit"))
  p1[i] <- summary(mod)$coefficient[4,4]
}
emp_power <- mean(p1 <= alpha)

### This is fairly closed to the true power 0.17






















