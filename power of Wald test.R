set.seed(100)
N <- 1000
G <- sample(c(0,1,2),size = N, replace = T, prob = c(0.16,0.48, 0.36))
E <- rnorm(N)
beta0 <- -1
betaG <- 0.3
betaE <- 1
betaGE <- 0.1
ylat <- beta0 + betaG*G + betaE*E + betaGE*G*E + rnorm(N)
y <- ifelse(ylat >=0, 1, 0)

mod <- glm(y~G+E+I(G*E), family = binomial(link = "probit"))


#### Get the design matrix:
X <- cbind(rep(1,N),mod$model[,-1])

### Compute the weight matrix W:
beta <- c(beta0,betaG,betaE,betaGE)
#beta <- as.numeric(mod$coefficients)

w <- c()
for (i in 1:N) {
  si <- as.numeric(as.numeric(X[i,]) %*% beta)
  w[i] <- (dnorm(si)^2)/(pnorm(si)*(1-pnorm(si)))
}

### The true information matrix
I <- as.matrix(t(X)) %*% diag(w,nrow = N,ncol = N) %*% as.matrix(X)



#### Invert to get the true covariance matrix 
V <- solve(I)
#### Compare with the covariance matrix estimated:
summary(mod)$cov.scaled


### Compute the power function 

### Assume d = beta0 - beta1, where beta0 = 0 is what we are testing as null:
delta <- sqrt(1/V[4,4])*(0-beta[4])
alpha <- 0.05
Power <- 1- pnorm(delta - qnorm(alpha/2)) + pnorm(delta + qnorm(alpha/2))


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

















######################## Repeat when the main effect of E is zero
set.seed(100)
G <- sample(c(0,1,2),size = N, replace = T, prob = c(0.16,0.48, 0.36))
E <- rnorm(N)
beta0 <- -1
betaG <- 0.3
betaE <- 0
betaGE <- 0.1
ylat <- beta0 + betaG*G + betaE*E + betaGE*G*E + rnorm(N)
y <- ifelse(ylat >=0, 1, 0)

mod2 <- glm(y~G+E+I(G*E), family = binomial("probit"))


#### Get the design matrix:
X <- cbind(rep(1,N),mod2$model[,-1])

### Compute the weight matrix W:
beta <- c(beta0,betaG,betaE,betaGE)


w <- c()
for (i in 1:N) {
  si <- as.numeric(as.numeric(X[i,]) %*% beta)
  w[i] <- (dnorm(si)^2)/(pnorm(si)*(1-pnorm(si)))
}

### The true information matrix
I2 <- as.matrix(t(X)) %*% diag(w,nrow = N,ncol = N) %*% as.matrix(X)
#### Invert to get the true covariance matrix 
V2 <- solve(I2)
#### Compare with the covariance matrix estimated:
summary(mod2)$cov.scaled


### Compute power
### Assume d = beta0 - beta1, where beta0 = 0 is what we are testing as null:
delta2 <- sqrt(1/V2[4,4])*(0-beta[4])
alpha <- 0.05
Power2 <- 1- pnorm(delta2 - qnorm(alpha/2)) + pnorm(delta2 + qnorm(alpha/2))


#### Illustrate that this power is correct:
set.seed(100)
p2 <- c()
for (i in 1:1000) {
  ylat <- beta0 + betaG*G + betaE*E + betaGE*G*E + rnorm(N)
  y <- ifelse(ylat >=0, 1, 0)
  mod <- glm(y~G+E+I(G*E), family = binomial(link = "probit"))
  p2[i] <- summary(mod)$coefficient[4,4]
}

emp_power2 <- mean(p2 <= alpha)

### Again, fairly closed to the theoretical power 0.35







