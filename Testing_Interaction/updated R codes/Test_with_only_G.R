library(tidyverse)

### Function to test one G's interaction:
Test_one_G <- function(data,case,G){
  mydata <- data.frame(y = data[[case]], G = data[[G]])
  n_0 <- sum(mydata$G == 0)
  n_1 <- sum(mydata$G == 1)
  n_2 <- sum(mydata$G == 2)
  p_0 <- sum(mydata$y[mydata$G == 0] == 1)/n_0
  p_1 <- sum(mydata$y[mydata$G == 1] == 1)/n_1
  p_2 <- sum(mydata$y[mydata$G == 2] == 1)/n_2
  v_0 <- (p_0*(1 - p_0)/n_0) * (1/dnorm(qnorm(p_0))^2)
  v_1 <- (p_1*(1 - p_1)/n_1) * (1/dnorm(qnorm(p_1))^2)
  v_2 <- (p_2*(1 - p_2)/n_2) * (1/dnorm(qnorm(p_2))^2)
  Teststats <- ((qnorm(p_2) - 2*qnorm(p_1) + qnorm(p_0))^2)/(v_0 + 4*v_1 + v_2)
  1 - pchisq(Teststats,df = 1)
}


#### Try it
set.seed(123)
n = 5000
p1 <- 0.2
q1 <- 1 - p1

G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
E <- rnorm(n, mean = 0, sd = 1)

beta0 <- -1.2
betaG <- 0.8
betaE <- 0.6
betaGE <- 0.6


### Without interaction:
latent_y <- beta0 + betaG*G + betaE*E  + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G = G)
p <- data %>% group_by(G) %>% summarise(p = mean(y))

Test_one_G(data = data, case = "y", G = "G")


### With interaction:
latent_y <- beta0 + betaG*G + betaE*E + betaGE * G*E + rnorm(n = n)
y <- ifelse(latent_y > 0,1,0)
data <- data.frame(y = y, G = G)
p <- data %>% group_by(G) %>% summarise(p = mean(y))

Test_one_G(data = data, case = "y", G = "G")





##### Take a look at the p-values:
Simulator_One_G <- function(beta0 = -1.2, betaG = 0.8, betaE = 0.6, betaGE = 0.6, muE = 0, sdE = 1, p1 = 0.2, size = 3000, num_trails = 1000){
  p_val <- numeric()
  q1 <- 1 - p1
  for (i in 1:num_trails) {
    n = size
    G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
    E <- rnorm(n, mean = muE, sd = sdE)
    latent_y <- beta0 + betaG*G + betaE*E + betaGE * G*E + rnorm(n = n)
    y <- ifelse(latent_y > 0,1,0)
    data <- data.frame(y = y, G = G)
    p_val[i] <- Test_one_G(data,case = 'y', G = 'G')
  }
  p_val
}

#### Without interaction:
set.seed(123)
p_val1 <- Simulator_One_G(betaGE = 0, size = 3000, num_trails = 5000)
hist(p_val1,breaks = 20, freq = F)

#### With large interaction:
set.seed(123)
p_val2 <- Simulator_One_G(betaGE = 1.2, size = 3000, num_trails = 5000)
hist(p_val2,breaks = 20, freq = F)


#### With medium interaction:
set.seed(123)
p_val2 <- Simulator_One_G(betaGE = 0.6, size = 3000, num_trails = 5000)
hist(p_val2,breaks = 20, freq = F)


#### With small interaction:
set.seed(123)
p_val2 <- Simulator_One_G(betaGE = 0.03, size = 3000, num_trails = 5000)
hist(p_val2,breaks = 20, freq = F)

