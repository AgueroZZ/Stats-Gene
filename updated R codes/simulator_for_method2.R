



###### Additive Model:

Additive_Simulator <- function(beta0 = -1.2, betaG = 0.8, betaE = 0.6, betaGE = 0.2, muE = 0, sdE = 3, p1 = 0.3, p2 = NULL, betaZ = NULL, size = 3000, num_trails = 1000, parallel = FALSE, nAGQ = 1L){
  n <- size
  q1 <- 1 - p1
  q2 <- 1 - p2
  if(parallel == F){
    p_val <- c()
    for (i in 1:num_trails) {
      G = apply(X = rmultinom(n,1,prob = c(p1^2,2*p1*q1,q1^2)) > 0, FUN = "which",MARGIN = 2) - 1
      E <- rnorm(n,muE,sdE)
      
      if(is.null(p2) && is.null(betaZ)){
        ylat <- beta0 + betaG*G +  betaE*E + betaGE*E*G + rnorm(n)
        y <- ifelse(ylat > 0,1,0)
        data <- data.frame(y = y, G = G)
        data$OLRE <- 1:nrow(data)
        model11 <- suppressMessages(glmer(y~ G + (-1+G|OLRE), family = binomial(link = "probit"), data = data,  nAGQ = nAGQ))
        model12 <- glm(y~ G, family = binomial(link = "probit"), data = data)
        p_val[i] <- as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
      else{
        Z <- apply(X = rmultinom(n,1,prob = c(p2^2,2*p2*q2,q2^2)) > 0, FUN = "which",MARGIN = 2) - 1
        ylat <- beta0 + betaG*G +  betaE*E + betaGE*E*G + betaZ*Z + rnorm(n) 
        y <- ifelse(ylat > 0,1,0)
        data <- data.frame(y = y, G = G, Z = Z)
        data$OLRE <- 1:nrow(data)
        model11 <- suppressMessages(glmer(y~ G + Z + (-1+G|OLRE), family = binomial(link = "probit"), data = data,  nAGQ = nAGQ))
        model12 <- glm(y~ G + Z, family = binomial(link = "probit"), data = data)
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
      E <- rnorm(n,muE,sdE)
      if(is.null(p2) && is.null(betaZ)){
        ylat <- beta0 + betaG*G +  betaE*E + betaGE*E*G + rnorm(n)
        y <- ifelse(ylat > 0,1,0)
        data <- data.frame(y = y, G = G)
        data$OLRE <- 1:nrow(data)
        model11 <- lme4::glmer(y~ G + (-1+G|OLRE), family = binomial(link = "probit"), data = data,  nAGQ = nAGQ)
        model12 <- glm(y~ G, family = binomial(link = "probit"), data = data)
        as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
      else{
        Z <- apply(X = rmultinom(n,1,prob = c(p2^2,2*p2*q2,q2^2)) > 0, FUN = "which",MARGIN = 2) - 1
        ylat <- beta0 + betaG*G +  betaE*E + betaGE*E*G + betaZ*Z + rnorm(n) 
        y <- ifelse(ylat > 0,1,0)
        data <- data.frame(y = y, G = G, Z = Z)
        data$OLRE <- 1:nrow(data)
        model11 <- lme4::glmer(y~ G + Z + (-1+G|OLRE), family = binomial(link = "probit"), data = data,  nAGQ = nAGQ)
        model12 <- glm(y~ G + Z, family = binomial(link = "probit"), data = data)
        as.numeric((1 - pchisq(2*(stats::logLik(model11) - stats::logLik(model12)), df = 1))/2)
      }
    }
    #stop cluster
    parallel::stopCluster(cl)
    p_val
  }
}


























###### Genotypic Model:
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















#### Check the robustness of additive model:
### Easier case with 8000 obs, without effect of E
### Without interaction
add_p1 <- Additive_Simulator(betaGE = 0, size = 8000, mu = 0, num_trails = 100, parallel = TRUE, betaE = 0, nAGQ = 8)
sum(add_p1 <= 0.05)/length(add_p1)
### with interaction
add_p2 <- Additive_Simulator(betaGE = 0.3, size = 8000, mu = 0, num_trails = 100, parallel = TRUE, betaE = 0, nAGQ = 8)
sum(add_p2 <= 0.05)/length(add_p2)


### Medium case with 8000 obs, with effect of E
### Without interaction
add_p1 <- Additive_Simulator(betaGE = 0, size = 8000, mu = 0, num_trails = 100, parallel = TRUE, betaE = 0.2, nAGQ = 8)
sum(add_p1 <= 0.05)/length(add_p1)
### with interaction
add_p2 <- Additive_Simulator(betaGE = 0.3, size = 8000, mu = 0, num_trails = 100, parallel = TRUE, betaE = 0.2, nAGQ = 8)
sum(add_p2 <= 0.05)/length(add_p2)



### Harder case with 8000 obs, with large effect of E
### Without interaction
add_p1 <- Additive_Simulator(betaG = 0.4, betaGE = 0, size = 8000, mu = 0, num_trails = 100, parallel = TRUE, betaE = 0.4, nAGQ = 8)
sum(add_p1 <= 0.05)/length(add_p1)
### with interaction
add_p2 <- Additive_Simulator(betaG = 0.4, betaGE = 0.3, size = 8000, mu = 0, num_trails = 100, parallel = TRUE, betaE = 0.4, nAGQ = 8)
sum(add_p2 <= 0.05)/length(add_p2)



### If we switch the sign, the robustness should still hold 
### Without interaction
add_p1 <- Additive_Simulator(betaG = 0.6, betaGE = 0, size = 8000, mu = 0, num_trails = 100, parallel = TRUE, betaE = 0, nAGQ = 8)
sum(add_p1 <= 0.05)/length(add_p1)
### with interaction
add_p2 <- Additive_Simulator(betaG = 0.6, betaGE = -0.6, size = 8000, mu = 0, num_trails = 100, parallel = TRUE, betaE = 0.1, nAGQ = 8)
sum(add_p2 <= 0.05)/length(add_p2)



### Seems not robust if both interaction and betaE are large













