library(tidyverse)

# set.seed(100,sample.kind = "Rounding")
set.seed(100)
N <- 15000


#### Note the two diseases are analyzed at two independent sets of individual
G <- sample(c(0,1,2),size = N, replace = T, prob = c(0.49,0.42,0.09))
E <- rnorm(N,sd = 3)
beta0 <- -0.5
betaG <- 0
betaE <- 0.95
ylat_A <- beta0 + betaG*G + betaE*E + rlogis(N)
Lung <- ifelse(ylat_A >=0, 1, 0)
LungData <- data.frame(G = G, E = E, Y = Lung)

G <- sample(c(0,1,2),size = N, replace = T, prob = c(0.49,0.42,0.09))
E <- rnorm(N,sd = 3)
beta0 <- -0.5
betaG <- 0.1
betaE <- 0
ylat_B <- beta0 + betaG*G + betaE*E + rlogis(N)
Breast <- ifelse(ylat_B >=0, 1, 0)
BreastData <- data.frame(G = G, E = E, Y = Breast)


##### Each fitted Model
modLung <- glm(Y~G+E, family = binomial('logit'), data = LungData)
summary(modLung)

modBreast <- glm(Y~G+E, family = binomial('logit'), data = BreastData)
summary(modBreast)


######## First model is not significant, but the second one is, what if we 
######## combine these P values using meta-analysis methods
pvec <- c(summary(modLung)$coefficients[2,4], summary(modBreast)$coefficients[2,4])

### Just take the minimum:
min(pvec)

### Fisher two sided test:
QU <- -2 * log(prod(pvec))
m = 2
pchisq(QU, df = 2*m, lower.tail = F)

### Weighted Z method: with sample size being weight:
weights <- c(N,N)
Zvec <- c(summary(modLung)$coefficients[2,3], summary(modBreast)$coefficients[2,3])
weighted_test <- sum(weights*Zvec/sqrt(sum(weights^2)))
1 - pnorm(weighted_test)

### Weighted Z method: with SE being weight:
weights <- c(1/summary(modLung)$coefficients[2,2],1/summary(modBreast)$coefficients[2,2])
Zvec <- c(summary(modLung)$coefficients[2,3], summary(modBreast)$coefficients[2,3])
weighted_test <- sum(weights*Zvec/sqrt(sum(weights^2)))
1 - pnorm(weighted_test)


### What if we use betaE as weight?
weights <- abs(c(1/summary(modLung)$coefficients[3,1],1/summary(modBreast)$coefficients[3,1]))
Zvec <- c(summary(modLung)$coefficients[2,3], summary(modBreast)$coefficients[2,3])
weighted_test <- sum(weights*Zvec/sqrt(sum(weights^2)))
1 - pnorm(weighted_test)


###### Aggregate our comparison:
compare_methods <- function(betaG = c(0.1,0.1), betaE = c(0.95,0.05), beta0 = c(-0.5, -0.5), p = 0.3, N = 15000, M = 100){
  agg <- data.frame()
  for (i in 1:M) {
    result <- numeric(5)
    q <- 1 - p
    G <- sample(c(0,1,2),size = N, replace = T, prob = c(p^2,2*p*q,q^2))
    E <- rnorm(N,sd = 3)
    ylat_A <- beta0[1] + betaG[1]*G + betaE[1]*E + rlogis(N)
    Lung <- ifelse(ylat_A >=0, 1, 0)
    LungData <- data.frame(G = G, E = E, Y = Lung)
    G <- sample(c(0,1,2),size = N, replace = T, prob = c(p^2,2*p*q,q^2))
    E <- rnorm(N,sd = 3)
    ylat_B <- beta0[2] + betaG[2]*G + betaE[2]*E + rlogis(N)
    Breast <- ifelse(ylat_B >=0, 1, 0)
    BreastData <- data.frame(G = G, E = E, Y = Breast)
    ##### Each fitted Model
    modLung <- glm(Y~G+E, family = binomial('logit'), data = LungData)
    modBreast <- glm(Y~G+E, family = binomial('logit'), data = BreastData)
    pvec <- c(summary(modLung)$coefficients[2,4], summary(modBreast)$coefficients[2,4])
    result[1] <- min(pvec)
    QU <- -2 * log(prod(pvec))
    result[2] <- pchisq(QU, df = 4, lower.tail = F)
    weights <- c(N,N)
    Zvec <- c(summary(modLung)$coefficients[2,3], summary(modBreast)$coefficients[2,3])
    weighted_test <- sum(weights*Zvec/sqrt(sum(weights^2)))
    result[3] <- 2 * (1 - pnorm(abs(weighted_test)))
    weights <- c(1/summary(modLung)$coefficients[2,2],1/summary(modBreast)$coefficients[2,2])
    Zvec <- c(summary(modLung)$coefficients[2,3], summary(modBreast)$coefficients[2,3])
    weighted_test <- sum(weights*Zvec/sqrt(sum(weights^2)))
    result[4] <- 2 * (1 - pnorm(abs(weighted_test)))
    weights <- abs(c(1/summary(modLung)$coefficients[3,1],1/summary(modBreast)$coefficients[3,1]))
    Zvec <- c(summary(modLung)$coefficients[2,3], summary(modBreast)$coefficients[2,3])
    weighted_test <- sum(weights*Zvec/sqrt(sum(weights^2)))
    result[5] <- 2 * (1 - pnorm(abs(weighted_test)))
    agg <- rbind(agg, result)
  }
  colnames(agg) <- c("Min","Fisher","WeightN", "WeightSE", "WeightE")
  agg
}
count_bumber_of_optimal <- function(agg){
  nameVec <- names(agg)
  ID <- agg %>% apply(., 1, which.min)
  table(nameVec[ID])
}
count_bumber_of_TypeIerror <- function(agg){
  error <- agg <= 0.05
  apply(error, 2, mean)
}

Aggregate_discordant <- compare_methods(M = 30, betaG = c(0.1,-0.1))
count_bumber_of_optimal(Aggregate_discordant[,-1])

Aggregate_One_True <- compare_methods(M = 30, betaG = c(0.1,0))
count_bumber_of_optimal(Aggregate_One_True[,-1])

Aggregate_One_True <- compare_methods(M = 30, betaG = c(0,0.1))
count_bumber_of_optimal(Aggregate_One_True[,-1])

Aggregate_Concordant <- compare_methods(M = 30, betaG = c(0.1,0.1))
count_bumber_of_optimal(Aggregate_Concordant[,-1])

Aggregate_Concordant2 <- compare_methods(M = 30, betaG = c(-0.1,-0.1))
count_bumber_of_optimal(Aggregate_Concordant2[,-1])

Aggregate_Concordant_noE <- compare_methods(M = 30, betaG = c(0.1,0), betaE = c(0,0))
count_bumber_of_optimal(Aggregate_Concordant_noE[,-1])

Aggregate_TypeIError <- compare_methods(M = 30, betaG = c(0,0), betaE = c(0,0))
count_bumber_of_TypeIerror(Aggregate_TypeIError[,-1])

