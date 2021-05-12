compute_effect_add <- function(beta0,betaG,betaGE,betaE, sig = 1){
  G <- c(0,1,2)
  diff(c((beta0+betaG*G)/sqrt((betaGE*G*sig)^2 + (betaE*sig)^2 + 1)))
}

compute_effect_add(beta0 = -2, betaG = 0.3, betaGE = 0.5, betaE = 0,sig = 3)



my_test_G <- G_using[,POS_using %in% POS_EFF]
my_test_E <- E
my_test_cases <- case_new

for (i in 1:ncol(my_test_G)) {
  modi <- glm(my_test_cases~factor(my_test_G[,i]), family = binomial(link = "probit"))
  print(as.numeric(wald.test(vcov(modi)[-1,-1], b = modi$coefficients[-1], H0 = matrix(0,nrow = 1,ncol = 1), L = matrix(c(2,-1),nrow = 1))$result$chi2[3]))
}