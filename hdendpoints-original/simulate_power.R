possible.ns <- seq(from=100, to=2000, by=50)
powers <- rep(NA, length(possible.ns))
powers.cov <- rep(NA, length(possible.ns))        # Need a second empty vector
alpha <- 0.05
sims <- 500
for (j in 1:length(possible.ns)){
  N <- possible.ns[j]
  
  significant.experiments <- rep(NA, sims)
  significant.experiments.cov <- rep(NA, sims)      # Need a second empty vector here too
  
  for (i in 1:sims){
    gender <- c(rep("F", N/2), rep("M", N/2))       # Generate "gender" covariate
    age <- sample(x=18:65, size=N, replace=TRUE)    # Generate "age" covariate
    effectofgender <- 10                            # Hypothesize the "effect" of gender on income
    effectofage <- 2                                # Hypothesize the "effect" of age on income
    
    ## Hypothesize Control Outcome as a function of gender, age, and error
    Y0 <- effectofgender*(gender=="M") + effectofage*age + rnorm(n=N, mean=100, sd=20)
    
    ## This is all the same ##
    tau <- 5
    Y1 <- Y0 + tau
    Z.sim <- rbinom(n=N, size=1, prob=.5)
    Y.sim <- Y1*Z.sim + Y0*(1-Z.sim)
    fit.sim <- lm(Y.sim ~ Z.sim)
    
    ## This is the novel analysis -- including two covariates to increase precision ##
    fit.sim.cov <- lm(Y.sim ~ Z.sim + (gender=="M") + age)
    
    ## extract p-values and calculate significance ##
    p.value <- summary(fit.sim)$coefficients[2,4]
    p.value.cov <- summary(fit.sim.cov)$coefficients[2,4]
    significant.experiments[i] <- (p.value <= alpha)
    significant.experiments.cov[i] <- (p.value.cov <= alpha)
  }
  
  powers[j] <- mean(significant.experiments)
  powers.cov[j] <- mean(significant.experiments.cov)
}

plot(possible.ns, powers, ylim=c(0,1))
points(possible.ns, powers.cov, col="red")



