## Effect of mis-specification of the baseline hazard

# Write simulate() function to repeatedly generate a sample
simulate <- function(){
  kappa <- 3
  gamma <- 1
  neglogU <- -log(runif(50))
  
  x1 <- rbinom(50, 1, 0.5)
  x2 <- runif(50, 0, 1)
  X <- cbind(x1, x2)
  beta_true = c(0.5, 0.5)
  eXtB = exp(X %*% beta_true)
  
  y <- gamma * (exp(neglogU/eXtB) - 1)^(1/kappa)
  event <- rep(1, 50)
  t <- y
  
  df = data.frame(t, event)
  df = cbind(df, X)
  return(df)
}

#Load flexsurv package for Weibull PH model
library(flexsurv)

#Run simulation
save <- matrix(0, nrow = 300, ncol = 12)
for(s in 1:300){
  df <- simulate()
  
  weib_fit <- flexsurvreg(Surv(t, event) ~ x1 + x2, data = df, dist = "weibullPH")
  cox_fit <- coxph(Surv(t, event) ~ x1 + x2, data = df)
  
  save[s,1] <- weib_fit$coefficients[3]
  save[s,2] <- sqrt(weib_fit$cov[3,3])
  save[s,3] <- as.numeric(weib_fit$coefficients[3] - 
                            1.96 * sqrt(weib_fit$cov[3,3]) < (0.5) &
                            weib_fit$coefficients[3] + 
                            1.96 * sqrt(weib_fit$cov[3,3]) > (0.5) )
  
  save[s,4] <- weib_fit$coefficients[4]
  save[s,5] <- sqrt(weib_fit$cov[4,4])
  save[s,6] <- as.numeric(weib_fit$coefficients[4] - 
                            1.96 * sqrt(weib_fit$cov[4,4]) < (0.5) &
                            weib_fit$coefficients[4] + 
                            1.96 * sqrt(weib_fit$cov[4,4]) > (0.5) )
  
  save[s,7] <- cox_fit$coefficients[1]
  save[s,8] <- sqrt(cox_fit$var[1,1])
  save[s,9] <- as.numeric(cox_fit$coefficients[1] - 
                            1.96 * sqrt(cox_fit$var[1,1]) < (0.5) &
                            cox_fit$coefficients[1] + 
                            1.96 * sqrt(cox_fit$var[1,1]) > (0.5) )
  
  save[s,10] <- cox_fit$coefficients[2]
  save[s,11] <- sqrt(cox_fit$var[2,2])
  save[s,12] <- as.numeric(cox_fit$coefficients[2] - 
                             1.96 * sqrt(cox_fit$var[2,2]) < (0.5) &
                             cox_fit$coefficients[2] + 
                             1.96 * sqrt(cox_fit$var[2,2]) > (0.5) )
}


results <- apply(save, 2, mean)
results <- round(results, 4)
results <- matrix(results, nrow = 4, byrow = TRUE)

