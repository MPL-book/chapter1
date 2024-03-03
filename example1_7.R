## Log normal survival model

#Generate event times
t <- rlnorm(n = 500, mean = 3, sd = 0.5)
# Generate censoring times
c <- rlnorm(n = 500, mean = 3.4, sd = 0.5)
# Identify censored individuals
delta <- as.numeric(t < c)
# Save minimum of event and censoring time
y <- t
y[which(delta == 0)] <- c[which(delta == 0)]
# Check censored proportion
sum(delta)/500

# 1. contribution to the estimation equation for mu
mu_MLE <- function(mu, sd, t, event) {
  u <- (log(t) - mu)/sd
  phi_u <- dnorm(u)
  sf = pnorm(u, lower.tail = FALSE)
  sf[sf<1e-30] = 1e-30
  contr <- sum(event*u + (1 - event)*phi_u/sf)
  return(contr)
}

# 2. contribution to the estimation equation for sd
sd_MLE <- function(sd, mu, t, event) {
  u <- (log(t) - mu)/sd
  phi_u <- dnorm(u)
  sf_u = pnorm(u, lower.tail = FALSE)
  sf_u[sf_u<1e-30] = 1e-20
  contr <- sum(event * (u^2 - 1) + (1-event) * u * phi_u/sf_u)
  return(contr)
}

# Solve two equations iteratively for mu and sd
#initial values
mu.old = 0
sd.old = 1
max.iter = 100

for (k in 1:max.iter) {
  #solve for mu
  solve_mu <- uniroot(mu_MLE, interval = c(-5, 5), 
                      sd = sd.old, t = y, event = delta
  )
  mu.new <- solve_mu$root
  
  # Solve for sd
  solve_sd <- uniroot(
    sd_MLE, interval = c(0.01, 5), 
    mu = mu.new, t = y, event = delta
  )
  sd.new <- solve_sd$root
  if (abs(mu.old-mu.new)<1e-5 & abs(sd.old-sd.new)<1e-5){
    break
  }
  else {
    mu.old = mu.new
    sd.old = sd.new
  }
  print(c(mu.new, sd.new))
}
mu_est = mu.new
sd_est = sd.new

# Estimated vs. true hazard and survival
# functions
v = seq(from = 0, to = max(t), by = 0.01)
u_est <- (log(v) - mu_est)/sd_est
u_true <- (log(v) - 3)/0.5
# Survival function
est_St <- pnorm(u_est, lower.tail = FALSE)
true_St <- pnorm(u_true, lower.tail = FALSE)
plot(est_St ~ v, type = "l", ylim = c(0, 1))
lines(true_St ~ v, type = "l", lty=2)
legend("topright", legend=c("Esti", "True"), lty=1:2, cex=0.9)
# Density function
est_f <- 1/(sd_est * v) * dnorm(u_est)
true_f <- 1/(0.5 * v) * dnorm(u_true)
# Hazard function
est_ht <- est_f/est_St
true_ht <- true_f/true_St
plot(est_ht ~ v, type = "l", ylim=c(0, 0.095))
lines(true_ht ~ v, type="l", lty=2) 
legend("topright", legend=c("Esti", "True"), lty=1:2, cex=0.9)

