## MPL estimate of hazard function

# Generate data with 30% censoring
kappa <- 5
gamma <- 1
neglogU <- -log(runif(500))
t <- gamma*(exp(neglogU) - 1)^(1/kappa)
c <- rexp(500, rate = 1/3)
event <- as.numeric(t < c)
y <- t
y[which(event == 0)] <- c[which(event == 0)]
sum(event)/500

# For each of m=8,20,40,60 and with \lambda = 0, 
#compute ML estimate of theta from (1.57)

estimate_theta <- function(m, y, event){
  w <- as.numeric(quantile(y, probs = seq(0, 1, 
                                          length.out = (m+1)))) 
  delta_v <- w[2:(m+1)] - w[1:m]
  theta = rep(0, m)
  
  for(v in 1:m){
    n_u <- sum(as.numeric(w[v] <= y & y < w[(v+1)]))
    N_u <- sum(as.numeric(w[v] <= t))
    theta[v] <- n_u/(N_u * delta_v[(v)])
  }
  out = list(theta = theta, w = w)
  return(out)
}

theta_m8 <- estimate_theta(m = 8, y, event)
theta_m20 <- estimate_theta(m = 20, y, event)
theta_m40 <- estimate_theta(m = 40, y, event)
theta_m60 <- estimate_theta(m = 60, y, event)

# Plot obtained theta i.e. discrete hazard function values

par(mfrow = c(2,2))
m=8
plot(c(theta_m8$theta, theta_m8$theta[m]) ~ theta_m8$w, 
     type = "s", col = "black", ylab = "hazard est", 
     xlab = "Time", main = "MLE of Theta; m = 8")
m=20
plot(c(theta_m20$theta, theta_m20$theta[m]) ~ theta_m20$w, 
     type = "s", col = "black", ylab = "hazard est", 
     xlab = "Time", main = "MLE of Theta; m = 20")
m=40
plot(c(theta_m40$theta, theta_m40$theta[m]) ~ theta_m40$w, 
     type = "s", col = "black", ylab = "hazard est", 
     xlab = "Time", main = "MLE of Theta; m = 40")
m=60
plot(c(theta_m60$theta, theta_m60$theta[m]) ~ theta_m60$w, 
     type = "s", col = "black", ylab = "hazard est", 
     xlab = "Time", main = "MLE of Theta; m = 60")

# compute MPL estimate with m = 60 and 
# smoothing values lambda = 1e-6, 1,10,50

# For m=60, create penalty matrix R
m <- 60
D <- matrix(0, (m), (m))
for(u in 1:(m-1)){
  D[u, (u+1)] = 1
  D[u, u] = -1
}
R <- t(D) %*% D

# Solve equation 1.63 for different smoothing values 
# using nloptr

library(nloptr)
eval_f0 <- function(theta, y, w, delta_v, R, lambda){
  m = length(theta)
  alld_u <- alln_u <- rep(0, m)
  for(v in 1:m){
    n_u <- sum(as.numeric(w[v] <= y & y < w[(v+1)]))
    d_u <- sum(as.numeric(w[v] <= t & t < w[(v+1)]))
    alln_u[v] = n_u
    alld_u[v] = d_u
  }
  theta[theta<1e-30] = 1e-30
  cumhaz = cumsum(theta * delta_v)
  p_log_lik <- sum(alln_u * log(theta)) 
  - sum(alld_u * cumhaz) 
  - lambda * t(theta) %*% R %*% theta
  np_log_lik <- -p_log_lik
  return(np_log_lik) 
}

w <- quantile(y, probs = seq(0, 1, length.out = (m+1)))
delta_v <- w[2:(m+1)] - w[1:m]

opts <- list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-6, "maxeval" = 10000 )



smtheta_1e6 <- nloptr(x0 = rep(1, m), eval_f = eval_f0,
                      lb = (rep(0, m)), ub = rep(Inf, m), 
                      opts = opts, y =y, w = w, delta_v = delta_v, 
                      R = R, lambda = 1e-6)$solution


smtheta_1 <- nloptr(x0 = rep(1, m), eval_f = eval_f0,
                    lb = (rep(0, m)), ub = rep(Inf, m), 
                    opts = opts, y =y, w = w, delta_v = delta_v, 
                    R = R, lambda = 1)$solution

smtheta_10 <- nloptr(x0 = rep(1, m), eval_f = eval_f0,
                     lb = (rep(0, m)), ub = rep(Inf, m), 
                     opts = opts, y =y, w = w, delta_v = delta_v, 
                     R = R, lambda = 10)$solution

smtheta_50 <- nloptr(x0 = rep(1, m), eval_f = eval_f0,
                     lb = (rep(0, m)), ub = rep(Inf, m), 
                     opts = opts, y =y, w = w, delta_v = delta_v, 
                     R = R, lambda = 50)$solution


par(mfrow = c(2,2))
plot(c(smtheta_1e6, smtheta_1e6[m]) ~ w,
     type = "s", col = "black",ylab = "hazard est", 
     xlab = "Time", main = "Theta, lambda = 1e-6")
plot(c(smtheta_1, smtheta_1[m]) ~ w,
     type = "s", col = "black",ylab = "hazard est", 
     xlab = "Time", main = "Theta, lambda = 1")
plot(c(smtheta_10, smtheta_10[m]) ~ w,
     type = "s", col = "black",ylab = "hazard est", 
     xlab = "Time", main = "Theta, lambda = 10")
plot(c(smtheta_50, smtheta_50[m]) ~ w,
     type = "s", col = "black",ylab = "hazard est", 
     xlab = "Time", main = "Theta, lambda = 50")
