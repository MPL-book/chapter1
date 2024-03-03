## KM method

library(survival)

# Find the KM estimate using the survfit() function
KM_est <- summary(survfit(Surv(time, status) ~1, data=lung))

# Plot the estimated survival function
plot(survfit(Surv(time, status) ~1, data=lung))

# Test for significant difference between two categories
survdiff(Surv(time, status) ~sex, data=lung)

# Derive and plot hazard function estimate from the KM estimate
diff <- -(log(KM_est$surv)[2:139] + log(KM_est$surv)[1:138])
t <- KM_est$time[1:138]
delta_t <- KM_est$time[2:139] - KM_est$time[1:138]
h_t <- diff/delta_t
plot(h_t ~ t, type = "l")
lo <- loess(h_t ~ t)
lines(lo$fitted ~ t, col='red', lwd=2)
