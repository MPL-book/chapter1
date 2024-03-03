## Nelson-Aalen method

# Compute Nelson-Aalen estimate of cumulative hazard via R
library(survival)
nelson_aalen = survfit(Surv(time, status) ~1, data=lung, ctype = 1)
nelson_aalen_cumhaz = nelson_aalen$cumhaz

# Find cumulative hazard using H_Y(t) = -logS_Y(t)
derived_cumhaz = -log(nelson_aalen$surv)

# Plot the above two functions and comment on similarities/differences
plot(nelson_aalen_cumhaz ~ nelson_aalen$time, type = "l")
lines(derived_cumhaz ~ nelson_aalen$time, type = "l", lty = 2)
