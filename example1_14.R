## Semiparametric AFT model with R

# Generate right-censored survival times from AFT model
epsilon_i <- rnorm(500)
sigma <- 1
x1 <- rbinom(500, 1, 0.5)
x2 <- runif(500, 0, 1)
X <- cbind(x1, x2)
beta_true = c(1, -0.5)

log_t <- X %*% beta_true + sigma * epsilon_i
t <- exp(log_t)
c <- rexp(500, 0.2)
delta <- as.numeric(t < c)
y <- t
y[which(delta == 0)] <- c[which(delta == 0)]

df <- data.frame(y, delta)
df <- cbind(df, X)

# Fit this model using R package aftgee
library(aftgee)
aft_fit <- aftgee(Surv(y, delta) ~ -1 + x1 + x2, data = df, corstr = "independence")
summary(aft_fit)
