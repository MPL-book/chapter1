## Partial likelihood example with R

library(survival)
summary(leukemia)
leukemia.sva <- coxph(Surv(time,status) ~ x, data=leukemia)
summary(leukemia.sva)