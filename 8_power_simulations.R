library(dglm)

# What is power of SLM and DGLM given true sample size, if parameter estimates were true?
set.seed(27599)

n <- 244
num_sims <- 1e4
sex <- rbinom(n = n, size = 1, prob = 0.5)
g <- rbinom(n = n, size = 2, prob = 0.5)

mean_mu <- 28.1915124
mean_sex <- -9.6256202
mean_het <- -0.6585349
mean_hom <- 3.8556717
var_mu <- 4.5160431
var_sex <- -0.2212854
var_het <- -0.1748949
var_hom <- -1.3465278

mean_lp <- mean_mu + sex*mean_sex + mean_het*(g == 1) + mean_hom*(g == 2)
var_lp <- var_mu + sex*var_sex + var_het*(g == 1) + var_hom*(g == 2)
g <- factor(g)

lm_LRT <- dglm_LRT <- rep(NA, num_sims)
for (sim_idx in 1:num_sims) {
  y <- rnorm(n = n,
             mean = mean_lp,
             sd = exp(0.5*var_lp))
  lm_nul <- lm(formula = y ~ sex)
  lm_alt <- lm(formula = y ~ sex + g)
  dglm_nul <- dglm(formula = y ~ sex, dformula = ~ sex + g)
  dglm_alt <- dglm(formula = y ~ sex + g, dformula = ~ sex + g)
  lm_LRT[sim_idx] <- 2*(logLik(object = lm_alt) - logLik(object = lm_nul))
  dglm_LRT[sim_idx] <- dglm_nul$m2loglik - dglm_alt$m2loglik
}

par(mfrow = c(3, 1))
plot(x = sort(lm_LRT), y = sort(dglm_LRT)); abline(0, 1)
hist(x = lm_ps <- pchisq(q = lm_LRT, df = 2, lower.tail = FALSE),
     xlim = c(0, 1), main = 'SLM')
hist(x = dglm_ps <- pchisq(q = dglm_LRT, df = 2, lower.tail = FALSE),
     xlim = c(0, 1), main = 'DGLM')

mean(lm_ps < 0.05/100)
mean(dglm_ps < 0.05/100)


# how many more mice would be required to get ~90% power from SLM?
set.seed(27599)

n <- 350
num_sims <- 1e4
sex <- rbinom(n = n, size = 1, prob = 0.5)
g <- rbinom(n = n, size = 2, prob = 0.5)

mean_mu <- 28.1915124
mean_sex <- -9.6256202
mean_het <- -0.6585349
mean_hom <- 3.8556717
var_mu <- 4.5160431
var_sex <- -0.2212854
var_het <- -0.1748949
var_hom <- -1.3465278

mean_lp <- mean_mu + sex*mean_sex + mean_het*(g == 1) + mean_hom*(g == 2)
var_lp <- var_mu + sex*var_sex + var_het*(g == 1) + var_hom*(g == 2)
g <- factor(g)

lm_LRT <- dglm_LRT <- rep(NA, num_sims)
for (sim_idx in 1:num_sims) {
  y <- rnorm(n = n,
             mean = mean_lp,
             sd = exp(0.5*var_lp))
  lm_nul <- lm(formula = y ~ sex)
  lm_alt <- lm(formula = y ~ sex + g)
  lm_LRT[sim_idx] <- 2*(logLik(object = lm_alt) - logLik(object = lm_nul))
}

lm_ps <- pchisq(q = lm_LRT, df = 2, lower.tail = FALSE)
mean(lm_ps < 0.05/100)
