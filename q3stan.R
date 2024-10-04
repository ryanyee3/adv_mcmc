y1 <- read_csv("Y1.csv", col_names = F)
x1 <- read_csv("x1.csv", col_names = F)
beta_star1 <- read_csv("beta_star1.csv", col_names = F)


y10 <- read_csv("Y10.csv", col_names = F)
x10 <- read_csv("x10.csv", col_names = F)
beta_star10 <- read_csv("beta_star10.csv", col_names = F)


y100 <- read_csv("Y100.csv", col_names = F)
x100 <- read_csv("x100.csv", col_names = F)
beta_star100 <- read_csv("beta_star100.csv", col_names = F)



library(rstan)

model = stan_model("log_stan.stan")


data1 <- list(N = nrow(y1),
              P = ncol(x1),
              x = as.matrix(x1),
              y = array(as.integer(y1$X1)))

fit1 = sampling(model,data=data1,iter=2000,chains=4)
print(fit1)
stan_rhat(fit1)
stan_ess(fit1)
summary(fit1)$summary


data10 <- list(N = nrow(y10),
              P = ncol(x10),
              x = as.matrix(x10),
              y = array(as.integer(y10$X1)))

fit10 = sampling(model,data=data10,iter=2000,chains=4)
print(fit10)
stan_rhat(fit10)
stan_ess(fit10)
summary(fit10)$summary



data100 <- list(N = nrow(y100),
               P = ncol(x100),
               x = as.matrix(x100),
               y = array(as.integer(y100$X1)))

fit100 = sampling(model,data=data100,iter=2000,chains=4)
print(fit100)
stan_rhat(fit100)
stan_ess(fit100)
summary(fit100)$summary
