params <-
list(prefix = "R_v_C", min.pomp.version = "2.0.3")


## ----precheck,include=FALSE----------------------------------------------
stopifnot(packageVersion("pomp") >= params$min.pomp.version)


## ----packages------------------------------------------------------------
library(pomp)
library(ggplot2)
library(magrittr)
library(microbenchmark)
library(tibble)
## ----seed,echo=FALSE-----------------------------------------------------
set.seed(56300069)

simulate(times=0:100,t0=0,
  params=c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=1),
  dmeasure=Csnippet("
    lik = dlnorm(Y,log(X),tau,give_log);"
  ),
  rmeasure=Csnippet("
    Y = rlnorm(log(X),tau);"
  ),
  rprocess=discrete_time(
    step.fun=Csnippet("
    double S = exp(-r*dt);
    double logeps = (sigma > 0.0) ? rnorm(0,sigma) : 0.0;
    X = pow(K,(1-S))*pow(X,S)*exp(logeps);"
    ),
    delta.t=1
  ),
  paramnames=c("r","K","sigma","tau"),
  obsnames="Y",
  statenames="X"
  ) -> Gompertz100

Gompertz100_out_50 <- pfilter(Gompertz100,Np=10000, time = 0:50)
Gompertz100_out_rest <- pfilter(Gompertz100_out_50, Np=10000, time = 51:100)

p2 <- Gompertz100@params

p2['sigma'] <- 0.3

pomp_r6 <- pomp_node$new(p2, 10000, 100, Gompertz100)

pomp_r6$run_pf(1:100)

pomp_r6_orig_params <- pomp_node$new(Gompertz100@params, 10000, 100, Gompertz100)

pomp_r6_orig_params$run_pf(1:100)
pfilter(Gompertz100, Np=10000)@loglik

simulate(times=0:10000,t0=0,
         params=c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=1),
         dmeasure=Csnippet("
    lik = dlnorm(Y,log(X),tau,give_log);"
    ),
    rmeasure=Csnippet("
    Y = rlnorm(log(X),tau);"
    ),
    rprocess=discrete_time(
      step.fun=Csnippet("
    double S = exp(-r*dt);
    double logeps = (sigma > 0.0) ? rnorm(0,sigma)  : 0.0;
    X = pow(K,(1-S))*pow(X,S)*exp(logeps);"
    ),
    delta.t=1
    ),
    paramnames=c("r","K","sigma","tau"),
    obsnames="Y",
    statenames="X"
    ) -> Gompertz10k

gompertz_bench_10k <- microbenchmark(pfilter(Gompertz10k,Np=10000), times = 10, unit = "s")
print(gompertz_bench_10k)
