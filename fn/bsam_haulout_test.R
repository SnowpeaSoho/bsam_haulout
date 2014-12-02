require(rjags)

tod <- TRUE
n.chains <- 2
n.iter <- 10000 #00
n.thin <- 10 #20
n.burn <- 20000 #00
seed = sample(1:1e+05, 1)
tstep_hr <- 6
tstep <- tstep_hr/24

#source("/Users/sophie_bes/BESTLEY/scripts/bsam/R/dat4jags.R")
source("fn/ssm3.R")

fit3c <- ssm(loc.list=dat1, model="haulout3", LM=FALSE, HO=TRUE,
            adapt=n.burn, samples=n.iter, thin=n.thin, chains=n.chains)
