require(rjags)
library(markovchain)
library(sp)
library(lattice)

load("data/Weddell_data.RData")
#unique(datafile1$argos$id)
#[1] wd04-880-11 wd04-896-11 wd04-897-11 wd04-899-11 wd04-911-11
#unique(datafile2$argos$id)
#[1] wd04-882-11 wd04-883-11 wd04-898-11 wd04-907-11

tod <- TRUE
n.chains <- 2
n.iter <- 1000 #00
n.thin <- 10 #20
n.burn <- 2000 #00
seed = sample(1:1e+05, 1)
tstep_hr <- 6
tstep <- tstep_hr/24

source("fn/check.haulouts.tstep.r") # function deals with the HOs
source("fn/dat4jags_reg.r") # adapted to regular start times
source("fn/plot.output.r")
source("fn/ssm3.r")

dat <- datafile2$argos
xyplot(lat~lon|id,data=dat,xlim=c(73,87),ylim=c(-70,-66),type="b",col="black")

use <- dat4jags(dat,tstep=tstep,tod=tod)
ho.check <- check.haulouts.tstep(use,hauldat)
ho.check$HO_cnt

# go directly to case study animal
id <- 2
dat1 <- ho.check[[1]][[id]]
dat1$id <- factor(dat1$id,levels=unique(dat1$id))
dat1 <- list(dat1)
names(dat1[[1]])[11] <- "ho"
xyplot(lat~lon|id,data=dat,xlim=c(73,87),ylim=c(-70,-66),
  type="b",col="black",subset=as.character(dat$id)==unique(dat1[[1]]$id))

st = proc.time()
fit <- ssm(loc.list=dat1, model="haulout3", LM=FALSE, HO=TRUE,
            adapt=n.burn, samples=n.iter, thin=n.thin, chains=n.chains)
cat("Elapsed time: ", round((proc.time() - st)[3]/60,2), "min \n")

save(fit,dat1,file=paste("output/",dat1[[1]]$id[1],"haulout4.RData",sep=""))

plot.output(fit,dat1,map_range=c(73,87,-70,-66),haulout4=TRUE)

table(dat1[[1]]$ho)
length(which(fit[[1]]$summary$b>2))
