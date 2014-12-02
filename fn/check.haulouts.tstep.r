# IDs SSM location dates encompassed entirely within a HO period
# NB usually very low for SES (and often at very tail)
# check.haulouts: inputs are 'data' from dat4bugCOV and 'HO'=direct from SMRU table
check.haulouts.tstep <- function(data,HO)  {
HO_cnt <- array(NA,dim=c(length(data),3))
for (i in 1:length(data)) {
tstep <- data[[i]]$tstep
dates <- data[[i]]$first.date + 86400*tstep*c(0:c(data[[i]]$RegN-1))

data[[i]]$haulouts <- matrix(0,nrow=data[[i]]$RegN,ncol=1)
ho <- HO[HO$ref==names(data)[i] ,]

for (ii in 1:length(dates) ) {
 # NB. in dat4bugsCOV the covariate info is averaged over the interval BETWEEN
 # regularised times (ie. not straddling/centred on)
 t1 <- dates[ii] ; t2 <- dates[ii]+c(86400*tstep-1)
 # June 2014 this has been modified so that the HO status allocated to each estimated 
 # position now STRADDLES the position (i.e. interval of interest is x-0.5*tstep : x+0.5*tstep
 t1 <- t1 -0.5*86400*tstep;  t2 <- t2 -0.5*86400*tstep;
 h1 <- ho$S_DATE ; h2 <- ho$E_DATE
 # 1. FULL HO - ***record 1 (can record row number of the haulout)*** 
 ff <- which(t1 >= h1  &  t2 <= h2);
 if (length(ff)) {data[[i]]$haulouts[ii] <- 1 #ff  
 } else {
 # 2. PARTIAL HO - ***records proportion of timespent in HO status***
 # A. Haulout straddles (starts before/at and ends during) timestep
 ff <- which(h1 <= t1 & h2 > t1 & h2 < t2)
 if (length(ff)) {x1=h2[ff];x2=t1 ;    #if (length(ff)>1) browser()
    data[[i]]$haulouts[ii] <- data[[i]]$haulouts[ii]+
    as.numeric(difftime(x1,x2,units="secs"))/c(86400*tstep)  }
 # B. Haulout straddles (starts during and ends after) timestep
 ff <- which(h1 >= t1 & h1 < t2 & h2 >= t2)
 if (length(ff)) {x1=t2; x2=h1[ff];        #if (length(ff)>1) browser()
    data[[i]]$haulouts[ii] <- data[[i]]$haulouts[ii]+
    as.numeric(difftime(x1,x2,units="secs"))/c(86400*tstep)  }
 # C. Haulout(s) start and end during timestep (may be multiple)
 ff <- which(t1 <= h1 & t2 > h2)
 if (length(ff)) {x1=h2[ff]; x2=h1[ff]; #if (length(ff)>1) browser()
    data[[i]]$haulouts[ii] <- data[[i]]$haulouts[ii]+
    sum(as.numeric(difftime(x1,x2,units="secs"))/c(86400*tstep))  }
} 
} # end date loop
data[[i]]$haulouts <- as.vector(data[[i]]$haulouts)

#data[[i]]$HO_cnt <- c( length(which(data[[i]]$haulouts>=1)) , # full HO
#  length(which(data[[i]]$haulouts>0 & data[[i]]$haulouts<1)), # partial HO
#  length(data[[i]]$haulouts) )   # number timesteps
HO_cnt[i,] <-  c( length(which(data[[i]]$haulouts>=1)) , # full HO
  length(which(data[[i]]$haulouts>0 & data[[i]]$haulouts<1)), # partial HO
  length(data[[i]]$haulouts) )   # number timesteps
} # end tagID loop 
list(data=data, HO_cnt=HO_cnt)
} # end function