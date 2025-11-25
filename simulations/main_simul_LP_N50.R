rm(list=ls())
setwd("D:\\GitHub\\panel-lp-replication\\simulations")

PackageList =c('lpirfs', 'Jmisc', 'dplyr','plm','sandwich','parallel',"ivreg")
NewPackages=PackageList[!(PackageList %in%
                            installed.packages()[,"Package"])]
if(length(NewPackages)) install.packages(NewPackages)
lapply(PackageList,require,character.only=TRUE)

source("LP_panel_fe.R")
source("simul_LP.R")
Niter = 1000
N = 50
tau = 0
tau2 = 0
lagY = 0


for (rho in c(0,0.2,0.5,0.8)){
  for (T0 in c(120)){
    file.name <- paste0("T",T0,"N",N,"rho",rho*10,"tau",tau*10,"Niter",Niter,"lagY",lagY,".RDS")
    out <- simul_LP(N,T0,rho,tau,tau2,Niter = Niter, lagY = lagY)
    saveRDS(out, file = file.name)
  }
}

