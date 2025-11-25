

simul_LP <- function(N,T0,rho,tau = 0, tau2 = 0, TE = F, diff = F,
                     Niter = 1000,H = 10,
                     beta = -0.6, 
                     delta = 0.2, eta = 0.2,
                     sigmaep = 1,sigmamu = 1, sigmachi = 1,
                     sig.level = 0.05,  
                     lagX = 1, lagY = 1, seed = 2023){
  set.seed(seed)
  
  if (lagY == 0){
    IRFs <- beta * rho^(0:H)
  }else{
    A.mat <- matrix(c(tau + beta*tau2,tau2,beta*rho,rho),2,2)
    IRFs = rep(NA,H+1)
    Ah = A.mat
    IRFs[1] = beta
    IRFs[2] = beta*rho
    for (hh in 2:H){
      Ah <- A.mat %*% Ah
      IRFs[hh+1] <- Ah[1,2]
    }
  }
  
  if (diff){
    IRFs = c(beta,cumsum(beta * rho^(1:H)) )
  }
  
  cover.FE <- matrix(NA,H+1,Niter)
  cover.jackknife <- matrix(NA,H+1,Niter)
  
  SE.FE <- matrix(NA,H+1,Niter)
  SE.jackknife <- matrix(NA,H+1,Niter)
  
  irf.FE <- matrix(NA,H+1,Niter)
  irf.jackknife <- matrix(NA,H+1,Niter)
  
  if (tau == 0){
    len.DB <- matrix(NA,H+1,Niter)
    cover.DB <- matrix(NA,H+1,Niter)
    SE.DB <- matrix(NA,H+1,Niter)
    irf.DB <- matrix(NA,H+1,Niter)
  }
  
  
  for (iter in 1:Niter){
    
    print(c(N,T0,rho,tau,iter))
    #Step 2: Data initialization & define random noise
    Y=matrix(0,T0,N)  
    X=matrix(0,T0,N)
    alpha=rep(0,N)
    Xbar=rep(0,1)
    epsilon=matrix(rnorm(N*T0,mean=0,sd=sigmaep),nrow=T0)
    mu=matrix(rnorm(N*T0,mean=0,sd=sigmamu),nrow=T0)
    chi=rnorm(n=1,mean=0,sd=sigmachi)
    X0=rnorm(N,mean=0,sd=1)
    Xm1=rnorm(N,mean=0,sd=1)
    # Step 3: Generate data using iteration
    
    for (i in 1:N) {
      
      ### AR(1) 
      if (lagY == 0){
        for(t in 1:T0){
          if(t==1){
            X[t,i]=delta+rho* X0[i]+mu[t,i]
          } else {
            X[t,i]=delta+rho* X[(t-1),i]+ mu[t,i] 
          }
        }
      }
      
      
      Xbar[i]=mean(X[,i])
      alpha[i]=eta*sqrt(T0)*Xbar[i]+chi #fixed effect
      
      for (t in 1:T0){
        
        if (t == 1){ 
          X[t,i]=delta+rho* X0[i]+mu[t,i]
        }
        
        te = TE*(0.025*t + 0.001 * t^2)
        
        if (lagY == 1 & t >= 2){
          X[t,i]=delta + tau2 * Y[t-1,i]+ rho* X[(t-1),i]+ mu[t,i] 
          Y[t,i]=alpha[i] + tau*Y[t-1,i]+beta*X[t,i]+epsilon[t,i] + te 
        }else{
          Y[t,i]= alpha[i]+ beta*X[t,i]+epsilon[t,i] + te
        }
        
      }
      
      if (diff){
        Y[,i] = cumsum(Y[,i])
      }
      
    }
    
    YX <- cbind(c(Y), c(X) )
    data <- as.data.frame(cbind(rep(1:N,each = T0),rep(1:T0,N),YX))
    
    colnames(data) <- c("id","time","Y","X")
    X.name <- c("X")
    fit.FE <- LP_panel(data,Y.name = "Y",X.name = X.name,
                       lagX = lagX, lagY = lagY, method = "FE", H = H, diff = diff)
    IRF.FE <-  fit.FE$IRF
    se.FE <-  fit.FE$se
    SE.FE[,iter] <- (IRF.FE - IRFs)^2
    cover.FE[,iter] <- (IRFs <= IRF.FE + qnorm(1-sig.level/2) *se.FE) &
      (IRFs >= IRF.FE - qnorm(1-sig.level/2) *se.FE)
    irf.FE[,iter] <- fit.FE$IRF
    
    if (tau == 0){
      x.f <- as.vector(X[2:(T0),])
      x.l <- as.vector(X[1:(T0-1),])
      AR.fit <- summary(lm(x.f~x.l))
      rho.hat <- AR.fit$coefficients[2]
      s.u <- sqrt(mean(AR.fit$residuals^2))
      
      y.vec=as.vector(demean(Y))
      x.vec=as.vector(demean(X))
      all.fit <- summary(lm(y.vec~x.vec-1))
      s.eps <- sqrt(mean(all.fit$residuals^2)) 
      b0.hat <- all.fit$coefficients[1]
      
      f.rho <- s.u^2*( (1-rho.hat^(0:H)) - (0:H/(T0-0:H)))/(1-rho.hat)^2  
      
      sx = sqrt(mean(x.vec^2))
      
      IRF.DB = IRF.FE + b0.hat * f.rho / ( (T0-0:H) * sx^2)
      SE.DB[,iter] <- (IRF.DB - IRFs)^2
      cover.DB[,iter] <- (IRFs <= IRF.DB + qnorm(1-sig.level/2) *se.FE) &
        (IRFs >= IRF.DB - qnorm(1-sig.level/2) *se.FE)
      irf.DB[,iter] <- IRF.FE + b0.hat * f.rho / ( (T0-0:H) * sx^2)
    }

    fit.jackknife <- LP_panel(data,Y.name = "Y",X.name = X.name,
                              lagX = lagX, lagY = lagY, method = "SPJ", H = H, diff = diff)
    IRF.jackknife <-  fit.jackknife$IRF
    se.jackknife <-  fit.jackknife$se
    SE.jackknife[,iter] <- (IRF.jackknife - IRFs)^2
    cover.jackknife[,iter] <- (IRFs <= IRF.jackknife + qnorm(1-sig.level/2) *se.jackknife) &
      (IRFs >= IRF.jackknife - qnorm(1-sig.level/2) *se.jackknife)
    irf.jackknife[,iter] <- fit.jackknife$IRF
  }
  
  cover = list()
  cover$FE = cover.FE
  cover$jackknife = cover.jackknife
  
  SE = list()
  SE$FE = SE.FE
  SE$jackknife = SE.jackknife
  
  irf = list()
  irf$FE = irf.FE
  irf$jackknife = irf.jackknife
  
  if (tau == 0){
    SE$DB = SE.DB
    cover$DB = cover.DB
    irf$DB = irf.DB
  }
 
  return(list(cover = cover,
              SE =  SE,
              irf = irf))
 
}










 
  






