LP_panel = function(data,Y.name,X.name,
                              c.name = NULL,
                              id.name = NULL,
                              time.name = NULL,
                              H = 10,
                              h.seq = 0:H,
                              lagY = 1, 
                              lagX = NULL,
                              method = "SPJ",
                              diff = 0) {
  #method: "FE","SPJ"
  
  px = length(X.name)
  IRF=matrix(NA,length(h.seq),px)
  se=matrix(NA,length(h.seq),px)
  weakiv = matrix(NA,length(h.seq),px)
  
  if (is.null(id.name)){
    id.name <- names(data)[1]
  }
  if (is.null(time.name)){
    time.name <- names(data)[2]
  }
  if (is.null(lagX)){
    lagX = max(lagY,1)
  }
  lag.max <- max(lagX,lagY)
  id.temp <- data[,id.name]
  N = length(unique(id.temp))
  for (nn in 1:N){
    data[,id.name][id.temp == unique(id.temp)[nn]] <- nn 
  }
  time.temp <- data[,time.name]
  T0 <- length( sort(unique(time.temp)) )
  for (tt in 1:T0){
    data[,time.name][time.temp == sort(unique(time.temp))[tt]] <- tt
  }
  
  id <- data[,id.name]
  time <- data[,time.name]
  
  data.new <- NULL 
  for (nn in 1:N){
    data_nn <- data[data[,id.name] == nn,]
    Missing <- !(1:T0 %in% data_nn[,time.name])
    if (sum(Missing) > 0){
      Missing.ind <- (1:T0)[Missing]
      add <- matrix(NA,nrow=length(Missing.ind),ncol = ncol(data) )
      add[,which(names(data) == id.name)] <- rep(nn,length(Missing.ind))
      add[,which(names(data) == time.name)] <- Missing.ind
      colnames(add) <- colnames(data_nn)
      data_nn <- rbind(data_nn,add)
      order <- sort(data_nn[,which(names(data) == time.name)],index.return = T)$ix 
      data_nn <- data_nn[order,] 
    }
    
    data.new <- rbind(data.new, data_nn)
  }
  
  data <- data.new 
  
  id <- data[,id.name]
  time <- data[,time.name]
  
  y.long <- as.matrix(data[,Y.name])
  X.long <- as.matrix(data[,X.name]) 
  
  
  if (is.null(c.name)){
    cc.ind <- which(!(colnames(data) %in% c(Y.name,X.name,id.name,time.name)))
    c.name <- names(data)[cc.ind]
  }
  pc = length(c.name)
  cc.long <- as.matrix(data[ , c.name])
  
  y <- matrix(y.long, ncol = N)
  
  x <- NULL
 
  for (ix in 1:px){
    x.ix.long <- X.long[,ix]
    x.ix <- matrix(x.ix.long, ncol = N)
    x = cbind(x,x.ix)
  }
  
  cc <- NULL
  if (pc > 0){
    for (ic in 1:pc){
      cc.ic.long <- cc.long[,ic]
      cc.ic <- matrix(cc.ic.long, ncol = N)
      cc = cbind(cc,cc.ic)
    }
  }
  
  
  
  for(ih in 1:length(h.seq)){
    
    h = h.seq[ih]
    
    if (lagY == 0){
      hh=h
    }else if (lagY > 0){
      hh=max(h,1)
    }
    
    lag_T=T0-hh
    
    
    if (!diff){
      y_h <- y[(hh+lag.max):(T0),]
    }else{
      if (hh+lag.max > 1){
        y_h = y[(hh+lag.max):(T0),] - y[ ((hh+lag.max):(T0) - max(hh,1)),]
      }else{
        y_h = y[(hh+lag.max):(T0),] - rbind( matrix(NA,max(hh,1),ncol(y)),y[1:(T0-max(hh,1)),])
      }
    
    }
    
 
  
    
    if (lagY >= 1){
      for (ilag in 0:(lagY-1)){
        
        if (!diff){
          
          if (h == 0){
            y.temp <- y[(hh+lag.max-ilag-1):(T0-ilag-1),]
          }else{
            y.temp <- y[(lag.max-ilag):(lag_T-ilag),]
          }
          
        }else{
          
          if (h == 0){
            
            if (hh+lag.max-ilag-1 > 1){
              y.temp <- y[(hh+lag.max-ilag-1):(T0-ilag-1),] - y[(hh+lag.max-ilag-1):(T0-ilag-1)-1,]
            }else{
              y.temp <- y[1:(T0-ilag-1),] - rbind(NA,y[1:(T0-ilag-2),])
            }
            
          }else{
            
            if (lag.max-ilag > 1){
              y.temp <- y[(lag.max-ilag):(lag_T-ilag),] - y[(lag.max-ilag):(lag_T-ilag)-1,]
            }else{
              y.temp <- y[1:(T0-ilag-1),] - rbind(NA,y[1:(T0-ilag-2),])
            } 
            
          }
        }
        
       
        
        y_h <- cbind(y_h,y.temp)
         
      }
    }
    
    
    x_h <- NULL 
    
    
    for (ipx in 1:px){
      x_h_ipx <- NULL
       
      if (h == 0){
        x_h_ipx <- x[(hh+lag.max):(T0),(ipx-1)*N+1:N]
      }else{
          for (ilag in 0:(lagX-1)){
            x.temp <- x[(lag.max-ilag):(lag_T-ilag),(ipx-1)*N+1:N]
            x_h_ipx <- cbind(x_h_ipx,x.temp)
          }
        }
      
      x_h <- cbind(x_h,x_h_ipx)
      
    }
    
    if (pc > 0){
      cc_h=cc[(hh+lag.max):T0,]
    }
    
    if(method == "FE") {
      
      yt <- c(y_h[,1:N])
      
      X <- matrix(x_h, nrow = length(yt) )
      
      CC = NULL 
      if (pc >0){
        CC <- matrix(cc_h, nrow = length(yt) )
      }
      
      ylag <- NULL
      if (lagY >= 1){
        ylag <- matrix(y_h[,-(1:N)], nrow = length(yt) )
      }
      
      yt.dm <- yt
      X.dm <- X
      CC.dm <- CC
      ylag.dm <- ylag 
      
      
      for (iN in 1:N){
        want <- 1:nrow(y_h) + (iN-1)*nrow(y_h)
        want.mean <- rowSums(is.na(cbind(yt[want],X[want,],CC[want,],ylag[want,]))) == 0
        yt.dm[want] <-  yt.dm[want] - mean(yt[want[want.mean]])
        for (ix in 1:ncol(X)){
          X.dm[want,ix] <-  X.dm[want,ix] - mean(X[want[want.mean],ix])
        }
        
        if (pc >0){
          for (ic in 1:ncol(CC)){
            CC.dm[want,ic] <-  CC.dm[want,ic] -  mean(CC[want[want.mean],ic])
          }
        }
        if (lagY >= 1){
          for (iylag in 1:ncol(ylag)){
          ylag.dm[want,iylag] <-  ylag.dm[want,iylag] - mean(ylag[want[want.mean],iylag]) 
          }
        }
      }
      
      indep_var <- cbind(X.dm,CC.dm,ylag.dm)
      fit_h=lm(yt.dm ~ indep_var + 0)
      
      if (h == 0){
        IRF[ih,] = fit_h$coefficients[1:px]
      }else{
        IRF[ih,] = fit_h$coefficients[(1:px-1)*lagX + 1]
      }
      
      res.vec = yt.dm - indep_var %*% fit_h$coefficients
      res <- matrix(res.vec,ncol = N)
      
      W = matrix(0,ncol(indep_var),ncol(indep_var))
      
      for (iN in 1:N){ 
        indep_iN <- as.matrix(indep_var[(iN-1)*nrow(y_h)+(1:nrow(y_h)),])
        want_iN <- !( rowSums(is.na(indep_iN)) > 0 | is.na(res[,iN]) )
        
        
        if (sum(want_iN) == 1){
          temp <- as.matrix(indep_iN[want_iN,])%*%res[want_iN,iN]
          W = W + temp %*% t(temp)
        }else{
          W = W + t(indep_iN[want_iN,])%*%res[want_iN,iN]%*%t(res[want_iN,iN])%*%indep_iN[want_iN,]
        }
       
      }
      want_NT <- !( rowSums(is.na(indep_var)) > 0 | is.na(res.vec) )
      smp <- length(res.vec[want_NT,])
      temp <- t(indep_var[want_NT,])%*%indep_var[want_NT,]
      var_mat <- solve(temp)%*%W%*%solve(temp) *(N/(N-1)*(smp-1)/(smp-ncol(indep_var)))
      if (h == 0){
        se[ih,] = sqrt(diag( as.matrix(var_mat) ))[1:px]
      }else{
        se[ih,] = sqrt(diag( as.matrix(var_mat) ))[(1:px-1)*lagX + 1]
      }
      
    } else if (method == "SPJ") {
      
      yt <- c(y_h[,1:N])
      
      X <- matrix(x_h, nrow = length(yt) )
      
      CC = NULL 
      if (pc >0){
        CC <- matrix(cc_h, nrow = length(yt) )
      }
      
      
      ylag <- NULL
      if (lagY >= 1){
        ylag <- matrix(y_h[,-(1:N)], nrow = length(yt) )
      }
      
      yt.dm <- yt
      X.dm <- X
      CC.dm <- CC
      ylag.dm <- ylag 
      
      for (iN in 1:N){
        want <- 1:nrow(y_h) + (iN-1)*nrow(y_h)
        want.mean <- rowSums(is.na(cbind(yt[want],X[want,],CC[want,],ylag[want,]))) == 0
        yt.dm[want] <-  yt.dm[want] - mean(yt[want[want.mean]])
        for (ix in 1:ncol(X)){
          X.dm[want,ix] <-  X.dm[want,ix] - mean(X[want[want.mean],ix])
        }
        
        if (pc >0){
          for (ic in 1:ncol(CC)){
            CC.dm[want,ic] <-  CC.dm[want,ic] -  mean(CC[want[want.mean],ic])
          }
        }
        if (lagY >= 1){
          for (iylag in 1:ncol(ylag)){
            ylag.dm[want,iylag] <-  ylag.dm[want,iylag] - mean(ylag[want[want.mean],iylag]) 
          }
        }
      }
      indep_var <- cbind(X.dm,CC.dm,ylag.dm)
      
      fit_h=lm(yt.dm ~ indep_var + 0)
      beta_all=fit_h$coefficients 
      
      ######################
      T0h <- nrow(y_h)

      smp=cbind(yt.dm,indep_var)
      smpna = complete.cases(smp)
      smpna.ti <- matrix(smpna, ncol = N)
      cut = matrix(NA,1,ncol(smpna.ti))
      for (iN in 1:ncol(cut)) {
        cut[1,iN] = floor(median(which(smpna.ti[,iN]==TRUE)))
      }
      
      ################## Sample a 
      y_h_a = matrix(NA,nrow(y_h),ncol(y_h))
      for (ixx in 1:ncol(y_h_a)) {
        if (ixx%%N !=0) {
          if (is.na(cut[1,ixx%%N])==0) {
            y_h_a[(1:cut[1,ixx%%N]),ixx]=y_h[(1:cut[1,ixx%%N]),ixx]
          }
        } else if (ixx%%N ==0) {
          if (is.na(cut[1,N])==0) {
            y_h_a[(1:cut[1,N]),ixx]=y_h[(1:cut[1,N]),ixx]
          }
        }
      }
      yt_a <- c(y_h_a[,1:N])
      
      x_h_a = matrix(NA,nrow(x_h),ncol(x_h))
      for (ixx in 1:ncol(x_h_a)) {
        if (ixx%%N !=0) {
          if (is.na(cut[1,ixx%%N])==0) {
            x_h_a[(1:cut[1,ixx%%N]),ixx]=x_h[(1:cut[1,ixx%%N]),ixx]
          }
        } else if (ixx%%N ==0) {
          if (is.na(cut[1,N])==0) {
            x_h_a[(1:cut[1,N]),ixx]=x_h[(1:cut[1,N]),ixx]
          }
        }
      }
      X_a <- matrix(x_h_a, nrow = length(yt_a))
      
      CC_a = NULL
      if (pc >0){
        cc_h_a = matrix(NA,nrow(cc_h),ncol(cc_h))
        for (ixx in 1:ncol(cc_h_a)) {
          if (ixx%%N !=0) {
            if (is.na(cut[1,ixx%%N])==0) {
              cc_h_a[(1:cut[1,ixx%%N]),ixx]=cc_h[(1:cut[1,ixx%%N]),ixx]
            }
          } else if (ixx%%N ==0) {
            if (is.na(cut[1,N])==0) {
              cc_h_a[(1:cut[1,N]),ixx]=cc_h[(1:cut[1,N]),ixx]
            }
          }
        }
        CC_a <- matrix(cc_h_a, nrow = length(yt_a) )
      }
      
      ylag_a <- NULL
      if (lagY >= 1){
        ylag_a <- matrix(y_h_a[,-(1:N)], nrow = length(yt_a) )
      }
      
      yt.dm_a <- yt_a
      X.dm_a <- X_a
      CC.dm_a <- CC_a
      ylag.dm_a <- ylag_a
     
      for (iN in 1:N){
        want <- 1:nrow(y_h_a) + (iN-1)*nrow(y_h_a)
        want.mean <- rowSums(is.na(cbind(yt_a[want],X_a[want,],CC_a[want,],ylag_a[want,]))) == 0
        yt.dm_a[want] <-  yt.dm_a[want] - mean(yt_a[want[want.mean]])
        for (ix in 1:ncol(X_a)){
          X.dm_a[want,ix] <-  X.dm_a[want,ix] - mean(X_a[want[want.mean],ix])
        }
        
        if (pc >0){
          for (ic in 1:ncol(CC_a)){
            CC.dm_a[want,ic] <-  CC.dm_a[want,ic] -  mean(CC_a[want[want.mean],ic])
          }
        }
        if (lagY >= 1){
          for (iylag in 1:ncol(ylag_a)){
            ylag.dm_a[want,iylag] <-  ylag.dm_a[want,iylag] - mean(ylag_a[want[want.mean],iylag]) 
          }
        }
      }
      indep_var_a <- cbind(X.dm_a,CC.dm_a,ylag.dm_a)
      
      fit_h_a=lm(yt.dm_a ~ indep_var_a + 0) 
      beta_a = fit_h_a$coefficients 
       
      ################## Sample b
      
      y_h_b = matrix(NA,nrow(y_h),ncol(y_h))
      for (ixx in 1:ncol(y_h_b)) {
        if (ixx%%N !=0) {
          if (is.na(cut[1,ixx%%N])==0) {
            y_h_b[-(1:cut[1,ixx%%N]),ixx]=y_h[-(1:cut[1,ixx%%N]),ixx]
          }
        } else if (ixx%%N ==0) {
          if (is.na(cut[1,N])==0) {
            y_h_b[-(1:cut[1,N]),ixx]=y_h[-(1:cut[1,N]),ixx]
          }
        }
      }
      yt_b <- c(y_h_b[,1:N])
      
      x_h_b = matrix(NA,nrow(x_h),ncol(x_h))
      for (ixx in 1:ncol(x_h_b)) {
        if (ixx%%N !=0) {
          if (is.na(cut[1,ixx%%N])==0) {
            x_h_b[-(1:cut[1,ixx%%N]),ixx]=x_h[-(1:cut[1,ixx%%N]),ixx]
          }
        } else if (ixx%%N ==0) {
          if (is.na(cut[1,N])==0) {
            x_h_b[-(1:cut[1,N]),ixx]=x_h[-(1:cut[1,N]),ixx]
          }
        }
      }
      X_b <- matrix(x_h_b, nrow = length(yt_b))
      
      CC_b = NULL
      if (pc >0){
        cc_h_b = matrix(NA,nrow(cc_h),ncol(cc_h))
        for (ixx in 1:ncol(cc_h_b)) {
          if (ixx%%N !=0) {
            if (is.na(cut[1,ixx%%N])==0) {
              cc_h_b[-(1:cut[1,ixx%%N]),ixx]=cc_h[-(1:cut[1,ixx%%N]),ixx]
            }
          } else if (ixx%%N ==0) {
            if (is.na(cut[1,N])==0) {
              cc_h_b[-(1:cut[1,N]),ixx]=cc_h[-(1:cut[1,N]),ixx]
            }
          }
        }
        CC_b <- matrix(cc_h_b, nrow = length(yt_b) )
      }
      
      ylag_b <- NULL
      if (lagY >= 1){
        ylag_b <- matrix(y_h_b[,-(1:N)], nrow = length(yt_b) )
      }
      
      yt.dm_b <- yt_b
      X.dm_b <- X_b
      CC.dm_b <- CC_b
      ylag.dm_b <- ylag_b
      
      for (iN in 1:N){
        want <- 1:nrow(y_h_b) + (iN-1)*nrow(y_h_b)
        want.mean <- rowSums(is.na(cbind(yt_b[want],X_b[want,],CC_b[want,],ylag_b[want,]))) == 0
        yt.dm_b[want] <-  yt.dm_b[want] - mean(yt_b[want[want.mean]])
        for (ix in 1:ncol(X_b)){
          X.dm_b[want,ix] <-  X.dm_b[want,ix] - mean(X_b[want[want.mean],ix])
        }
        
        if (pc >0){
          for (ic in 1:ncol(CC_b)){
            CC.dm_b[want,ic] <-  CC.dm_b[want,ic] -  mean(CC_b[want[want.mean],ic])
          }
        }
        if (lagY >= 1){
          for (iylag in 1:ncol(ylag_b)){
            ylag.dm_b[want,iylag] <-  ylag.dm_b[want,iylag] - mean(ylag_b[want[want.mean],iylag]) 
          }
        }
      }
      indep_var_b <- cbind(X.dm_b,CC.dm_b,ylag.dm_b)
      
      fit_h_b=lm(yt.dm_b ~ indep_var_b + 0) 
      beta_b = fit_h_b$coefficients 
      
      ######### estimate IRF 
      beta.hat = 2*beta_all-0.5*(beta_a+beta_b)
      
      if (h == 0){
        IRF[ih,] = beta.hat[1:px]
      }else{
        IRF[ih,] = beta.hat[(1:px-1)*lagX + 1]
      }
      ################### s.e. 
      
      res.vec <- yt.dm - indep_var %*% beta.hat
      res <- matrix(res.vec, ncol = N) 
      
      XX.mat <- matrix(indep_var, nrow = T0h)
      XX.mat.a <- matrix(indep_var_a, nrow = T0h)
      XX.mat.b <- matrix(indep_var_b, nrow = T0h)
      
      XX.mat.sub = matrix(NA,nrow(XX.mat),ncol(XX.mat))
      for (ixx in 1:ncol(XX.mat.sub)) {
        if (ixx%%N !=0) {
          if (is.na(cut[1,ixx%%N])==0) {
            XX.mat.sub[(1:cut[1,ixx%%N]),ixx]=XX.mat.a[(1:cut[1,ixx%%N]),ixx]
            XX.mat.sub[-(1:cut[1,ixx%%N]),ixx]=XX.mat.b[-(1:cut[1,ixx%%N]),ixx]
          }
        } else if (ixx%%N ==0) {
          if (is.na(cut[1,N])==0) {
            XX.mat.sub[(1:cut[1,N]),ixx]=XX.mat.a[(1:cut[1,N]),ixx]
            XX.mat.sub[-(1:cut[1,N]),ixx]=XX.mat.b[-(1:cut[1,N]),ixx]
          }
        }
      }
      
      dd.mat <- 2*XX.mat - XX.mat.sub
      dd.mat <- matrix(dd.mat, ncol = ncol(indep_var) )
      
      R.hat = matrix(0,ncol(indep_var),ncol(indep_var))
      for (iN in 1:N){ 
        dd_iN <- as.matrix(dd.mat[(iN-1)*nrow(y_h)+(1:nrow(y_h)),])
        want_iN <- !( rowSums(is.na(dd_iN)) > 0 | is.na(res[,iN]) )
        if (sum(want_iN) == 1){
          temp <- as.matrix(dd_iN[want_iN,])%*%res[want_iN,iN]
          R.hat = R.hat + temp %*% t(temp)
        }else{
        R.hat = R.hat+t(dd_iN[want_iN,])%*%res[want_iN,iN]%*%t(res[want_iN,iN])%*%dd_iN[want_iN,]
        }
      }
      
      want_NT <- !( rowSums(is.na(indep_var)) > 0 | is.na(res.vec) )
      smp <- length(res.vec[want_NT,])
      R.hat <- R.hat / smp
      Q.hat = t(indep_var[want_NT,]) %*% indep_var[want_NT,] / smp
      var.hat <- solve(Q.hat) %*% R.hat %*% solve(Q.hat) *(N/(N-1)*(smp-1)/(smp-ncol(indep_var)))
      if (h == 0){
        se[ih,] = sqrt( diag(var.hat / smp ))[1:px]
      }else{
        se[ih,] = sqrt( diag(var.hat / smp ))[(1:px-1)*lagX + 1]
      }
      
    } else{
      cat("pls speficy the method")
    }
  }
  
  results=list(IRF,se,X.name,h.seq,weakiv)
  names(results)=c("IRF","se","X.name","h.seq","weakiv")
  return(results)
}