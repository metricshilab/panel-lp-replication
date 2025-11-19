LP_panel = function(data,Y.name,X.name,
                    c.name = NULL,
                    id.name = NULL,
                    time.name = NULL,
                    H = 10,
                    h.seq = 0:H,
                    lagY = 1, 
                    lagX = NULL,
                    method = "SPJ",
                    te = F,
                    two_way_cluster = F,
                    robust = F,
                    eigens = F) {
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
  
  #generate time dummy
  data[ , time.name] <- as.factor(data[ , time.name])
  td.long <- model.matrix(~data[ , time.name]-1,data)
  td <- NULL
  for (itd in 1:T0){
    td.itd.long <- td.long[,itd]
    td.itd <- matrix(td.itd.long, ncol = N)
    td = cbind(td,td.itd)
  }
  
  
  for(ih in 1:length(h.seq)){
    
    h = h.seq[ih]
    
    # When h=0, it is necessary to reduce one period of samples to ensure that lagy does not report an error
    if (lagY == 0){
      hh=h
    }else if (lagY > 0){
      hh=max(h,1)
    }
    
    lag_T=T0-hh
    
    y_h <- y[(hh+lag.max):(T0),]
    
    
    if (lagY >= 1){
      for (ilag in 0:(lagY-1)){
        
        if (h == 0){
          y.temp <- y[(hh+lag.max-ilag-1):(T0-ilag-1),]
        }else{
          y.temp <- y[(lag.max-ilag):(lag_T-ilag),]
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
    
    td_h=td[(hh+lag.max):T0,]
    
    
    if(method == "FE") {
      
      #demean_t
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
      
      if (te == F) {
        
        dep_var = yt.dm
        indep_var <- cbind(X.dm,CC.dm,ylag.dm)
        
      } else if (te == T) {
        
        #demean_n
        yt.dm.nt <- t(matrix(yt.dm, ncol = N))
        yt.dm.ntl <- c(yt.dm.nt[1:N,])
        
        X.dm.ntl = NULL 
        for (ix in 1:ncol(X.dm)){
          x.dm.ix <- X.dm[,ix]
          x.dm.ix.nt <- t(matrix(x.dm.ix, ncol = N))
          x.dm.ix.ntl <- matrix(x.dm.ix.nt, nrow = length(yt) )
          X.dm.ntl = cbind(X.dm.ntl,x.dm.ix.ntl)
        }
        
        CC.dm.ntl = NULL 
        if (pc >0){
          for (ic in 1:ncol(CC.dm)){
            cc.dm.ic <- CC.dm[,ic]
            cc.dm.ic.nt <- t(matrix(cc.dm.ic, ncol = N))
            cc.dm.ic.ntl <- matrix(cc.dm.ic.nt, nrow = length(yt) )
            CC.dm.ntl = cbind(CC.dm.ntl,cc.dm.ic.ntl)
          }
        }
        
        ylag.dm.ntl <- NULL
        if (lagY >= 1){
          ylag.dm.nt <- t(matrix(ylag.dm, ncol = N))
          ylag.dm.ntl <- matrix(ylag.dm.nt, nrow = length(yt) )
        }
        
        yt.dmdm.ntl <- yt.dm.ntl
        X.dmdm.ntl <- X.dm.ntl
        CC.dmdm.ntl <- CC.dm.ntl
        ylag.dmdm.ntl <- ylag.dm.ntl
        
        
        for (iT in 1:lag_T){
          want <- 1:N + (iT-1)*N
          want.mean <- rowSums(is.na(cbind(yt.dm.ntl[want],X.dm.ntl[want,],CC.dm.ntl[want,],ylag.dm.ntl[want,]))) == 0
          yt.dmdm.ntl[want] <-  yt.dmdm.ntl[want] - mean(yt.dm.ntl[want[want.mean]])
          for (ix in 1:ncol(X.dm.ntl)){
            X.dmdm.ntl[want,ix] <-  X.dmdm.ntl[want,ix] - mean(X.dm.ntl[want[want.mean],ix])
          }
          
          if (pc >0){
            for (ic in 1:ncol(CC.dm.ntl)){
              CC.dmdm.ntl[want,ic] <-  CC.dmdm.ntl[want,ic] -  mean(CC.dm.ntl[want[want.mean],ic])
            }
          }
          if (lagY >= 1){
            for (iylag in 1:ncol(ylag.dm.ntl)){
              ylag.dmdm.ntl[want,iylag] <-  ylag.dmdm.ntl[want,iylag] - mean(ylag.dm.ntl[want[want.mean],iylag]) 
            }
          }
        }
        
        #transfer to T*N
        yt.dmdm.tn <- t(matrix(yt.dmdm.ntl, nrow = N))
        yt.dmdm <- c(yt.dmdm.tn[,1:N])
        
        X.dmdm = NULL 
        for (ix in 1:ncol(X.dmdm.ntl)){
          x.dmdm.ix.ntl <- X.dmdm.ntl[,ix]
          x.dmdm.ix.tn <- t(matrix(x.dmdm.ix.ntl, nrow = N))
          x.dmdm.ix <- matrix(x.dmdm.ix.tn, nrow = length(yt) )
          X.dmdm = cbind(X.dmdm,x.dmdm.ix)
        }
        
        CC.dmdm = NULL 
        if (pc >0){
          for (ic in 1:ncol(CC.dmdm.ntl)){
            cc.dmdm.ic.ntl <- CC.dmdm.ntl[,ic]
            cc.dmdm.ic.tn <- t(matrix(cc.dmdm.ic.ntl, nrow = N))
            cc.dmdm.ic <- matrix(cc.dmdm.ic.tn, nrow = length(yt) )
            CC.dmdm = cbind(CC.dmdm,cc.dmdm.ic)
          }
        }
        
        ylag.dmdm <- NULL
        if (lagY >= 1){
          ylag.dmdm.tn <- t(matrix(ylag.dmdm.ntl, nrow = N))
          ylag.dmdm <- matrix(ylag.dmdm.tn, nrow = length(yt) )
        }
        
        dep_var = yt.dmdm
        indep_var <- cbind(X.dmdm,CC.dmdm,ylag.dmdm)
        
      }
      
      fit_h=lm(dep_var ~ indep_var + 0)
      res.vec = dep_var - indep_var %*% fit_h$coefficients

      if (h == 0){
        IRF[ih,] = fit_h$coefficients[1:px]
      }else{
        IRF[ih,] = fit_h$coefficients[(1:px-1)*lagX + 1]
      }

      W_N = matrix(0,ncol(indep_var),ncol(indep_var))
      res_N <- matrix(res.vec,ncol = N)
      for (iN in 1:N){ 
        indep_iN <- as.matrix(indep_var[(iN-1)*nrow(y_h)+(1:nrow(y_h)),])
        want_iN <- !( rowSums(is.na(indep_iN)) > 0 | is.na(res_N[,iN]) )
        
        if (sum(want_iN) == 1){
          temp_N <- as.matrix(indep_iN[want_iN,])%*%res_N[want_iN,iN]
          W_N = W_N + temp_N %*% t(temp_N)
        }else{
          W_N = W_N + t(indep_iN[want_iN,])%*%res_N[want_iN,iN]%*%t(res_N[want_iN,iN])%*%indep_iN[want_iN,]
        }
      }
      
      want_NT <- !( rowSums(is.na(indep_var)) > 0 | is.na(res.vec) )
      smp <- length(res.vec[want_NT,])
      
      if (two_way_cluster == F){
        
        if (robust == F){
          W = W_N
        } else if (robust == T){
          W = W_N *(N/(N-1)*(smp-1)/(smp-ncol(indep_var)))
        }
        
      } else if (two_way_cluster == T){
        
        #transfer to N*T
        indep_var_NT = NULL 
        for (iv in 1:ncol(indep_var)){
          indep_var.iv.tnl <- indep_var[,iv]
          indep_var.iv.nt <- t(matrix(indep_var.iv.tnl, ncol = N))
          indep_var.iv <- matrix(indep_var.iv.nt, nrow = length(yt) )
          indep_var_NT = cbind(indep_var_NT,indep_var.iv)
        }
        
        W_T = matrix(0,ncol(indep_var),ncol(indep_var))
        res_T <- t(matrix(res.vec,ncol = N))
        for (iT in 1:T0){
          indep_iT <- as.matrix(indep_var_NT[(iT-1)*N+(1:N),])
          want_iT <- !( rowSums(is.na(indep_iT)) > 0 | is.na(res_T[,iT]) )
          
          if (sum(want_iT) == 1){
            temp_T <- as.matrix(indep_iT[want_iT,])%*%res_T[want_iT,iT]
            W_T = W_T + temp_T %*% t(temp_T)
          }else{
            W_T = W_T + t(indep_iT[want_iT,])%*%res_T[want_iT,iT]%*%t(res_T[want_iT,iT])%*%indep_iT[want_iT,]
          }
        }
        
        W_NT = t(indep_var[want_NT,]*res.vec[want_NT,])%*%(indep_var[want_NT,]*res.vec[want_NT,])
        
        if (robust == F){
          W = W_N+W_T-W_NT
        } else if (robust == T){
          W = (W_N  + W_T  - W_NT ) *(min(N,T0)/(min(N,T0)-1)*(smp-1)/(smp-N-ncol(indep_var)))
        }
      }
      
      temp <- t(indep_var[want_NT,])%*%indep_var[want_NT,]
      var0_mat <- solve(temp)%*%W%*%solve(temp)
      
      if (eigens == F) {
        var_mat <- var0_mat
      } else if (eigens == T) {
        #eigens
        eigens <- eigen(var0_mat)
        U = eigens$vectors
        Lambda = diag(eigens$values * (eigens$values >= 0))
        var_mat = U %*% Lambda %*% solve(U)
      }
      
      if (h == 0){
        se[ih,] = sqrt(diag( as.matrix(var_mat) )[1:px])
      }else{
        se[ih,] = sqrt(diag( as.matrix(var_mat) )[(1:px-1)*lagX + 1])
      }
      
    } else if (method == "SPJ") {
      
      #demean_t
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
      
      if (te ==F) {
        
        dep_var = yt.dm
        indep_var <- cbind(X.dm,CC.dm,ylag.dm)
        
      } else if (te ==T) {
        
        #demean_n
        yt.dm.nt <- t(matrix(yt.dm, ncol = N))
        yt.dm.ntl <- c(yt.dm.nt[1:N,])
        
        X.dm.ntl = NULL 
        for (ix in 1:ncol(X.dm)){
          x.dm.ix <- X.dm[,ix]
          x.dm.ix.nt <- t(matrix(x.dm.ix, ncol = N))
          x.dm.ix.ntl <- matrix(x.dm.ix.nt, nrow = length(yt) )
          X.dm.ntl = cbind(X.dm.ntl,x.dm.ix.ntl)
        }
        
        CC.dm.ntl = NULL 
        if (pc >0){
          for (ic in 1:ncol(CC.dm)){
            cc.dm.ic <- CC.dm[,ic]
            cc.dm.ic.nt <- t(matrix(cc.dm.ic, ncol = N))
            cc.dm.ic.ntl <- matrix(cc.dm.ic.nt, nrow = length(yt) )
            CC.dm.ntl = cbind(CC.dm.ntl,cc.dm.ic.ntl)
          }
        }
        
        ylag.dm.ntl <- NULL
        if (lagY >= 1){
          ylag.dm.nt <- t(matrix(ylag.dm, ncol = N))
          ylag.dm.ntl <- matrix(ylag.dm.nt, nrow = length(yt) )
        }
        
        yt.dmdm.ntl <- yt.dm.ntl
        X.dmdm.ntl <- X.dm.ntl
        CC.dmdm.ntl <- CC.dm.ntl
        ylag.dmdm.ntl <- ylag.dm.ntl
        
        
        for (iT in 1:lag_T){
          want <- 1:N + (iT-1)*N
          want.mean <- rowSums(is.na(cbind(yt.dm.ntl[want],X.dm.ntl[want,],CC.dm.ntl[want,],ylag.dm.ntl[want,]))) == 0
          yt.dmdm.ntl[want] <-  yt.dmdm.ntl[want] - mean(yt.dm.ntl[want[want.mean]])
          for (ix in 1:ncol(X.dm.ntl)){
            X.dmdm.ntl[want,ix] <-  X.dmdm.ntl[want,ix] - mean(X.dm.ntl[want[want.mean],ix])
          }
          
          if (pc >0){
            for (ic in 1:ncol(CC.dm.ntl)){
              CC.dmdm.ntl[want,ic] <-  CC.dmdm.ntl[want,ic] -  mean(CC.dm.ntl[want[want.mean],ic])
            }
          }
          if (lagY >= 1){
            for (iylag in 1:ncol(ylag.dm.ntl)){
              ylag.dmdm.ntl[want,iylag] <-  ylag.dmdm.ntl[want,iylag] - mean(ylag.dm.ntl[want[want.mean],iylag]) 
            }
          }
        }
        
        #transfer to T*N
        yt.dmdm.tn <- t(matrix(yt.dmdm.ntl, nrow = N))
        yt.dmdm <- c(yt.dmdm.tn[,1:N])
        
        X.dmdm = NULL 
        for (ix in 1:ncol(X.dmdm.ntl)){
          x.dmdm.ix.ntl <- X.dmdm.ntl[,ix]
          x.dmdm.ix.tn <- t(matrix(x.dmdm.ix.ntl, nrow = N))
          x.dmdm.ix <- matrix(x.dmdm.ix.tn, nrow = length(yt) )
          X.dmdm = cbind(X.dmdm,x.dmdm.ix)
        }
        
        CC.dmdm = NULL 
        if (pc >0){
          for (ic in 1:ncol(CC.dmdm.ntl)){
            cc.dmdm.ic.ntl <- CC.dmdm.ntl[,ic]
            cc.dmdm.ic.tn <- t(matrix(cc.dmdm.ic.ntl, nrow = N))
            cc.dmdm.ic <- matrix(cc.dmdm.ic.tn, nrow = length(yt) )
            CC.dmdm = cbind(CC.dmdm,cc.dmdm.ic)
          }
        }
        
        ylag.dmdm <- NULL
        if (lagY >= 1){
          ylag.dmdm.tn <- t(matrix(ylag.dmdm.ntl, nrow = N))
          ylag.dmdm <- matrix(ylag.dmdm.tn, nrow = length(yt) )
        }
        
        dep_var = yt.dmdm
        indep_var <- cbind(X.dmdm,CC.dmdm,ylag.dmdm)
        
      }
      
      fit_h=lm(dep_var ~ indep_var + 0)
      beta_all=fit_h$coefficients 
      
      ######################
      T0h <- nrow(y_h)
      
      smp=cbind(dep_var,indep_var)
      smpna = complete.cases(smp)
      smpna.ti <- matrix(smpna, ncol = N)
      cut = matrix(NA,1,ncol(smpna.ti))
      for (iN in 1:ncol(cut)) {
        cut[1,iN] = floor(median(which(smpna.ti[,iN]==TRUE)))
      }
      
      ################## Sample a 
      
      #demean_t
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
      
      if (te == F) {
        
        dep_var_a = yt.dm_a
        indep_var_a <- cbind(X.dm_a,CC.dm_a,ylag.dm_a)
        
      } else if (te == T) {
        
        #demean_n
        yt.dm_a.nt <- t(matrix(yt.dm_a, ncol = N))
        yt.dm_a.ntl <- c(yt.dm_a.nt[1:N,])
        
        X.dm_a.ntl = NULL 
        for (ix in 1:ncol(X.dm_a)){
          x.dm_a.ix <- X.dm_a[,ix]
          x.dm_a.ix.nt <- t(matrix(x.dm_a.ix, ncol = N))
          x.dm_a.ix.ntl <- matrix(x.dm_a.ix.nt, nrow = length(yt_a) )
          X.dm_a.ntl = cbind(X.dm_a.ntl,x.dm_a.ix.ntl)
        }
        
        CC.dm_a.ntl = NULL 
        if (pc >0){
          for (ic in 1:ncol(CC.dm_a)){
            cc.dm_a.ic <- CC.dm_a[,ic]
            cc.dm_a.ic.nt <- t(matrix(cc.dm_a.ic, ncol = N))
            cc.dm_a.ic.ntl <- matrix(cc.dm_a.ic.nt, nrow = length(yt_a) )
            CC.dm_a.ntl = cbind(CC.dm_a.ntl,cc.dm_a.ic.ntl)
          }
        }
        
        ylag.dm_a.ntl <- NULL
        if (lagY >= 1){
          ylag.dm_a.nt <- t(matrix(ylag.dm_a, ncol = N))
          ylag.dm_a.ntl <- matrix(ylag.dm_a.nt, nrow = length(yt_a) )
        }
        
        yt.dmdm_a.ntl <- yt.dm_a.ntl
        X.dmdm_a.ntl <- X.dm_a.ntl
        CC.dmdm_a.ntl <- CC.dm_a.ntl
        ylag.dmdm_a.ntl <- ylag.dm_a.ntl
        
        
        for (iT in 1:lag_T){
          want <- 1:N + (iT-1)*N
          want.mean <- rowSums(is.na(cbind(yt.dm_a.ntl[want],X.dm_a.ntl[want,],CC.dm_a.ntl[want,],ylag.dm_a.ntl[want,]))) == 0
          yt.dmdm_a.ntl[want] <-  yt.dmdm_a.ntl[want] - mean(yt.dm_a.ntl[want[want.mean]])
          for (ix in 1:ncol(X.dm_a.ntl)){
            X.dmdm_a.ntl[want,ix] <-  X.dmdm_a.ntl[want,ix] - mean(X.dm_a.ntl[want[want.mean],ix])
          }
          
          if (pc >0){
            for (ic in 1:ncol(CC.dm_a.ntl)){
              CC.dmdm_a.ntl[want,ic] <-  CC.dmdm_a.ntl[want,ic] -  mean(CC.dm_a.ntl[want[want.mean],ic])
            }
          }
          if (lagY >= 1){
            for (iylag in 1:ncol(ylag.dm_a.ntl)){
              ylag.dmdm_a.ntl[want,iylag] <-  ylag.dmdm_a.ntl[want,iylag] - mean(ylag.dm_a.ntl[want[want.mean],iylag]) 
            }
          }
        }
        
        #transfer to T*N
        yt.dmdm_a.tn <- t(matrix(yt.dmdm_a.ntl, nrow = N))
        yt.dmdm_a <- c(yt.dmdm_a.tn[,1:N])
        
        X.dmdm_a = NULL 
        for (ix in 1:ncol(X.dmdm_a.ntl)){
          x.dmdm_a.ix.ntl <- X.dmdm_a.ntl[,ix]
          x.dmdm_a.ix.tn <- t(matrix(x.dmdm_a.ix.ntl, nrow = N))
          x.dmdm_a.ix <- matrix(x.dmdm_a.ix.tn, nrow = length(yt_a) )
          X.dmdm_a = cbind(X.dmdm_a,x.dmdm_a.ix)
        }
        
        CC.dmdm_a = NULL 
        if (pc >0){
          for (ic in 1:ncol(CC.dmdm_a.ntl)){
            cc.dmdm_a.ic.ntl <- CC.dmdm_a.ntl[,ic]
            cc.dmdm_a.ic.tn <- t(matrix(cc.dmdm_a.ic.ntl, nrow = N))
            cc.dmdm_a.ic <- matrix(cc.dmdm_a.ic.tn, nrow = length(yt_a) )
            CC.dmdm_a = cbind(CC.dmdm_a,cc.dmdm_a.ic)
          }
        }
        
        ylag.dmdm_a <- NULL
        if (lagY >= 1){
          ylag.dmdm_a.tn <- t(matrix(ylag.dmdm_a.ntl, nrow = N))
          ylag.dmdm_a <- matrix(ylag.dmdm_a.tn, nrow = length(yt_a) )
        }
        
        dep_var_a = yt.dmdm_a
        indep_var_a <- cbind(X.dmdm_a,CC.dmdm_a,ylag.dmdm_a)
        
      }
      
      fit_h_a=lm(dep_var_a ~ indep_var_a + 0)
      beta_a = fit_h_a$coefficients 
      
      ################## Sample b
      
      #demean_t
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
      
      if (te == F) {
        
        dep_var_b = yt.dm_b
        indep_var_b <- cbind(X.dm_b,CC.dm_b,ylag.dm_b)
        
      } else if (te == T) {
        
        #demean_n
        yt.dm_b.nt <- t(matrix(yt.dm_b, ncol = N))
        yt.dm_b.ntl <- c(yt.dm_b.nt[1:N,])
        
        X.dm_b.ntl = NULL 
        for (ix in 1:ncol(X.dm_b)){
          x.dm_b.ix <- X.dm_b[,ix]
          x.dm_b.ix.nt <- t(matrix(x.dm_b.ix, ncol = N))
          x.dm_b.ix.ntl <- matrix(x.dm_b.ix.nt, nrow = length(yt_b) )
          X.dm_b.ntl = cbind(X.dm_b.ntl,x.dm_b.ix.ntl)
        }
        
        CC.dm_b.ntl = NULL 
        if (pc >0){
          for (ic in 1:ncol(CC.dm_b)){
            cc.dm_b.ic <- CC.dm_b[,ic]
            cc.dm_b.ic.nt <- t(matrix(cc.dm_b.ic, ncol = N))
            cc.dm_b.ic.ntl <- matrix(cc.dm_b.ic.nt, nrow = length(yt_b) )
            CC.dm_b.ntl = cbind(CC.dm_b.ntl,cc.dm_b.ic.ntl)
          }
        }
        
        ylag.dm_b.ntl <- NULL
        if (lagY >= 1){
          ylag.dm_b.nt <- t(matrix(ylag.dm_b, ncol = N))
          ylag.dm_b.ntl <- matrix(ylag.dm_b.nt, nrow = length(yt_b) )
        }
        
        yt.dmdm_b.ntl <- yt.dm_b.ntl
        X.dmdm_b.ntl <- X.dm_b.ntl
        CC.dmdm_b.ntl <- CC.dm_b.ntl
        ylag.dmdm_b.ntl <- ylag.dm_b.ntl
        
        
        for (iT in 1:lag_T){
          want <- 1:N + (iT-1)*N
          want.mean <- rowSums(is.na(cbind(yt.dm_b.ntl[want],X.dm_b.ntl[want,],CC.dm_b.ntl[want,],ylag.dm_b.ntl[want,]))) == 0
          yt.dmdm_b.ntl[want] <-  yt.dmdm_b.ntl[want] - mean(yt.dm_b.ntl[want[want.mean]])
          for (ix in 1:ncol(X.dm_b.ntl)){
            X.dmdm_b.ntl[want,ix] <-  X.dmdm_b.ntl[want,ix] - mean(X.dm_b.ntl[want[want.mean],ix])
          }
          
          if (pc >0){
            for (ic in 1:ncol(CC.dm_b.ntl)){
              CC.dmdm_b.ntl[want,ic] <-  CC.dmdm_b.ntl[want,ic] -  mean(CC.dm_b.ntl[want[want.mean],ic])
            }
          }
          if (lagY >= 1){
            for (iylag in 1:ncol(ylag.dm_b.ntl)){
              ylag.dmdm_b.ntl[want,iylag] <-  ylag.dmdm_b.ntl[want,iylag] - mean(ylag.dm_b.ntl[want[want.mean],iylag]) 
            }
          }
        }
        
        #transfer to T*N
        yt.dmdm_b.tn <- t(matrix(yt.dmdm_b.ntl, nrow = N))
        yt.dmdm_b <- c(yt.dmdm_b.tn[,1:N])
        
        X.dmdm_b = NULL 
        for (ix in 1:ncol(X.dmdm_b.ntl)){
          x.dmdm_b.ix.ntl <- X.dmdm_b.ntl[,ix]
          x.dmdm_b.ix.tn <- t(matrix(x.dmdm_b.ix.ntl, nrow = N))
          x.dmdm_b.ix <- matrix(x.dmdm_b.ix.tn, nrow = length(yt_b) )
          X.dmdm_b = cbind(X.dmdm_b,x.dmdm_b.ix)
        }
        
        CC.dmdm_b = NULL 
        if (pc >0){
          for (ic in 1:ncol(CC.dmdm_b.ntl)){
            cc.dmdm_b.ic.ntl <- CC.dmdm_b.ntl[,ic]
            cc.dmdm_b.ic.tn <- t(matrix(cc.dmdm_b.ic.ntl, nrow = N))
            cc.dmdm_b.ic <- matrix(cc.dmdm_b.ic.tn, nrow = length(yt_b) )
            CC.dmdm_b = cbind(CC.dmdm_b,cc.dmdm_b.ic)
          }
        }
        
        ylag.dmdm_b <- NULL
        if (lagY >= 1){
          ylag.dmdm_b.tn <- t(matrix(ylag.dmdm_b.ntl, nrow = N))
          ylag.dmdm_b <- matrix(ylag.dmdm_b.tn, nrow = length(yt_b) )
        }
        
        dep_var_b = yt.dmdm_b
        indep_var_b <- cbind(X.dmdm_b,CC.dmdm_b,ylag.dmdm_b)
        
      }
      
      fit_h_b=lm(dep_var_b ~ indep_var_b + 0)
      beta_b = fit_h_b$coefficients 
      
      ######### estimate IRF 
      beta.hat = 2*beta_all-0.5*(beta_a+beta_b)
      
      if (h == 0){
        IRF[ih,] = beta.hat[1:px]
      }else{
        IRF[ih,] = beta.hat[(1:px-1)*lagX + 1]
      }
      ################### s.e. 
      
      res.vec <- dep_var - indep_var %*% beta.hat
      
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
      
      R_N.hat = matrix(0,ncol(indep_var),ncol(indep_var))
      res_N <- matrix(res.vec, ncol = N)
      for (iN in 1:N){ 
        dd_iN <- as.matrix(dd.mat[(iN-1)*nrow(y_h)+(1:nrow(y_h)),])
        want_iN <- !( rowSums(is.na(dd_iN)) > 0 | is.na(res_N[,iN]) )
        
        if (sum(want_iN) == 1){
          temp <- as.matrix(dd_iN[want_iN,])%*%res_N[want_iN,iN]
          R_N.hat = R_N.hat + temp %*% t(temp)
        }else{
          R_N.hat = R_N.hat+t(dd_iN[want_iN,])%*%res_N[want_iN,iN]%*%t(res_N[want_iN,iN])%*%dd_iN[want_iN,]
        }
      }
      
      want_NT <- !( rowSums(is.na(indep_var)) > 0 | is.na(res.vec) )
      smp <- length(res.vec[want_NT,])
      
      if (two_way_cluster == F){
        
        if (robust == F){
          R.hat = R_N.hat
        } else if (robust == T){
          R.hat = R_N.hat *(N/(N-1)*(smp-1)/(smp-ncol(indep_var)))
        }
        
      } else if (two_way_cluster == T){
        
        #transfer to N*T
        dd.mat_NT = NULL 
        for (iv in 1:ncol(dd.mat)){
          dd.mat.iv.tnl <- dd.mat[,iv]
          dd.mat.iv.nt <- t(matrix(dd.mat.iv.tnl, ncol = N))
          dd.mat.iv <- matrix(dd.mat.iv.nt, nrow = length(yt) )
          dd.mat_NT = cbind(dd.mat_NT,dd.mat.iv)
        }
        
        R_T.hat = matrix(0,ncol(indep_var),ncol(indep_var))
        res_T <- t(matrix(res.vec,ncol = N))
        for (iT in 1:T0){
          dd_iT <- as.matrix(dd.mat_NT[(iT-1)*N+(1:N),])
          want_iT <- !( rowSums(is.na(dd_iT)) > 0 | is.na(res_T[,iT]) )
          
          if (sum(want_iT) == 1){
            temp_T <- as.matrix(dd_iT[want_iT,])%*%res_T[want_iT,iT]
            R_T.hat = R_T.hat + temp_T %*% t(temp_T)
          }else{
            R_T.hat = R_T.hat + t(dd_iT[want_iT,])%*%res_T[want_iT,iT]%*%t(res_T[want_iT,iT])%*%dd_iT[want_iT,]
          }
        }
        
        R_NT.hat = t(dd.mat[want_NT,]*res.vec[want_NT,])%*%(dd.mat[want_NT,]*res.vec[want_NT,])
        
        if (robust == F){
          R.hat = R_N.hat  + R_T.hat  - R_NT.hat
        } else if (robust == T){
          R.hat = (R_N.hat  + R_T.hat  - R_NT.hat ) *(min(N,T0)/(min(N,T0)-1)*(smp-1)/(smp-N-ncol(indep_var)))
        }
        
      }
      
      R.hat <- R.hat
      Q.hat = t(indep_var[want_NT,]) %*% indep_var[want_NT,]
      var0.hat <- solve(Q.hat) %*% R.hat %*% solve(Q.hat)
      
      if (eigens == F) {
        var.hat <- var0.hat
      } else if (eigens == T) {
        #eigens
        eigens <- eigen(var0.hat)
        U = eigens$vectors
        Lambda = diag(eigens$values * (eigens$values >= 0))
        var.hat = U %*% Lambda %*% solve(U)
      }
      
      if (h == 0){
        se[ih,] = sqrt( diag(var.hat )[1:px])
      }else{
        se[ih,] = sqrt( diag(var.hat )[(1:px-1)*lagX + 1])
      }
      
    } else{
      cat("pls speficy the method")
    }
  }
  
  results=list(IRF,se,X.name,h.seq,weakiv)
  names(results)=c("IRF","se","X.name","h.seq","weakiv")
  return(results)
}