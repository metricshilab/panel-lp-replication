rm(list=ls())
setwd("D:\\GitHub\\panel-lp-replication\\simulations")

if(!require(writexl)) install.packages("writexl")
library(writexl)

H=10
beta <- -0.6

COVER <- NULL
RMSE <- NULL 
IRFm <- NULL 
All_Results <- NULL

Labels <- NULL
col_names <- c()
for(h in 0:H){
  col_names <- c(col_names, paste0("h", h, "_FE"), paste0("h", h, "_SPJ"), paste0("h", h, "_DB"))
}

for (rho in c(0,0.2,0.5,0.8)){
  for (T0 in c(60,120)){
    for (N in c(30,50)){
      if (N == 50 && T0 == 60) next
      file.name <- paste0("T",T0,"N",N,"rho",rho*10,"tau",0,"Niter1000lagY0.RDS")
      out <- readRDS(file = file.name)
      
      cover1 <- rowMeans(out$cover$FE)
      cover2 <- rowMeans(out$cover$jackknife)
      cover3 <- rowMeans(out$cover$DB)
      cover <- as.vector(rbind(cover1,cover2,cover3))
      COVER <- rbind(COVER,cover)
      
      rmse1 <- sqrt(rowMeans(out$SE$FE))
      rmse2 <- sqrt(rowMeans(out$SE$jackknife))
      rmse3 <- sqrt(rowMeans(out$SE$DB))
      rmse <- as.vector(rbind(rmse1,rmse2,rmse3))
      RMSE <- rbind(RMSE,rmse)
      
      irfm1 <- rowMeans(out$irf$FE)
      irfm2 <- rowMeans(out$irf$jackknife)
      irfm3 <- rowMeans(out$irf$DB)
      irfm <- as.vector(rbind(irfm1,irfm2,irfm3))
      IRFm <- rbind(IRFm,irfm)
      
      Labels <- rbind(Labels, data.frame(rho = rho, T = T0, N = N))
    }
  }
}

df_cover <- cbind(Labels, as.data.frame(COVER))
colnames(df_cover)[4:ncol(df_cover)] <- col_names

df_rmse <- cbind(Labels, as.data.frame(RMSE))
colnames(df_rmse)[4:ncol(df_rmse)] <- col_names

df_irfm <- cbind(Labels, as.data.frame(IRFm))
colnames(df_irfm)[4:ncol(df_irfm)] <- col_names

sheets_list <- list(
  "COVER"   = df_cover,
  "RMSE"    = df_rmse,
  "IRFmean" = df_irfm
)
file_xlsx <- paste0("toy_", beta, ".xlsx")
write_xlsx(sheets_list, path = file_xlsx)
