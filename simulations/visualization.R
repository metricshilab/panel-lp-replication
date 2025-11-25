rm(list=ls())
setwd("D:\\GitHub\\panel-lp-replication\\simulations")

beta = -0.6
dgp = "toy"

rho_seq_vals = c(0, 0.2, 0.5, 0.8)
N_seq = c(30, 30, 50)
T_seq = c(60, 120, 120)
col_vec = c("blue", "firebrick3", "black")
lwd = 1.5
cex = 2
cex1 = 1.5
cex2 = 1.5

library("readxl")

rho_labels <- list(
  expression(rho == 0),
  expression(rho == 0.2),
  expression(rho == 0.5),
  expression(rho == 0.8)
)

#-------------------------------------------------------------------------------
# 1. Figure 1
#-------------------------------------------------------------------------------
file_xlsx <- paste0(dgp, "_", beta, ".xlsx")
toy_results <- as.data.frame(read_excel(file_xlsx, sheet = 3)) 

ylim_irf <- c(-0.7, 0.1)
file_eps <- paste0(dgp, "_trueIRF_beta", beta, ".eps")

setEPS()
postscript(file_eps, width = 12, height = 15, family = "Times")
par(oma = c(7, 7, 7, 7), mfrow = c(4, 3), mar = c(1, 2, 3, 2), family = "Times")

for (irho in 1:4) {
  rho <- rho_seq_vals[irho]
  
  for (iTN in 1:3) {
    irow <- (irho - 1) * 3 + iTN
    
    y_fe  <- as.numeric(toy_results[irow, seq(4, 36, by = 3)])
    y_spj <- as.numeric(toy_results[irow, seq(5, 36, by = 3)])
    y_db  <- as.numeric(toy_results[irow, seq(6, 36, by = 3)])
    
    main_title <- ""
    if (irho == 1) {
      main_title <- paste0("N=", N_seq[iTN], ", T=", T_seq[iTN])
    }
    
    matplot(0:10, y_fe, type = "o", pch = 21,
            col = "blue", lty = 1, lwd = lwd,
            main = main_title, ylab = "", xlab = "horizon",
            ylim = ylim_irf, cex.main = cex, cex.axis = cex1, family = "Times")
    
    lines(0:10, y_db, type = "o", col = "black", lwd = lwd, pch = 23)
    
    IRFs <- beta * rho^(0:10)
    lines(0:10, IRFs, col = "#009900", type = "o", lwd = lwd, pch = 22)
    
    lines(0:10, y_spj, type = "o", col = "firebrick3", lwd = lwd, pch = 24)
    lines(0:10, rep(0, 11), lty = 2)
    
    if (iTN == 1) {
      mtext(rho_labels[[irho]], side = 2, line = 3, las = 0, cex = cex2)
    }
  }
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x = 'bottom',
       c("FE", "SPJ", "DB", "true IRF"),
       col = c(col_vec, "#009900"), lty = 1, lwd = lwd, cex = cex, 
       horiz = TRUE, bty = "n", pch = c(21, 24, 23, 22))
dev.off()


#-------------------------------------------------------------------------------
# 2. Figure 2
#-------------------------------------------------------------------------------
toy_results <- as.data.frame(read_excel(file_xlsx, sheet = 2))
ylim_rmse <- c(0, 0.1)
file_eps <- paste0(dgp, "_RMSE_beta", beta, ".eps")

setEPS()
postscript(file_eps, width = 12, height = 15, family = "Times")
par(oma = c(7, 7, 7, 7), mfrow = c(4, 3), mar = c(1, 2, 3, 2), family = "Times")

for (irho in 1:4) {
  for (iTN in 1:3) {
    irow <- (irho - 1) * 3 + iTN
    
    y_fe  <- as.numeric(toy_results[irow, seq(4, 36, by = 3)])
    y_spj <- as.numeric(toy_results[irow, seq(5, 36, by = 3)])
    y_db  <- as.numeric(toy_results[irow, seq(6, 36, by = 3)])
    
    main_title <- ifelse(irho == 1, paste0("N=", N_seq[iTN], ", T=", T_seq[iTN]), "")
    
    matplot(0:10, y_fe, type = "o", pch = 21,
            col = "blue", lty = 1, lwd = lwd,
            main = main_title, ylab = "", xlab = "horizon",
            ylim = ylim_rmse, cex.main = cex, cex.axis = cex1, family = "Times")
    
    lines(0:10, y_db, type = "o", col = "black", lwd = lwd, pch = 23)
    lines(0:10, y_spj, type = "o", col = "firebrick3", lwd = lwd, pch = 24)
    
    if (iTN == 1) {
      mtext(rho_labels[[irho]], side = 2, line = 3, las = 0, cex = cex2)
    }
  }
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x = 'bottom',
       c("FE", "SPJ", "DB"),
       col = col_vec, lty = 1, lwd = lwd, cex = cex, 
       horiz = TRUE, bty = "n", pch = c(21, 24, 23))
dev.off()


#-------------------------------------------------------------------------------
# 3. Figure 3
#-------------------------------------------------------------------------------
toy_results <- as.data.frame(read_excel(file_xlsx, sheet = 1))
ylim_cov <- c(0, 1)
file_eps <- paste0(dgp, "_coverage_beta", beta, ".eps")

setEPS()
postscript(file_eps, width = 12, height = 15, family = "Times")
par(oma = c(7, 7, 7, 7), mfrow = c(4, 3), mar = c(1, 2, 3, 2), family = "Times")

for (irho in 1:4) {
  for (iTN in 1:3) {
    irow <- (irho - 1) * 3 + iTN
    
    y_fe  <- as.numeric(toy_results[irow, seq(4, 36, by = 3)])
    y_spj <- as.numeric(toy_results[irow, seq(5, 36, by = 3)])
    y_db  <- as.numeric(toy_results[irow, seq(6, 36, by = 3)])
    
    main_title <- ifelse(irho == 1, paste0("N=", N_seq[iTN], ", T=", T_seq[iTN]), "")
    
    matplot(0:10, y_fe, type = "o", pch = 21,
            col = "blue", lty = 1, lwd = lwd,
            main = main_title, ylab = "", xlab = "horizon",
            ylim = ylim_cov, cex.main = cex, cex.axis = cex1, family = "Times")
    
    lines(0:10, y_db, type = "o", col = "black", lwd = lwd, pch = 23)
    lines(0:10, rep(0.95, 11), lty = 2)
    lines(0:10, y_spj, type = "o", col = "firebrick3", lwd = lwd, pch = 24)
    
    if (iTN == 1) {
      mtext(rho_labels[[irho]], side = 2, line = 3, las = 0, cex = cex2)
    }
  }
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x = 'bottom',
       c("FE", "SPJ", "DB"),
       col = col_vec, lty = 1, lwd = lwd, cex = cex, 
       horiz = TRUE, bty = "n", pch = c(21, 24, 23))
dev.off()