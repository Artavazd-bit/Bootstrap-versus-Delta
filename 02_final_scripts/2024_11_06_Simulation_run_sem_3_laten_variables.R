library(MASS)
library(cSEM)
library(cSEM.DGP)
library(lavaan)
library(foreach)
library(doParallel)
library(diptest)

model_base <- "
# Structural model
Y ~ a*X1 + b*X2
X1 ~~ c*X2
# Measurement model
# Ladungen zwischen -1 und 1 
# Werte zwischen 0.6 und 0.9
X1 <~ 0.6*x11 + 0.81*x12 + 0.778*x13
x11 ~~ 0.3*x12 + 0.4*x13
x12 ~~d*x13
X2 <~ 0.8*x21 + 0.9*x22 + 0.3*x23
x21 ~~ 0.5*x22 + 0.3*x23
x22 ~~e*x23
Y <~ 0.9*y1 + 0.6*y2 + 0.4*y3
y1 ~~ 0.2*y2 + 0.3*y3
y2 ~~ 0.1*y3
"

model_est <- "
# Structural model
Y ~ X1 + X2
# Measurement model
# Ladungen zwischen -1 und 1 
# Werte zwischen 0.6 und 0.9
X1 <~ x11 + x12 + x13
X2 <~ x21 + x22 + x23
Y <~ y1 + y2 + y3
"

model_est_sem <- "
# Structural model
Y ~ X1 + X2
# Measurement model
# Ladungen zwischen -1 und 1 
# Werte zwischen 0.6 und 0.9
X1 =~ x11 + x12 + x13
X2 =~ x21 + x22 + x23
Y =~ y1 + y2 + y3
"

Anzahl_path_estimate <- 2
set.seed(123)
sim_data <- generateData(model_base, 
                         .N= 5000, 
                         .empirical = TRUE, 
                         a = c(0, 0.1, 0.2, 0.3), 
                         b = c(0.5),
                         c = c(0, 0.1, 0.2, 0.3),
                         d = c(0.4),
                         e = c(0.5),
                         .return_type = "cor")

cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

o_table <- foreach(jj = 1: 1, .packages = c("cSEM", "MASS"), .combine = "rbind") %:%
  foreach(n = c(500), .combine = "rbind") %:%
  foreach(sim_runs = 1:3, .combine = "rbind") %dopar% 
  {
    set.seed(50+jj+sim_runs+n)
    # Ziehe Daten aus einer Normalverteilung mit Varianz-Covarianz-Matrix aus sim_data  
    data_sim <- MASS::mvrnorm(n = n, 
                              mu= rep(0, nrow(sim_data$dgp[[jj]])), 
                              Sigma =  sim_data$dgp[[jj]],
                              empirical = F
    )
    res_sem <- sem(model = model_est_sem, 
                   data = data_sim
    )
    paramees <- parameterEstimates(res_sem)
    paramees2 <- paramees[paramees$lhs == "Y" & paramees$op == "~" & paramees$rhs == "X1",]
  
    # Schätze mit den aus den oben gezogenen Daten die Parameter
    res <- csem(.data = data_sim, 
                .model = model_est,
                .resample_method = 'bootstrap',
                .R = 500,   
                .PLS_weight_scheme_inner = 'factorial',
                .tolerance = 1e-06
    )
    # Berechne die Konfidenzintervalle aus den Bootstrap-Werten
    infer_res <- infer(res, .alpha = 0.05, .quantity = "CI_percentile")
    # Delta-Methode 
    # Für die Delta-Methode werden die Estimates aus obiger Schätzung als Referenzwert verwendet.
    bt = res$Estimates$Path_estimates["Y", 1:Anzahl_path_estimate]
    ####################################################################################
    ## Berechnung der Varianz-Covarianz-Matrix der Korrelationskoeffizienten gilt nur asymptotisch Dykstra(2013) S.11
    ## nach Dykstra(2013) Gleichung 27/28 und Isserlis(2019) Gleichung 21 + Anmerkungen von Flo
    ####################################################################################
    size_n_ <- ((ncol(data_sim)*ncol(data_sim)) - ncol(data_sim))/2
    vc_r <- matrix(0 , nrow = size_n_ , ncol = size_n_)
    cor_data_sim <- cor(data_sim)
    cnter_y <- 1
    cnter_t <- 1
    for(x in 1:(ncol(data_sim)-1)){
      bar_x <- mean(data_sim[,x])
      sigma_x <- sd(data_sim[,x])
      for(y in (x+1):ncol(data_sim)){
        bar_y <- mean(data_sim[,y])
        sigma_y <- sd(data_sim[,y])
        for(z in 1:(ncol(data_sim)-1)){
          bar_z <- mean(data_sim[,z])
          sigma_z <- sd(data_sim[,z])
          for(t in (z+1):ncol(data_sim)){
            bar_t <- mean(data_sim[,t])
            sigma_t <- sd(data_sim[,t])
            
            mu_xyzt_temp <- 0
            mu_xxzt_temp <- 0 
            mu_yyzt_temp <- 0
            mu_xyzz_temp <- 0 
            mu_xytt_temp <- 0
            
            mu_xxtt_temp <- 0
            mu_xxzz_temp <- 0
            
            mu_yytt_temp <- 0
            mu_yyzz_temp <- 0
            
            for(i in 1:nrow(data_sim)){
              mu_xyzt_temp <- mu_xyzt_temp + (data_sim[i,x] - bar_x) * (data_sim[i,y] - bar_y) * (data_sim[i,z] - bar_z) * (data_sim[i,t] - bar_t)
              
              mu_xxzt_temp <- mu_xxzt_temp + (data_sim[i,x] - bar_x) * (data_sim[i,x] - bar_x) * (data_sim[i,z] - bar_z) * (data_sim[i,t] - bar_t)
              mu_yyzt_temp <- mu_yyzt_temp + (data_sim[i,y] - bar_y) * (data_sim[i,y] - bar_y) * (data_sim[i,z] - bar_z) * (data_sim[i,t] - bar_t)
              mu_xyzz_temp <- mu_xyzz_temp + (data_sim[i,x] - bar_x) * (data_sim[i,y] - bar_y) * (data_sim[i,z] - bar_z) * (data_sim[i,z] - bar_z)
              mu_xytt_temp <- mu_xytt_temp + (data_sim[i,x] - bar_x) * (data_sim[i,y] - bar_y) * (data_sim[i,t] - bar_t) * (data_sim[i,t] - bar_t)
              
              mu_xxtt_temp <- mu_xxtt_temp + (data_sim[i,x] - bar_x) * (data_sim[i,x] - bar_x) * (data_sim[i,t] - bar_t) * (data_sim[i,t] - bar_t)
              mu_xxzz_temp <- mu_xxzz_temp + (data_sim[i,x] - bar_x) * (data_sim[i,x] - bar_x) * (data_sim[i,z] - bar_z) * (data_sim[i,z] - bar_z)
              
              mu_yytt_temp <- mu_yytt_temp + (data_sim[i,y] - bar_y) * (data_sim[i,y] - bar_y) * (data_sim[i,t] - bar_t) * (data_sim[i,t] - bar_t)
              mu_yyzz_temp <- mu_yyzz_temp + (data_sim[i,y] - bar_y) * (data_sim[i,y] - bar_y) * (data_sim[i,z] - bar_z) * (data_sim[i,z] - bar_z)
            }
            
            mu_xyzt <- 1/nrow(data_sim) * mu_xyzt_temp
            
            mu_xxzt <- 1/nrow(data_sim) * mu_xxzt_temp
            
            mu_yyzt <- 1/nrow(data_sim) * mu_yyzt_temp
            
            mu_xyzz <- 1/nrow(data_sim) * mu_xyzz_temp
            
            mu_xytt <- 1/nrow(data_sim) * mu_xytt_temp
            
            
            mu_xxtt <- 1/nrow(data_sim) * mu_xxtt_temp
            mu_xxzz <- 1/nrow(data_sim) * mu_xxzz_temp
            
            mu_yytt <- 1/nrow(data_sim) * mu_yytt_temp
            mu_yyzz <- 1/nrow(data_sim) * mu_yyzz_temp
            
            vc_r[cnter_y, cnter_t]<- (mu_xyzt - 1/2*cor_data_sim[x,y]*(mu_xxzt + mu_yyzt) - 1/2*cor_data_sim[z,t]*(mu_xyzz + mu_xytt) 
                                      + 1/4*cor_data_sim[x,y]*cor_data_sim[z,t]*(mu_xxzz+mu_xxtt+mu_yyzz+mu_yytt))/nrow(data_sim)
            cnter_t <- cnter_t + 1
          }
        }
        cnter_y <- cnter_y + 1
        cnter_t <- 1
      }
    }
    ################################################################################
    # Ein paar Werte werden initialisiert
    # für den Gradienten der pfadestimates
    gr = matrix(0 , nrow = size_n_ , ncol = 1)
    gdr = matrix(0 , nrow = size_n_ , ncol = 1)
    # Path estimates
    Gbt = matrix(0 , nrow = Anzahl_path_estimate , ncol = size_n_)
    Gbt2 = matrix(0 , nrow = Anzahl_path_estimate, ncol = size_n_)
    # für die h - Methode
    list_dh_all <- c(0.00001)
    #Übersicht für verschiedene dh
    Gbt_overview_delta_x1 <- data.frame(matrix(0, nrow = size_n_, ncol = 1))
    colnames(Gbt_overview_delta_x1) <- list_dh_all
    
    Gbt_overview_delta_x2 <-data.frame(matrix(0, nrow = size_n_, ncol = 1))
    colnames(Gbt_overview_delta_x2) <- list_dh_all
    
    
    # Counter, sodass jeder Korrelationskoeffizient durchgegangen wird
    cnter = 1
    for(dh in list_dh_all){
      #print(dh)
      for (i in 1:(ncol(data_sim)-1)) {
        for (j in (i+1):ncol(data_sim)) {
          cor_sim1 <- cor(data_sim)
          cor_sim1[i,j] = cor_sim1[i,j] + dh
          cor_sim1[j,i] = cor_sim1[i,j]
          set.seed(50+jj+sim_runs+n)
          data_after_dh <- MASS::mvrnorm(n = n, 
                                         mu = rep(0,nrow(sim_data$dgp[[jj]])), 
                                         Sigma = cor_sim1, 
                                         empirical = TRUE
          )
          out1 <- csem(.data = data_after_dh, 
                       .model = model_est,
                       .PLS_weight_scheme_inner = 'factorial',
                       .tolerance = 1e-06
          )
          # in h-Methode: f(x+h):
          bt1 = out1$Estimates$Path_estimates["Y",1:Anzahl_path_estimate]
          # f(x+h) - f(x) / h
          dlta_bt = (bt1 - bt)/dh
          Gbt[,cnter] = dlta_bt[1]
          Gbt2[,cnter] = dlta_bt[2]
          cnter = cnter + 1
        }
      }
      cnter <- 1
      Gbt_overview_delta_x1[, paste0(dh)]  <- Gbt[1,]
      Gbt_overview_delta_x2[, paste0(dh)]  <- Gbt2[1,]
      Gbt = matrix(0 , nrow = Anzahl_path_estimate , ncol = size_n_)
      Gbt2 = matrix(0 , nrow = Anzahl_path_estimate , ncol = size_n_)
    }
    
    Gbt_x1 <- t(as.matrix(Gbt_overview_delta_x1[,"1e-05"]))
    Gbt_x2 <- t(as.matrix(Gbt_overview_delta_x2[,"1e-05"]))
    
    dip <- diptest::dip.test(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,'Y ~ X1'])
    
    csv_write <- data.frame(a = sim_data$a[jj],
                            b = sim_data$b[jj],
                            c = sim_data$c[jj],
                            n = n,
                            path_estimate_y_x1 = res$Estimates$Estimates_resample$Estimates1$Path_estimates$Original['Y ~ X1'], 
                            path_estimate_y_x2 = res$Estimates$Estimates_resample$Estimates1$Path_estimates$Original['Y ~ X2'],
                            sd_bootstrap_y_x1 = sd(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,'Y ~ X1']), 
                            Z_test_boot_y_x1 = abs(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Original['Y ~ X1']/ sd(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,'Y ~ X1'])) > qnorm(0.975),
                            Z_test_boot_y_x2 = abs(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Original['Y ~ X2']/ sd(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,'Y ~ X2'])) > qnorm(0.975),
                            perc_ci_lb_boot_y_x1 = infer_res$Path_estimates$CI_percentile["95%L", "Y ~ X1"], 
                            perc_ci_ub_boot_y_x1 = infer_res$Path_estimates$CI_percentile["95%U", "Y ~ X1"], 
                            perc_ci_sgn_same_y_x1 = sign(infer_res$Path_estimates$CI_percentile["95%L", "Y ~ X1"]) == sign(infer_res$Path_estimates$CI_percentile["95%U", "Y ~ X1"]), 
                            sd_delta_y_x1 =  sqrt(diag( Gbt_x1 %*% vc_r %*% t(Gbt_x1))), 
                            Z_test_delta_y_x1 = abs(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Original['Y ~ X1']/sqrt(diag( Gbt_x1 %*% vc_r %*% t(Gbt_x1) ))) > qnorm(0.975), 
                            sd_delta_y_x2 = sqrt(diag( Gbt_x2 %*% vc_r %*% t(Gbt_x2) )), 
                            Z_test_delta_y_x2 = abs(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Original['Y ~ X2']/sqrt(diag( Gbt_x2 %*% vc_r %*% t(Gbt_x2) ))) > qnorm(0.975),
                            Sim_run = sim_runs, 
                            simulation_run = jj, 
                            c_estimate = res$Estimates$Construct_VCV["X1", "X2"],
                            dip_test_p_value = dip$p.value,
                            ML_estimate  = paramees2$est,
                            ML_se = paramees2$se,
                            ML_z = paramees2$z,
                            ML_z_p_value = paramees2$pvalue,
                            ML_z_test  = paramees2$pvalue < 0.05
    )
    csv_write$list_with_bootstrap <- list(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,'Y ~ X1'])
    csv_write
  } 
closeAllConnections()
saveRDS(o_table, file = "./Data/2024_11_06_sem_3_latent.rds")
 
