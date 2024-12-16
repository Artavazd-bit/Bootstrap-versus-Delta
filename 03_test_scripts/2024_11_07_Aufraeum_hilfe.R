set.seed(50+2+1+100)
data_sim <- MASS::mvrnorm(n = 100, 
                          mu= rep(0, nrow(sim_data$dgp[[2]])), 
                          Sigma =  sim_data$dgp[[2]],
                          empirical = F
)
# Schätze mit den aus den oben gezogenen Daten die Parameter
res <- csem(.data = data_sim, 
            .model = model_est,
            .resample_method = 'bootstrap',
            .R = 500,   
            .PLS_weight_scheme_inner = 'factorial',
            .tolerance = 1e-06
)
# Delta-Methode 
# Für die Delta-Methode werden die Estimates aus obiger Schätzung als Referenzwert verwendet.
bt = res$Estimates$Path_estimates["X3", 1:Anzahl_path_estimate]
size_n_ <- ((ncol(data_sim)*ncol(data_sim)) - ncol(data_sim))/2
# für die h - Methode
dh <- 0.00001
Gbt = matrix(0 , nrow = 1 , ncol = size_n_)

# Counter, sodass jeder Korrelationskoeffizient durchgegangen wird
cnter = 1
  for (i in 1:(ncol(data_sim)-1)) {
    for (j in (i+1):ncol(data_sim)) {
      cor_sim1 <- cor(data_sim)
      cor_sim1[i,j] = cor_sim1[i,j] + dh
      cor_sim1[j,i] = cor_sim1[i,j]
      set.seed(50+2+1+100)
      data_after_dh <- MASS::mvrnorm(n = 100, 
                                     mu = rep(0,nrow(sim_data$dgp[[2]])), 
                                     Sigma = cor_sim1, 
                                     empirical = TRUE
      )
      out1 <- csem(.data = data_after_dh, 
                   .model = model_est,
                   .PLS_weight_scheme_inner = 'factorial',
                   .tolerance = 1e-06
      )
      # in h-Methode: f(x+h):
      bt1 = out1$Estimates$Path_estimates["X3",1:Anzahl_path_estimate]
      # f(x+h) - f(x) / h
      dlta_bt = (bt1 - bt)/dh
      Gbt[,cnter] = dlta_bt["X2"]
      cnter = cnter + 1
    }
  }
