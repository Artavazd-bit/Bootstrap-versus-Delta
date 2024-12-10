# Die VC für den Vektor der 45 Korrelationen wird
# durch Simulation und auch durch Bootstrap geschätzt

# Das u aus y = u %*% L ist Grundlage für alle
# Berechnungen, cor(y) wird partiell verändert. 

# Zudem werden die Bootstrap Resultate 
# aus csem analysiert


# requires cSem
# library(cSem)

setwd('cSEM')
stndrdz=function(q){
     n = nrow(as.matrix(q))
     cn = sqrt( (n-1) / n )
     # cn = sqrt( (n) / n ) # Ergebnisse bleiben auch mit dieser Zeile gleich !!!
     return( (q-mean(q)) / ( sd(q) * cn) )
}
Russett <- as.data.frame(readxl::read_excel("Russett.xlsx"))
model_Russett = ' # Specify composite models
              AgrIneq <~ gini + farm +rent
              IndDev  <~ gnpr + labo
              PolInst <~ inst + ecks + deat 
              + stab + dict 
              
              # Specify the relation among 
              # the emergent variables
              PolInst ~ AgrIneq + IndDev
              '
corRussett = cor(Russett)
cholcorRu = chol(corRussett)

r_sim = matrix(0,nrow = 4999, ncol = 45)

for (i in 1:4999) {
  u_sim = matrix(rnorm(470), nrow = 47 ,ncol = 10)
  y_sim = u_sim %*% (cholcorRu)
  cor_sim = cor(y_sim)
  r_sim[i,] = c(cor_sim[1,2:10],cor_sim[2,3:10],
                cor_sim[3,4:10],cor_sim[4,5:10],
                cor_sim[5,6:10],cor_sim[6,7:10],
                cor_sim[7,8:10],cor_sim[8,9:10],
                cor_sim[9,10])
}
vc_r_sim = cov(r_sim)
rm(r_sim)

# Chapter 3 (see 3.4.2)
out <- csem(.data = Russett, 
            .model = model_Russett,
            # To reproduce the Adanco results
            .PLS_weight_scheme_inner = 'factorial',
            #.PLS_weight_scheme_inner = 'factorial',
              .resample_method = 'bootstrap',
            .tolerance = 1e-06,
            .seed = -552280054
)
summarize(out)

out_2 <- csem(.data = Russett, 
            .model = model_Russett,
            # To reproduce the Adanco results
            .PLS_weight_scheme_inner = 'factorial',
            #.PLS_weight_scheme_inner = 'factorial',
            .resample_method = 'bootstrap',
            .tolerance = 1e-06,
            .dominant_indicators = c('IndDev'='gnpr','AgrIneq'='gini','PolInst'='inst'),
            .seed = -552280054
)

 summarize(out_2)

 # Analyse der Resampling Resulate
 Pfad_boot = out$Estimates$Estimates_resample$Estimates1
plot(density(Pfad_boot$Path_estimates$Resampled[,1]), col = 'red')
plot(density(Pfad_boot$Path_estimates$Resampled[,2]), col = 'red')

#colMeans(Pfad_boot$Path_estimates$Resampled)
#sd(abs(Pfad_boot$Path_estimates$Resampled[,1]))
#sd(abs(Pfad_boot$Path_estimates$Resampled[,2]))
colMeans(Pfad_boot$Path_estimates$Resampled)
apply(Pfad_boot$Path_estimates$Resampled, 2, FUN = sd)

colMeans(abs(Pfad_boot$Path_estimates$Resampled))
apply(abs(Pfad_boot$Path_estimates$Resampled), 2, FUN = sd)

b2_boot= Pfad_boot$Path_estimates$Resampled[,2]
w4_boot= Pfad_boot$Weight_estimates$Resampled[,4]
w5_boot= Pfad_boot$Weight_estimates$Resampled[,5]
sum(b2_boot<0)
sum(b2_boot<0&w4_boot>0)
sum(b2_boot<0&w5_boot<0)
sum(b2_boot<0&w4_boot>0&w5_boot<0)
sum(b2_boot>0&w4_boot<0&w5_boot>0)


 # Ende der Resampling Analyse

bt = out$Estimates$Path_estimates[3,1:2]
wt = c(out$Estimates$Loading_estimates[1,1:3],
       out$Estimates$Loading_estimates[2,4:5],
       out$Estimates$Loading_estimates[3,6:10])
lt = c(out$Estimates$Weight_estimates[1,1:3],
       out$Estimates$Weight_estimates[2,4:5],
       out$Estimates$Weight_estimates[3,6:10])
Russett1 = Russett
for (i in 1:10) {
  Russett1[,i] = stndrdz(Russett1[,i])   
}
corRussett = cor(Russett)
cholcorRu = chol(corRussett)
# sum(abs(t(cholcorRu) %*% (cholcorRu) -corRussett ))    # = 0
u_Russ = as.matrix(Russett1) %*% solve(cholcorRu)

gr = matrix(0 , nrow = 45 , ncol = 1)
gdr = matrix(0 , nrow = 45 , ncol = 1)
Gbt = matrix(0 , nrow = 2 , ncol = 45)
Gwt = matrix(0 , nrow = 10 , ncol = 45)
Glt = matrix(0 , nrow = 10 , ncol = 45)
dh = 1e-05
cnter = 1
for (i in 1:9) {
 for (j in (i+1):10) {
  corRu = corRussett
  gr[cnter,1] = corRu[i,j]
  gdr[cnter,1] = (1 - corRu[i,j]^2 )^2
  corRu[i,j] = corRu[i,j] + dh
  corRu[j,i] = corRu[i,j] 
  
  chlcrRu = chol(corRu)
  yin = u_Russ %*% chlcrRu
  out1 <- csem(.data = yin, 
            .model = model_Russett,
            # To reproduce the Adanco results
            .PLS_weight_scheme_inner = 'factorial',
            #.PLS_weight_scheme_inner = 'factorial',
            #  .resample_method = 'bootstrap',
            .tolerance = 1e-06
  )
  bt1 = out1$Estimates$Path_estimates[3,1:2]
  dlta_bt = (bt1 - bt)/dh
  Gbt[,cnter] = dlta_bt 
  wt1 = c(out1$Estimates$Loading_estimates[1,1:3],
         out1$Estimates$Loading_estimates[2,4:5],
         out1$Estimates$Loading_estimates[3,6:10])
  dlta_wt = (wt1 - wt)/dh
  Gwt[,cnter] = dlta_wt 
  lt1 = c(out1$Estimates$Weight_estimates[1,1:3],
         out1$Estimates$Weight_estimates[2,4:5],
         out1$Estimates$Weight_estimates[3,6:10])
  dlta_lt = (lt1 - lt)/dh
  Glt[,cnter] = dlta_lt 
  
    cnter = cnter + 1
 }
}


sqrt(diag( Gbt %*% vc_r_sim %*% t(Gbt) /(4 - 3) ) )
sqrt(diag( Gwt %*% vc_r_sim %*% t(Gwt) /(4 - 3) ) )
sqrt(diag( Glt %*% vc_r_sim %*% t(Glt) /(4 - 3) ) )

r_boot = matrix(0,nrow = 4999, ncol = 45)

for (i in 1:4999) {
  cor_sim = cor(Russett1[
    floor(runif(47,min = 1,max = 48)) ,])
  r_boot[i,] = c(cor_sim[1,2:10],cor_sim[2,3:10],
                 cor_sim[3,4:10],cor_sim[4,5:10],
                 cor_sim[5,6:10],cor_sim[6,7:10],
                 cor_sim[7,8:10],cor_sim[8,9:10],
                 cor_sim[9,10])
  
}
vc_r_boot = cov(r_boot)
rm(r_boot)

sqrt(diag( Gbt %*% vc_r_boot %*% t(Gbt) /(4 - 3) ) )



agr = Russett1[,c('gini','farm','rent')]
ind = Russett1[,c('gnpr','labo')]
pol = Russett1[,c('inst' , 'ecks' , 'deat' , 'stab' , 'dict' )]



