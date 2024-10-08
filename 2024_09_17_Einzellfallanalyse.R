# This code is to analyze some cases
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
iteration <- 2
set.seed(123)
sim_data_einzel <- generateData(model_base, 
                                .N= 5000, 
                                .empirical = TRUE, 
                                a = 0.02, 
                                b = c(0.5),
                                c = 0.01,
                                d = c(0.4),
                                e = c(0.5),
                                .return_type = "cor" )

sim_data_einzel <- generateData(model_base, 
                                .N= 5000, 
                                .empirical = TRUE, 
                                a = 0.5, 
                                b = 0.5,
                                c = 0.6,
                                d = 0.4,
                                e = 0.5,
                                .return_type = "cor" )

set.seed(50123141)
n = 100
data_sim_einzelanalyse <- MASS::mvrnorm(n = n, 
                          mu= rep(0, nrow(sim_data_einzel$dgp[[1]])), 
                          Sigma =  sim_data_einzel$dgp[[1]],
                          empirical = F)
set.seed(50123141)
res_einzelanalyse <- csem(.data = data_sim_einzelanalyse, 
                          .model = model_est,
                          .resample_method = 'bootstrap',
                          .R = 500,   
                          .PLS_weight_scheme_inner = 'factorial',
                          .tolerance = 1e-06)

infer_res <- infer(res_einzelanalyse, .alpha = 0.05, .quantity = "CI_percentile")
infer_res$Path_estimates$CI_percentile["95%L", "Y ~ X1"]

plot(density(res_einzelanalyse$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,'Y ~ X1']), main = "Empirical distribution of parameter estimates")

wt1 = c(res_einzelanalyse$Estimates$Weight_estimates[1,1:3],
        res_einzelanalyse$Estimates$Weight_estimates[2,4:6],
        res_einzelanalyse$Estimates$Weight_estimates[3,7:9])
lt1 = c(res_einzelanalyse$Estimates$Loading_estimates[1,1:3],
        res_einzelanalyse$Estimates$Loading_estimates[2,4:6],
        res_einzelanalyse$Estimates$Loading_estimates[3,7:9])
##############################################################################################################

overview_table$difference_sd <- overview_table$sd_delta_y_x1 - overview_table$sd_bootstrap_y_x1
overview_table$a_sq <- overview_table$a * overview_table$a
overview_table$c_sq <- 
lm_model <- lm(difference_sd ~ a + c + a:c + a*a + c*c, data = overview_table)

###############################################################################################################


o_table2 <- read.csv("C:/Users/jab49wd/iCloudDrive/Geteilt/Data/overview_table")


o_table_df <- as.data.frame(o_table2)



o_table_df$`Z_test_boot_y_x1.Y ~ X1`<- 1-o_table_df$`Z_test_boot_y_x1.Y ~ X1`
o_table_df$`Z_test_delta_y_x1.Y ~ X1` <- 1-o_table_df$`Z_test_delta_y_x1.Y ~ X1`

o_table_df_bildchen <- o_table_df[o_table_df$a != 0,]

o_table_df_bildchen <- o_table_df_bildchen %>% select(a, c, n, `path_estimate_y_x1.Y ~ X1`, sd_bootstrap_y_x1, `Z_test_boot_y_x1.Y ~ X1`, sd_delta_y_x1, `Z_test_delta_y_x1.Y ~ X1`)

o_table_df2 <- o_table_df %>% 
                group_by(a, c, n) %>%
                summarize(Rejection_rate_boot = mean(`Z_test_boot_y_x1.Y...X1`), 
                          Rejection_rate_delta = mean(`Z_test_delta_y_x1.Y...X1`))

o_table_df2 <- o_table_df2[o_table_df2$a != 0, ]
o_table_df2$delta_better <- as.factor(o_table_df2$Rejection_rate_delta > o_table_df2$Rejection_rate_boot)

tree_model_rpart <- rpart::rpart(delta_better ~ a + c + n, 
                    data = o_table_df2, 
                    method = "class" 
                    )

rpart.plot(tree_model_rpart)
printcp(tree_model_rpart) 
optimal_cp <- tree_model_rpart$cptable[which.min(tree_model_rpart$cptable[,"xerror"]),"CP"]
pruned_tree <- prune(tree_model_rpart, cp = optimal_cp)
rpart.plot(pruned_tree)
text(pruned_tree, use.n = TRUE) 

############################################################################################################

o_table_df2_1 <- o_table_df2
o_table_df2_2 <- o_table_df2


o_table_df2_1$treatment <- 0
o_table_df2_1$Rejection_rate <- o_table_df2_1$Rejection_rate_boot
o_table_df2_1$Rejection_rate_boot <- NULL
o_table_df2_1$Rejection_rate_delta <- NULL

o_table_df2_2$treatment <- 1
o_table_df2_2$Rejection_rate <- o_table_df2_2$Rejection_rate_delta
o_table_df2_2$Rejection_rate_boot <- NULL
o_table_df2_2$Rejection_rate_delta <- NULL

o_table_optimal_policy_learning <- rbind(o_table_df2_2, o_table_df2_1)




library(policytree)    # library for running a tree after double machine learning that provides optimal policies for covariate-determined subgroups
library(DiagrammeR)    # library for plotting 

y <- o_table_optimal_policy_learning$Rejection_rate
d <- o_table_optimal_policy_learning$treatment
x <- cbind(o_table_optimal_policy_learning$a, o_table_optimal_policy_learning$c, o_table_optimal_policy_learning$n)

d=factor(d)                                         # redefine treatment to be a factor
forest=multi_arm_causal_forest(X = x, Y = y, W = d) # estimate propensity scores and conditional means by random forests
scores=double_robust_scores(forest)                 # obtain doubly robust scores
xpol=cbind(a = o_table_optimal_policy_learning$a, c = o_table_optimal_policy_learning$c, n =o_table_optimal_policy_learning$n)   # define x variables based on which subgroups are to be defined
tree=policy_tree(X=xpol, Gamma=scores, depth = 3)   # define 8 subgroups based on 7 splits

# plot the tree defining the subgroups (in terms of X) and showing the optimal policies for each subgroup
plot(tree, leaf.labels = levels(d))


