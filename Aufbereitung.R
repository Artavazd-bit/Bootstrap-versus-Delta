library(dplyr)
library(rpart)
library(rpart.plot)
overview_table <- read.csv("./Data/2024_10_14_overview_table.csv")

overview_table_2 <- as.data.frame(overview_table)

overview_table_2$Z_test_boot_y_x1.Y...X1<- 1-overview_table_2$Z_test_boot_y_x1.Y...X1
overview_table_2$Z_test_delta_y_x1.Y...X1 <- 1-overview_table_2$Z_test_delta_y_x1.Y...X1

overview_table_3 <- overview_table_2 %>% 
                                    group_by(a, c, n) %>%
                                    summarize(Rejection_rate_boot = mean(`Z_test_boot_y_x1.Y...X1`), 
                                              Rejection_rate_delta = mean(`Z_test_delta_y_x1.Y...X1`),
                                              Rejection_rate_boot_sd = sd(Z_test_boot_y_x1.Y...X1),
                                              Rejection_rate_delta_sd = sd(`Z_test_delta_y_x1.Y...X1`))

overview_table_3$delta_besser <- 0
overview_table_3$delta_besser[(overview_table_3$a !=0 & overview_table_3$Rejection_rate_boot< overview_table_3$Rejection_rate_delta) | (overview_table_3$a ==0 & overview_table_3$Rejection_rate_boot> overview_table_3$Rejection_rate_delta) ] <-1 

overview_table_3$delta_besser <- as.factor(overview_table_3$delta_besser)

overview_table_4 <- overview_table_3[overview_table_3$a ==0,]

tree_model_rpart <- rpart::rpart(delta_besser ~ a + c + n, 
                                 data = overview_table_3, 
                                 method = "class" 
)

rpart.plot(tree_model_rpart)
printcp(tree_model_rpart) 
optimal_cp <- tree_model_rpart$cptable[which.min(tree_model_rpart$cptable[,"xerror"]),"CP"]
pruned_tree <- prune(tree_model_rpart, cp = optimal_cp)
rpart.plot(pruned_tree)
text(pruned_tree, use.n = TRUE) 

overview_table_3$Rejection_rate_boot_sd <- NULL
overview_table_3$Rejection_rate_delta_sd <- NULL
