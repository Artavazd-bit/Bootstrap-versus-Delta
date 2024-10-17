library(dplyr)
library(rpart)
library(rpart.plot)
library(ggplot2)
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


overview_table_a_0_c_0.4_n_500 <- overview_table_2[overview_table_2$a == 0 & overview_table_2$c == 0.4 & overview_table_2$n == 500,]
plot(density(overview_table_a_0_c_0.4_n_500$path_estimate_y_x1.Y...X1))

overview_table_a_0_c_0_n_500 <- overview_table_2[overview_table_2$a == 0 & overview_table_2$c == 0 & overview_table_2$n == 500,]
plot(density(overview_table_a_0_c_0_n_500$path_estimate_y_x1.Y...X1))

data_for_plot <- overview_table_a_0_c_0.4_n_500
plot(data_for_plot$path_estimate_y_x1.Y...X1, data_for_plot$sd_bootstrap_y_x1, type = "p", col = "blue")
points(data_for_plot$path_estimate_y_x1.Y...X1, data_for_plot$sd_delta_y_x1, col="red")

overview_table_5 <- overview_table_2
overview_table_5$sd_difference  <- overview_table_5$sd_delta_y_x1 - overview_table_5$sd_bootstrap_y_x1


ggplot(overview_table_5, aes(x = path_estimate_y_x1.Y...X1, y = sd_difference, color = a)) +
  geom_point(size = 3) +  # Größere Punkte für bessere Sichtbarkeit
  scale_color_gradient(low = "red", high = "blue") +  # Farbverlauf von rot (kleine Werte) bis blau (große Werte)
  labs(
    x = "Pfad-Schätzung (Y...X1)", 
    y = "Differenz: SD Delta vs SD Bootstrap",
    color = "Wert von a"
  ) +
  theme_minimal()



ggplot(overview_table_5, aes(x = path_estimate_y_x1.Y...X1, y = sd_difference, color = c)) +
  geom_point(size = 3) +  # Größere Punkte für bessere Sichtbarkeit
  scale_color_gradient(low = "red", high = "blue") +  # Farbverlauf von rot (kleine Werte) bis blau (große Werte)
  labs(
    x = "Pfad-Schätzung (Y...X1)", 
    y = "Differenz: SD Delta vs SD Bootstrap",
    color = "Wert von c"
  ) +
  theme_minimal()

overview_table_a_0_c_0.1_n_50 <- overview_table_2[overview_table_2$a == 0 & overview_table_2$c == 0.1 & overview_table_2$n == 50,]
plot(density(overview_table_a_0_c_0.1_n_50$path_estimate_y_x1.Y...X1))
lines(x = seq(-1, 1, length = 100), y = dnorm(x = seq(-1, 1, length = 100),mean= 0.00363605, sd = 0.2371697), col = "blue")

qnorm(p = c(0.025, 0.975),mean= 0.00363605, sd = 0.2371697)

quantile(overview_table_a_0_c_0.1_n_50$path_estimate_y_x1.Y...X1, probs = c(0.025, 0.975))
mean(overview_table_a_0_c_0.1_n_50$path_estimate_y_x1.Y...X1)



