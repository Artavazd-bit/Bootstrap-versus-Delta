sigma1 <- matrix(c(1, 0.3, 0.2, 0.3, 1, 0.4, 0.2, 0.4, 1), nrow = 3, ncol = 3)
test_1 <- MASS::mvrnorm(n = 100, mu = rep(0,3), Sigma = sigma1, empirical = TRUE)
cov(test_1)- sigma1


library(dplyr)

o_table <- read.csv(file = "C:/Users/jab49wd/iCloudDrive/Geteilt/overview.csv")

o_table_df <- as.data.frame(o_table)

o_table_df$Z_test_boot_y_x1.Y...X1 <- 1-o_table_df$Z_test_boot_y_x1.Y...X1
o_table_df$Z_test_delta_y_x1.Y...X1 <- 1-o_table_df$Z_test_delta_y_x1.Y...X1

o_table_df2 <- o_table_df %>% 
                        group_by(a, c, n) %>%
                        summarize(Rejection_rate_boot = mean(`Z_test_boot_y_x1.Y...X1`), 
                                  Rejection_rate_delta = mean(`Z_test_delta_y_x1.Y...X1`))

o_table_df2_a0 <- o_table_df2[o_table_df2$a == 0, ]

o_table_df2


