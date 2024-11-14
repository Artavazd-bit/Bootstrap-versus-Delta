library(tidyverse)
library(tidyr)
library(ggplot2)
library(dplyr)

o_table <- readRDS("./01_Data/2024_11_06_sem_3_latent.rds")


o_table_2 <- o_table %>% 
  group_by(a, c, n) %>%
  summarize(Rejection_rate_boot = mean(`Z_test_boot_y_x1`), 
            Rejection_rate_delta = mean(`Z_test_delta_y_x1`),
            dip_test_mean = mean(dip_test_p_value))


################################################################################
o_table_lr <- readRDS("./01_Data/2024_11_06_linear_regression.rds")

o_table_lr_2 <- o_table_lr %>% 
  group_by(a, b, n) %>%
  summarize(Rejection_rate_boot = mean(`Z_test_boot_y_x1`), 
            Rejection_rate_delta = mean(`Z_test_delta_y_x1`),
            dip_test_mean = mean(dip_test_p_value)
            )

################################################################################
o_table_4_latent <- readRDS("./01_Data/2024_11_11_sem_4_latent.rds")

o_table_4_latent_2 <- o_table_4_latent %>% 
  group_by(a, b, c, d, e, n) %>%
  summarize(Rejection_rate_boot = mean(`Z_test_boot_x3_x2`), 
            Rejection_rate_delta = mean(`Z_test_delta_x3_x2`),
            dip_test_mean = mean(dip_test_p_value)
  )


################################################################################
ggplot(o_table, aes(x = path_estimate_y_x1, color = interaction(a, c, n))) + 
  geom_density() + 
  facet_wrap(~ interaction(a, c, n), scales = "free") + 
  labs(title = "Dichte der Werte von path_estimate_y_x1 für jede Kombination von a, c und n",
       x = "path_estimate_y_x1",
       y = "Dichte") +
  theme_minimal() +
  theme(legend.position = "none") 



###############################################################################
# Dichteplots
df_long <- o_table %>%
  unnest(col = c(list_with_bootstrap)) %>%
  rename(Bootstrap = list_with_bootstrap)  # Umbenennen für Klarheit

# Dichte-Plot erstellen
ggplot() +
  # Helle graue Schattierung für Bootstrap-Werte
  geom_density(data = df_long, aes(x = Bootstrap, group = Sim_run), color = "lightgray", alpha = 0.5) +
  # Dichte-Linie für path_estimate_y_x1.Y...X1 für jede Population
  geom_density(data = o_table, aes(x = path_estimate_y_x1, color = as.factor(a)), size = 1) +
  labs(title = "Dichteplot der Pfadschätzung mit Bootstrap-Verteilung",
       x = "Pfad-Koeffizient",
       y = "Dichte",
       color = "Population") +
  theme_minimal() +
  facet_wrap(a ~ n + c) 
