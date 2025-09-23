# Get coefficient for csi_past
data = mod_df_scaled
sem_model = best_model
b_csi_past <- coef(sem_model)["csi_past"]

# Sequence of csi_past values for plotting
csi_past_seq <- seq(min(mod_df_scaled$csi_past), max(mod_df_scaled$csi_past), length.out = 100)
    
# Predicted response holding other predictors at their mean
partial_mod_vars <- mod_vars[!mod_vars %in% c("csi_past", "Ho")]
partial_mod_vars[partial_mod_vars == "glacier_factor"] <- "glacier_factor1"  
data$glacier_factor1 <- data$glacier
other_effect <- as.matrix(data[, partial_mod_vars]) %*% coef(sem_model)[partial_mod_vars]
y_pred <- b_csi_past * csi_past_seq + mean(other_effect)

# Plot
library(ggplot2)
ggplot(data, aes(x = csi_past, y = residuals(sem_model) + fitted(sem_model))) +
  geom_point(alpha = 0.5) +
  geom_line(data = data.frame(csi_past = csi_past_seq, y = y_pred), aes(x = csi_past, y = y), color = "blue", size = 1) +
  labs(y = "Partial effect of csi_past on response")



pdf(here("DELETE.pdf"), height = 6, width = 8)

partial_resid <- residuals(sem_model) + coef(sem_model)["csi_past"] * data$csi_past
ggplot(data, aes(x = csi_past, y = partial_resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue") +
  labs(y = "Partial residuals (adjusted for spatial error)") +
  ggpubr::stat_cor()


partial_resid <- residuals(sem_model) + coef(sem_model)["tmean_dif"] * data$tmean_dif
ggplot(data, aes(x = tmean_dif, y = partial_resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue") +
  labs(y = "Partial residuals (adjusted for spatial error)") +
  ggpubr::stat_cor()

partial_resid <- residuals(sem_model) + coef(sem_model)["Q"] * data$Q
ggplot(data, aes(x = Q , y = partial_resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue") +
  labs(y = "Partial residuals (adjusted for spatial error)") +
  ggpubr::stat_cor()


partial_resid <- residuals(sem_model) + coef(sem_model)["bio1"] * data$bio1
ggplot(data, aes(x = bio1, y = partial_resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue") +
  labs(y = "Partial residuals (adjusted for spatial error)") +
  ggpubr::stat_cor()
dev.off()