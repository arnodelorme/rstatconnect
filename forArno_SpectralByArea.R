#Code for Arno to run GAM on spectral data, for delta band only, by area
#Aug 31 2022

#Setup: 
setwd('/home/arno/nemar/trance_bids_2021/connectivity')


library(R.matlab)
library(dplyr)
library(tidyr)
library(mgcv)
library(gratia)
library(ggplot2)

source("utils.R")


# Load data
DATA <- read.csv('trance_area_power_fixed.csv')


#Update: Aug 16 2022
# Switch area number for xyz coordinates
xyz <- read.csv('area_coordinates.csv', header = F, col.names = c('area_x', 'area_y', 'area_z'))

DATA$area_x <- xyz$area_x[DATA$area]
DATA$area_y <- xyz$area_y[DATA$area]
DATA$area_z <- xyz$area_z[DATA$area]




# Separate DATA by band
DATA <- DATA %>% 
  mutate(band = case_when(
    frequency >= 1 & frequency <= 3 ~ "delta",
    frequency >= 4 & frequency <= 6 ~ "theta", 
    frequency >= 8 & frequency <= 12 ~ "alpha", 
    frequency >= 18 & frequency <= 22 ~ "beta", 
    frequency >= 30 & frequency <= 50 ~ "gamma"
  ))

bands_DATA <- DATA %>% 
  filter(band != "NA") %>% 
  group_by(subject, session, trial, area, area_x, area_y, area_z, condition, band) %>% summarise(power = mean(power))


# Set all as.factor
bands_DATA$subject <- as.factor(bands_DATA$subject)
bands_DATA$session <- as.factor(bands_DATA$session)
bands_DATA$trial <- as.factor(bands_DATA$trial)
bands_DATA$condition <- as.factor(bands_DATA$condition)
bands_DATA$band <- as.factor(bands_DATA$band)

options(bitmapType='cairo')

#Run models within band
#Updated Aug 16 2022: shifting area from numeric to xyz
model_delta <- bam(log(power) ~ condition + te(area_x, area_y, area_z, bs = "tp", by = condition, m = 2, k = c(15, 15, 15)) + s(subject, bs = "re") + s(session, bs = "re") + s(trial, bs = "re"), 
                   data = bands_DATA %>% filter(band == "delta"))
model_theta <- bam(log(power) ~ condition + te(area_x, area_y, area_z, bs = "tp", by = condition, m = 2, k = c(15, 15, 15)) + s(subject, bs = "re") + s(session, bs = "re") + s(trial, bs = "re"), 
                   data = bands_DATA %>% filter(band == "theta"))
model_alpha <- bam(log(power) ~ condition + te(area_x, area_y, area_z, bs = "tp", by = condition, m = 2, k = c(15, 15, 15)) + s(subject, bs = "re") + s(session, bs = "re") + s(trial, bs = "re"), 
                   data = bands_DATA %>% filter(band == "alpha"))
model_beta <- bam(log(power) ~ condition + te(area_x, area_y, area_z, bs = "tp", by = condition, m = 2, k = c(15, 15, 15)) + s(subject, bs = "re") + s(session, bs = "re") + s(trial, bs = "re"), 
                  data = bands_DATA %>% filter(band == "beta"))
model_gamma <- bam(log(power) ~ condition + te(area_x, area_y, area_z, bs = "tp", by = condition, m = 2, k = c(15, 15, 15)) + s(subject, bs = "re") + s(session, bs = "re") + s(trial, bs = "re"), 
                   data = bands_DATA %>% filter(band == "gamma"))


# Arno, if you want to print results into another folder, create that folder, and reassign the working direction. 
# setwd('<path to new working directory>')
save(model_delta, file = "model_delta.rda")
save(model_theta, file = "model_theta.rda")
save(model_alpha, file = "model_alpha.rda")
save(model_beta , file = "model_beta.rda")
save(model_gamma, file = "model_gamma.rda")

# DELTA BAND
summary(model_delta)
#png("DELTA_RESIDUAL_PLOT.png") 
par(mfrow=c(2,2))
gam.check(model_delta)
dev.off()

sink("gam_check_delta.txt")
gam.check(model_delta)
sink()

# # THETA BAND
# summary(model_theta)
# jpeg('THETA_RESIDUAL_PLOT.jpg')
# par(mfrow=c(2,2))
# gam.check(model_theta)
# dev.off()
# 
# # ALPHA BAND
# summary(model_alpha)
# jpeg('ALPHA_RESIDUAL_PLOT.jpg')
# par(mfrow=c(2,2))
# gam.check(model_alpha)
# dev.off()
# 
# # BETA BAND
# summary(model_beta)
# jpeg('BETA_RESIDUAL_PLOT.jpg')
# par(mfrow=c(2,2))
# gam.check(model_beta)
# dev.off()
# 
# # GAMMA BAND
# summary(model_gamma)
# jpeg('GAMMA_RESIDUAL_PLOT.jpg')
# par(mfrow=c(2,2))
# gam.check(model_gamma)
# dev.off()





# Compute and test difference scores.

pdat <- expand.grid(area = seq(min(bands_DATA$area), max(bands_DATA$area), 1),
                    condition = c("channeling", "mindwandering"),
                    subject = 1, session = 1, trial = 1) 


#Updated August 16 2022: Shift area to xyz
pdat$area_x <- xyz$area_x[pdat$area]
pdat$area_y <- xyz$area_y[pdat$area]
pdat$area_z <- xyz$area_z[pdat$area]

pdat <- pdat[, c("area_x", "area_y", "area_z", "condition", "subject", "session", "trial")]


delta_pairwise_table <- smooth_diff(model_delta, pdat, 'channeling', 'mindwandering', 'condition')
# theta_pairwise_table <- smooth_diff(model_theta, pdat, 'channeling', 'mindwandering', 'condition')
# alpha_pairwise_table <- smooth_diff(model_alpha, pdat, 'channeling', 'mindwandering', 'condition')
# beta_pairwise_table <- smooth_diff(model_beta, pdat, 'channeling', 'mindwandering', 'condition')
# gamma_pairwise_table <- smooth_diff(model_gamma, pdat, 'channeling', 'mindwandering', 'condition')


# alpha_pairwise_table <- cbind(unique(bands_DATA$area), xyz, alpha_pairwise_table)
# theta_pairwise_table <- cbind(unique(bands_DATA$area), xyz, theta_pairwise_table)
delta_pairwise_table <- cbind(unique(bands_DATA$area), xyz, delta_pairwise_table)
# beta_pairwise_table  <- cbind(unique(bands_DATA$area), xyz, beta_pairwise_table)
# gamma_pairwise_table <- cbind(unique(bands_DATA$area), xyz, gamma_pairwise_table)


# colnames(alpha_pairwise_table)[1] <- "area"
# colnames(theta_pairwise_table)[1] <- "area"
colnames(delta_pairwise_table)[1] <- "area"
# colnames(beta_pairwise_table)[1] <- "area"
# colnames(gamma_pairwise_table)[1] <- "area"


# alpha_pairwise_table$Significant <- ifelse(alpha_pairwise_table$fdr_adj_p < 0.05, "yes", "no")
# theta_pairwise_table$Significant <- ifelse(theta_pairwise_table$fdr_adj_p < 0.05, "yes", "no")
delta_pairwise_table$Significant <- ifelse(delta_pairwise_table$fdr_adj_p < 0.05, "yes", "no")
# beta_pairwise_table$Significant <- ifelse(beta_pairwise_table$fdr_adj_p < 0.05, "yes", "no")
# gamma_pairwise_table$Significant <- ifelse(gamma_pairwise_table$fdr_adj_p < 0.05, "yes", "no")


# write.csv(alpha_pairwise_table, "alpha_pairwise_table.csv", row.names = FALSE)
# write.csv(theta_pairwise_table, "theta_pairwise_table.csv", row.names = FALSE)
write.csv(delta_pairwise_table, "delta_pairwise_table.csv", row.names = FALSE)
# write.csv(beta_pairwise_table, "beta_pairwise_table.csv", row.names = FALSE)
# write.csv(gamma_pairwise_table, "gamma_pairwise_table.csv", row.names = FALSE)


# save(model_alpha, file="model_alpha.Rdata")
# save(model_theta, file="model_theta.Rdata")
save(model_delta, file="model_delta.Rdata")
# save(model_beta, file="model_beta.Rdata")
# save(model_gamma, file="model_gamma.Rdata")



# Plotting differences in each area

# In order to determine 1 fixed y-axis, first plot as is, then assess min and max for y across all plots. 
# Assign diff_ylim with these values, and re-plot with the "coord_cartesian()" code active.


#diff_ylim <- c(<min>, <max>) 

# png('Alpha diff plot.png', width=400, height=450)
# ggplot(alpha_pairwise_table, aes(x = area, y = diff, col = diff)) + 
#   theme_bw() + 
#   geom_point() +
#   facet_wrap(~ pair, ncol = 2) +
#   # coord_cartesian(ylim = diff_ylim) +
#   ggtitle("Alpha smooth differences (log scale)")
# dev.off()
# 
# png('Theta diff plot.png', width=400, height=450)
# ggplot(theta_pairwise_table, aes(x = area, y = diff, col = diff)) + 
#   theme_bw() + 
#   geom_point() +
#   facet_wrap(~ pair, ncol = 2) +
#   # coord_cartesian(ylim = diff_ylim) +
#   ggtitle("Theta smooth differences (log scale)")
# dev.off()

#png('Delta diff plot.png', width=400, height=450)
#ggplot(delta_pairwise_table, aes(x = area, y = diff, col = diff)) + 
#  theme_bw() + 
#  geom_point() +
#  facet_wrap(~ pair, ncol = 2) +
  # coord_cartesian(ylim = diff_ylim) +
#  ggtitle("Delta smooth differences (log scale)")
#dev.off()

# png('Beta diff plot.png', width=400, height=450)
# ggplot(beta_pairwise_table, aes(x = area, y = diff, col = diff)) + 
#   theme_bw() + 
#   geom_point() +
#   facet_wrap(~ pair, ncol = 2) +
#   # coord_cartesian(ylim = diff_ylim) +
#   ggtitle("Beta smooth differences (log scale)")
# dev.off()
# 
# png('Gamma diff plot.png', width=400, height=450)
# ggplot(gamma_pairwise_table, aes(x = area, y = diff, col = diff)) + 
#   theme_bw() + 
#   geom_point() +
#   facet_wrap(~ pair, ncol = 2) +
#   # coord_cartesian(ylim = diff_ylim) +
#   ggtitle("Gamma smooth differences (log scale)")
# dev.off()


