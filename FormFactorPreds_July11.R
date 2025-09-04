library(RobustGaSP)
library(ggplot2)
library(dplyr)
library(lhs)
library(plot3D)
library(ggplot2)
library(colorspace)
library(patchwork)
library(gganimate)
library(reshape2)

# load('Datasets_Jun2.Rdata')

source('Phase_Functions.R')


# LOAD ALL DATASETS ###########
X_train <- as.matrix(read.csv('All_Datasets/X_train.csv'))
G_entries_train <- as.matrix(read.csv('All_Datasets/G_entries_train.csv'))

X_test_2d <- as.matrix(read.csv('All_Datasets/X_test_2d.csv'))
G_entries_test_2d <- as.matrix(read.csv('All_Datasets/G_entries_test_2d.csv'))

X_test <- as.matrix(read.csv('All_Datasets/X_test.csv'))
y_test <- as.vector(read.csv('All_Datasets/y_test.csv'))$x

X_test_grid <- read.csv('All_Datasets/X_test_grid.csv')
phase_record_grid <- read.csv('All_Datasets/phase_record_grid.csv')

X_train_paired <- read.csv('All_Datasets/X_train_paired.csv')
y_train_paired <- as.vector(read.csv('All_Datasets/y_train_paired.csv'))$x

X_train_nonpaired <- read.csv('All_Datasets/X_train_nonpaired.csv')
y_train_nonpaired <- as.vector(read.csv('All_Datasets/y_train_nonpaired.csv'))$x

all_X_train <- array(0, dim = c(500, 2, 10))
all_X_train_nonpaired <- array(0, dim = c(500, 13, 10))
all_G_entries_train <- array(0, dim = c(500, 600, 10))
all_phases <- matrix(0, nrow = 500, ncol = 10)

for(i in 1:10){
  filename_for_2d_training_data <- paste0('All_Datasets/all_X_train_',i,'.csv')
  filename_for_13d_training_data <- paste0('All_Datasets/all_X_train_nonpaired_',i,'.csv')
  filename_for_2d_training_data_response <- paste0('All_Datasets/all_G_entries_train_',i,'.csv')
  filename_for_13d_training_data_response <- paste0('All_Datasets/all_phases_',i,'.csv')
  
  all_X_train[,,i] <- as.matrix(read.csv(filename_for_2d_training_data))
  all_X_train_nonpaired[,,i] <- as.matrix(read.csv(filename_for_13d_training_data))
  all_G_entries_train[,,i] <- as.matrix(read.csv(filename_for_2d_training_data_response))
  all_phases[,i] <- as.vector(read.csv(filename_for_13d_training_data_response))$x
  
}


n_train <- 50
n_train_big <- 100
ppgp <- ppgasp(X_train[1:n_train,], log(G_entries_train[1:n_train,]),
               isotropic = F, nugget.est = T, optimization = 'nelder-mead', method = 'mle')
ppgp_big <- ppgasp(X_train[1:n_train_big,], log(G_entries_train[1:n_train_big,]),
                   isotropic = F, nugget.est = T, optimization = 'nelder-mead', method = 'mle')
preds <- predict(ppgp, X_test_2d)
preds_big <- predict(ppgp_big, X_test_2d)

log_G11_df <- data.frame(
  Nc = X_test_2d[,1],
  Nu = X_test_2d[,2],
  G11_true = log(G_entries_test_2d[,1]),
  G11_pred = preds$mean[,1],
  G11_pred_big = preds_big$mean[,1],
  G11_CI = preds$upper95[,1] - preds$lower95[,1],
  G11_CI_big = preds_big$upper95[,1] - preds_big$lower95[,1]
)


log_G12_df <- data.frame(
  Nc = X_test_2d[,1],
  Nu = X_test_2d[,2],
  G12_true = log(G_entries_test_2d[,201]),
  G12_pred = preds$mean[,201],
  G12_pred_big = preds_big$mean[,201],
  G12_CI = preds$upper95[,201] - preds$lower95[,201],
  G12_CI_big = preds_big$upper95[,201] - preds$lower95[,201]
)

log_G22_df <- data.frame(
  Nc = X_test_2d[,1],
  Nu = X_test_2d[,2],
  G22_true = log(G_entries_test_2d[,401]),
  G22_pred = preds$mean[,401],
  G22_pred_big = preds_big$mean[,401],
  G22_CI = preds$upper95[,401] - preds$lower95[,401],
  G22_CI_big = preds_big$upper95[,401] - preds_big$lower95[,401]
)

near_white <- '#EAF0F3'
high_val <- '#74B0E4'
low_val <- '#F18344'

# clim_G11 <- c(min(c(log_G11_df$G11_true, log_G11_df$G11_pred, log_G11_df$G11_pred_big)), max(c(log_G11_df$G11_true, log_G11_df$G11_pred, log_G11_df$G11_pred_big)))
clim_G11 <- c(-2,11)
clim_G11_CI <- c(min(c(log_G11_df$G11_CI,log_G11_df$G11_CI_big)),max(c(log_G11_df$G11_CI,log_G11_df$G11_CI_big)))

G11_true_plt <- ggplot(log_G11_df) + 
  geom_raster(aes(Nc, Nu, fill = G11_true)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  # scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G11) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G11, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G11_pred_plt <- ggplot(log_G11_df) + 
  geom_raster(aes(Nc, Nu, fill = G11_pred)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G11) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G11, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G11_Residual_plt <- ggplot(log_G11_df) + 
  geom_raster(aes(Nc, Nu, fill = G11_true - G11_pred)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G11) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G11, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G11_pred_plt_big <- ggplot(log_G11_df) +
  geom_raster(aes(Nc, Nu, fill = G11_pred_big)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G11) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G11, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))


G11_CI_plt <- ggplot(log_G11_df) +
  geom_raster(aes(Nc, Nu, fill = G11_CI)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G11, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G11_CI_plt_big <- ggplot(log_G11_df) +
  geom_raster(aes(Nc, Nu, fill = G11_CI_big)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G11, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))






#clim_G12 <- c(min(c(log_G12_df$G12_true, log_G12_df$G12_pred, log_G12_df$G12_pred_big)), max(c(log_G12_df$G12_true, log_G12_df$G12_pred, log_G12_df$G12_pred_big)))
clim_G12 <- c(-2,11)
clim_G12_CI <- c(min(c(log_G12_df$G12_CI,log_G12_df$G12_CI_big)),max(c(log_G12_df$G12_CI,log_G12_df$G12_CI_big)))

G12_true_plt <- ggplot(log_G12_df) + 
  geom_raster(aes(Nc, Nu, fill = G12_true)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G12) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G12, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G12_pred_plt <- ggplot(log_G12_df) + 
  geom_raster(aes(Nc, Nu, fill = G12_pred)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G12) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G12, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G12_Residual_plt <- ggplot(log_G12_df) + 
  geom_raster(aes(Nc, Nu, fill = G12_true - G12_pred)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G11) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G12, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G12_pred_plt_big <- ggplot(log_G12_df) + 
  geom_raster(aes(Nc, Nu, fill = G12_pred_big)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G12) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G12, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G12_CI_plt <- ggplot(log_G12_df) + 
  geom_raster(aes(Nc, Nu, fill = G12_CI)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G12, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G12_CI_plt_big <- ggplot(log_G12_df) + 
  geom_raster(aes(Nc, Nu, fill = G12_CI_big)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  scale_fill_gradient2(name='',low = '#7FE6F0',
                       mid = '#74B0E4', high = '#005A91', 
                       midpoint = 0.5*sum(clim_G12_CI), limits = clim_G12_CI) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

clim_G22 <- c(-2,11)
clim_G22_CI <- c(min(c(log_G12_df$G12_CI,log_G12_df$G12_CI_big)),max(c(log_G12_df$G12_CI,log_G12_df$G12_CI_big)))

G22_true_plt <- ggplot(log_G22_df) + 
  geom_raster(aes(Nc, Nu, fill = G22_true)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G12) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G22, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G22_pred_plt <- ggplot(log_G22_df) + 
  geom_raster(aes(Nc, Nu, fill = G22_pred)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G12) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G22, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G22_Residual_plt <- ggplot(log_G22_df) + 
  geom_raster(aes(Nc, Nu, fill = G22_true - G22_pred)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G11) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G22, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G22_pred_plt_big <- ggplot(log_G22_df) + 
  geom_raster(aes(Nc, Nu, fill = G22_pred_big)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  #scale_fill_distiller(palette = 'RdYlBu', name = "", limits = clim_G12) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G22, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G22_CI_plt <- ggplot(log_G22_df) + 
  geom_raster(aes(Nc, Nu, fill = G22_CI)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  scale_fill_gradient2(name = '', low=low_val, mid = near_white, high = high_val,
                       limits = clim_G12, midpoint = 0) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))

G22_CI_plt_big <- ggplot(log_G22_df) + 
  geom_raster(aes(Nc, Nu, fill = G22_CI_big)) +
  labs(x = expression(N[c]), y = expression(N[u])) +
  scale_fill_gradient2(name='',low = '#7FE6F0',
                       mid = '#74B0E4', high = '#005A91', 
                       midpoint = 0.5*sum(clim_G22_CI), limits = clim_G22_CI) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA))


ggsave('FF_plots_July/true_log_G11.png', plot = G11_true_plt, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/pred_log_G11.png', plot = G11_pred_plt, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/pred_log_G11_big.png', plot = G11_pred_plt_big, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/CI_log_G11.png', plot = G11_CI_plt, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/CI_log_G11_big.png', plot = G11_CI_plt_big, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/log_G11_residuals.png', plot = G11_Residual_plt, width = 3.75, height = 2.75, dpi = 300)

ggsave('FF_plots_July/true_log_G12.png', plot = G12_true_plt, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/pred_log_G12.png', plot = G12_pred_plt, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/pred_log_G12_big.png', plot = G12_pred_plt_big, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/CI_log_G12.png', plot = G12_CI_plt, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/CI_log_G12_big.png', plot = G12_CI_plt_big, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/log_G12_residuals.png', plot = G12_Residual_plt, width = 3.75, height = 2.75, dpi = 300)

ggsave('FF_plots_July/true_log_G22.png', plot = G22_true_plt, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/pred_log_G22.png', plot = G22_pred_plt, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/pred_log_G22_big.png', plot = G22_pred_plt_big, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/CI_log_G22.png', plot = G22_CI_plt, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/CI_log_G22_big.png', plot = G22_CI_plt_big, width = 3.75, height = 2.75, dpi = 300)
ggsave('FF_plots_July/log_G22_residuals.png', plot = G22_Residual_plt, width = 3.75, height = 2.75, dpi = 300)



NRMSE <- function(preds, truth, g_1 = 0, g_2 = 0){
  if(g_1 == 1 & g_2 == 1){
    idx <- 1:200
  } else if((g_1 == 1 & g_2 == 2) | (g_1 == 2 & g_1 == 1)){
    idx <- 201:400
  } else if(g_1 == 2 & g_2 == 2){
    idx <- 401:600
  } else {
    idx <- 1:600
  }
  
  
  # Numerator
  NRMSE_Numerator <- sqrt(mean((truth[,idx] - exp(preds$mean[,idx]))^2))
  # Denominator
  sample_means <- apply(truth[,idx], 2, mean)
  NRMSE_Denominator <- sd(truth[,idx])
  
  return(NRMSE_Numerator/NRMSE_Denominator)
  
}

EmpiricalCoverage <- function(preds, truth, idx, x_idx = NULL){
  if(is.null(x_idx)) x_idx <- 1:nrow(truth)
  return(
    mean(exp(preds$lower95[x_idx,idx]) < truth[x_idx,idx] & exp(preds$upper95[x_idx,idx]) > truth[x_idx,idx])
  )
}

# n = 50
nrmse_11 <- NRMSE(preds, G_entries_test_2d, 1, 1)
nrmse_12 <- NRMSE(preds, G_entries_test_2d, 1, 2)
nrmse_22 <- NRMSE(preds, G_entries_test_2d, 2, 2)
nrmse_all <- NRMSE(preds, G_entries_test_2d)

nrmse_11 # .020
nrmse_12 # .068
nrmse_22 # .006
nrmse_all # .022

# n = 100
nrmse_11_big <- NRMSE(preds_big, G_entries_test_2d, 1, 1)
nrmse_12_big <- NRMSE(preds_big, G_entries_test_2d, 1, 2)
nrmse_22_big <- NRMSE(preds_big, G_entries_test_2d, 2, 2)
nrmse_all_big <- NRMSE(preds_big, G_entries_test_2d)

nrmse_11_big # .007
nrmse_12_big # .028
nrmse_22_big # .002
nrmse_all_big # .008

# Coverage_G11 <- EmpiricalCoverage(preds, G_entries_test_2d, 1:200)
# Coverage_G12 <- EmpiricalCoverage(preds, G_entries_test_2d, 201:400)
# Coverage_G22 <- EmpiricalCoverage(preds, G_entries_test_2d, 401:600)
# Coverage_Overall <- EmpiricalCoverage(preds, G_entries_test_2d, 1:600)
# 
# Coverage_G11 # 0.830
# Coverage_G12 # 0.735
# Coverage_G22 # 0.820
# Coverage_Overall # 0.795
# 
# # Coverage without boundary
# no_bdry_idx <- which(X_test_2d[,1] > 10 & X_test_2d[,2] > 10)
# Coverage_G11_no_bdry <- EmpiricalCoverage(preds, G_entries_test_2d, 1:200, no_bdry_idx)
# Coverage_G12_no_bdry <- EmpiricalCoverage(preds, G_entries_test_2d, 201:400, no_bdry_idx)
# Coverage_G22_no_bdry <- EmpiricalCoverage(preds, G_entries_test_2d, 401:600, no_bdry_idx)
# Coverage_Overall_no_bdry <- EmpiricalCoverage(preds, G_entries_test_2d, 1:600, no_bdry_idx)
# 
# Coverage_G11_no_bdry # 0.88
# Coverage_G12_no_bdry # 0.848
# Coverage_G22_no_bdry # 0.876
# Coverage_Overall_no_bdry # 0.870
# 
# coverages <- rep(NA, nrow(X_test_2d))
# for(i_test in 1:nrow(X_test_2d)) coverages[i_test] <- mean(exp(preds$lower95[i_test,1:200]) < truth[i_test,1:200] & exp(preds$upper95[i_test,1:200]) > truth[i_test,1:200])
# bad_idx <- sort(coverages, index.return = T)$ix
# 
# 
# i_test <- bad_idx[1]
# plot(G_entries_test_2d[i_test,1:200],type='l')
# lines(exp(preds$mean[i_test,1:200]),col='blue')
# lines(exp(preds$upper95[i_test,1:200]),col='blue',lty=2)
# lines(exp(preds$lower95[i_test,1:200]),col='blue',lty=2)
# print(X_test_2d[i_test,])
# print(i_test)
# i_test <- i_test + 1
# 
# X_test_2d[bad_idx[1:100],]


# NRMSE v.s. n
n_vec <- seq(20, 250, by = 10)
nrmse_11_record <- nrmse_12_record <- nrmse_22_record <- rep(NA, length(n_vec))
for(i in 1:length(n_vec)){
  n_train <- n_vec[i]

  ppgp <- ppgasp(X_train[1:n_train,], log(G_entries_train[1:n_train,]),
                 isotropic = F, nugget.est = T, optimization = 'nelder-mead', method = 'mle')
  preds <- predict(ppgp, X_test_2d)
  nrmse_11 <- NRMSE(preds, G_entries_test_2d, 1, 1)
  nrmse_12 <- NRMSE(preds, G_entries_test_2d, 1, 2)
  nrmse_22 <- NRMSE(preds, G_entries_test_2d, 2, 2)
  print(n_train)
  print(nrmse_11)
  print(nrmse_12)
  print(nrmse_22)
  nrmse_11_record[i] <- nrmse_11
  nrmse_12_record[i] <- nrmse_12
  nrmse_22_record[i] <- nrmse_22
}


nrmse_df <- data.frame(
  n = n_vec,
  G11 = nrmse_11_record,
  G12 = nrmse_12_record,
  G22 = nrmse_22_record
) %>% melt(id.vars = c('n'))
label_vec <- c(expression(G[11]),expression(G[12]),expression(G[22]))
nrmse_plt <- ggplot(nrmse_df, aes(x = n, y = value, color = variable)) +
  geom_point(aes(shape = variable), size = 3) +
  geom_line(aes(linetype = variable), linewidth = .8) +
  scale_color_manual(values = c('G11' = '#FE938C', 'G12' = '#76C173', 'G22' = '#4281A4'),
                     labels = label_vec) +
  scale_shape_manual(values = 21:23, labels = label_vec) +
  scale_linetype_manual(values = c(1,2,4), labels = label_vec) +
  labs(x = expression(n), y = 'NRMSE',
       color = 'Entry', shape = 'Entry', shape = 'Entry', linetype = 'Entry') +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.key.spacing.y = unit(.01, 'cm'))

nrmse_plt
ggsave(filename = 'figures/nrmse_plt.png', plot = nrmse_plt, width = 7, height = 3.25, dpi = 300)


# save.image('FormFactorPreds_July11.Rdata')
