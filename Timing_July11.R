library(RobustGaSP)
library(dplyr)
library(ggplot2)
library(colorspace)
library(patchwork)
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
ppgp <- ppgasp(X_train[1:n_train,], log(G_entries_train[1:n_train,]),
               isotropic = F, nugget.est = T, method = 'mle', optimization = 'nelder-mead')


# FORM FACTORS ONLY

set.seed(1)
n_non_truncated_compare <- 8500
N_charge_compare <- 1+199*runif(n_non_truncated_compare)
N_uncharge_compare <- 1+199*runif(n_non_truncated_compare)
X_compare <- round(cbind(N_charge_compare, N_uncharge_compare))
X_compare <- X_compare[which(X_compare[,1] + X_compare[,2] <= 200 & X_compare[,1] + X_compare[,2] >= 50),]

report_interval_form_factor <- 50
total_points_form_factor <- dim(X_compare)[1]
num_loops_form_factor <- floor(total_points_form_factor / report_interval_form_factor)

all_times_form_factor <- rep(NA, num_loops_form_factor)
sim_times_form_factor <- pred_times_form_factor <- rep(NA, num_loops_form_factor)

for(i in 1:num_loops_form_factor){
  
  idx_here <- ((i-1)*report_interval_form_factor + 1):(i*report_interval_form_factor)
  print(c(min(idx_here), max(idx_here)))
  
  loop_sim_time_form_factor <- system.time(
    for(j in idx_here) G <- calculateG(X_compare[j,1], X_compare[j,2])
  )
  
  loop_pred_time_form_factor <- system.time(
    for(j in 1:1){
      
      entries <- exp(predict(ppgp, X_compare[1:(i*report_interval_form_factor),])$mean)
      det_inv <- 1 / (entries[,1:200]*entries[,401:600] - entries[,201:400]^2)
      
    }
  )
  
  print('Simulation Time:')
  print(loop_sim_time_form_factor)
  print('Prediction Time:')
  print(loop_pred_time_form_factor)
  
  if(i == 1){
    sim_times_form_factor[i] <- loop_sim_time_form_factor[[3]]
    pred_times_form_factor[i] <- loop_pred_time_form_factor[[3]]
  } else {
    sim_times_form_factor[i] <- sim_times_form_factor[i-1] + loop_sim_time_form_factor[[3]]
    pred_times_form_factor[i] <- loop_pred_time_form_factor[[3]]
  }
  
}


# FULL SIMULATION

lb <- c(0.5,0.02,0.02,0.01,0.01,0.0,0.0,0.0,0.0,1.0,0.0,0.1,0.1)
ub <- c(2.0,0.98,0.98,0.99,0.99,1.0,100.0,100.0,2.0,5.0,10.0,10.0,10.0)

report_interval_full <- 50
num_loops_full <- 80

sim_times_full <- pred_times_full <- rep(NA, num_loops_full)
set.seed(1)
X_full <- matrix(runif(13*report_interval_full*num_loops_full), ncol=13)
for(i in 1:13) X_full[,i] <- X_full[,i] * (ub[i]-lb[i]) + lb[i]
# X_full <- as.data.frame(X_full)
# colnames(X_full) <- c('r','alpha_A','alpha_B','f_A','f_p','rho_s','chiABN','chiASN','chiRatio',
#                       'epsRatio','lb','a_plus','a_minus')


for(i in 1:num_loops_full){
  
  idx_here <- ((i-1)*report_interval_full + 1):(i*report_interval_full)
  
  sim_time_loop_full_here <- system.time(
    for(j in 1:1){
      out <- full_sim_process_w_fast(X_full[idx_here,])
    }
  )
  
  pred_time_loop_full_here <- system.time(
    for(j in 1:1){
      out_pred <- predict_phase_w_correction(X_full[1:(i*report_interval_full),], ppgp)
    }
  )
  
  print('Simulation Time:')
  print(sim_time_loop_full_here)
  print('Prediction Time:')
  print(pred_time_loop_full_here)
  
  if(i == 1){
    sim_times_full[i] <- sim_time_loop_full_here[[3]]
    pred_times_full[i] <- pred_time_loop_full_here[[3]]
  } else {
    sim_times_full[i] <- sim_time_loop_full_here[[3]] + sim_times_full[i-1]
    pred_times_full[i] <- pred_time_loop_full_here[[3]]
  }
  
}

time_df <- data.frame(
  n = seq(50, 4000, by = 50),
  G_sim = sim_times_form_factor,
  G_pred = pred_times_form_factor,
  All_sim = sim_times_full,
  All_pred = pred_times_full
)

melt_time_df <- melt(time_df, id.vars = 'n')

melt_time_df$setting <- ifelse(melt_time_df$variable == 'All_pred', 'Phase, Prediction',
                               ifelse(melt_time_df$variable == 'All_sim', 'Phase, Truth',
                                      ifelse(melt_time_df$variable == 'G_pred', 'Form Factor, Prediction', 'Form Factor, Truth')))


options(scipen=999)
time_plt <- ggplot(melt_time_df, aes(x = n, y = value, color = setting, linetype = setting)) +
  geom_line() + 
  # geom_point(size=.3) +
  #scale_color_manual(values = c('Truth' = 'red', 'Prediction' = 'blue')) +
  scale_color_manual(values = c('Form Factor, Prediction' = '#4f9da6', 'Form Factor, Truth' = 'black',
                                'Phase, Prediction' = '#4f9da6', 'Phase, Truth' = 'black')) +
  scale_linetype_manual(values = c('Form Factor, Prediction' = 1, 'Form Factor, Truth' = 1,
                                   'Phase, Prediction' = 2, 'Phase, Truth' = 2)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = 'Number of Generated Samples', y = 'Time (s)', color = 'Setting', linetype = 'Setting') +
  theme_minimal()
time_plt

ggsave('TimeComparison.png', plot = time_plt, width = 7, height = 2.75, dpi = 320)

# save.image('Timing_July11.Rdata')

#load('Timing_July11.Rdata')







