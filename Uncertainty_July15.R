library(RobustGaSP)
library(reshape2)
library(ggplot2)
library(patchwork)

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


phase_with_raw_predictions <- function(ppgp_obj, X_test){
  
  all_U_and_pre <- gen_all_U_and_prefactors_vec(X_test, nrow(X_test))
  
  N_A <- as.numeric(100)
  X_test_A <- cbind(round(N_A * X_test[,2]), round(N_A - N_A * X_test[,2]))
  
  N_B <- round(N_A * X_test[,1])
  X_test_B <- cbind(round(N_B * X_test[,3]), round(N_B - N_B * X_test[,3]))
  
  preds_A <- exp(predict(ppgp_obj, X_test_A)$mean)
  preds_B <- exp(predict(ppgp_obj, X_test_B)$mean)
  
  out <- phase_from_G_U_pre_fast(preds_A[,1:200], preds_A[,201:400], preds_A[,401:600],
                                 preds_B[,1:200], preds_B[,201:400], preds_B[,401:600],
                                 all_U_and_pre$all_U, all_U_and_pre$all_pre)
  
  return(out$phase_record)
  
}

uncertainty_at_test_point <- function(ppgp_obj, x_test, sample_size = 100){
  
  # ppgp_obj is a trained ppgasp model
  # x_test is a 13 length vector of parameters
  
  U_and_pre <- get_U_and_prefactors_fast(x_test[1], x_test[2], x_test[3],
                                         x_test[4], x_test[5], x_test[6],
                                         x_test[7], x_test[8], x_test[9],
                                         x_test[10], x_test[11], x_test[12],
                                         x_test[13])
  
  N_A <- as.numeric(100)
  N_Ac <- round(x_test[2]*N_A)
  N_Au <- N_A-N_Ac
  x_test_A <- matrix(c(N_Ac,N_Au),nrow=1)
  
  N_B <- round(x_test[1]*N_A)
  N_Bc <- round(N_B*x_test[3])
  N_Bu <- N_B-N_Bc
  x_test_B <- matrix(c(N_Bc,N_Bu),nrow=1)
  
  preds_A <- predict(ppgp_obj, x_test_A)
  preds_B <- predict(ppgp_obj, x_test_B)
  
  samples <- get_det_S_inv_samples(
    log_G11_A = preds_A$mean[,1:200],
    log_G12_A = preds_A$mean[,201:400],
    log_G22_A = preds_A$mean[,401:600],
    log_G11_B = preds_B$mean[,1:200],
    log_G12_B = preds_B$mean[,201:400],
    log_G22_B = preds_B$mean[,401:600],
    log_sd_G11_A = preds_A$sd[,1:200],
    log_sd_G12_A = preds_A$sd[,201:400],
    log_sd_G22_A = preds_A$sd[,401:600],
    log_sd_G11_B = preds_B$sd[,1:200],
    log_sd_G12_B = preds_B$sd[,201:400],
    log_sd_G22_B = preds_B$sd[,401:600],
    U = U_and_pre$U_mat,
    prefactors = U_and_pre$prefactors,
    sample_size = sample_size
  )
  
  return(c(
    sum(samples$phase_sample_record == 0) / sum(samples$phase_sample_record != -1), # Proportion of valid predictions that are homo
    sum(samples$phase_sample_record == 1) / sum(samples$phase_sample_record != -1), # ^ that are macro
    sum(samples$phase_sample_record == 2) / sum(samples$phase_sample_record != -1), # ^^ that are micro
    sum(samples$phase_sample_record == -1) # Number invalid for our recordkeeping
  ))
  
}

check_unique_valids <- function(row){
  return(length(unique(row[which(row != -1)])))
}

n_vec <- seq(20, 40, by = 1)
# small_det_correct_num <- # Corrected, accurate predictions
#   small_det_incorrect_num <- # Corrected, inaccurate predictions
#   small_det_uncertain_num <- # Corrected, uncertain predictions
#   small_det_unpredictable_num <-  # Corrected, unpredictable samples
#   large_det_correct_num <- # Uncorrected, accurate predictions
#   large_det_incorrect_num <- # Uncorrected, inaccurate predictions
#   rep(NA, length(n_vec))



sample_size <- 100
set.seed(1)

uncertainty_array <- array(0, dim = c(length(n_vec), nrow(X_test), 8))
dimnames(uncertainty_array) <- list(paste0('n=',n_vec), NULL, 
                                    c('Large_Det_Correct','Large_Det_Incorrect',
                                      'Small_Det_Correct','Small_Det_Incorrect',
                                      'Small_Det_Unpredictable','Small_Det_Uncertain',
                                      'Raw_Correct', 'Raw_Incorrect'))


for(i in 1:length(n_vec)){
  
  n_train <- n_vec[i]
  print(n_train)
  
  ppgp <- ppgasp(design = X_train[1:n_train,], response = log(G_entries_train[1:n_train,]),
                 method = 'mle', optimization = 'nelder-mead', isotropic = F, nugget.est = T)
  
  preds <- predict_phase_w_correction(X_test, ppgp, correction_sample_size = sample_size)
  
  need_correction_idx <- which(preds$pred_phase == -1)
  no_correction_idx <- (1:2000)[-need_correction_idx]
  
  # Uncorrected, accurate predictions
  large_det_correct_idx <- which(preds$pred_phase == y_test & preds$pred_phase != -1)
  uncertainty_array[i,large_det_correct_idx,1] <- 1
  
  # Uncorrected, inaccurate predictions
  large_det_incorrect_idx <- which(preds$pred_phase != y_test & preds$pred_phase != -1)
  uncertainty_array[i,large_det_incorrect_idx,2] <- 1
  
  # Corrected, accurate predictions
  small_det_correct_idx <- which(preds$pred_phase_corrected == y_test & preds$pred_phase == -1)
  uncertainty_array[i,small_det_correct_idx,3] <- 1
  
  # Corrected, inaccurate predictions (But still output the phase)
  small_det_incorrect_idx <- which(preds$pred_phase_corrected != y_test & preds$pred_phase == -1 & preds$pred_phase_corrected != -1)
  uncertainty_array[i,small_det_incorrect_idx,4] <- 1
  
  # Corrected, unpredictable predictions
  small_det_unpredictable_idx <- which(preds$pred_phase_corrected == -1)
  uncertainty_array[i,small_det_unpredictable_idx,5] <- 1
  
  # The next category overlaps with above columns 3, 4 and only exists for record keeping
  
  # Corrected, uncertain predictions (Output more than one sampled phase)
  sampled_phases <- preds$pred_phase_sample_record
  num_unique_samples <- apply(sampled_phases, 1, check_unique_valids)
  small_det_uncertain_idx <- which(num_unique_samples > 1)
  uncertainty_array[i,small_det_uncertain_idx,6] <- 1
  
  #### RAW PREDICTIONS ####
  raw_phase_preds <- phase_with_raw_predictions(ppgp, X_test)
  raw_phase_correct_idx <- which(raw_phase_preds == y_test)
  raw_phase_incorrect_idx <- which(raw_phase_preds != y_test)
  
  uncertainty_array[i,raw_phase_correct_idx,7] <- 1
  uncertainty_array[i,raw_phase_incorrect_idx,8] <- 1
  
  
  print((length(small_det_correct_idx) + length(large_det_correct_idx)) / 2000)
  print(length(raw_phase_correct_idx)/2000)
  
}

# num_record <- matrix(NA, length(n_vec), 6)
# for(i in 1:length(n_vec)){
#   num_record[i,] <- apply(!is.na(uncertainty_array[i,,]),2,sum)
# }
#colnames(num_record) <- c('Large_Det_Correct','Large_Det_Incorrect','Small_Det_Correct','Small_Det_Incorrect','Small_Det_Unpredictable','Small_Det_Uncertain')

num_record <- apply(uncertainty_array,c(1,3),function(x) sum(x))
apply(num_record[,1:5],1,sum) # check - the sum should be 2000 (n_test)
num_record <- as.data.frame(num_record)

# more columns
## samples need to be simulated
num_record$Need_Sim = num_record$Small_Det_Unpredictable + num_record$Small_Det_Uncertain

## after correction, but before simulation
num_record$Correct_before_Sim = num_record$Large_Det_Correct + num_record$Small_Det_Correct 

## after correction and simulation
Correct_after_Sim = apply(uncertainty_array[,,c('Large_Det_Correct','Small_Det_Correct','Small_Det_Unpredictable','Small_Det_Uncertain')],
                          c(1,2), function(x) any(x == 1,na.rm = T))
num_record$Correct_after_Sim = apply(Correct_after_Sim,1,sum)

# add column for n:
num_record$n <- as.numeric(gsub("n=", "", rownames(num_record)))
rownames(num_record) <- NULL

# plot
p1 <- ggplot(num_record,aes(x=n, y = Need_Sim/2000))+
  geom_bar(stat = "identity",fill='#a6bcd0')+
  ylim(0,.08)+ylab('Proportion')+
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.key.height = unit(32,'pt'),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
p1

p2 <- ggplot(num_record)+
  geom_line(aes(x=n, y = Correct_before_Sim/2000, color = 'Sample-Corrected Predictions'))+#, linetype = 'Sample-Corrected Predictions'))+
  geom_point(aes(x=n, y = Correct_before_Sim/2000, color = 'Sample-Corrected Predictions', shape = 'Sample-Corrected Predictions'), size = 2.5)+
  geom_line(aes(x=n, y = Correct_after_Sim/2000, color = 'Uncertainty-Driven Simulations'))+#, linetype = 'Uncertainty-Driven Simulation'))+
  geom_point(aes(x=n, y = Correct_after_Sim/2000, color = 'Uncertainty-Driven Simulations', shape = 'Uncertainty-Driven Simulations'), size = 2.5)+ 
  geom_line(aes(x=n, y = Raw_Correct/2000, color = 'Raw Predictions'))+#, linetype = 'Raw Predictions')) +
  geom_point(aes(x=n, y = Raw_Correct/2000, color = 'Raw Predictions', shape = 'Raw Predictions'), size = 3.5) +
  geom_hline(yintercept = 1.0, linetype = 2, color = 'purple') +
  ylim(0.82,1)+
  ylab('Accuracy')+
  scale_color_manual(
    name = NULL,
    values = c('Sample-Corrected Predictions' = "#D04D33", 'Uncertainty-Driven Simulations' = "#56c6f0", 'Raw Predictions' = '#44D361'),
    breaks = c('Raw Predictions', 'Sample-Corrected Predictions', 'Uncertainty-Driven Simulations')
  ) +
  scale_shape_manual(
    name = NULL,
    values = c('Sample-Corrected Predictions' = 16, 'Uncertainty-Driven Simulations' = 17, 'Raw Predictions' = 18)
  ) +
  # scale_linetype_manual(
  #   name = NULL,
  #   values = c('Sample-Corrected Predictions' = 1, 'Uncertainty-Driven Simulation' = 2, 'Raw Predictions' = 4)
  # ) +
  theme(panel.background = element_rect(fill = 'white', color = NA),
        panel.grid = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.key.height = unit(32,'pt'),
        legend.key.width = unit(24, 'pt'),
        legend.position = "inside",
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
p2

pdf('uncertainty.pdf', width = 12, height = 4)
p2 + p1 #+ plot_layout(heights = c(1, 1))
dev.off()
#save.image('Uncertainty_July15.Rdata')
#load('Uncertainty_July15.Rdata')


















