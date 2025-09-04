# load('Datasets_Jun2.Rdata')

library(RobustGaSP)
library(ggplot2)
library(reshape)
source('Phase_Functions.R')
# source('Phase_Correction_Clayton_Testing_July10.R')
library(xgboost)
library(randomForest)
library(zeallot)
library(reticulate)
library(tensorflow)
library(dplyr)
library(tidyr)

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



virtualenv_create("r-tensorflow")
virtualenv_install("r-tensorflow", packages = c("tensorflow==2.16.1"))
use_virtualenv("r-tensorflow", required = TRUE)

do_method_comparisons <- function(X_train, G_entries_train,
                                  X_train_nonpaired, y_train_nonpaired,
                                  X_test, n_vec, phase_true){
  # all_U_and_pre <- gen_all_U_and_prefactors_vec(X_test, batch_size = nrow(X_test), num_cores = 8)
  # all_U <- all_U_and_pre$all_U
  # all_pre <- all_U_and_pre$all_pre
  
  X_train_nonpaired <- as.matrix(X_train_nonpaired)
  #X_test <- as.matrix(X_test)
  
  N_A <- 100
  N_B <- round(N_A * X_test[,1])
  N_Ac <- round(N_A * X_test[,2])
  N_Au <- N_A - N_Ac
  N_Bc <- round(N_B * X_test[,3])
  N_Bu <- N_B - N_Bc
  
  X_test_A <- cbind(N_Ac, N_Au)
  X_test_B <- cbind(N_Bc, N_Bu)
  
  colnames(X_test) <- colnames(X_train_nonpaired) <-
    c('r','alpha_a','alpha_b','f_A','f_p','rho_s','chiABN','chiASN','chiRatio','epsRatio','lb','a_plus','a_minus')
  train_ind <- 1:500 #sample(nrow(X_train_nonpaired), size = max(n_vec), replace = F)
  ppgp_accs <- ppgp_accs_no_nug <- nn_accs <- rf_accs <- gb_accs <- rep(NA, length(n_vec))
  
  for(i in 1:length(n_vec)){
    print(paste0('Beginning Iteration ', i,' with ', n_vec[i], ' training data'))
    Sys.sleep(2)
    n_train <- n_vec[i]
    
    X_train_nonpaired_here <- as.matrix(X_train_nonpaired[train_ind[1:n_train],])
    y_train_nonpaired_here <- y_train_nonpaired[train_ind[1:n_train]]
    y_train_nonpaired_here_ohe <- tf$keras$utils$to_categorical(y_train_nonpaired_here, num_classes = as.integer(3))
    print(1)
    # ppgp models
    ppgp_entries <- ppgasp(design = X_train[1:n_train,], response = log(G_entries_train[1:n_train,]),
                           isotropic = F, nugget.est = T, optimization = 'nelder-mead', method = 'mle')
    # ppgp_det_inv <- ppgasp(design = X_train[1:n_train,], response = log(G_det_inv_train[1:n_train,]),
    #                        isotropic = F, nugget.est = T, optimization = 'nelder-mead', method = 'mle')
    # 
    #ppgp_entries_no_nug <- ppgasp(design = X_train[1:n_train,], response = log(G_entries_train[1:n_train,]),
    #                       isotropic = F, nugget.est = F, optimization = 'nelder-mead')
    #ppgp_det_inv_no_nug <- ppgasp(design = X_train[1:n_train,], response = log(G_det_inv_train[1:n_train,]),
    #                       isotropic = F, nugget.est = F, optimization = 'nelder-mead')
    
    # NN
    nn_classifier <- tf$keras$Sequential(list(
      tf$keras$layers$Input(list(as.integer(13))),
      tf$keras$layers$Dense(units = as.integer(512), activation = 'relu'),
      tf$keras$layers$Dense(units = as.integer(512), activation = 'relu'),
      tf$keras$layers$Dense(units = as.integer(512), activation = 'relu'),
      tf$keras$layers$Dense(units = as.integer(3), activation = 'softmax')
    ))
    nn_classifier$compile(
      optimizer = tf$keras$optimizers$RMSprop(),
      loss = 'categorical_crossentropy',
      metrics = list('accuracy')
    )
    
    nn_classifier$fit(
      x = X_train_nonpaired_here,
      y = y_train_nonpaired_here_ohe,
      validation_split = 0.2,
      epochs = as.integer(200),
      batch_size = min(as.integer(64), as.integer(round(dim(X_train_nonpaired_here)[1]/2)))
    )
    
    rf <- randomForest(x = X_train_nonpaired_here, y = factor(y_train_nonpaired_here))
    D_gb <- xgb.DMatrix(data = X_train_nonpaired_here, label = y_train_nonpaired_here)
    gb <- xgboost(data = D_gb, params = list(objective = 'multi:softmax'),
                  nrounds = 100, num_class = 3)
    
    
    ##### Preds
    
    # PP-GP
    preds_ppgp <- predict_phase_w_correction(X_test, ppgp_entries)
    phase_gp <- preds_ppgp$pred_phase_corrected
    ppgp_accs[i] <- sum(phase_gp == phase_true) / length(phase_true)
    print(paste0('PP-GP Accuracy: ', ppgp_accs[i]))

    # NN
    nn_preds <- nn_classifier$predict(X_test) %>% apply(1, which.max) - 1
    nn_accs[i] <- sum(nn_preds == phase_true) / length(phase_true)
    print(paste0('NN Accuracy: ', nn_accs[i]))
    # rf
    print(names(X_train_nonpaired_here))
    print(names(X_test))
    rf_preds <- predict(rf, X_test)
    rf_accs[i] <- sum(rf_preds == phase_true) / length(phase_true)
    print(paste0('RF Accuracy: ', rf_accs[i]))
    # gb
    gb_preds <- predict(gb, X_test)
    gb_accs[i] <- sum(gb_preds == phase_true) / length(phase_true)
    print(paste0('GB Accuracy: ', gb_accs[i]))
    
  }
  
  accs_df <- data.frame(
    n = n_vec,
    PPGP = ppgp_accs,
    #PPGP_no_nug = ppgp_accs_no_nug,
    NN = nn_accs,
    RF = rf_accs,
    GB = gb_accs
  )
  
  return(accs_df)
}



set.seed(1)
n_vec <- c(10,seq(50, 500, by = 50))

all_ppgp_accs <- all_nn_accs <- all_rf_accs <- all_gb_accs <- matrix(NA, 5, length(n_vec))
for(i in 1:5){
  print(paste0('Beginning Dataset ', i))
  Sys.sleep(5)
  accs_i<-do_method_comparisons(all_X_train[,,i],all_G_entries_train[,,i],
                                all_X_train_nonpaired[,,i],all_phases[,i],
                                X_test,n_vec,y_test)
  all_ppgp_accs[i,] <- accs_i$PPGP
  all_nn_accs[i,] <- accs_i$NN
  all_rf_accs[i,] <- accs_i$RF
  all_gb_accs[i,] <- accs_i$GB
  
}

all_ppgp_accs_grid <- all_nn_accs_grid <- all_rf_accs_grid <- all_gb_accs_grid <- matrix(NA, 5, length(n_vec))
for(i in 1:5){
  print(paste0('Beginning Dataset ', i))
  Sys.sleep(5)
  accs_i<-do_method_comparisons(all_X_train[,,i],all_G_entries_train[,,i],
                                all_X_train_nonpaired[,,i],all_phases[,i],
                                as.matrix(X_test_grid),n_vec,phase_record_grid)
  all_ppgp_accs_grid[i,] <- accs_i$PPGP
  all_nn_accs_grid[i,] <- accs_i$NN
  all_rf_accs_grid[i,] <- accs_i$RF
  all_gb_accs_grid[i,] <- accs_i$GB
  
}

# Uniform Testing
mean_ppgp_accs <- apply(all_ppgp_accs, 2, mean)
mean_nn_accs <- apply(all_nn_accs, 2, mean)
mean_rf_accs <- apply(all_rf_accs, 2, mean)
mean_gb_accs <- apply(all_gb_accs, 2, mean)

min_ppgp_accs <- apply(all_ppgp_accs, 2, min)
min_nn_accs <- apply(all_nn_accs, 2, min)
min_rf_accs <- apply(all_rf_accs, 2, min)
min_gb_accs <- apply(all_gb_accs, 2, min)

max_ppgp_accs <- apply(all_ppgp_accs, 2, max)
max_nn_accs <- apply(all_nn_accs, 2, max)
max_rf_accs <- apply(all_rf_accs, 2, max)
max_gb_accs <- apply(all_gb_accs, 2, max)

uniform_accuracy_data <- data.frame(
  n = n_vec,
  PPGP_mean = mean_ppgp_accs,
  NN_mean = mean_nn_accs,
  RF_mean = mean_rf_accs,
  GB_mean = mean_gb_accs,
  PPGP_min = min_ppgp_accs,
  NN_min = min_nn_accs,
  RF_min = min_rf_accs,
  GB_min = min_gb_accs,
  PPGP_max = max_ppgp_accs,
  NN_max = max_nn_accs,
  RF_max = max_rf_accs,
  GB_max = max_gb_accs
)

# Grid Testing
mean_ppgp_accs_grid <- apply(all_ppgp_accs_grid, 2, mean)
mean_nn_accs_grid <- apply(all_nn_accs_grid, 2, mean)
mean_rf_accs_grid <- apply(all_rf_accs_grid, 2, mean)
mean_gb_accs_grid <- apply(all_gb_accs_grid, 2, mean)

min_ppgp_accs_grid <- apply(all_ppgp_accs_grid, 2, min)
min_nn_accs_grid <- apply(all_nn_accs_grid, 2, min)
min_rf_accs_grid <- apply(all_rf_accs_grid, 2, min)
min_gb_accs_grid <- apply(all_gb_accs_grid, 2, min)

max_ppgp_accs_grid <- apply(all_ppgp_accs_grid, 2, max)
max_nn_accs_grid <- apply(all_nn_accs_grid, 2, max)
max_rf_accs_grid <- apply(all_rf_accs_grid, 2, max)
max_gb_accs_grid <- apply(all_gb_accs_grid, 2, max)

grid_accuracy_data <- data.frame(
  n = n_vec,
  PPGP_mean = mean_ppgp_accs_grid,
  NN_mean = mean_nn_accs_grid,
  RF_mean = mean_rf_accs_grid,
  GB_mean = mean_gb_accs_grid,
  PPGP_min = min_ppgp_accs_grid,
  NN_min = min_nn_accs_grid,
  RF_min = min_rf_accs_grid,
  GB_min = min_gb_accs_grid,
  PPGP_max = max_ppgp_accs_grid,
  NN_max = max_nn_accs_grid,
  RF_max = max_rf_accs_grid,
  GB_max = max_gb_accs_grid
)

uniform_accuracy_data_longer <- pivot_longer(uniform_accuracy_data, cols = -n, 
                                             names_to = c('Model','.value'), 
                                             names_pattern = '(.+)_(.+)')

grid_accuracy_data_longer <- pivot_longer(grid_accuracy_data, cols = -n,
                                          names_to = c('Model','.value'),
                                          names_pattern = '(.+)_(.+)')

labels_map = c(
  'PPGP' = 'PPGP',
  'NN' = 'NN',
  'RF' = 'RF',
  'GB' = 'GB'
)

levels(uniform_accuracy_data_longer$Model) <- c('PPGP', 'NN', 'RF', 'GB')
levels(grid_accuracy_data_longer$Model) <- c('PPGP', 'NN', 'RF', 'GB')

uniform_plt <- ggplot(uniform_accuracy_data_longer, aes(x = n, y = mean, color = Model)) +
  geom_line(aes(linetype = Model), linewidth = 0.8) +
  geom_errorbar(aes(ymax = max, ymin = min)) +
  geom_point(aes(shape = Model)) +
  labs(x = expression(n[train]), y = 'Accuracy') +
  scale_color_manual(values = c('PPGP' = '#4f9da6', 'NN' = '#1a0841', 'RF' = '#ffad5a', 'GB' = '#ff5959'), labels = labels_map,
                     breaks = c('PPGP', 'NN', 'RF', 'GB')) +
  scale_shape_manual(values = 21:24, labels = labels_map, breaks = c('PPGP', 'NN', 'RF', 'GB')) +
  scale_linetype_manual(values = 1:4, labels = labels_map, breaks = c('PPGP', 'NN', 'RF', 'GB')) +
  ylim(0.45, 1.0) +
  theme_classic() + 
  theme(legend.title = element_text(size=12),legend.text=element_text(size=10),
        axis.title=element_text(size=13),axis.text=element_text(size=11),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA)) 

grid_plt <- ggplot(grid_accuracy_data_longer, aes(x = n, y = mean, color = Model)) +
  geom_line(aes(linetype = Model), linewidth = 0.8) +
  geom_errorbar(aes(ymax = max, ymin = min)) +
  geom_point(aes(shape = Model)) +
  labs(x = expression(n[train]), y = 'Accuracy') +
  scale_color_manual(values = c('PPGP' = '#4f9da6', 'NN' = '#1a0841', 'RF' = '#ffad5a', 'GB' = '#ff5959'), labels = labels_map,
                     breaks = c('PPGP', 'NN', 'RF', 'GB')) +
  scale_shape_manual(values = 21:24, labels = labels_map, breaks = c('PPGP', 'NN', 'RF', 'GB')) +
  scale_linetype_manual(values = 1:4, labels = labels_map, breaks = c('PPGP', 'NN', 'RF', 'GB')) +
  ylim(0.45, 1.0) +
  theme_classic() + 
  theme(legend.title = element_text(size=12),legend.text=element_text(size=10),
        axis.title=element_text(size=13),axis.text=element_text(size=11),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA)) 


pdf('uniform_accuracy_plt.pdf', width = 4, height = 2.7)
uniform_plt
dev.off()

pdf('grid_accuracy_plt.pdf', width = 4, height = 2.7)
grid_plt
dev.off()

# save.image('MethodComparison_July10.Rdata')
