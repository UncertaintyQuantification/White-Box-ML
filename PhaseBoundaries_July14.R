library(RobustGaSP)
library(ggplot2)
library(dplyr)
library(parallel)
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
ppgp_entries <- ppgasp(X_train[1:n_train,], log(G_entries_train[1:n_train,]),
               isotropic = F, nugget.est = T, optimization = 'nelder-mead')

vars <- c('r','alpha_A','alpha_B','f_A','f_p','rho_s','chiABN',
          'chiASN','chiRatio','epsRatio','lb','a_plus','a_minus')

sweep_bounds <- list(
  r = c(0.5,2.0),
  alpha_A = c(0.02,0.5),
  alpha_B = c(0.02,0.5),
  f_A = c(0.01,0.5),
  f_p = c(0.01,1.0),
  rho_s = c(0.01,0.99),
  chiABN = c(0.0,250.0),
  chiASN = c(0.0,100.0),
  chiRatio = c(0.0,2.0),
  epsRatio = c(1.0,5.0),
  lb = c(0.0,0.2),
  a_plus = c(0.1,10.0),
  a_minus = c(0.1,10.0)
)

all_diagrams <- list()
grid_col1_arr <- c('f_A', 'alpha_A', 'f_A')
grid_col2_arr <- c('chiABN','chiABN','chiABN')
change_vars <- c('r','lb','alpha_A')

lims <- list(
  r = c(0.5,2.0),
  lb = c(1e-6,0.15),
  f_p = c(0.01,1.0),
  f_A = c(0.01,0.49),
  alpha_A = c(0.03,0.49),
  chiABN = c(0,150),
  epsRatio = c(1,5)
)

# changes <- list(
#   #lb = c(0.00001,0.0001,0.001),
#   lb = c(1e-6,1e-5,1e-4),
#   alpha_A = c(0.1,0.2,0.3)
# )

changes <- matrix(c(c(1.25,1.5,2.0),c(1e-6,1e-5,1e-4),c(0.1,0.2,0.3)),ncol=3)
empty_idx <- c()
do_sim = F


# Generating on a coarse grid

for(l in 1:length(grid_col1_arr)){
  #if(l != 1) next
  print('-----')
  print(l/length(grid_col1_arr))
  print('-----')
  print(changes[,l])
  phase_diagram <- gen_phase_diagram_set(ppgp_entries, N=50,
                                         grid_col1 = grid_col1_arr[l],
                                         grid_col2 = grid_col2_arr[l],
                                         grid_col1_lim = lims[[grid_col1_arr[l]]],
                                         grid_col2_lim = lims[[grid_col2_arr[l]]],
                                         alpha_A = 0.25,
                                         change_var = change_vars[l],
                                         change_vals = changes[,l],
                                         do_sim = do_sim, match_alpha_A_B = T,
                                         salt_lb_match = F)
  
  all_diagrams[[l]] <- phase_diagram
  print(unique(phase_diagram$predicted_phase_diagrams[[1]]$phase))
  print(unique(phase_diagram$predicted_phase_diagrams[[2]]$phase))
  print(unique(phase_diagram$predicted_phase_diagrams[[3]]$phase))
  print(phase_diagram$phase_boundary_plt)
  
}

ggplot(all_diagrams[[2]]$predicted_phase_diagrams[[1]],
       aes(x = alpha_A, y = chiABN, color = phase)) + geom_point()

filenames <- c('phase_bdry_fA_chiABN_r.png','phase_bdry_alphaA_chiABN_lb.png',
               'phase_bdry_fA_chiABN_alphaA.png')

options(scipen=999)
library(scales)
for(l in 1:3){
  #if(l != 1) next
  plt <- all_diagrams[[l]]$phase_boundary_plt
  plt <- plt +
    scale_color_manual(values = c('#FF6256','#5DE789','#68CAFF')) +
    guides(linetype = guide_legend(override.aes = list(color = 'black'), order = 2),
           color = guide_legend(order = 1)) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.position = 'inside',
          legend.position.inside = c(.9,.7),
          panel.background = element_blank())
  plt
  all_diagrams[[l]]$phase_boundary_plt <- plt
  # ggsave(paste0('PhaseBoundariesJun18/',filenames[l]), plot = plt,
  #        width = 5, height = 4, dpi = 300)
  # pdf(paste0('PhaseBoundariesJun18/',filenames[l]), width = 4, height = 2.6)
  # print(plt)
  # dev.off()
}

all_fine_datas <- list()
all_fine_datas_ppgp <- list()

run_sim <- T

n_grid_fine <- 100
n_grid_fine_ppgp <- 500

thr_sim <- 0.025
thr_pred <- 0.015

lb <- c(0.01,0.03)
ub <- c(0.5,150.0)

mem.maxVSize(10^10)

for(l in 1:3){
  #if(l != 1) next
  print(paste0('Plot ', l))
  plt <- all_diagrams[[l]]$phase_boundary_plt
  
  change_var <- change_vars[[l]]
  color_scale <- c('#FF6256' = changes[1,l], '#5DE789' = changes[2,l], '#68CAFF' = changes[3,l])
  
  layer <- layer_data(plt, 1)
  contours <- layer[,c('x','y')] %>% as.matrix()
  
  if(run_sim) fine_grid <- cbind(expand.grid(x = seq(lb[1],ub[1],length.out=n_grid_fine), y = seq(lb[2],ub[2],length.out=n_grid_fine)), rep(NA, n_grid_fine^2))
  fine_grid_ppgp <- expand.grid(x = seq(lb[1],ub[1],length.out=n_grid_fine_ppgp), y = seq(lb[2],ub[2],length.out=n_grid_fine_ppgp))
  fine_grid_ppgp[,3] <- NA
  if(run_sim){
    colnames(fine_grid) <- c('x','y','color')
    keep <- rep(F, n_grid_fine)
  }
  colnames(fine_grid_ppgp) <- c('x','y','color')
  keep_ppgp <- rep(F, n_grid_fine_ppgp)
  
  if(run_sim){
    pb <- txtProgressBar(min = 1, max = n_grid_fine, style = 3)
    for(i in 1:n_grid_fine){
      setTxtProgressBar(pb, i)
      for(j in 1:n_grid_fine){
        dists <- ((contours[,1] - fine_grid[(i-1)*n_grid_fine + j,1]))^2 + ((contours[,2] - fine_grid[(i-1)*n_grid_fine + j,2])/150.0)^2
        if(min(dists) < thr_sim){
          idx <- which.min(dists)
          fine_grid[(i-1)*n_grid_fine+j,3] <- color_scale[layer$colour[idx]]
          keep[(i-1)*n_grid_fine+j] <- T
        }
      }
    }
    n_new <- length(which(keep==T))
    print(n_new)
    close(pb)
  }
  
  pb <- txtProgressBar(min = 1, max = n_grid_fine_ppgp, style = 3)
  for(i in 1:n_grid_fine_ppgp){
    setTxtProgressBar(pb, i)
    for(j in 1:n_grid_fine_ppgp){
      dists <- ((contours[,1] - fine_grid_ppgp[(i-1)*n_grid_fine_ppgp + j,1]))^2 + ((contours[,2] - fine_grid_ppgp[(i-1)*n_grid_fine_ppgp + j,2])/150.0)^2
      if(min(dists) < thr_pred){
        idx <- which.min(dists)
        fine_grid_ppgp[(i-1)*n_grid_fine_ppgp+j,3] <- color_scale[layer$colour[idx]]
        keep_ppgp[(i-1)*n_grid_fine_ppgp+j] <- T
      }
    }
  }
  close(pb)
  
  n_new_ppgp <- length(which(keep_ppgp==T))
  if(run_sim) new_data <- fine_grid[which(keep),]
  new_data_ppgp <- fine_grid_ppgp[which(keep_ppgp),]
  #plot(new_data[,1],new_data[,2])
  print('')
  if(run_sim) print(n_new)
  if(run_sim){ data_set <- data.frame(
    r = rep(1.0, n_new), alpha_A = rep(0.25, n_new), alpha_B = rep(0.25, n_new),
    f_A = rep(0.5, n_new), f_p = rep(1.0, n_new), rho_s = rep(0.0, n_new),
    chiABN = rep(50.0, n_new), chiASN = rep(0.0, n_new), chiRatio = rep(0.0, n_new),
    epsRatio = rep(1.0, n_new), lb = rep(2.0, n_new), a_plus = rep(1.0, n_new), a_minus = rep(1.0, n_new)
  )
  }
  
  data_set_ppgp <- data.frame(
    r = rep(1.0, n_new_ppgp), alpha_A = rep(0.25, n_new_ppgp), alpha_B = rep(0.25, n_new_ppgp),
    f_A = rep(0.5, n_new_ppgp), f_p = rep(1.0, n_new_ppgp), rho_s = rep(0.0, n_new_ppgp),
    chiABN = rep(50.0, n_new_ppgp), chiASN = rep(0.0, n_new_ppgp), chiRatio = rep(0.0, n_new_ppgp),
    epsRatio = rep(1.0, n_new_ppgp), lb = rep(2.0, n_new_ppgp), a_plus = rep(1.0, n_new_ppgp), a_minus = rep(1.0, n_new_ppgp)
  )
  
  
  if(run_sim){
    data_set[,grid_col1_arr[l]] <- new_data$x
    data_set[,grid_col2_arr[l]] <- new_data$y
    data_set[,change_var] <- new_data[,3]
    data_set$alpha_B <- data_set$alpha_A
  }
  
  data_set_ppgp[,grid_col1_arr[l]] <- new_data_ppgp$x
  data_set_ppgp[,grid_col2_arr[l]] <- new_data_ppgp$y
  data_set_ppgp[,change_var] <- new_data_ppgp[,3]
  data_set_ppgp$alpha_B <- data_set_ppgp$alpha_A
  if(run_sim) full_data_set <- data_set
  full_data_set_ppgp <- data_set_ppgp
  # rm(data_set, data_set_ppgp, fine_grid, fine_grid_ppgp, keep, keep_ppgp)
  # for(j in 1:3){
  #   data_here <- data_set
  #   data_here[,change_var] <- changes[[change_var]][j]
  #   if(j == 1) full_data_set <- data_here
  #   else full_data_set <- rbind(full_data_set, data_here)
  # }
  if(run_sim) print(nrow(full_data_set))
  print(nrow(full_data_set_ppgp))
  
  pred_output <- predict_phase_w_correction(full_data_set_ppgp, ppgp_entries)
  phase_pred <- ifelse(pred_output$pred_phase_corrected != 1, 0, 1)
  
  print(unique(phase_pred))
  # Simulation
  if(run_sim){
    print(nrow(full_data_set))
    full_data_set <- as.matrix(full_data_set)
    sim_output <- full_sim_process_w_fast(full_data_set)
    phase_sim <- sim_output$phase
    print(unique(phase_sim))
    full_data_set <- as.data.frame(full_data_set)
    full_data_set$phase <- ifelse(phase_sim != 1, 0, 1)
    # full_data_set$phase_sim <- phase_sim
    
    all_fine_datas[[l]] <- full_data_set
  }
  full_data_set_ppgp$phase <- phase_pred
  all_fine_datas_ppgp[[l]] <- full_data_set_ppgp
}

expressions_list <- list(
  r = expression(r),
  alpha_A = expression(alpha),
  alpha_B = expression(alpha[B]),
  f_A = expression(f[A]),
  f_p = expression(f[p]),
  rho_s = expression(rho[s]),
  chiABN = expression(N[A] ~ chi[AB]),
  chiASN = expression(N[A] ~ chi[AS]),
  chiRatio = expression(r[chi]),
  epsRatio = expression(epsilon[A] ~ '/' ~ epsilon[B]),
  lb = expression(l[B]~'/'~b),
  a_plus = expression(a['+']),
  a_minus = expression(a['-'])
)

for(l in 1:3){
  #if(l != 1) next
  sim_df_here <- all_fine_datas[[l]]
  pred_df_here <- all_fine_datas_ppgp[[l]]
  
  sim_df_here$phase_sim <- sim_df_here$phase
  pred_df_here$phase_pred <- pred_df_here$phase
  
  sim_df_here <- sim_df_here[,c(grid_col1_arr[l], grid_col2_arr[l], change_vars[l], 'phase_sim')]
  pred_df_here <- pred_df_here[,c(grid_col1_arr[l], grid_col2_arr[l], change_vars[l], 'phase_pred')]
  
  sim_df_here_melt <- melt(sim_df_here, id.vars = c(grid_col1_arr[l], grid_col2_arr[l], change_vars[l]))
  pred_df_here_melt <- melt(pred_df_here, id.vars = c(grid_col1_arr[l], grid_col2_arr[l], change_vars[l]))
  
  all_df_here_melt <- rbind(sim_df_here_melt, pred_df_here_melt)
  all_df_here_melt$variable <- ifelse(all_df_here_melt$variable == 'phase_pred', 'Prediction', 'Truth')
  
  bdry_plt <- ggplot(all_df_here_melt, aes(x = !!sym(grid_col1_arr[l]), y = !!sym(grid_col2_arr[l]),
                                           color = factor(!!sym(change_vars[l])))) +
    geom_contour(aes(z = value, linetype = variable), breaks = 0.5) +
    xlim(0.0, 0.5) + ylim(0.0, 155.0) +
    labs(x = expressions_list[[grid_col1_arr[l]]], y = expressions_list[[grid_col2_arr[l]]],
         color = expressions_list[[change_vars[l]]], linetype = 'Source') +
    scale_color_manual(values = c('#F48FB1', '#AB47BC', '#6A1B9A')) +
    guides(linetype = guide_legend(override.aes = list(color = 'black'), order = 1),
           color = guide_legend(order = 2)) +
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          panel.background = element_rect(fill = 'white', color = NA),
          panel.grid = element_blank(),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8))
  
  ggsave(paste0('figures/',filenames[l]), plot = bdry_plt,
         width = 4.8, height = 3, dpi = 300)
}




#save.image('PhaseBoundaries_July14.Rdata')
#load('PhaseBoundaries_July14.Rdata')






