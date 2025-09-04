source('Phase_Functions.R')
library(RobustGaSP)
library(lhs)
library(plot3D)
library(ggplot2)
library(colorspace)
library(patchwork)
library(gganimate)
library(reshape2)
library(dplyr)



var_names = c('r','alpha_A','alpha_B','f_A','f_p','rho_s','chiABN','chiASN','chiRatio',
              'epsRatio','lb','a_plus','a_minus')

# 13d paired training data
set.seed(1)
print('Starting Paired Data Generation...')
n_archs_paired_train <- 100 # Use this many architectures
n_paired_per_arch_train <- 750 # Pair this many points with each architecture
X_train_paired <- matrix(NA, n_archs_paired_train * n_paired_per_arch_train, 13)
y_train_paired <- matrix(NA, n_archs_paired_train * n_paired_per_arch_train, 1)

archs_lhs_train <- maximinLHS(n_archs_paired_train, 3) # Simulate architectures LHS
archs_lhs_train[,1] <- archs_lhs_train[,1]*1.5 + 0.5
archs_lhs_train[,2] <- archs_lhs_train[,2] * 0.95 + 0.025 # Slightly smaller than [0,1] since boundaries are fast special cases
archs_lhs_train[,3] <- archs_lhs_train[,3] * 0.95 + 0.025

# Do batched LHS samples, all together is too slow
others <- runif(10*n_archs_paired_train*n_paired_per_arch_train) %>% matrix(nrow = n_archs_paired_train*n_paired_per_arch_train)
others <- matrix(NA, nrow = n_archs_paired_train*n_paired_per_arch_train, 10)
print('Starting LHS Generation of Non-Architecture Data')
for(i in 1:n_archs_paired_train){
  if(i %% 5 == 0) print(paste0('Sampling ', i, 'th LHS Sample out of ', n_archs_paired_train, ' of size ', n_paired_per_arch_train))
  # Pair with ith architecture
  others[(1+(i-1)*n_paired_per_arch_train):(i*n_paired_per_arch_train),] <- maximinLHS(n_paired_per_arch_train, 10)
}
# Rescale
others[,4] <- others[,4] * 100 # chiABN
others[,5] <- others[,5] * 100 # chiASN
others[,6] <- others[,6] * 2 # chiRatio
others[,7] <- 1 + others[,7] * 4 # epsRatio
others[,8] <- others[,8] * 10 # lb
others[,9] <- others[,9] * 9.9 + 0.1 # a_plus
others[,10] <- others[,10] * 9.9 + 0.1 # a_minus

# Combine
X_train_paired[,4:13] <- others
for(i in 1:n_archs_paired_train){
  # Plug ith architecture into X_train_paired
  X_train_paired[(1+(i-1)*n_paired_per_arch_train):(i*n_paired_per_arch_train),1] <- 
    rep(archs_lhs_train[i,1],n_paired_per_arch_train)
  X_train_paired[(1+(i-1)*n_paired_per_arch_train):(i*n_paired_per_arch_train),2] <- 
    rep(archs_lhs_train[i,2],n_paired_per_arch_train)
  X_train_paired[(1+(i-1)*n_paired_per_arch_train):(i*n_paired_per_arch_train),3] <- 
    rep(archs_lhs_train[i,3],n_paired_per_arch_train)
}

# Simulating determinants and phase
det_record_train_paired <- matrix(NA, dim(X_train_paired)[1], 200)
for(i in 1:dim(X_train_paired)[1]){
  if(i %% 50 == 0) print(paste0('Simulation ', i, ' out of ', dim(X_train_paired)[1]))
  
  # Avoid resimulating G.  Otherwise, we're looking at ~6 hours for 56k pts
  if(i > 1){
    if(sum(X_train_paired[i,1:3] == X_train_paired[i-1,1:3]) == 3){
      G_A_input = sim_out$G_A
      G_B_input = sim_out$G_B
      
    } else {
      print('New Architecture')
      G_A_input = NULL
      G_B_input = NULL
    }
  } else {
    G_A_input = NULL
    G_B_input = NULL
  }
  
  
  sim_out <- full_sim_process(r=X_train_paired[i,1],alpha_a=X_train_paired[i,2],alpha_b=X_train_paired[i,3],
                              f_A=X_train_paired[i,4],f_p=X_train_paired[i,5],rho_s=X_train_paired[i,6],
                              chiABN=X_train_paired[i,7],chiASN=X_train_paired[i,8],chiRatio=X_train_paired[i,9],
                              lb=X_train_paired[i,10],epsilonRatio=X_train_paired[i,11],a_plus=X_train_paired[i,12],a_minus=X_train_paired[i,13],
                              G_A_input = G_A_input, G_B_input = G_B_input)
  
  det_record_train_paired[i,] <- sim_out$detSinv
  y_train_paired[i] <- sim_out$phase
  
}


########### Non-Paired Data (13 pars)
set.seed(1)
print('Starting Generation of Non-Paired Data...')
# Initialize datasets
n_train_nonpaired <- 5000
X_train_nonpaired <- runif(n_train_nonpaired*13) %>% matrix(n_train_nonpaired, 13) # Very Slow for LHS
y_train_nonpaired <- rep(NA, n_train_nonpaired)
det_train_nonpaired <- matrix(NA, n_train_nonpaired, 200)

# Rescale nonpaired train data to proper ranges.  Slightly shrink alpha ranges
# since the case of 0 charged or 0 uncharged beads are special cases with simple and fast conditionals
X_train_nonpaired[,1] <- X_train_nonpaired[,1] * 1.5 + 0.5 # r
X_train_nonpaired[,2] <- X_train_nonpaired[,2] * 0.95 + 0.025 # alpha_a
X_train_nonpaired[,3] <- X_train_nonpaired[,3] * 0.95 + 0.025 # alpha_b
X_train_nonpaired[,7] <- X_train_nonpaired[,7] * 100 # chiABN
X_train_nonpaired[,8] <- X_train_nonpaired[,8] * 100 # chiASn
X_train_nonpaired[,9] <- X_train_nonpaired[,9] * 2 # chiRatio
X_train_nonpaired[,10] <- X_train_nonpaired[,10] * 4 + 1 # epsilonRatio
X_train_nonpaired[,11] <- X_train_nonpaired[,11] * 10 # lb
X_train_nonpaired[,12] <- X_train_nonpaired[,12] * 9.9 + 0.1 # a_plus
X_train_nonpaired[,13] <- X_train_nonpaired[,13] * 9.9 + 0.1 # a_minus

det_record_train_nonpaired <- matrix(NA, n_train_nonpaired, 200)
for(i in 1:n_train_nonpaired){
  if(i %% 10 == 0) print(paste0('Simulating Iteration ', i, ' out of ', n_train_nonpaired, '...'))
  # Simulate nonpaired data
  sim_out <- full_sim_process(r=X_train_nonpaired[i,1],alpha_a=X_train_nonpaired[i,2],alpha_b=X_train_nonpaired[i,3],
                              f_A=X_train_nonpaired[i,4],f_p=X_train_nonpaired[i,5],rho_s=X_train_nonpaired[i,6],
                              chiABN=X_train_nonpaired[i,7],chiASN=X_train_nonpaired[i,8],chiRatio=X_train_nonpaired[i,9],
                              lb=X_train_nonpaired[i,10],epsilonRatio=X_train_nonpaired[i,11],a_plus=X_train_nonpaired[i,12],a_minus=X_train_nonpaired[i,13],
                              G_A_input = NULL, G_B_input = NULL)
  
  det_record_train_nonpaired[i,] <- sim_out$detSinv # update nonpaired det(S^-1)
  y_train_nonpaired[i] <- sim_out$phase # update nonpaired phase record
}


########### Testing Data for 13d space
set.seed(1)
print('Starting Generation of 13D Testing Data...')
# Same process as nonpaired data.  This data is just heldout for testing
n_test_13d <- 2000
X_test <- runif(13*n_test_13d) %>% matrix(nrow = n_test_13d, ncol = 13)
G11_A_test <- G12_A_test <- G22_A_test <- G_A_det_inv_test <- 
  G11_B_test <- G12_B_test <- G22_B_test <- G_B_det_inv_test <- matrix(NA, n_test_13d, 200)
y_test <- rep(NA, n_test_13d)
det_Sinv_record_test <- matrix(NA, n_test_13d, 200)

X_test[,1] <- X_test[,1] * 1.5 + 0.5 # r
X_test[,2] <- X_test[,2] * 0.95 + 0.025 # alpha_a
X_test[,3] <- X_test[,3] * 0.95 + 0.025 # alpha_b
X_test[,7] <- X_test[,7] * 100 # chiABN
X_test[,8] <- X_test[,8] * 100 # chiASN
X_test[,9] <- X_test[,9] * 2 # chiRatio
X_test[,10] <- X_test[,10] * 4 + 1 # epsilonRatio
X_test[,11] <- X_test[,11] * 10 # lb
X_test[,12] <- X_test[,12] * 9.9 + 0.1 # a_plus
X_test[,13] <- X_test[,13] * 9.9 + 0.1 # a_minus

for(i in 1:n_test_13d){
  if(i %% 10 == 0) print(paste0('Simulating Iteration ', i, ' out of ', n_test_13d, '...'))
  sim_out <- full_sim_process(r=X_test[i,1], alpha_a=X_test[i,2], alpha_b=X_test[i,3],
                              f_A=X_test[i,4], f_p=X_test[i,5], rho_s=X_test[i,6],
                              chiABN=X_test[i,7], chiASN=X_test[i,8], chiRatio=X_test[i,9],
                              epsilonRatio=X_test[i,10], lb=X_test[i,11], a_plus=X_test[i,12], a_minus=X_test[i,13])
  G_A_here <- sim_out$G_A
  G_B_here <- sim_out$G_B
  det_S_inv_here <- sim_out$detSinv
  phase_here <- sim_out$phase
  
  G11_A_test[i,] <- G_A_here[,1,1]
  G12_A_test[i,] <- G_A_here[,1,2]
  G22_A_test[i,] <- G_A_here[,2,2]
  G_A_det_inv_test[i,] <- 1/(G_A_here[,1,1]*G_A_here[,2,2]-G_A_here[,1,2]^2)
  
  G11_B_test[i,] <- G_B_here[,1,1]
  G12_B_test[i,] <- G_B_here[,1,2]
  G22_B_test[i,] <- G_B_here[,2,2]
  G_B_det_inv_test[i,] <- 1/(G_B_here[,1,1]*G_B_here[,2,2]-G_B_here[,1,2]^2)
  
  y_test[i] <- phase_here
  det_Sinv_record_test[i,] <- det_S_inv_here
}


# Testing

print('Starting Generation of 2D Testing Data')
# Grid for 2d data
test_grid_length <- 40

all_pairs <- matrix(0,nrow=test_grid_length^2,ncol=2)
n_range <- round(seq(1, 200, length.out = test_grid_length))
for(i in 1:test_grid_length){
  for(j in 1:test_grid_length){
    all_pairs[(i-1)*test_grid_length+j,] <- c(n_range[i],n_range[j])
  }
}
# Truncate the grid
valid <- which(all_pairs[,1] + all_pairs[,2] >= 50 & all_pairs[,1] + all_pairs[,2] <= 200)
X_test_2d <- all_pairs[valid,]

n_test_2d <- dim(X_test_2d)[1]

G11_test_2d <- matrix(NA, n_test_2d, 200)
G12_test_2d <- matrix(NA, n_test_2d, 200)
G22_test_2d <- matrix(NA, n_test_2d, 200)
G_det_inv_test_2d <- matrix(NA, n_test_2d, 200)
# Simulate
for(i in 1:n_test_2d){
  if(i %% 20 == 0) print(paste0('Simulating ', i, 'th 2d testing data out of ', n_test_2d))
  G_data_here <- G_sim_process(X_test_2d[i,1], X_test_2d[i,2])
  G11_test_2d[i,] <- G_data_here$G11
  G12_test_2d[i,] <- G_data_here$G12
  G22_test_2d[i,] <- G_data_here$G22
  G_det_inv_test_2d[i,] <- G_data_here$det_G_inv
}

G_entries_test_2d <- cbind(G11_test_2d, G12_test_2d, G22_test_2d)


n_datas <- 10

N <- 500
n0 <- 20
n_check <- 75

all_X_train <- array(NA, dim = c(N, 2, n_datas))
all_X_train_nonpaired <- array(NA, dim = c(N, 13, n_datas))
all_phases <- matrix(NA, N, n_datas)
all_G_entries_train <- array(NA, dim = c(N, 600, n_datas))
all_G_det_inv_train <- array(NA, dim = c(N, 200, n_datas))


set.seed(1)
for(i in 1:n_datas){
  print('-------------------------------------------')
  print(paste0('----------Generating ', i, 'th data set----------'))
  print('-------------------------------------------')
  data_i <- gen_training_data(N, n_check, n0, plot_it = T)
  X_train_i <- data_i$input
  G_entries_train_i <- data_i$G_entries_train
  G_det_inv_train_i <- data_i$G_det_inv_train
  
  X_nonpaired_i <- runif(N*13) %>% matrix(N, 13)
  lb <- c(0.5,0.025,0.025,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.1,0.1)
  ub <- c(2.0,0.975,0.975,1.0,1.0,1.0,100.0,100.0,2.0,5.0,10.0,10.0,10.0)
  
  for(l in 1:13) X_nonpaired_i[,l] <- lb[l] + X_nonpaired_i[,l] * (ub[l] - lb[l])
  X_nonpaired_i_df <- as.data.frame(X_nonpaired_i)
  colnames(X_nonpaired_i_df) <- c('r','alpha_a','alpha_b','f_A','f_p','rho_s',
                                  'chiABN','chiASN','chiRatio',
                                  'epsRatio','lb','a_plus','a_minus')
  all_phases_i <- rep(NA, N) %>% as.matrix()
  for(l in 1:N){
    if(l %% 10 == 0) print(paste0('Simulating for iteration ', l, '...'))
    sim_output <- full_sim_process(X_nonpaired_i_df[l,1],X_nonpaired_i_df[l,2],X_nonpaired_i_df[l,3],
                                   X_nonpaired_i_df[l,4],X_nonpaired_i_df[l,5],X_nonpaired_i_df[l,6],
                                   X_nonpaired_i_df[l,7],X_nonpaired_i_df[l,8],X_nonpaired_i_df[l,9],
                                   X_nonpaired_i_df[l,10],X_nonpaired_i_df[l,11],X_nonpaired_i_df[l,12],
                                   X_nonpaired_i_df[l,13])
    all_phases_i[l] <- sim_output$phase
  }
  
  all_X_train[,,i] <- X_train_i
  all_G_entries_train[,,i] <- G_entries_train_i
  all_G_det_inv_train[,,i] <- G_det_inv_train_i
  all_X_train_nonpaired[,,i] <- X_nonpaired_i
  all_phases[,i] <- all_phases_i
  
}

X_train <- all_X_train[,,1]
G_entries_train <- all_G_entries_train[,,1]
G_det_inv_train <- all_G_det_inv_train[,,1]

grid_num <- 10
n_match <- 10
X_test_grid <- runif(grid_num^3 * n_match * 13) %>% matrix(ncol=13)
lb <- c(0.5,0.02,0.02,0.1,0.1,0.1,0.0,0.0,0.0,1.0,0.0,0.1,0.1)
ub <- c(2.0,0.98,0.98,1.0,1.0,1.0,100.0,100.0,2.0,5.0,10.0,10.0,10.0)
r_grid <- seq(0.5,2.0,length.out=grid_num)
alpha_a_grid <- alpha_b_grid <- seq(0.02,0.98,length.out=grid_num)
for(i in 4:13){
  X_test_grid[,i] <- X_test_grid[,i] * (ub[i]-lb[i]) + lb[i]
}
colnames(X_test_grid) <- c('r','alpha_a','alpha_b','f_A','f_p','rho_s','chiABN','chiASN','chiRatio','epsilonRatio',
                           'lb','a_plus','a_minus')
arch_grid <- expand.grid(r_grid, alpha_a_grid, alpha_b_grid)
for(i in 1:(nrow(X_test_grid)/n_match)){
  idx_here <- (1+(i-1)*n_match):(i*n_match)
  X_test_grid[idx_here,1] <- arch_grid[i,1]
  X_test_grid[idx_here,2] <- arch_grid[i,2]
  X_test_grid[idx_here,3] <- arch_grid[i,3]
}

phase_record_grid <- rep(NA, nrow(X_test_grid))
det_Sinv_record_grid <- matrix(NA, nrow(X_test_grid), 200)
G_A_grid <- G_B_grid <- array(NA, dim = c(nrow(X_test_grid), 2, 2, 200))

print(paste0('Simulating ', nrow(X_test_grid), ' points'))
pb <- txtProgressBar(min=1, max = nrow(X_test_grid), style = 3)
for(i in 1:nrow(X_test_grid)){
  setTxtProgressBar(pb, i)
  sim_out_here <- full_sim_process(X_test_grid[i,1], X_test_grid[i,2], X_test_grid[i,3],
                                   X_test_grid[i,4], X_test_grid[i,5], X_test_grid[i,6],
                                   X_test_grid[i,7], X_test_grid[i,8], X_test_grid[i,9],
                                   X_test_grid[i,10], X_test_grid[i,11], X_test_grid[i,12], X_test_grid[i,13])
  phase_record_grid[i] <- sim_out_here$phase
  det_Sinv_record_grid[i,] <- sim_out_here$detSinv
  G_A_grid[i,,,] <- sim_out_here$G_A
  G_B_grid[i,,,] <- sim_out_here$G_B
}


# Write CSVs
write.csv(X_train, 'All_Datasets/X_train.csv', row.names = F) # 2d training data
write.csv(G_entries_train, 'All_Datasets/G_entries_train.csv', row.names = F)

write.csv(X_test_2d, 'All_Datasets/X_test_2d.csv', row.names = F) # 2d testing data
write.csv(G_entries_test_2d, 'All_Datasets/G_entries_test_2d.csv', row.names = F)

colnames(X_test) <- var_names
write.csv(X_test, 'All_Datasets/X_test.csv', row.names = F, col.names = var_names) # nongrid testing data
write.csv(y_test, 'All_Datasets/y_test.csv', row.names = F)

colnames(X_test_grid) <- var_names
write.csv(X_test_grid, 'All_Datasets/X_test_grid.csv', row.names = F) # grid testing data
write.csv(phase_record_grid, 'All_Datasets/phase_record_grid.csv', row.names = F)
# write.csv(G_A_grid, 'All_Datasets/G_A_grid.csv', row.names = F)

colnames(X_train_paired) <- var_names
write.csv(X_train_paired, 'All_Datasets/X_train_paired.csv', row.names = F)
write.csv(y_train_paired, 'All_Datasets/y_train_paired.csv', row.names = F)

write.csv(X_train_nonpaired, 'All_Datasets/X_train_nonpaired.csv', row.names = F, col.names = var_names)
write.csv(y_train_nonpaired, 'All_Datasets/y_train_nonpaired.csv', row.names = F)

for(i in 1:n_datas){
  filename_for_2d_training_data <- paste0('All_Datasets/all_X_train_',i,'.csv')
  filename_for_13d_training_data <- paste0('All_Datasets/all_X_train_nonpaired_',i,'.csv')
  filename_for_2d_training_data_response <- paste0('All_Datasets/all_G_entries_train_',i,'.csv')
  filename_for_13d_training_data_response <- paste0('All_Datasets/all_phases_',i,'.csv')
  
  write.csv(all_X_train[,,i], filename_for_2d_training_data, row.names = F)
  all_X_train_nonpaired_here <- all_X_train_nonpaired[,,i]
  colnames(all_X_train_nonpaired_here) <- var_names
  write.csv(all_X_train_nonpaired_here, filename_for_13d_training_data, row.names = F)
  write.csv(all_G_entries_train[,,i], filename_for_2d_training_data_response, row.names = F)
  write.csv(all_phases[,i], filename_for_13d_training_data_response, row.names = F)
}




# save.image('Datasets_Jun2.Rdata')

