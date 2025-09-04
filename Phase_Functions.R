library(parallel)



compute_phi <- function(Na, fA, r, alpha_a, alpha_b, phi_p){ # Volume fractions for blend
  if(fA == 0){
    charge <- 1
    A <- matrix(c(-round(alpha_b*Na*r), 1,1,1), 2, 2)
    b <- matrix(c(0,phi_p), 2)
    sola <- solve(A, t(b))
    sol <- rep(0,3)
    sol[2] <- sola[1]
    sol[3] <- sola[2]
  } else if(fA == 1){
    charge <- -1
    A <- matrix(c(round(alpha_a*Na),1,-1,1),2,2)
    b <- matrix(c(0,phi_p),2)
    sola <- solve(A, t(b))
    sol <- rep(0,3)
    sol[1] <- sola[1]
    sol[3] <- sola[2]
  } else {
    if(fA*round(alpha_a*Na) > (1-fA)*round(alpha_b*Na*r)){
      charge <- -1
      A <- matrix(c(round(alpha_a*Na),1,1,-round(alpha_b*Na*r),1,0,-1,1,fA),3,3)
      b <- matrix(c(0,phi_p,fA*phi_p),1,3)
    } else {
      charge <- 1
      A <- matrix(c(round(alpha_a*Na),1,1,-round(alpha_b*Na*r),1,0,1,1,fA),3,3)
      b <- matrix(c(0,phi_p,fA*phi_p),1,3)
    }
    sol <- solve(A, t(b))
  }
  return(list(sol=sol, charge=charge))
}

chain_tracker_return_obj <- function(phi_list, olist){ # Track Chain elements
  element <- list()
  index <- c()
  count <- 0
  for(i in 1:length(phi_list)){
    phi <- phi_list[i]
    ele <- olist[i]
    if(phi != 0){
      element[[i]] <- ele
      index <- c(index, count)
    }
    count <- count + 1
  }
  return(list(element = element, index = index))
}

chain_tracker_update_phi <- function(phi_list, Glist){ # Update Volume Fractions
  return_obj_out <- chain_tracker_return_obj(phi_list, Glist)
  Glist_new <- return_obj_out$element
  chain_list <- return_obj_out$index
  new_phi_list <- rep(NA, length(chain_list))
  #print(chain_list)
  for(i in chain_list){
    new_phi_list[i+1] <- phi_list[i+1] 
  }
  
  return(
    list(phi_list = new_phi_list, Glist = Glist_new)
  )
  
}

composition_solver_solve <- function(Na, r, alpha_a, alpha_b,
                                     phi_p, phi_c, fA){ # Determine Blend Composition
  phi_return <- rep(0,5)
  solved_compute_phi <- compute_phi(Na, fA, r, alpha_a, alpha_b, phi_p)
  composition <- solved_compute_phi$sol
  charge <- solved_compute_phi$charge
  phi_return[1] <- composition[1] #A
  phi_return[2] <- composition[2] #B
  if(phi_c > 0){
    phi_return[3] <- phi_c/2
    phi_return[4] <- phi_c/2
  }
  if(charge > 0) phi_return[3] <- phi_return[3] + composition[3]
  else if(charge < 0) phi_return[4] <- phi_return[4] + composition[3]
  
  for(i in 1:4){
    if(phi_return[i] < 1e-10) phi_return[i] <- 0
  }
  
  phi_return[5] <- 1 - phi_p - phi_c
  
  return(phi_return)
}

create_chain_details <- function(N, r, alpha_a, alpha_b, blist, Nref){ # Records
  return(
    list(Na = N, r = r, alpha_a = alpha_a, alpha_b = alpha_b,
         Charge_Pattern = 'even', blist = blist, Nref = Nref)
  )
}


read_chain_details <- function(chain_details){ # Read Records
  pattern <- chain_details$Charge_Pattern
  Na <- as.numeric(round(chain_details$Na))
  r <- chain_details$r
  alpha_a <- chain_details$alpha_a
  alpha_b <- chain_details$alpha_b
  Nb <- round(Na*r)
  Na_c <- round(alpha_a*Na)
  Nb_c <- round(alpha_b*Nb)
  blist <- chain_details$blist
  Nref <- chain_details$Nref
  
  Na_seq <- sequence(Na,Na_c)
  Nb_seq <- sequence(Nb,Nb_c)
  
  return(list(
    Na = Na, Na_seq = Na_seq, Nb = Nb, Nb_seq = Nb_seq,
    blist = blist, Nref = Nref
  ))
}


sequence <- function(N, Ncharge){ # Evenly spaced charged beads
  if(Ncharge == 0) return(list(charged = c(), uncharged = 1:N))
  else if(Ncharge == N) return(list(charged = 1:N, uncharged = c()))
  else {
    charged <- round(seq(1,N,N/Ncharge) - 1 + floor(N/Ncharge/2), 0)
    uncharged <- c()
    for(i in 0:(N-1)){
      if(!(i %in% charged)) uncharged <- c(uncharged, i)
    }
    return(list(charged=charged, uncharged=uncharged))
  }
}

Phi <- function(k, c=1) exp(-c * k^2/6) # Linker Function

doubleSum <- function(k, seqA, seqB){ # Entries @ specific k
  Phik <- Phi(k)
  sumTemp <- 0
  for(i in seqA){
    for(j in seqB){
      #arg <- sqrt(abs(i-j)) * k
      sumTemp <- sumTemp + Phik^abs(i-j)
    }
  }
  return(sumTemp)
}


computeGij <- function(krange, seqA, seqB){ # Compute double sum across valid k
  tempVec <- rep(NA, length(krange))
  for(i in 1:length(krange)){
    k <- krange[i]
    tempVec[i] <- doubleSum(k, seqA, seqB)
  }
  return(tempVec)
}

logspace <- function(a, b, length.out=200, base = 10){
  # Space data on logscale from base^a to base^b.  Their code uses base 10
  return(base^seq(a,b,length.out=length.out))
}

calculateG <- function(Ncharge, Nuncharge){ # Determine G for A/B block
  kpoints <- 200
  krange <- logspace(-2, 1, length.out=kpoints)
  
  #Ncharge <- round(N * alpha)
  sequences <- sequence(Ncharge + Nuncharge, Ncharge)
  
  G <- array(0, dim = c(kpoints, 2, 2))
  
  if(Ncharge == 0) { 
    G[,1,1] <- computeGij(krange, sequences$uncharged, sequences$uncharged)
  } else if(Ncharge == Ncharge + Nuncharge) {
    G[,2,2] <- computeGij(krange, sequences$charged, sequences$charged)
  } else {
    G[,1,1] <- computeGij(krange, sequences$uncharged, sequences$uncharged)
    G[,2,2] <- computeGij(krange, sequences$charged, sequences$charged)
    G[,1,2] <- G[,2,1] <- computeGij(krange, sequences$charged, sequences$uncharged)
  }
  return(G)
}


# Make U

# For Record Purposes (Line 200-228)
charge_to_index <- list(
  zero = 0, negone = 1, one = 2
)

species_to_index <- list(
  A = 0, B = 1, S_plus = 2, S_minus = 3, S = 4
)

charge_str_to_num <- list(
  zero = 0, negone = -1, one = 1
)

index_to_charge <- c('zero','negone','one')

index_to_species <- c('A','B','S_+','S_-','S')

poss.charges <- c(0,-1,1)

getlenbeadtypes <- function(BeadTypes) length(BeadTypes[,1])

getspecies <- function(BeadTypes, i) index_to_species[[BeadTypes[i,1]+1]]

getcharge <- function(BeadTypes, i){
  
  charge_i_here <- index_to_charge[BeadTypes[i,2]+1]
  #print(c(BeadTypes[i,2]+1, charge_i_here, BeadTypes[i,]))
  
  return(charge_str_to_num[[charge_i_here]])
}

# For Calculating Interactions
Gamma2 <- function(k2, a) exp(-k2*(a^2))^2

BornSolvation <- function(interact) -interact[['lb']]/interact[['a']]/2*(interact[['epsAepsB']]-1)

# For Recording Interactions
create_interactions <- function(chiab=0.0, chias=0.0, chi_contrast=0.0, epsilon_ratio=0.0,
                                lb=0.0, a_plus=0.0, a_minus=0.0, zeta=0.0, a_i=0.0){
  chibs <- chias * chi_contrast
  interactions <- list( list1 = list(type = 'chi', zero = 'A', one = 'B', value = chiab, smear_length = a_i),
                        list2 = list(type = 'born', zero = 'A', one = 'S_+', value = list(epsAepsB = epsilon_ratio,
                                                                                  lb = lb, a = a_plus),
                             smear_length = a_i),
                        list3 = list(type = 'born', zero = 'A', one = 'S_-', value = list(epsAepsB = epsilon_ratio,
                                                                                  lb = lb, a = a_minus),
                             smear_length = a_i),
                        list4 = list(type = 'chi', zero = 'A', one = 'S', value = chias, smear_length = a_i),
                        list5 = list(type = 'chi', zero = 'B', one = 'S', value = chibs, smear_length = a_i),
                        list6 = list(type = 'zeta', value = zeta, smear_length = a_i),
                        list7 = list(type = 'electrostatic', value = lb, smear_length = a_i)
                       )
  return(interactions)
}

# Creating a homopolymer w/ N beads where Nc are charged and Z = \pm 1, species \in {A,B}
create_charged_homopolymer <- function(N, Nc, Z, species){ # Nc is the charge sequence
  G <- list()
  for(i in 1:N){
    if(i %in% Nc){
      G[[i]] <-  list(species = species, charge = Z)
    }
    else{
      G[[i]] <- list(species = species, charge = 'zero')
    }
  }
  return(G)
}

# Same as above except single-bead salt/solvent
create_salt_solvent <- function(Z, species){
  G <- list()
  G[[1]] <- list(species = species, charge = Z)
  return(G)
}

InitializeChains <- function(Na, Na_seq, Nb, Nb_seq){
  Glist <- list(
    A = create_charged_homopolymer(Na, Na_seq, 'one', 'A'),
    B = create_charged_homopolymer(Nb, Nb_seq, 'negone', 'B'),
    saltPlus = create_salt_solvent('one', 'S_plus'),
    saltMinus = create_salt_solvent('negone', 'S_minus'),
    solvent = create_salt_solvent('zero', 'S')
  )
  return(Glist)
}

# Record beads present in the mix of A, B, cationic salt, anionic salt, solvent
sort_charge <- function(Glist, species_to_index, charge_to_index){
  BeadTypes <- matrix(0, nrow=1, ncol=2)
  for(i in 1:length(Glist)){
    G_here <- Glist[[i]]
    for(j in 1:length(G_here)){
      G_here_here <- G_here[[j]]
      species <- G_here_here$species
      charge <- (G_here_here$charge)
      #print(charge)
      #print(c(species, charge))
      #print(species)
      #print(species_to_index[[species]])
     
      if(!is.null(species) & !is.null(charge)){
        BeadTypes <- rbind(BeadTypes, matrix(c(species_to_index[[species]], charge_to_index[[charge]]), nrow = 1, ncol=2 ))
      }#print(species_to_index[G[[i]][['species']]])
      #print(charge_to_index[G[[i]][['charge']]])
      #BeadTypes <- rbind(BeadTypes, 
      #                    matrix(unlist(c(species_to_index[G[[i]][['species']]], charge_to_index[G[[i]][['charge']]])),
      #                           nrow = 1, ncol=2))
    }
  }
  BeadTypes <- BeadTypes[-1,]
  #print(BeadTypes)
  BeadTypes <- unique(BeadTypes)
  #print(BeadTypes)
  sort.idx <- sort(BeadTypes[,2], index.return=T)$ix
  #print(sort.idx)
  return(BeadTypes[sort.idx,])
}


# Calculate U(k_1),...,U(k_m)
return_interaction_matrix <- function(interactions, karray, BeadTypes){
  k2 <- karray^2
  interaction_matrix <- array(0, dim=c(length(karray), getlenbeadtypes(BeadTypes), getlenbeadtypes(BeadTypes)))
  for(i in 1:getlenbeadtypes(BeadTypes)){
    for(j in 1:getlenbeadtypes(BeadTypes)){
      #print(paste0('New Interaction', i, j))
      species_i <- getspecies(BeadTypes, i)
      species_j <- getspecies(BeadTypes, j)
      charge_i <- as.numeric(getcharge(BeadTypes, i))
      charge_j <- as.numeric(getcharge(BeadTypes, j))
      #print(c(species_i,species_j,charge_i,charge_j))
      for(k in 1:length(interactions)){
        #print(k)
        #print((species_j == interactions[[k]][['zero']]))
        G2 <- try(Gamma2(k2, interactions[[k]][['smear_length']]))
        if(class(G2) == 'try-error') G2 <- 1
        if(interactions[[k]][['type']] == 'chi'){
          if(species_i == interactions[[k]][['zero']] & species_j == interactions[[k]][['one']]){
            #print(1)
            interaction_matrix[,i,j] <- interaction_matrix[,i,j] + (interactions[[k]][['value']] * G2)
          } else if(species_j == interactions[[k]][['zero']] & species_i == interactions[[k]][['one']]){
            #print(2)
            interaction_matrix[,i,j] <- interaction_matrix[,i,j] + (interactions[[k]][['value']] * G2)
          }
        } else if(interactions[[k]][['type']] == 'zeta'){
          #print(3)
          interaction_matrix[,i,j] <- interaction_matrix[,i,j] + interactions[[k]][['value']] * G2
        } else if(interactions[[k]][['type']] == 'electrostatic' & charge_i != 0 & charge_j != 0){
          #print(4)
          #print(c(charge_i, charge_j))
          e_ <- ((4*pi*interactions[[k]][['value']])/k2)
          netcharge <- charge_i * charge_j
          interaction_matrix[,i,j] <- interaction_matrix[,i,j] + e_ * netcharge * G2
        } else if(interactions[[k]][['type']] == 'born'){
          if(species_i == interactions[[k]][['zero']] & species_j == interactions[[k]][['one']]){
            #print(5)
            interaction_matrix[,i,j] <- interaction_matrix[,i,j] + (BornSolvation(interactions[[k]][['value']]) * G2)
          } else if(species_j == interactions[[k]][['zero']] & species_i == interactions[[k]][['one']]){
            #print(6)
            interaction_matrix[,i,j] <- interaction_matrix[,i,j] + (BornSolvation(interactions[[k]][['value']]) * G2)
          }
        }
      }
    }
  }
  return(interaction_matrix)
}



# Eliminate Zeros rows/columns from G to ensure invertibility
slice_zeros <- function(G){
  keep_rows <- c()
  keep_cols <- c()
  
  G_sum <- matrix(0, 7, 7)
  for(i in 1:7){
    for(j in 1:7){
      G_sum[i,j] <- sum(G[,i,j])
    }
  }
  for(i in 1:7){
    #print(sum(G_sum[i,]))
    #print(sum(G_sum[,i]))
    if(sum(G_sum[i,]) != 0 || sum(is.na(G_sum[i,]))>0 || sum(is.nan(G_sum[i,]))>0) keep_rows <- c(keep_rows, i)
    if(sum(G_sum[,i]) != 0 || sum(is.na(G_sum[i,]))>0 || sum(is.nan(G_sum[i,]))>0) keep_cols <- c(keep_cols, i)
  }
  #print(keep_rows)
  #print(keep_cols)
  
  #print(dim(G[,keep_rows,keep_cols]))
  return(G[,keep_rows,keep_cols])
}

# Permute G to match U
permute_matrix <- function(A){
  #print(dim(A))
  # Does the weird row-reordering stuff to make G be in the non-block form
  perm_order <- c(1,3,7,4,6,2,5)
  P <- matrix(0, 7, 7)
  for(i in 1:7){
    j <- perm_order[i]
    P[i,j] <- 1
  }
  #print(P)
  return(P %*% A %*% t(P))
}

# Run the simulation
simulate <- function(r, alpha_a, alpha_b, f_A, f_p, rho_s, chiABN, chiASN,
                     chiRatio, lb, epsilonRatio, a_plus, a_minus,
                     calc_G = T, use_model = F, models_list = NULL){
  
  N <- 100
  blist <- rep(1,5)
  Nref <- 1
  zeta <- 0.0
  a_i <- 0.0
  kmin <- -2
  kmax <- 1
  kpoints <- 200
  kinfo <- list(kmin = kmin, kmax = kmax, kpoints = kpoints)
  project_down <- F
  
  chain_details <- create_chain_details(N, r, alpha_a, alpha_b, blist, Nref)
  
  phi_list <- rep(0.2, 5) # default phi_list
  karray <- logspace(kmin, kmax, kpoints)
  
  read_chain_details_out <- read_chain_details(chain_details)
  
  # Access read_chain_details results
  Na <- read_chain_details_out$Na
  Na_seq <- read_chain_details_out$Na_seq
  Nb <- read_chain_details_out$Nb
  Nb_seq <- read_chain_details_out$Nb_seq
  blist_sys <- read_chain_details_out$blist
  Nref_sys <- read_chain_details_out$Nref
  
  Glist <- InitializeChains(Na, Na_seq$charged, Nb, Nb_seq$charged)
  
  Nlist <- rep(0, length(Glist))
  for(i in 1:length(Glist)) Nlist[i] <- length(Glist[[i]])
  #print(length(Glist))
  BeadTypes <- sort_charge(Glist, species_to_index, charge_to_index)
  #print(BeadTypes)
  
  rho_p <- (1-rho_s)*f_p
  rho_c <- (1-rho_s)*(1-f_p)
  phi_list <- composition_solver_solve(Na, r, alpha_a, alpha_b, rho_p, rho_c, f_A)
  
  prefactors <- rep(NA, 5)
  for(i in 1:5) prefactors[i] <- Nlist[i]*phi_list[i]/Nref_sys/length(Glist[[i]])^2
  
  G <- array(0, dim=c(200,7,7))
  NA_charge <- round(Na * alpha_a)
  NA_uncharge <- Na - NA_charge
  NB_charge <- round(Nb * alpha_b)
  NB_uncharge <- Nb - NB_charge
  #print(c(NA_charge, NA_uncharge, NB_charge, NB_uncharge))

  G[,1:2,1:2] <- calculateG(NA_charge, NA_uncharge) #* prefactors[1]
  G[,3:4,3:4] <- calculateG(NB_charge, NB_uncharge) #* prefactors[2]
  G[,5,5] <- 1 #*prefactors[3]
  G[,6,6] <- 1 #*prefactors[4]
  G[,7,7] <- 1 #*prefactors[5]

  
  updated_phi <- chain_tracker_update_phi(phi_list, Glist)
  phi_list <- updated_phi$phi_list
  Glist <- list()
  Glist$A <- updated_phi$Glist[[1]]$A
  Glist$B <- updated_phi$Glist[[2]]$B
  Glist$saltPlus <- updated_phi$Glist[[3]]$saltPlus
  Glist$saltMinus <- updated_phi$Glist[[4]]$saltMinus
  Glist$solvent <- updated_phi$Glist[[5]]$solvent
  remove_ind <- c()
  for(i in 1:5){
    if(is.null(Glist[i][[1]])) remove_ind <- c(remove_ind, i)
  }
  #print(remove_ind)
  if(!is.null(remove_ind)) Glist <- Glist[-c(remove_ind)]
  BeadTypes <- sort_charge(Glist, species_to_index, charge_to_index)
  # phi_list = ChainTracker.return_list(phi_list)
  
  interactions_list <- create_interactions(chiABN/N, chiASN/N, chiRatio, epsilonRatio,
                                           lb, a_plus, a_minus, zeta, a_i)
  
  
  U <- return_interaction_matrix(interactions_list, karray, BeadTypes)
  
  return(list(G_mat = G, U_mat = U, prefactors = prefactors))
}


# Update G with the prefactors
add_prefactors <- function(G, prefactors){
  G[,1:2,1:2] <- G[,1:2,1:2] * prefactors[1]
  G[,3:4,3:4] <- G[,3:4,3:4] * prefactors[2]
  G[,5,5] <- prefactors[3]
  G[,6,6] <- prefactors[4]
  G[,7,7] <- prefactors[5]
  return(G)
}

# Calculate G inverse manually w/ block structure, mindful of zero rows/columns
invert_G <- function(G){
  Ginv <- array(0, dim=c(200,7,7))
  for(i in 1:200){
    det_GA <- G[i,1,1]*G[i,2,2] - G[i,1,2]^2
    Ginv[i,1,1] <- G[i,2,2] / det_GA
    Ginv[i,1,2] <- Ginv[i,2,1] <- -G[i,1,2] / det_GA
    Ginv[i,2,2] <- G[i,1,1] / det_GA
    det_GB <- G[i,3,3]*G[i,4,4] - G[i,3,4]^2
    Ginv[i,3,3] <- G[i,4,4] / det_GB
    Ginv[i,3,4] <- Ginv[i,4,3] <- -G[i,3,4] / det_GB
    Ginv[i,4,4] <- G[i,3,3] / det_GB
    if(G[i,5,5] != 0) Ginv[i,5,5] <- 1/G[i,5,5]
    if(G[i,6,6] != 0) Ginv[i,6,6] <- 1/G[i,6,6]
    if(G[i,7,7] != 0) Ginv[i,7,7] <- 1/G[i,7,7]
  }
  return(Ginv)
}

# permuting G(k_1),...,G(k_m)
permute_G_array <- function(Garr){
  for(i in 1:200){
    Garr[i,,] <- permute_matrix(Garr[i,,])
  }
  return(Garr)
}

# Calculate det S^-1(k) = det( G^-1(k) + U(k) ) for k = k_1,...,k_m
determinant_Sinv <- function(Ginv, U){
  Sinv <- Ginv + U
  det_Sinv <- rep(NA, dim(Ginv)[1])
  for(i in 1:length(det_Sinv)) det_Sinv[i] <- det(Ginv[i,,] + U[i,,])
  return(det_Sinv)
}

det_S <- function(S){
  dets <- rep(NA, dim(S)[1])
  for(k in 1:dim(S)[1]){
    det_here <- determinant(S[k,,], logarithm=T)
    dets[k] <- det_here$sign * exp(det_here$modulus)
  }
  return(dets)
}

# Determine phase homo/macro/micro
checkphase <- function(detcurve, n_macro = 1){
  if(sum(detcurve > 0) == length(detcurve)) return(0)
  else{
    if(n_macro > 1){
      if(sum(detcurve[1:n_macro] < 0)/n_macro > 0.5 & which.min(detcurve) %in% 1:n_macro){
        # If majority of first n_macro points are negative, macro
        return(1)
      }
    } else {
      # Only checking first point
      if(detcurve[1] < 0 & which.min(detcurve) == 1) return(1)
    }
  }
  
  return(2)
}

# U equivalent for removing zero rows/columns. This time NA
slice_na <- function(arr){
  keep_rows <- c()
  keep_cols <- c()
  sum_arr <- matrix(0,dim(arr)[2],dim(arr)[3])
  for(i in 1:dim(arr)[2]){
    for(j in 1:dim(arr)[3]){
      sum_arr[i,j] <- sum(arr[,i,j])
    }
  }
  for(i in 1:dim(sum_arr)[1]){
    if(sum(is.na(sum_arr[i,])) != dim(sum_arr)[1]) keep_rows <- c(keep_rows, i)
    if(sum(is.na(sum_arr[,i])) != dim(sum_arr)[2]) keep_cols <- c(keep_cols, i)
  }
  return(list(rows=keep_rows, cols=keep_cols))
}

# Function to get only U and the prefactors:
get_U_and_prefactors <- function(r, alpha_a, alpha_b, f_A, f_p, rho_s, chiABN, chiASN,
                                 chiRatio, lb, epsilonRatio, a_plus, a_minus){
  
  N <- 100
  blist <- rep(1,5)
  Nref <- 1
  zeta <- 0.0
  a_i <- 0.0
  kmin <- -2
  kmax <- 1
  kpoints <- 200
  kinfo <- list(kmin = kmin, kmax = kmax, kpoints = kpoints)
  project_down <- F
  
  chain_details <- create_chain_details(N, r, alpha_a, alpha_b, blist, Nref)
  
  phi_list <- rep(0.2, 5) # default phi_list
  karray <- logspace(kmin, kmax, kpoints)
  
  read_chain_details_out <- read_chain_details(chain_details)
  
  # Access read_chain_details results
  Na <- read_chain_details_out$Na
  Na_seq <- read_chain_details_out$Na_seq
  Nb <- read_chain_details_out$Nb
  Nb_seq <- read_chain_details_out$Nb_seq
  blist_sys <- read_chain_details_out$blist
  Nref_sys <- read_chain_details_out$Nref
  
  Glist <- InitializeChains(Na, Na_seq$charged, Nb, Nb_seq$charged)
  
  Nlist <- rep(0, length(Glist))
  for(i in 1:length(Glist)) Nlist[i] <- length(Glist[[i]])
  #print(length(Glist))
  BeadTypes <- sort_charge(Glist, species_to_index, charge_to_index)
  #print(BeadTypes)
  
  rho_p <- (1-rho_s)*f_p
  rho_c <- (1-rho_s)*(1-f_p)
  phi_list <- composition_solver_solve(Na, r, alpha_a, alpha_b, rho_p, rho_c, f_A)
  
  prefactors <- rep(NA, 5)
  for(i in 1:5) prefactors[i] <- Nlist[i]*phi_list[i]/Nref_sys/length(Glist[[i]])^2
  
  #G <- array(0, dim=c(200,7,7))
  NA_charge <- round(Na * alpha_a)
  NA_uncharge <- Na - NA_charge
  NB_charge <- round(Nb * alpha_b)
  NB_uncharge <- Nb - NB_charge
  #print(c(NA_charge, NA_uncharge, NB_charge, NB_uncharge))
  
  updated_phi <- chain_tracker_update_phi(phi_list, Glist)
  phi_list <- updated_phi$phi_list
  Glist <- list()
  if(length(updated_phi$Glist) >= 1) Glist$A <- updated_phi$Glist[[1]]$A
  if(length(updated_phi$Glist) >= 2) Glist$B <- updated_phi$Glist[[2]]$B
  if(length(updated_phi$Glist) >= 3) Glist$saltPlus <- updated_phi$Glist[[3]]$saltPlus
  if(length(updated_phi$Glist) >= 4) Glist$saltMinus <- updated_phi$Glist[[4]]$saltMinus
  if(length(updated_phi$Glist) >= 5) Glist$solvent <- updated_phi$Glist[[5]]$solvent
  # Glist$A <- updated_phi$Glist[[1]]$A
  # Glist$B <- updated_phi$Glist[[2]]$B
  # Glist$saltPlus <- updated_phi$Glist[[3]]$saltPlus
  # Glist$saltMinus <- updated_phi$Glist[[4]]$saltMinus
  # Glist$solvent <- updated_phi$Glist[[5]]$solvent
  remove_ind <- c()
  for(i in 1:5){
    if(is.null(Glist[i][[1]])) remove_ind <- c(remove_ind, i)
  }
  #print(remove_ind)
  if(!is.null(remove_ind)) Glist <- Glist[-c(remove_ind)]
  BeadTypes <- sort_charge(Glist, species_to_index, charge_to_index)
  # phi_list = ChainTracker.return_list(phi_list)
  
  interactions_list <- create_interactions(chiABN/N, chiASN/N, chiRatio, epsilonRatio,
                                           lb, a_plus, a_minus, zeta, a_i)
  
  
  U <- return_interaction_matrix(interactions_list, karray, BeadTypes)
  
  return(list(U_mat = U, prefactors = prefactors))
}


# return list of BeadCounts, G, Ginv, U, Sinv, detSinv, phase s.t.
# BeadCounts == c(N_A,N_Ac,N_Au,N_B,N_Bc,N_Bu)
full_sim_process <- function(r, alpha_a, alpha_b, f_A, f_p, rho_s,
                             chiABN, chiASN, chiRatio, 
                             epsilonRatio, lb, a_plus, a_minus, G_A_input = NULL,
                             G_B_input = NULL){
  N <- 100
  blist <- rep(1,5)
  Nref <- 1
  zeta <- 0.0
  a_i <- 0.0
  kmin <- -2
  kmax <- 1
  kpoints <- 200
  kinfo <- list(kmin = kmin, kmax = kmax, kpoints = kpoints)
  project_down <- T
  
  N_A <- 100
  N_B <- round(r * N_A)
  N_Ac <- round(N_A*alpha_a)
  N_Bc <- round(N_B*alpha_b)
  N_Au <- N_A - N_Ac
  N_Bu <- N_B - N_Bc
  
  if(is.null(G_A_input) || is.null(G_B_input)){
    G_A <- calculateG(N_Ac, N_Au)
    G_B <- calculateG(N_Bc, N_Bu)
  } else {
    G_A <- G_A_input
    G_B <- G_B_input
  }
  
  
  U_and_prefactors <- get_U_and_prefactors(r, alpha_a, alpha_b, 
                                           f_A, f_p, rho_s,
                                           chiABN, chiASN, chiRatio, 
                                           lb, epsilonRatio, a_plus, a_minus)
  
  prefactors <- U_and_prefactors$prefactors
  U <- U_and_prefactors$U_mat
  
  G <- array(0, dim = c(kpoints, 7, 7))
  G[,1:2,1:2] <- G_A * prefactors[1]
  G[,3:4,3:4] <- G_B * prefactors[2]
  G[,5,5] <- prefactors[3]
  G[,6,6] <- prefactors[4]
  G[,7,7] <- prefactors[5]
  
  Ginv <- invert_G(G)
  Ginv <- slice_zeros(permute_G_array(Ginv))
  
  Sinv <- Ginv + U
  det_curve <- rep(NA, kpoints)
  for(k in 1:kpoints){
    det_curve[k] <- det(Sinv[k,,])
  }
  
  phase <- checkphase(det_curve)
  
  return(list(
    BeadCounts = c(N_A, N_Ac, N_Au, N_B, N_Bc, N_Bu),
    G_A = G_A,
    G_B = G_B,
    Ginv = Ginv,
    prefactors = prefactors,
    #phi_list = phi_list,
    inputs = list(r=r,alpha_a=alpha_a,alpha_b=alpha_b,
                  f_A=f_A,f_p=f_p,rho_s=rho_s,
                  chiABN=chiABN,chiASN=chiASN,chiRatio=chiRatio,
                  epsilonRatio=epsilonRatio,lb=lb,a_plus=a_plus,a_minus=a_minus),
    U = U,
    Sinv = Sinv,
    detSinv = det_curve,
    phase = phase
  ))
  
}



# From a dataset, generate every U and prefactors together
gen_all_U_and_prefactors <- function(data){
  n_data <- nrow(data)
  all_U <- array(NA, dim=c(n_data,200,7,7))
  all_pre <- matrix(NA, n_data, 5)
  for(i in 1:n_data){
    #print(i)
    if(i %% 20 == 0) print(paste0(i,'/',n_data))
    U_and_pre_here <- get_U_and_prefactors(r=data[i,1], alpha_a=data[i,2], alpha_b=data[i,3],
                                           f_A=data[i,4], f_p=data[i,5], rho_s=data[i,6],
                                           chiABN=data[i,7], chiASN=data[i,8], chiRatio=data[i,9],
                                           epsilonRatio=data[i,10], lb=data[i,11], a_plus=data[i,12], a_minus=data[i,13])
    U_here <- U_and_pre_here$U_mat
    U_size <- dim(U_here)[2]
    pre_here <- U_and_pre_here$prefactors
    
    all_U[i,,1:U_size,1:U_size] <- U_here
    all_pre[i,] <- pre_here
    
  }
  return(list(
    all_U = all_U,
    all_pre = all_pre
  ))
}





Ginv_from_preds_and_prefactors <- function(G11_A, G12_A, G22_A, G_det_inv_A, 
                                           G11_B, G12_B, G22_B, G_det_inv_B,
                                           prefactors){
  Ginv <- array(0, dim = c(200, 7, 7))
  
  Ginv[,1,1] <- G22_A * G_det_inv_A / prefactors[1]
  Ginv[,1,2] <- Ginv[,2,1] <- -G12_A * G_det_inv_A / prefactors[1]
  Ginv[,2,2] <- G11_A * G_det_inv_A / prefactors[1]
  
  Ginv[,3,3] <- G22_B * G_det_inv_B / prefactors[2]
  Ginv[,3,4] <- Ginv[,4,3] <- -G12_B * G_det_inv_B / prefactors[2]
  Ginv[,4,4] <- G11_B * G_det_inv_B / prefactors[2]
  
  if(prefactors[3] != 0) Ginv[,5,5] <- 1/prefactors[3]
  if(prefactors[4] != 0) Ginv[,6,6] <- 1/prefactors[4]
  if(prefactors[5] != 0) Ginv[,7,7] <- 1/prefactors[5]
  
  return(Ginv)
}
# Compile a group of G predictions/entries, U, prefactors into phase vector
phase_from_G_U_pre <- function(all_G11_A, all_G12_A, all_G22_A, all_G_det_inv_A,
                               all_G11_B, all_G12_B, all_G22_B, all_G_det_inv_B,
                               all_U, all_prefactors, silent = T, do_exp = T){
  N <- nrow(all_G11_A)
  if(is.null(N)) N <- 1
  phase_record_out <- rep(NA, N)
  det_S_inv_record_out <- matrix(NA, N, 200)
  for(i in 1:N){
    
    if(i > 4 & silent == F) if(i %% round(N/4) == 0) print(paste0('Checking ', i, 'th phase of ', N, '...'))
    Ginv_here <- Ginv_from_preds_and_prefactors(all_G11_A[i,], all_G12_A[i,], all_G22_A[i,], all_G_det_inv_A[i,],
                                                all_G11_B[i,], all_G12_B[i,], all_G22_B[i,], all_G_det_inv_B[i,],
                                                all_prefactors[i,])
    
    Ginv_here <- slice_zeros(permute_G_array(Ginv_here))
    U_here <- all_U[i,,,]
    
    keep <- slice_na(U_here)
    U_here_no_NA <- U_here[,keep$rows,keep$cols]
    
    det_S_inv_here <- determinant_Sinv(Ginv_here, U_here_no_NA)
    
    det_S_inv_record_out[i,] <- det_S_inv_here
    phase_record_out[i] <- checkphase(det_S_inv_here)
  }
  
  return(
    list(phase_record = phase_record_out, det_S_inv_record = det_S_inv_record_out)
  )
}

# matrix of determinants -> phase.  Good for the Data -> det(Sinv) maps
phase_from_det_matrix <- function(det_matrix){
  phase_out <- apply(det_matrix, 1, checkphase)
  return(phase_out)
}


case_study_alpha_vary <- function(ppgp, 
                                  r = 1.0, f_A = 0.5, f_p = 1.0, rho_s = 0.0,
                                  chiABN = 50.0, chiASN = 0.0, chiRatio = 1.0,
                                  epsRatio = 1.0, lb = 2.0, a_plus = 1.0, a_minus = 1.0,
                                  do_sim = T, grid_size = 26){
  
  alpha_grid = expand.grid(x = seq(0.02, 0.98, length.out = grid_size),
                           y = seq(0.02, 0.98, length.out = grid_size))
  
  case_study_data <- matrix(NA, nrow = grid_size^2, ncol = 13)
  
  case_study_data[,1] <- r
  case_study_data[,2] <- alpha_grid[,1]
  case_study_data[,3] <- alpha_grid[,2]
  
  case_study_data[,4] <- f_A
  case_study_data[,5] <- f_p
  case_study_data[,6] <- rho_s
  
  case_study_data[,7] <- chiABN
  case_study_data[,8] <- chiASN
  case_study_data[,9] <- chiRatio
  
  case_study_data[,10] <- epsRatio
  case_study_data[,11] <- lb
  case_study_data[,12] <- a_plus
  case_study_data[,13] <- a_minus
  
  
  
  pred_time <- system.time(
    for(j in 1:1) pred_output <- predict_phase_w_correction(case_study_data, ppgp)
  )
  
  
  if(do_sim){
    
    sim_time <- system.time(
      for(j in 1:1) sim_output <- full_sim_process_w_fast(case_study_data)
    )
    
    case_study_data <- as.data.frame(case_study_data)
    colnames(case_study_data) <- c('r','alpha_A','alpha_B','f_A','f_p','rho_s',
                                   'chiABN','chiASN','chiRatio','epsRatio',
                                   'lb','a_plus','a_minus')
    
    case_study_data$pred_phase <- pred_output$pred_phase_corrected
    case_study_data$true_phase <- sim_output$phase
    
    return(list(
      sim_time = sim_time,
      pred_time = pred_time,
      case_study_data = case_study_data,
      do_sim = do_sim
    ))
    
  } else {
    
    case_study_data$pred_phase <- pred_output$pred_phase_corrected
    
    case_study_data <- as.data.frame(case_study_data)
    colnames(case_study_data) <- c('r','alpha_A','alpha_B','f_A','f_p','rho_s',
                                   'chiABN','chiASN','chiRatio','epsRatio',
                                   'lb','a_plus','a_minus')
    
    return(list(
      pred_time = pred_time,
      case_study_data = case_study_data,
      do_sim = do_sim
    ))
    
  }
  
}

plot_case_output <- function(case_study, title = ''){
  case_study_data <- case_study$case_study_data[,c('alpha_A','alpha_B','pred_phase','true_phase')]
  
  case_study_data$pred_phase <- factor(case_study_data$pred_phase, levels = c(0,1,2), labels = c('Homo','Macro','Micro'))
  case_study_data$true_phase <- factor(case_study_data$true_phase, levels = c(0,1,2), labels = c('Homo','Macro','Micro'))
  
  dummy_data <- data.frame(
    alpha_A = rep(0, 3),
    alpha_B = rep(0, 3),
    pred_phase = factor(c(0,1,2), levels = c(0,1,2), labels = c('Homo','Macro','Micro')),
    true_phase = factor(c(0,1,2), levels = c(0,1,2), labels = c('Homo','Macro','Micro'))
  )
  
  color_scale <- scale_color_manual(
    name = "Phase",
    #values = c("Homo" = "#369acc", "Macro" = "#f8e16f", "Micro" = "#95cf92"),
    #values = c("Homo" = "#369acc", "Macro" = "#f4b63a", "Micro" = "#5cbf66"),
    #values = c("Homo" = "#1a96e3", "Macro" = "#eec164", "Micro" = "#46ada7"),
    values = c("Homo" = "#369acc", "Macro" = "#f8c100", "Micro" = "#66bb77"),
    
    breaks = c('Homo','Macro','Micro'),
    drop = F
  )
  
  shape_scale <- scale_shape_manual(
    name = "Phase Source",
    values = c("Truth" = 15, "Prediction" = 16)
  )
  
  data <- rbind(case_study_data, dummy_data)
  
  # truth_size <- 4
  # pred_size <- 4/5.5 * truth_size
  
  plt <- ggplot(data, aes(x = alpha_A, y = alpha_B)) +
    geom_point(aes(color = true_phase, shape = 'Truth'), alpha = 0.5, size = 4.5) +
    geom_point(aes(color = pred_phase, shape = 'Prediction'), alpha = 1.0, size = 3.5) +
    color_scale + shape_scale +
    labs(x = expression(alpha[A]), y = expression(alpha[B]), title = title,
         subtitle = paste0('Accuracy: ', round(mean(data$true_phase == data$pred_phase), 3) * 100, '%')) +
    theme_minimal() + coord_fixed() +
    xlim(0.02,0.98) + ylim(0.02,0.98)
  
  return(plt)
}


# Training Data generation w/ space-filling method

fill_dist <- function(x,X) min(apply(
  X, 1, function(y) sqrt(sum((matrix(x,nrow=1) - y)^2))
))

gen_training_data <- function(n, n_check, n_init = NULL, plot_it = F){
  
  full_square <- matrix(NA, nrow=200^2,ncol=2)
  for(i in 1:200){
    for(j in 1:200){
      full_square[(i-1)*200+j,] <- c(i,j)
    }
  }
  valid_ind <- which(
    full_square[,1]+full_square[,2] >= 50 & full_square[,1] + full_square[,2] <= 200 & full_square[,1] > 1 & full_square[,2] > 1
  )
  bdry_ind <- which((full_square[,1] + full_square[,2] == 49 | full_square[,1] + full_square[,2] == 201 | full_square[,1] == 1 | full_square[,2] == 1) &
                      full_square[,1] + full_square[,2] >= 49)
  
  #bdry_ind <- which(
  #  (full_square[,1] == 1 | full_square[,2] == 1) &
  #    (full_square[,1] + full_square[,2] > 49) &
  #    (full_square[,2] + full_square[,2] < 201)
  #  )
  
  full_region <- full_square[valid_ind,]
  region_bdry <- full_square[bdry_ind,]
  
  n_train_full <- n
  if(is.null(n_init)) n_init <- round(0.1*n_train_full)
  training_index <- sample(1:nrow(full_region), n_init, replace = F)
  print('Generating Input Space w/ Space-Filling Design...')
  pb = txtProgressBar(min = 0, max = n_train_full-n_init, style = 3,
                      initial = 0)
  for(j in 1:(n_train_full-n_init)){
    setTxtProgressBar(pb, j)
    check_idx <- setdiff(sample(1:nrow(full_region), n_check), training_index)
    fill_dists <- rep(NA, n_check)
    for(i in 1:n_check){
      if(check_idx[i] %in% training_index){
        fill_dists[i] <- -1e6
      } else {
        fill_dists[i] <- #(1-alpha)*fill_dist(full_region[check_idx[i],], (full_region[training_index,])) +
          fill_dist(full_region[check_idx[i],], rbind(full_region[training_index,],region_bdry))
      }
    }
    to_add <- check_idx[which.max(fill_dists)]
    training_index <- c(training_index, to_add)
    if(plot_it){
      plot(full_region[training_index[1:n_init],],xlim=c(0,200),ylim=c(0,200))
      points(full_region[training_index[(n_init+1):length(training_index)],],pch=20)
      abline(100,-1,lty=2,col='red')
      abline(50,-1,col='blue')
      abline(200,-1,col='blue')
    }
  }
  close(pb)
  #const_ind <- which(X_train_2d[,1] < 2 | X_train_2d[,2] < 2)
  #training_index <- setdiff(training_index, const_ind)
  X_train <- full_region[training_index,]
  
  
  G_entries_train <- matrix(NA, nrow(X_train), 600)
  print('Simulating Training Inputs...')
  pb = txtProgressBar(min=0,max=nrow(X_train),style=3,initial=0)
  for(i in 1:nrow(X_train)){
    setTxtProgressBar(pb,i)
    #print(c(i, X_train[i,]))
    G <- calculateG(X_train[i,1],X_train[i,2])
    G_entries_train[i,1:200] <- G[,1,1]
    G_entries_train[i,201:400] <- G[,1,2]
    G_entries_train[i,401:600] <- G[,2,2]
  }
  close(pb)
  G_det_inv_train <- 1/(G_entries_train[,1:200]*G_entries_train[,401:600]-G_entries_train[,201:400]^2)
  
  return(list(
    input = X_train,
    G_entries_train = G_entries_train,
    G_det_inv_train = G_det_inv_train
  ))
  
}


# Vectorized solving of blend composition
composition_solver_solve_fast <- function(Na, r, alpha_a, alpha_b,
                                          phi_p, phi_c, fA){
  phi_return <- rep(0,5)
  solved_compute_phi <- compute_phi(Na, fA, r, alpha_a, alpha_b, phi_p)
  composition <- solved_compute_phi$sol
  charge <- solved_compute_phi$charge
  phi_return[1] <- composition[1] #A
  phi_return[2] <- composition[2] #B
  if(phi_c > 0){
    phi_return[3] <- phi_c/2
    phi_return[4] <- phi_c/2
  }
  if(charge > 0) phi_return[3] <- phi_return[3] + composition[3]
  else if(charge < 0) phi_return[4] <- phi_return[4] + composition[3]
  
  small_idx <- phi_return[1:4] < 1e-10
  if(any(small_idx)){
    phi_return[small_idx] <- 0
  }
  
  phi_return[5] <- 1 - phi_p - phi_c
  
  return(phi_return)
}

# Vectorized calculation of U(k)
return_interaction_matrix_fast <- function(interactions, karray, BeadTypes){
  k2 <- karray^2
  n_beads <- nrow(BeadTypes)
  species_all <- sapply(1:n_beads, function(i) getspecies(BeadTypes, i))
  charge_all <- sapply(1:n_beads, function(i) getcharge(BeadTypes, i))
  
  G2_lookup <- matrix(NA, length(interactions),length(k2))
  int_types <- rep(NA, length(interactions))
  chi_lookup <- born_lookup <- elec_lookup <- zeta_lookup <- list()
  
  for(i in 1:length(interactions)){
    G2 <- tryCatch(Gamma2(k2, interactions[[i]][['smear_length']]))
    if(class(G2) == 'try-error') G2<-1
    G2_lookup[i,] <- G2
    
    int_here <- interactions[[i]]
    int_type_here <- int_here$type
    interaction <- list()
    
    if(int_type_here == 'chi'){
      int_types[i] <- 'chi'
      interaction$zero <- int_here$zero
      interaction$one <- int_here$one
      interaction$value <- int_here$value * G2
      chi_lookup[[length(chi_lookup) + 1]] <- interaction
      
    } else if(int_type_here == 'born'){
      int_types[i] <- 'born'
      interaction$zero <- int_here$zero
      interaction$one <- int_here$one
      interaction$value <- BornSolvation(int_here$value) * G2
      born_lookup[[length(born_lookup) + 1]] <- interaction
      
    } else if(int_type_here == 'electrostatic'){
      int_types[i] <- 'electrostatic'
      e_base <- 4*pi*int_here$value / k2
      elec_lookup[[length(elec_lookup) + 1]] <- list(e_base = e_base, G2 = G2)
      
    } else if(int_type_here == 'zeta'){
      int_types[i] <- 'zeta'
      zeta_lookup[[length(zeta_lookup) + 1]] <- list(value = int_here$value * G2)
    }
  }
  
  interaction_matrix <- array(0, dim = c(length(karray),n_beads,n_beads))
  for(i in 1:n_beads){
    species_i <- species_all[i]
    charge_i <- charge_all[i]
    for(j in 1:n_beads){
      species_j <- species_all[j]
      charge_j <- charge_all[j]
      
      for(chi in chi_lookup){
        if((species_i == chi$zero & species_j == chi$one) ||
           (species_i == chi$one & species_j == chi$zero)){
          interaction_matrix[,i,j] <- interaction_matrix[,i,j] + chi$value
        }
      }
      
      for(born in born_lookup){
        if((species_i == born$zero & species_j == born$one) ||
           (species_i == born$one & species_j == born$zero)){
          interaction_matrix[,i,j] <- interaction_matrix[,i,j] + born$value
        }
      }
      
      for(elec in elec_lookup){
        if(charge_i != 0 & charge_j != 0){
          netcharge <- charge_i * charge_j
          interaction_matrix[,i,j] <- interaction_matrix[,i,j] + elec$e_base * netcharge * elec$G2
        }
      }
      
      for(zeta in zeta_lookup){
        interaction_matrix[,i,j] <- interaction_matrix[,i,j] + zeta$value
      }
      
    }
  }
  return(interaction_matrix)
}


get_U_and_prefactors_fast <- function(r, alpha_a, alpha_b, f_A, f_p, rho_s, chiABN, chiASN,
                                      chiRatio, lb, epsilonRatio, a_plus, a_minus){
  
  N <- 100
  blist <- rep(1,5)
  Nref <- 1
  zeta <- 0.0
  a_i <- 0.0
  kmin <- -2
  kmax <- 1
  kpoints <- 200
  kinfo <- list(kmin = kmin, kmax = kmax, kpoints = kpoints)
  project_down <- F
  
  chain_details <- create_chain_details(N, r, alpha_a, alpha_b, blist, Nref)
  
  phi_list <- rep(0.2, 5) # default phi_list
  karray <- logspace(kmin, kmax, kpoints)
  
  read_chain_details_out <- read_chain_details(chain_details)
  
  # Access read_chain_details results
  Na <- read_chain_details_out$Na
  Na_seq <- read_chain_details_out$Na_seq
  Nb <- read_chain_details_out$Nb
  Nb_seq <- read_chain_details_out$Nb_seq
  blist_sys <- read_chain_details_out$blist
  Nref_sys <- read_chain_details_out$Nref
  
  Glist <- InitializeChains(Na, Na_seq$charged, Nb, Nb_seq$charged)
  
  Nlist <- rep(0, length(Glist))
  for(i in 1:length(Glist)) Nlist[i] <- length(Glist[[i]])
  #print(length(Glist))
  BeadTypes <- sort_charge(Glist, species_to_index, charge_to_index)
  #print(BeadTypes)
  
  rho_p <- (1-rho_s)*f_p
  rho_c <- (1-rho_s)*(1-f_p)
  phi_list <- composition_solver_solve_fast(Na, r, alpha_a, alpha_b, rho_p, rho_c, f_A)
  
  prefactors <- rep(NA, 5)
  for(i in 1:5) prefactors[i] <- Nlist[i]*phi_list[i]/Nref_sys/length(Glist[[i]])^2
  
  #G <- array(0, dim=c(200,7,7))
  NA_charge <- round(Na * alpha_a)
  NA_uncharge <- Na - NA_charge
  NB_charge <- round(Nb * alpha_b)
  NB_uncharge <- Nb - NB_charge
  #print(c(NA_charge, NA_uncharge, NB_charge, NB_uncharge))
  
  updated_phi <- chain_tracker_update_phi(phi_list, Glist)
  phi_list <- updated_phi$phi_list
  Glist <- list()
  if(length(updated_phi$Glist) >= 1) Glist$A <- updated_phi$Glist[[1]]$A
  if(length(updated_phi$Glist) >= 2) Glist$B <- updated_phi$Glist[[2]]$B
  if(length(updated_phi$Glist) >= 3) Glist$saltPlus <- updated_phi$Glist[[3]]$saltPlus
  if(length(updated_phi$Glist) >= 4) Glist$saltMinus <- updated_phi$Glist[[4]]$saltMinus
  if(length(updated_phi$Glist) >= 5) Glist$solvent <- updated_phi$Glist[[5]]$solvent
  # Glist$A <- updated_phi$Glist[[1]]$A
  # Glist$B <- updated_phi$Glist[[2]]$B
  # Glist$saltPlus <- updated_phi$Glist[[3]]$saltPlus
  # Glist$saltMinus <- updated_phi$Glist[[4]]$saltMinus
  # Glist$solvent <- updated_phi$Glist[[5]]$solvent
  remove_ind <- c()
  for(i in 1:5){
    if(is.null(Glist[i][[1]])) remove_ind <- c(remove_ind, i)
  }
  #print(remove_ind)
  if(!is.null(remove_ind)) Glist <- Glist[-c(remove_ind)]
  BeadTypes <- sort_charge(Glist, species_to_index, charge_to_index)
  # phi_list = ChainTracker.return_list(phi_list)
  
  interactions_list <- create_interactions(chiABN/N, chiASN/N, chiRatio, epsilonRatio,
                                           lb, a_plus, a_minus, zeta, a_i)
  
  
  U <- return_interaction_matrix_fast(interactions_list, karray, BeadTypes)
  
  return(list(U_mat = U, prefactors = prefactors))
}


gen_all_U_and_prefactors_vec <- function(data, batch_size = 50){
  n_data <- nrow(data)
  
  # if(n_data %% batch_size != 0) stop('batch_size should evenly divide nrow(data)')
  data <- matrix(data, ncol = 13)
  if(nrow(data) == 1){
    results <- get_U_and_prefactors_fast(
      r = data[1],
      alpha_a = data[2],
      alpha_b = data[3],
      f_A = data[4],
      f_p = data[5],
      rho_s = data[6],
      chiABN = data[7],
      chiASN = data[8],
      chiRatio = data[9],
      epsilonRatio = data[10],
      lb = data[11],
      a_plus = data[12],
      a_minus = data[13]
    )
    
    all_pre <- results$prefactors
    U_size <- dim(results$U_mat)[2]
    all_U <- array(NA, dim = c(200, 7, 7))
    all_U[,1:U_size,1:U_size] <- results$U_mat
    
    
  } else {
    
    all_U <- array(NA, dim=c(n_data,200,7,7))
    all_pre <- matrix(NA, n_data, 5)
    
    # pb <- txtProgressBar(min=0, max = num_batches, style = 3)
    
    for(i in 1:n_data){
     
      U_and_pre_here <- get_U_and_prefactors_fast(
                              r = data[i,1],
                              alpha_a = data[i,2],
                              alpha_b = data[i,3],
                              f_A = data[i,4],
                              f_p = data[i,5],
                              rho_s = data[i,6],
                              chiABN = data[i,7],
                              chiASN = data[i,8],
                              chiRatio = data[i,9],
                              epsilonRatio = data[i,10],
                              lb = data[i,11],
                              a_plus = data[i,12],
                              a_minus = data[i,13]
                            )
      
      # U_here <- results[[j]]$U_mat
      pre_here <- U_and_pre_here$prefactors
      U_size <- dim(U_and_pre_here$U_mat)[2]
      #print('made U_here, pre_here, U_size')
      all_pre[i, ] <- pre_here
      #print('updated all_pre')
      if(U_size == 7) {
        all_U[i, , , ] <- U_and_pre_here$U_mat
      } else {
        #print('time for for loop')
        for(k in 1:200) {
          all_U[i, k, 1:U_size, 1:U_size] <- U_and_pre_here$U_mat[k, , ]
        }
        # temp <- results[[j]]$U_mat[k,,]
        # all_U[batch_idx[j], , 1:U_size, 1:U_size] <- temp
      }
      
      
      #print(Sys.time() - t)
      
      
      # rm(results)
      # gc(verbose=F)
      # setTxtProgressBar(pb, i)
      
    }
    # close(pb)
  }
  return(list(
    all_U = all_U,
    all_pre = all_pre
  ))
  
}

phase_from_G_U_pre_fast <- function(all_G11_A, all_G12_A, all_G22_A, #all_G_det_inv_A, ##no need to have all_G_det_inv_A predict
                                  all_G11_B, all_G12_B, all_G22_B, #all_G_det_inv_B,
                                  all_U, all_prefactors){
  all_G11_A <- matrix(all_G11_A, ncol = 200)
  all_G12_A <- matrix(all_G12_A, ncol = 200)
  all_G22_A <- matrix(all_G22_A, ncol = 200)
  
  all_G11_B <- matrix(all_G11_B, ncol = 200)
  all_G12_B <- matrix(all_G12_B, ncol = 200)
  all_G22_B <- matrix(all_G22_B, ncol = 200)
  
  all_prefactors <- matrix(all_prefactors, ncol = 5)
  
  N <- if(is.null(nrow(all_G11_A))) 1 else nrow(all_G11_A)
  
  if(length(dim(all_U)) == 3){
    all_U <- array(all_U, dim = c(1, 200, 7, 7))
  }
  
  phase_record_out <- rep(-1, N) # -1 will correspond to the phase not being computable. Not important v.s. NA, just easier to deal with accuracy. can just do mean(y_test == pred_phase)
                                 # or mean(y_test == pred_phase_corrected) instead of needing to deal with removing NA values.  Can also find in the vector using which(pred_phase == -1)
  det_S_inv_record_out <- matrix(NA, N, 200)
  
  ##MG: suppose we use inverse det here
  #all_G_det_inv_A_here=1/(all_G11_A*all_G22_A-all_G12_A^2)
  #all_G_det_inv_B_here=1/(all_G11_B*all_G22_B-all_G12_B^2)
  all_G_det_inv_A=(1/(all_G11_A*all_G22_A-all_G12_A^2)) ##when this is smaller than zero, this could be problematic
  all_G_det_inv_B=(1/(all_G11_B*all_G22_B-all_G12_B^2))
  
  index_det_small_zero_A=which(apply(all_G_det_inv_A < 0, 1, FUN = sum) > 0) ##only look at the first one 
  index_det_small_zero_B=which(apply(all_G_det_inv_B < 0, 1, FUN = sum) > 0) ##only look at the first one 
  
  index_det_small_zero=union(index_det_small_zero_A,index_det_small_zero_B)
  index_det_not_small <- if(length(index_det_small_zero) > 0) (1:N)[-index_det_small_zero] else 1:N
  
  for(i in index_det_not_small){
    Ginv_here <- Ginv_from_preds_and_prefactors(all_G11_A[i,], all_G12_A[i,], all_G22_A[i,], all_G_det_inv_A[i,],
                                                all_G11_B[i,], all_G12_B[i,], all_G22_B[i,], all_G_det_inv_B[i,],
                                                all_prefactors[i,])
    Ginv_here <- slice_zeros(permute_G_array(Ginv_here))
    U_here <- all_U[i,,,]
    #print(dim(U_here))
    keep <- slice_na(U_here)
    U_here_no_NA <- U_here[,keep$rows,keep$cols]
    # print(dim(U_here_no_NA))
    # print(dim(Ginv_here))
    # print(all_prefactors[i,])
    # print('-----------------')
    det_S_inv_here <- determinant_Sinv(Ginv_here, U_here_no_NA)

    #det_S_inv_record_out[i,] <- det_S_inv_here
    phase <- checkphase(det_S_inv_here)
    phase_record_out[i] <- phase
    det_S_inv_record_out[i,] <- det_S_inv_here
  }
  
  return(
    list(phase_record = phase_record_out,
         det_S_inv_record = det_S_inv_record_out,
         index_det_small_zero=index_det_small_zero)
  )
}


get_det_S_inv_samples <- function(log_G11_A, log_G12_A, log_G22_A, #log_G_det_inv_A,
                                  log_G11_B, log_G12_B, log_G22_B, #log_G_det_inv_B,
                                  log_sd_G11_A, log_sd_G12_A, log_sd_G22_A, #log_sd_G_det_inv_A,
                                  log_sd_G11_B, log_sd_G12_B, log_sd_G22_B, #log_sd_G_det_inv_B,
                                  U, prefactors, #batch_size = 50,
                                  sample_size=100,k=200){
  
  normal_samples <- matrix(rnorm(sample_size), sample_size, k) # Perfectly Correlated Samples
  
  G11_A_sample <- exp(t(matrix(log_G11_A,k,sample_size)) + t(matrix(log_sd_G11_A,k,sample_size)) * normal_samples)
  G12_A_sample <- exp(t(matrix(log_G12_A,k,sample_size)) + t(matrix(log_sd_G12_A,k,sample_size)) * normal_samples)
  G22_A_sample <- exp(t(matrix(log_G22_A,k,sample_size)) + t(matrix(log_sd_G22_A,k,sample_size)) * normal_samples)
  G_det_inv_A_sample <- 1 / (G11_A_sample*G22_A_sample - G12_A_sample^2)
  
  G11_B_sample <- exp(t(matrix(log_G11_B,k,sample_size)) + t(matrix(log_sd_G11_B,k,sample_size)) * normal_samples)
  G12_B_sample <- exp(t(matrix(log_G12_B,k,sample_size)) + t(matrix(log_sd_G12_B,k,sample_size)) * normal_samples)
  G22_B_sample <- exp(t(matrix(log_G22_B,k,sample_size)) + t(matrix(log_sd_G22_B,k,sample_size)) * normal_samples)
  G_det_inv_B_sample <- 1 / (G11_B_sample*G22_B_sample - G12_B_sample^2)
  
  invalid_samples_A <- which(apply(G_det_inv_A_sample < 0, MARGIN = 1, FUN = sum) > 0)
  invalid_samples_B <- which(apply(G_det_inv_B_sample < 0, MARGIN = 1, FUN = sum) > 0)
  
  invalid_samples <- union(invalid_samples_A, invalid_samples_B)
  
  valid_samples <- if(length(invalid_samples) > 0) (1:sample_size)[-invalid_samples] else 1:sample_size
  
  phase_sample_record <- rep(-1, sample_size)
  det_S_inv_sample_record <- matrix(-1, sample_size, k)
  
  keep <- slice_na(U)
  U_here_no_NA <- U[,keep$rows,keep$cols]
  
  for(valid_idx in valid_samples){
    Ginv_here <- Ginv_from_preds_and_prefactors(G11_A_sample[valid_idx,], G12_A_sample[valid_idx,], G22_A_sample[valid_idx,], G_det_inv_A_sample[valid_idx,],
                                                G11_B_sample[valid_idx,], G12_B_sample[valid_idx,], G22_B_sample[valid_idx,], G_det_inv_B_sample[valid_idx,],
                                                prefactors)
    
    Ginv_here <- slice_zeros(permute_G_array(Ginv_here))
    
    det_S_inv_sample_record[valid_idx,] <- determinant_Sinv(Ginv_here, U_here_no_NA)
    phase_sample_record[valid_idx] <- checkphase(det_S_inv_sample_record[valid_idx,])
    
  }
  
  return(list(det_S_inv_sample_record=det_S_inv_sample_record,
              phase_sample_record=phase_sample_record))
  
}



predict_phase_w_correction <- function(X_test, ppgp_obj, correction_sample_size = 100){
  
  X_test <- as.matrix(X_test)
  # print(X_test[,1])
  
  N_A <- as.numeric(100)
  N_B <- round(as.numeric(X_test[,1]) * N_A)
  
  N_Ac <- round(X_test[,2] * N_A)
  N_Au <- N_A - N_Ac
  N_Bc <- round(X_test[,3] * N_B)
  N_Bu <- N_B - N_Bc
  
  # There are some edge cases we cannot handle:
  No_Charge_A <- which(N_Ac == 0)
  All_Charge_A <- which(N_Au == 0)
  No_Charge_B <- which(N_Bc == 0)
  All_Charge_B <- which(N_Bu == 0)
  
  if(length(No_Charge_A) > 0 | length(No_Charge_B) > 0 | length(All_Charge_A) > 0 | length(All_Charge_B) > 0){
    stop('Ensure that each Polymer A and B have >0 charged beads and >0 uncharged beads')
  }
  
  X_test_A <- cbind(N_Ac, N_Au)
  X_test_B <- cbind(N_Bc, N_Bu)
  
  preds_entries_A <- predict(ppgp_obj, X_test_A)
  preds_entries_B <- predict(ppgp_obj, X_test_B)
  
  preds_entries_A_mean <- exp(preds_entries_A$mean)
  preds_entries_B_mean <- exp(preds_entries_B$mean)
  
  all_U_and_pre <- gen_all_U_and_prefactors_vec(X_test)
  
  out <- phase_from_G_U_pre_fast(preds_entries_A_mean[,1:200], preds_entries_A_mean[,201:400], preds_entries_A_mean[,401:600],
                                 preds_entries_B_mean[,1:200], preds_entries_B_mean[,201:400], preds_entries_B_mean[,401:600],
                                 all_U_and_pre$all_U, all_U_and_pre$all_pre)
  
  pred_det_S_inv <- out$det_S_inv_record
  
  pred_phase <- out$phase_record
  pred_phase_corrected <- out$phase_record
  
  pred_phase_sample_corrected <- matrix(-1, nrow(X_test), correction_sample_size)
  pred_det_S_inv_sample_corrected <- array(NA, dim = c(nrow(X_test), correction_sample_size, 200))
  
  if(length(out$index_det_small_zero) > 0){ # These are invalid samples and need to be corrected
    
    # Every time det(1/\Gamma_A) or det(1/\Gamma_B) have at least one negative value,
    # we sample to correct the prediction by using only valid samples
    for(invalid_idx in out$index_det_small_zero){
      
      det_S_inv_samples_list <- get_det_S_inv_samples(
        log_G11_A = preds_entries_A$mean[invalid_idx,1:200],
        log_G12_A = preds_entries_A$mean[invalid_idx,201:400],
        log_G22_A = preds_entries_A$mean[invalid_idx,401:600],
        log_G11_B = preds_entries_B$mean[invalid_idx,1:200],
        log_G12_B = preds_entries_B$mean[invalid_idx,201:400],
        log_G22_B = preds_entries_B$mean[invalid_idx,401:600],
        log_sd_G11_A = preds_entries_A$sd[invalid_idx,1:200],
        log_sd_G12_A = preds_entries_A$sd[invalid_idx,201:400],
        log_sd_G22_A = preds_entries_A$sd[invalid_idx,401:600],
        log_sd_G11_B = preds_entries_B$sd[invalid_idx,1:200],
        log_sd_G12_B = preds_entries_B$sd[invalid_idx,201:400],
        log_sd_G22_B = preds_entries_B$sd[invalid_idx,401:600],
        U = all_U_and_pre$all_U[invalid_idx,,,],
        prefactors = all_U_and_pre$all_pre[invalid_idx,],
        sample_size = correction_sample_size,
        k = 200
      )
      
      pred_phase_sample_corrected[invalid_idx,] <- det_S_inv_samples_list$phase_sample_record
      pred_det_S_inv_sample_corrected[invalid_idx,,] <- det_S_inv_samples_list$det_S_inv_sample_record
      
      num_phases <- c(
        length(which(pred_phase_sample_corrected[invalid_idx,] == 0)),
        length(which(pred_phase_sample_corrected[invalid_idx,] == 1)),
        length(which(pred_phase_sample_corrected[invalid_idx,] == 2))
      )
      
      num_cant_predict <- length(which(pred_phase_sample_corrected == -1))
      
      pred_phase_corrected[invalid_idx] <- if(max(num_phases) > 0) which.max(num_phases) - 1 else -1 # 1,2,3 -> 0,1,2
    }
    
  }
  
  return(list(
    pred_det_S_inv_record = pred_det_S_inv,
    pred_det_S_inv_sample_corrected = pred_det_S_inv_sample_corrected,
    pred_phase = pred_phase,
    pred_phase_corrected = pred_phase_corrected,
    pred_phase_sample_record = pred_phase_sample_corrected
  ))
  
}


full_sim_process_w_fast <- function(data){
  r <- data[,1]
  alpha_a <- data[,2]
  alpha_b <- data[,3]
  f_A <- data[,4]
  f_p <- data[,5]
  rho_s <- data[,6]
  chiABN <- data[,7]
  chiASN <- data[,8]
  chiRatio <- data[,9]
  epsilonRatio <- data[,10]
  lb <- data[,11]
  a_plus <- data[,12]
  a_minus <- data[,13]
  
  n_data <- nrow(data)
  
  N <- 100
  blist <- rep(1,5)
  Nref <- 1
  zeta <- 0.0
  a_i <- 0.0
  kmin <- -2
  kmax <- 1
  kpoints <- 200
  kinfo <- list(kmin = kmin, kmax = kmax, kpoints = kpoints)
  project_down <- T
  
  N_A <- 100
  N_B <- round(r * N_A)
  N_Ac <- round(N_A*alpha_a)
  N_Bc <- round(N_B*alpha_b)
  N_Au <- N_A - N_Ac
  N_Bu <- N_B - N_Bc
  
  all_G_A <- all_G_B <- array(NA, dim = c(n_data, 200, 2, 2))
  print('Generating G')
  pb <- txtProgressBar(min=0,max=n_data,style=3)
  for(i in 1:n_data){
    setTxtProgressBar(pb, i)
    G_A <- calculateG(N_Ac[i], N_Au[i])
    G_B <- calculateG(N_Bc[i], N_Bu[i])
    all_G_A[i,,,] <- G_A
    all_G_B[i,,,] <- G_B
  }
  close(pb)
  print('Generating U')
  all_U_and_prefactors <- gen_all_U_and_prefactors_vec(data)
  
  all_prefactors <- all_U_and_prefactors$all_pre
  all_U <- all_U_and_prefactors$all_U
  
  all_G_A_det_inv <- 1 / (all_G_A[,,1,1]*all_G_A[,,2,2] - all_G_A[,,1,2]^2)
  all_G_B_det_inv <- 1 / (all_G_B[,,1,1]*all_G_B[,,2,2] - all_G_B[,,1,2]^2)
  
  sim_out <- phase_from_G_U_pre_fast(all_G_A[,,1,1], all_G_A[,,1,2],all_G_A[,,2,2],
                                all_G_B[,,1,1], all_G_B[,,1,2],all_G_B[,,2,2],
                                all_U,all_prefactors)
  
  return(list(
    BeadCounts = c(N_A, N_Ac, N_Au, N_B, N_Bc, N_Bu),
    G_A = all_G_A,
    G_B = all_G_B,
    #Ginv = Ginv,
    prefactors = all_prefactors,
    #phi_list = phi_list,
    inputs = data,
    U = all_U,
    #Sinv = Sinv,
    detSinv = sim_out$det_S_inv_record,
    phase = sim_out$phase_record
  ))
}



gen_phase_diagram_general <- function(ppgp_entries, N,
                                      grid_col1 = 'alpha_A', grid_col2 = 'alpha_B',
                                      grid_col1_lim = NULL, grid_col2_lim = NULL,
                                      match_alpha_A_B = F,
                                      salt_lb_match = F, mask_microphase = F,
                                      r = 1.0, alpha_A = 0.5, alpha_B = 0.5,
                                      f_A = 0.5, f_p = 1.0, rho_s = 0.0,
                                      chiABN = 50.0, chiASN = 0.0, chiRatio = 0.0,
                                      epsRatio = 1.0, lb = 2.0, a_plus = 1.0, a_minus = 1.0){
  
  # grid_col1, grid_col2 are the values to vary for the plot
  # Can change the bounds, they are only relevant when the corresponding parameter
  # is being varied.  The odd lb value is to be consistent with gamma range
  
  bounds <- list(
    r = c(0.5,2.0), alpha_A = c(0.02,0.98), alpha_B = c(0.02,0.98),
    f_A = c(0.01,0.99), f_p = c(0.01,1.0), rho_s = c(0.01,0.99),
    chiABN = c(0.0,250.0), chiASN = c(0.0,100.0), chiRatio = c(0.0,2.0),
    epsRatio = c(1.0,5.0), lb = c(0.001,5.0), a_plus = c(0.1,10), a_minus = c(0.1,10)
  )
  
  if(!is.null(grid_col1_lim) & length(grid_col1_lim) == 2) bounds[[grid_col1]] <- grid_col1_lim
  if(!is.null(grid_col2_lim) & length(grid_col2_lim) == 2) bounds[[grid_col2]] <- grid_col2_lim
  
  phase_diagram_data <- data.frame(
    r = rep(r, N^2), alpha_A = rep(alpha_A, N^2), alpha_B = rep(alpha_B, N^2),
    f_A = rep(f_A, N^2), f_p = rep(f_p, N^2), rho_s = rep(rho_s, N^2),
    chiABN = rep(chiABN, N^2), chiASN = rep(chiASN, N^2), chiRatio = rep(chiRatio, N^2),
    epsRatio = rep(epsRatio, N^2), lb = rep(lb, N^2), a_plus = rep(a_plus, N^2), a_minus = rep(a_minus, N^2)
  )
  
  grid <- expand.grid(
    x = seq(bounds[[grid_col1]][1], bounds[[grid_col1]][2], length.out = N),
    y = seq(bounds[[grid_col2]][1], bounds[[grid_col2]][2], length.out = N)
  )
  #print(c(unique(grid$x),unique(grid$y)))
  if(salt_lb_match) phase_diagram_data$a_plus <- phase_diagram_data$a_minus <- phase_diagram_data$lb / 6
  
  phase_diagram_data[,c(grid_col1,grid_col2)] <- grid
  if(match_alpha_A_B & length(unique(alpha_A)) > 0){ 
    phase_diagram_data$alpha_B <- phase_diagram_data$alpha_A
  } else if(match_alpha_A_B & length(unique(alpha_B)) > 0) phase_diagram_data$alpha_A <- phase_diagram_data$alpha_B
  #else phase_diagram_data$alpha_B <- phase_diagram_data$alpha_A
  N_A <- 100
  
  N_B <- round(N_A*phase_diagram_data$r)
  phase_diagram_X_A <- cbind(round(N_A*phase_diagram_data$alpha_A), N_A - round(N_A*phase_diagram_data$alpha_A))
  phase_diagram_X_B <- cbind(round(N_B*phase_diagram_data$alpha_B), N_B - round(N_B*phase_diagram_data$alpha_B))
  
  preds_entries_A_phase_diagram <- predict(ppgp_entries, phase_diagram_X_A)$mean %>% exp()
  preds_entries_B_phase_diagram <- predict(ppgp_entries, phase_diagram_X_B)$mean %>% exp()
  
  # U_and_pre_phase_diagram <- gen_all_U_and_prefactors_vec(phase_diagram_data, N^2, num_cores = 1)
  outputs <- predict_phase_w_correction(phase_diagram_data, ppgp_entries)
  phase_diagram_data$phase <- outputs$pred_phase_corrected
  
  centroids <- data.frame(x = rep(NA, length(unique(phase_diagram_data$phase))),
                          y = rep(NA, length(unique(phase_diagram_data$phase))),
                          phase = rep(NA, length(unique(phase_diagram_data$phase))))
  
  idx <- 1
  for(i in c(0,1,2)){
    if(!(i %in% unique(phase_diagram_data$phase))) next
    look_at <- phase_diagram_data[which(phase_diagram_data$phase == i),c(grid_col1,grid_col2)]
    centroid <- apply(look_at, 2, mean)
    
    centroids$x[idx] = centroid[1]
    centroids$y[idx] = centroid[2]
    centroids$phase[idx] = i
    idx <- idx + 1
    
  }
  
  phases <- c('Not Macro','Macro')
  if(mask_microphase) phase_diagram_data$phase[which(phase_diagram_data$phase == 2)] <- 0
  p <- ggplot(phase_diagram_data, aes(x = !!sym(grid_col1), y = !!sym(grid_col2))) +
    geom_point(aes(color = factor(phase))) +
    #geom_text(aes(x = x, y = y, label = phases[phase+1]), data = centroids) +
    xlim(min(bounds[[grid_col1]]),max(bounds[[grid_col1]])) +
    ylim(min(bounds[[grid_col2]]),max(bounds[[grid_col2]]))
  
  if(mask_microphase){ 
    p <- p + geom_contour(aes(z=(phase)),breaks = c(0.5))
  }else p <- p + geom_contour(aes(z = phase), breaks = c(0.5,1.5))
  
  
  return(list(phase_diagram = phase_diagram_data, plt = p))
  
}


gen_phase_diagram_set <- function(ppgp_entries, N,
                                  change_var, change_vals,
                                  grid_col1 = 'alpha_A', grid_col2 = 'alpha_B',
                                  grid_col1_lim = NULL, grid_col2_lim = NULL,
                                  match_alpha_A_B = F,
                                  salt_lb_match = F, do_sim = F, 
                                  num_cores = 1, mask_microphase = T,
                                  r = 1.0, alpha_A = 0.5, alpha_B = 0.5,
                                  f_A = 0.5, f_p = 1.0, rho_s = 0.0,
                                  chiABN = 50.0, chiASN = 0.0, chiRatio = 0.0,
                                  epsRatio = 1.0, lb = 2.0, a_plus = 1.0, a_minus = 1.0){
  
  
  n_diagrams <- length(change_vals)
  
  
  predicted_phase_diagrams <- list()
  simulated_phase_diagrams <- list()
  predicted_plots <- list()
  simulated_plots <- list()
  
  # Generating phase diagrams
  for(i in 1:n_diagrams){
    
    # Update Variable
    args <- list(
      ppgp_entries = ppgp_entries,
      N = N,
      grid_col1 = grid_col1,
      grid_col2 = grid_col2,
      grid_col1_lim = grid_col1_lim,
      grid_col2_lim = grid_col2_lim,
      match_alpha_A_B = match_alpha_A_B,
      salt_lb_match = salt_lb_match,
      mask_microphase = mask_microphase,
      r = r,
      alpha_A = alpha_A,
      alpha_B = alpha_B,
      rho_s = rho_s,
      f_A = f_A,
      f_p = f_p,
      chiABN = chiABN,
      chiASN = chiASN,
      chiRatio = chiRatio,
      epsRatio = epsRatio,
      lb = lb,
      a_plus = a_plus,
      a_minus = a_minus
    )
    args[[change_var]] <- change_vals[i]
    
    # Predicted Phase Diagram
    phase_diagram_i <- do.call(gen_phase_diagram_general, args)
    predicted_phase_diagrams[[i]] <- phase_diagram_i$phase_diagram
    
    # Simulated Phase Diagram if Desired
    if(do_sim){
      
      sim_data_in_i <- as.matrix(phase_diagram_i$phase_diagram[,-14])
      sim_data_i_allInfo <- full_sim_process_w_fast(sim_data_in_i)
      
      sim_data_i <- as.data.frame(sim_data_in_i)
      sim_data_i$phase <- sim_data_i_allInfo$phase
      if(mask_microphase) sim_data_i$phase[which(sim_data_i$phase == 2)] <- 0
      
      simulated_phase_diagrams[[i]] <- sim_data_i
    }
    
  }
  
  expressions_list <- list(
    r = expression(r),
    alpha_A = expression(alpha[A]),
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
  
  bounds <- list(
    r = c(0.5,2.0), alpha_A = c(0.02,0.98), alpha_B = c(0.02,0.98),
    f_A = c(0.01,0.99), f_p = c(0.01,1.0), rho_s = c(0.01,0.99),
    chiABN = c(0.0,250.0), chiASN = c(0.0,100.0), chiRatio = c(0.0,2.0),
    epsRatio = c(1.0,5.0), lb = c(0.001,5.0), a_plus = c(0.1,10), a_minus = c(0.1,10)
  )
  
  if(is.null(grid_col1_lim)){
    x0 <- bounds[[grid_col1]][1]
    x1 <- bounds[[grid_col1]][2]
  } else {
    x0 <- grid_col1_lim[1]
    x1 <- grid_col1_lim[2]
  }
  
  if(is.null(grid_col2_lim)){
    y0 <- bounds[[grid_col2]][1]
    y1 <- bounds[[grid_col2]][2]
  } else {
    y0 <- grid_col2_lim[1]
    y1 <- grid_col2_lim[2]
  }
  
  for(i in 1:n_diagrams){
    if(i == 1){
      if(do_sim) all_sim_data <- simulated_phase_diagrams[[i]]
      all_pred_data <- predicted_phase_diagrams[[i]]
    } else {
      if(do_sim) all_sim_data <- rbind(all_sim_data, simulated_phase_diagrams[[i]])
      all_pred_data <- rbind(all_pred_data, predicted_phase_diagrams[[i]])
    }
  }
  
  
  
  
  if(do_sim){
    all_pred_data$label <- 'Prediction'
    all_sim_data$label <- 'Simulation'
    
    phase_bdry_plt <- ggplot(rbind(all_pred_data,all_sim_data), aes(x = !!sym(grid_col1), y = !!sym(grid_col2))) +
      stat_contour(aes(z = phase, color = factor(!!sym(change_var)), linetype = label), breaks = c(0.5)) +
      labs(x = expressions_list[[grid_col1]], y = expressions_list[[grid_col2]],
           color = expressions_list[[change_var]], linetype = 'Source') +
      xlim(x0,x1) + ylim(y0,y1)
    
    return(list(
      phase_boundary_plt=phase_bdry_plt,
      all_pred_data = all_pred_data,
      all_sim_data = all_sim_data,
      predicted_phase_diagrams=predicted_phase_diagrams,
      #predicted_plots=predicted_plots,
      simulated_phase_diagrams=simulated_phase_diagrams
      #simulated_plots=simulated_plots#,
      #contours=contours
    ))
    
  } else {
    
    phase_bdry_plt <- ggplot(all_pred_data, aes(x = !!sym(grid_col1), y = !!sym(grid_col2))) +
      stat_contour(aes(z = phase, color = factor(!!sym(change_var))), breaks = c(0.5)) +
      labs(x = expressions_list[[grid_col1]], y = expressions_list[[grid_col2]],
           color = expressions_list[[change_var]]) +
      xlim(x0,x1) + ylim(y0,y1)
    
    return(list(
      phase_boundary_plt=phase_bdry_plt,
      all_pred_data = all_pred_data,
      predicted_phase_diagrams=predicted_phase_diagrams
      #predicted_plots=predicted_plots
      #contours=contours
    ))
  }
  
}














