# Generating the correlation structure.
generate_sigma <- function(m, type.index = "En"){
  sigma <- matrix(0, m, m)
  if (type.index != "En"){
    for (i in 1:m){
      for(j in i:m){
        if(type.index == "AR1"){
          sigma[i,j] = 0.25^abs(i-j)
        }else if(type.index == "AR2"){
          sigma[i,j] = 0.75^abs(i-j)
        }else if(type.index == "B2"){
          if(abs(i-j) ==1){
            sigma[i,j]= 0.6
          }else if (abs(i-j)==2){
            sigma[i,j] = 0.3
          }
        }else if(type.index == "B1"){
          if(abs(i-j) ==1){
            sigma[i,j]= 0.3
          }
        }else if(type.index == "AR_E"){
          sigma[i,j]=0.2
        }
        sigma[j,i] = sigma[i,j]
      }
    }
  }
  diag(sigma)= 1
  return(sigma)
}


# input: p --> X dim, q --> Z dim
#         n --> sample size
#         covX --> X covariance type; covZ --> Z covariance type
# output: X, Z, data(combination of X & Z) n*(p+q)
generate_data_beta_alpha <- function(n, p, q, cotype_x = "En", cotype_z = "En"){
  sigma_x <- generate_sigma(p, cotype_x)
  sigma_z <- generate_sigma(q, cotype_z)
  # sigma_main <- as.matrix(Matrix::bdiag(sigma_x, sigma_z))
  X <- MASS::mvrnorm(n, rep(0,p), sigma_x)
  Z <- MASS::mvrnorm(n, rep(0,q), sigma_z)
  return(list("X" = X, "Z" = Z, "data_full" = cbind(X, Z)))
}

# generate true coef
# input: beta_nonzero, beta_vlen --> beta = first beta_vlen row (beta_nonzero) else 0
generate_coef <- function(p, q, beta_nonzero, alpha_nonzero, beta_vlen, alpha_vlen,
                          reverse = FALSE){
  K_up <- length(beta_nonzero)
  # beta
  coef_beta <- matrix(0, nrow = p, ncol = K_up)
  coef_beta[1:beta_vlen,] <- t(matrix(rep(beta_nonzero, beta_vlen), nrow = K_up))
  # alpha
  coef_alpha <- matrix(0, nrow = q, ncol = K_up)
  coef_alpha[1:alpha_vlen,] <- t(matrix(rep(alpha_nonzero, alpha_vlen), nrow = K_up))
  if(reverse){
    coef_beta <- coef_beta[nrow(coef_beta):1,]
    coef_alpha <- coef_alpha[nrow(coef_alpha):1,]
  }
  return(list("coef_beta" = coef_beta, 
              "coef_alpha" = coef_alpha, 
              "coef_full" = rbind(coef_beta, coef_alpha)))
}

# Generating the hierarchical subgroup structure.
ci_generate <- function(n, pr_sub = pr_sub, hier_struc = hier_struc){
  ci_sub_list=list()
  n_index = sample(n, n, replace = F)
  ci_sub_list[[1]] = sort(n_index[1:ceiling(n*pr_sub[1])])
  for (g in 2:(length(pr_sub))) {
    ci_sub_list[[g]] = sort(n_index[(ceiling(n*sum(pr_sub[1:(g-1)]))+1):ceiling(n*sum(pr_sub[1:g]))])
  }
  ci_main_list=list()
  for (g in 1:(length(hier_struc))) {
    ci_set = unlist(hier_struc[g])
    ci_sub = NULL
    for (gg in 1:length(ci_set)) {
      ci_sub = c(ci_sub,ci_sub_list[[ci_set[gg]]])
    }
    ci_main_list[[g]] = sort(ci_sub)
  }
  ci_sim <- rep(0, n)
  for(i in 1:length(pr_sub)){
    ci_sim[ci_sub_list[[i]]] <- i
  }
  ci_sim_main <- rep(0, n)
  ci_sim_sub <- rep(0, n)
  for(k in 1:length(hier_struc)){
    ci_sim_main[ci_main_list[[k]]] <- k
  }
  for(k in 1:length(pr_sub)){
    ci_sim_sub[ci_sub_list[[k]]] <- k
  }
  return(list("ci_main" = ci_main_list, 
              "ci_sub" = ci_sub_list,
              "ci_sim" = ci_sim,
              "ci_sim_main" = ci_sim_main,
              "ci_sim_sub" = ci_sim_sub))
}

# generate_all_data
generate_all_data <- function(data_seed, n, p, q, prob_sub, hier_struc, 
                              beta_nonzero, alpha_nonzero, 
                              beta_vlen, alpha_vlen,
                              cotype_x, cotype_z, epsilon_sd,
                              reverse = FALSE){
  set.seed(data_seed)
  ci <- ci_generate(n, prob_sub, hier_struc)
  data <- generate_data_beta_alpha(n, p, q, cotype_x, cotype_z)
  coef <- generate_coef(p, q, beta_nonzero, alpha_nonzero, beta_vlen, alpha_vlen, reverse)
  ys <- data$data_full %*% coef$coef_full
  y <- NULL
  for(i in 1:n){
    y[i] <- ys[i,ci_sim[i]]
  }
  y <- y + rnorm(length(y), mean = 0, sd = epsilon_sd)
  # View(data.frame(cbind(ys, ci_sim, y))) # 验证正确性
  return(list("data" = data,
              "coef" = coef,
              "y" = y,
              "ci_sim" = ci$ci_sim,
              "ci_sim_main" = ci$ci_sim_main,
              "ci_sim_sub" = ci$ci_sim_sub))
}

# get true q_c matrix
get_q_c_matrix <- function(n, K_up, ci_sim){
  q_c_matrix <- matrix(0, nrow = n, ncol = K_up)
  for(i in 1:n){
    q_c_matrix[i,ci_sim[i]] <- 1
  }
  return(q_c_matrix)
}











