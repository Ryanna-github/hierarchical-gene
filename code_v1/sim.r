# indx=c(1, 2), then return column vector c(1,-1,0,...0)
dMatrixFun <- function(indx){
  e.vec<-matrix(0,K,1)
  e.vec[indx[1],]<-1
  e.vec[indx[2],]<-(-1)
  return(e.vec)
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
  return(list(ci_main=ci_main_list, ci_sub=ci_sub_list))
}

# Generating the hierarchical regression coefficients
generate_coef_true <-function(n, p, q, ci_main, ci_sub, coef_main_value, coef_sub_value, incl_inter = TRUE,
                              len_beta = 2, len_alpha = 2, len_gamma = 4){
  
  beta <-  matrix(0, nrow = p, ncol = length(ci_sub))
  alpha <-  matrix(0, nrow = q, ncol = length(ci_sub))
  gamma <-  matrix(0, nrow = pq, ncol = length(ci_sub))
  eta <- matrix(0, nrow = pq, ncol = length(ci_sub))
  sub_num <- length(ci$ci_sub)/length(ci$ci_main) # subgroup num within a given group
  
  for (i in 1:length(ci_main)) {
    # 可以多列同时填充，注意不是整列填充，只填充非0项
    beta[1:len_beta, (1+(i-1)*sub_num):(i*sub_num)] = rep(coef_main_value[i], len_beta)
    alpha[1:len_alpha, (1+(i-1)*sub_num):(i*sub_num)] = rep(coef_main_value[i], len_alpha)
  }
  if(incl_inter){
    for (i in 1:length(ci_sub)) {
      gamma[1:len_gamma, i] = rep(coef_sub_value[i], len_gamma)
      eta[, i] = kronecker(matrix(1, nrow=q, ncol=1), beta[, i]) * gamma[, i]
    }
  }
  
  
  beta_full <- matrix(0, nrow = p, ncol = n)
  alpha_full <- matrix(0, nrow = q, ncol = n)
  gamma_full <- matrix(0, nrow = pq, ncol = n)
  eta_full <- matrix(0, nrow = pq, ncol = n)
  
  # 多行填充时会默认按列填充
  for(i in 1:length(ci_main)){
    beta_full[,unlist(ci_main[[i]])] <- beta[,i*sub_num]
    alpha_full[,unlist(ci_main[[i]])] <- alpha[,i*sub_num]
  }
  if(incl_inter){
    for(i in 1:length(ci_sub)){
      gamma_full[,unlist(ci_sub[[i]])] <- gamma[,i]
      eta_full[,unlist(ci_sub[[i]])] <- eta[,i]
    }
  }
  
  coef_true_K <- do.call(rbind, list(beta, alpha, gamma))
  rownames(coef_true_K) <- paste0("X.", 1:(p+q+pq))
  colnames(coef_true_K) <- paste0("Comp.", 1:length(ci_sub))
  return(list(coef_true_K = coef_true_K,
              beta_true_K = beta,
              alpha_true_K = alpha, 
              gamma_true_K = gamma,
              eta_true_K = eta,
              beta_true_mat = beta_full,
              alpha_true_mat = alpha_full,
              gamma_true_mat = gamma_full,
              eta_true_mat = eta_full
              ))
}


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

# Generating the simulated data
generate_sim_data <- function(n, p, q, coef_true, cotype_x = "AR_E", cotype_z = "AR_E", 
                              corZX = F, sd_error=0.5, discrete_ratio = 0){
  
  sigma_x <- generate_sigma(p, cotype_x)
  sigma_z <- generate_sigma(q, cotype_z)
  sigma_main <- as.matrix(Matrix::bdiag(sigma_x, sigma_z))
  
  if(corZX){
    sigma_main[sigma_main == 0] <- min(abs(sigma_main[sigma_main>0]))/20
  }
  
  XZ <- MASS::mvrnorm(n, rep(0,p+q), sigma_main)
  X <- XZ[,1:p]
  Z <- XZ[,(p+1):(p+q)]
  
  if(discrete_ratio > 0){
    discrete_x <- ceiling(discrete_ratio * p)
    discrete_z <- ceiling(discrete_ratio * q)
    X[,1:discrete_x] <- as.numeric(X[,1:discrete_x] > 0 )
    Z[,1:discrete_z] <- as.numeric(Z[,1:discrete_z] > 0 )
  }
  W <- matrix(0, n, pq)
  for(i in 1:n){
    W[i,] <- kronecker(Z[i,], X[i,])
  }
  
  # error <- rnorm(n, 0, sd_error)
  Iq <- matrix(1, q, 1)
  data_X <- bdiag(lapply(split(X, 1:n), function(x) t(x)))
  data_Z <- bdiag(lapply(split(Z, 1:n), function(x) t(x)))
  data_W <- bdiag(lapply(split(W, 1:n), function(x) t(x)))
  coef_beta <- as.vector(coef_true$beta_true_mat)
  coef_alpha <- as.vector(coef_true$alpha_true_mat)
  coef_gamma <- as.vector(coef_true$gamma_true_mat)
  # coef_eta <- kronecker(Iq, coef_beta) * coef_gamma
  coef_eta <- as.vector(coef_true$eta_true_mat)
  
  # data_y <- data_X %*% coef_beta + data_Z %*% coef_alpha + data_W %*% coef_eta + error
  data_y <- data_X %*% coef_beta + data_Z %*% coef_alpha + data_W %*% coef_eta
  # data_y <- data_X %*% coef_beta + data_Z %*% coef_alpha
  data_y <- matrix(data_y)
  return(list(data_y = data_y,
              X_all = cbind(X, Z, W),
              X_main = cbind(X, Z),
              X = X,
              Z = Z,
              W = W,
              data_X = data_X,
              data_Z = data_Z,
              data_W = data_W,
              # data_e = error,
              data_X_all = cbind(data_X, data_Z, data_W),
              coef_true = coef_true$coef_true_K))
  
}