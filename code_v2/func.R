
flexmix_init <- function(q_c_seed){
  set.seed(q_c_seed)
  pi_init <- rep(1/K_up, K_up)
  rho_init <- c(1, 1, 1, 1)/0.5
  q_c_matrix <- abs(t(kronecker(pi_init, matrix(1, ncol = n))) + rnorm(n*K_up, mean = 0, sd = .1))
  q_c_matrix <- q_c_matrix / apply(q_c_matrix, 1, sum)
  
  result <- tryCatch({
    m_glm <- flexmix(y~cbind(X, Z)-1, cluster = q_c_matrix,
            model = FLXMRglm(),
            control = list(minprior = 0))
  }, error = function(err) {
    cat("Error occurred:", conditionMessage(err), "\n")
    m_glm <- flexmix(y~cbind(X, Z)-1, k = K_up,
                     model = FLXMRglm(),
                     control = list(minprior = 0))
  })
  coef_est <- parameters(m_glm)[1:(p+q),] *
    t(kronecker(rho_init, matrix(1, ncol = p+q)))
  row.names(coef_est) <- NULL
  cdist <- coef_dist(coef_est, coef$coef_full)
  sc_score <- sc(m_glm@cluster, ci_sim)
  return(list(cdist = cdist, 
              ci_prob_mean = 999,
              coef_full_ori = coef_est,
              mse = 999,
              sc_score = sc_score))
}

random_init <- function(q_c_seed){
  set.seed(q_c_seed)
  pi_init <- rep(1/K_up, K_up)
  rho_init <- c(1, 1, 1, 1)/0.5
  q_c_matrix <- abs(t(kronecker(pi_init, matrix(1, ncol = n))) + rnorm(n*K_up, mean = 0, sd = .1))
  q_c_matrix <- q_c_matrix / apply(q_c_matrix, 1, sum)
  coef_est <- coef$coef_full - coef$coef_full + rnorm(K_up*(p+q), 0, 1)
  cdist <- coef_dist(coef_est, coef$coef_full)
  print(cdist)
  ci_est <- apply(q_c_matrix, 1, which.max)
  ci_prob_mean <- mean(apply(q_c_matrix, 1, max))
  sc_score <- sc(ci_est, ci_sim)
  return(list(cdist = cdist,
              ci_prob_mean = ci_prob_mean,
              coef_full_ori = coef_est,
              mse = 111,
              sc_score = sc_score))
}

ADMM_trail <- function(aa, tau, lambda_1, lambda_2, lambda_3, q_c_seed, 
                       coef_full_init, eps =1e-7){
  tau <- ifelse(tau == 0, 1e-4, tau) # 防止 tau == 0 导致分母为0情况
  rho_init <- c(1, 1, 1, 1)/0.5
  
  kj <- function(dim) {return((k-1)*dim+j)}
  ks <- function(dim) {return((k-1)*dim+s)}
  rho_list <- list(rho_init)
  pi_init <- rep(1/K_up, K_up)
  pi_list <- list(pi_init)
  set.seed(q_c_seed)
  # q_c_matrix <- get_q_c_matrix(n, K_up, ci_sim)
  # 真值附近
  # q_c_matrix <- abs(get_q_c_matrix(n, K_up, ci_sim) + rnorm(n*K_up, mean = 1, sd = .1))
  # 全随机
  q_c_matrix <- abs(t(kronecker(pi_init, matrix(1, ncol = n))) + rnorm(n*K_up, mean = 0, sd = .1))
  q_c_matrix <- q_c_matrix / apply(q_c_matrix, 1, sum)
  
  # coef full
  # coef_full_init <- coef$coef_full - coef$coef_full
  
  # coef_full_init[(p+1):(p+q),] <- coef$coef_alpha
  # coef_full_init[1:p,] <- coef$coef_beta
  # coef_full_init <- coef$coef_full * t(kronecker(rho_init, matrix(1, ncol = p+q)))
  coef_full_list <- list(coef_full_init)
  coef_full_ori_list <- list(coef_full_init/t(kronecker(rho_init, 
                                                        matrix(1, ncol = p+q))))
  
  
  
  # coef diff (v,w):=u
  diff_v_init <- matrix(H_p %*% as.vector(coef_full_init[1:p,]), ncol = 1)
  # diff_v_init <- matrix(H_p %*% as.vector(coef$coef_beta), ncol = 1)
  diff_v_list <- list(diff_v_init)
  diff_w_init <- matrix(H_q %*% as.vector(coef_full_init[(p+1):(p+q),]), ncol = 1)
  # diff_w_init <- matrix(H_q %*% as.vector(coef$coef_alpha), ncol = 1)
  diff_w_list <- list(diff_w_init)
  
  # dual variable (xi,zeta):=eta
  dual_xi_init <- matrix(0, nrow = (p)*K_up*(K_up-1)/2, ncol = 1)
  dual_xi_list <- list(dual_xi_init)
  dual_zeta_init <- matrix(0, nrow = (q)*K_up*(K_up-1)/2, ncol = 1)
  dual_zeta_list <- list(dual_zeta_init)
  
  case <- NULL
  cu_tracker <- NULL
  cv_tracker <- NULL
  
  for(iter in 2:100){
    rho_est <- rho_list[[iter-1]]
    coef_full_est <- coef_full_list[[iter-1]]
    coef_beta_est <- coef_full_est[1:p,]
    coef_alpha_est <- coef_full_est[(p+1):(p+q),]
    
    diff_v_est <- diff_v_list[[iter-1]]
    diff_w_est <- diff_w_list[[iter-1]]
    
    dual_xi_est <- dual_xi_list[[iter-1]]
    dual_zeta_est <- dual_zeta_list[[iter-1]]
    
    # =========================== beta, alpha update =========================
    # beta,alpha （相比矩阵形式，只修改 beta,alpha 更新）
    for(k in 1:K_up){
      W_k <- diag(q_c_matrix[,k])
      for(j in 1:p){
        lower_1 <- t(X[,j]) %*% W_k %*% X[,j]
        lower_2 <- (t(H_p)%*%H_p)[kj(p), kj(p)]
        # lower_2 <- 0
        Lower <- lower_1/n + tau*lower_2
        upper_1 <- t(X[,j]) %*% W_k %*% (rho_est[k]*y - X[,-j]%*%(coef_beta_est[-j,k]) - Z%*%coef_alpha_est[,k])
        upper_2 <- (t(H_p)%*%(diff_v_est-dual_xi_est))[kj(p)] - 
          as.vector(coef_beta_est)[-kj(p)]%*%((t(H_p)%*%H_p)[kj(p),-kj(p)])
        # upper_2 <- 0
        Upper <- upper_1/n + tau*upper_2
        # coef_beta_kj <- Upper/Lower
        coef_beta_kj <- mcp_solution(Lower, Upper/Lower, aa, lambda_1)
        coef_beta_est[j,k] <- coef_beta_kj
      }
      for(s in 1:q){
        lower_1 <- t(Z[,s]) %*% W_k %*% Z[,s]
        lower_2 <- (t(H_q)%*%H_q)[ks(q), ks(q)]
        # lower_2 <- 0
        Lower <- lower_1/n + tau*lower_2
        upper_1 <- t(rho_est[k]*y - Z[,-s]%*%(coef_alpha_est[-s,k]) - X%*%coef_beta_est[,k]) %*% W_k %*% Z[,s]
        upper_2 <- (t(H_q)%*%(diff_w_est-dual_zeta_est))[ks(q)] - 
          as.vector(coef_alpha_est)[-ks(q)]%*%((t(H_q)%*%H_q)[ks(q),-ks(q)])
        # upper_2 <- 0
        Upper <- upper_1/n + tau*upper_2
        # coef_alpha_ks <- Upper/Lower
        coef_alpha_ks <- mcp_solution(Lower, Upper/Lower, aa, lambda_1)
        coef_alpha_est[s,k] <- coef_alpha_ks
      }
    }
    coef_full_est <- rbind(coef_beta_est, coef_alpha_est)
    
    
    # rho 的更新
    # rho_list[[iter]] <- rho_list[[iter-1]]
    # coef_full_list[[iter]] <- coef_full_est
    # coef_full_ori_list[[iter]] <- coef_full_est/extend_x_to_row(rho_est,p+q)
    rho_list[[iter]] <- rho_est
    coef_full_list[[iter]] <- coef_full_est

    coef_full_ori_list[[iter]] <- coef_full_est/extend_x_to_row(rho_est,p+q)
    
    diff_v_est <- matrix(H_p %*% as.vector(coef_beta_est), ncol = 1)
    diff_w_est <- matrix(H_q %*% as.vector(coef_alpha_est), ncol = 1)
    diff_v_last <- diff_v_list[[iter-1]]
    diff_w_last <- diff_w_list[[iter-1]]
    
    # v,w 的更新 2 means prime
    diff_v2_est <- H_p%*%as.vector(coef_beta_est) + dual_xi_est
    diff_w2_est <- H_q%*%as.vector(coef_alpha_est) + dual_zeta_est
    
    for(k1 in 1:(K_up-1)){
      for(k2 in (k1+1):K_up){
        ikk <- which(apply(comb_pair, 2, function(x){all(x==c(k1,k2))}))
        
        v_kk <- diff_v_est[((ikk-1)*p+1):(ikk*p)]
        w_kk <- diff_w_est[((ikk-1)*q+1):(ikk*q)]
        u_kk <- c(v_kk, w_kk)
        
        v_kk_last <- diff_v_last[((ikk-1)*p+1):(ikk*p)]
        w_kk_last <- diff_w_last[((ikk-1)*q+1):(ikk*q)]
        u_kk_last <- c(v_kk_last, w_kk_last)
        
        v_kk.F <- norm(matrix(v_kk_last), type = "2")
        w_kk.F <- norm(matrix(w_kk_last), type = "2")
        u_kk.F <- sqrt(v_kk.F**2 + w_kk.F**2)
        
        v_kk2 <- diff_v2_est[((ikk-1)*p+1):(ikk*p)]
        w_kk2 <- diff_w2_est[((ikk-1)*q+1):(ikk*q)]
        
        v_kk2.F <- norm(matrix(v_kk2), type = "2")
        w_kk2.F <- norm(matrix(w_kk2), type = "2")
        u_kk2.F <- sqrt(v_kk2.F**2 + w_kk2.F**2)
        
        cu <- positive_part(1-lambda_2/tau/u_kk2.F)/(1-1/aa/tau)
        cv <- positive_part(1-lambda_3/tau/v_kk2.F)/(1-1/aa/tau)
        if(is.na(cu)){ print(paste(iter, "cu na of k1 k2", k1, k2)); cu <- 0 }
        if(is.na(cv)){ print(paste(iter, "cv na of k1 k2", k1, k2)); cv <- 0 }
        # 判断压缩类别
        case_kk <- case_when(
          (u_kk2.F > aa*lambda_2) & (v_kk2.F > aa*lambda_3) ~ 1,
          (u_kk2.F <= aa*lambda_2) & (cu*v_kk2.F > aa*lambda_3) ~ 2,
          (w_kk2.F**2+(cv*v_kk2.F)**2 > (aa*lambda_2)**2) &
            (v_kk2.F <= aa*lambda_3) ~ 3,
          TRUE ~ 4
        )
        # tracker
        case <- c(case, case_kk)
        cu_tracker <- c(cu_tracker, cu)
        cv_tracker <- c(cv_tracker, cv)
      
        
        if(case_kk == 1){
          diff_v_est[((ikk-1)*p+1):(ikk*p)] <- v_kk2
          diff_w_est[((ikk-1)*q+1):(ikk*q)] <- w_kk2
        }else if(case_kk == 2){
          diff_v_est[((ikk-1)*p+1):(ikk*p)] <- cu*v_kk2
          diff_w_est[((ikk-1)*q+1):(ikk*q)] <- cu*w_kk2
        }else if(case_kk == 3){
          diff_v_est[((ikk-1)*p+1):(ikk*p)] <- cv*v_kk2
          diff_w_est[((ikk-1)*q+1):(ikk*q)] <- w_kk2
        }else{
          diff_v_est[((ikk-1)*p+1):(ikk*p)] <- v_kk2/
            (1+mcp_d(u_kk.F,aa,lambda_2,TRUE)/tau/(u_kk.F+eps) +
               mcp_d(v_kk.F,aa,lambda_3,TRUE)/tau/(v_kk.F+eps))
          diff_w_est[((ikk-1)*q+1):(ikk*q)] <- w_kk2/
            (1+mcp_d(u_kk.F,aa,lambda_2,TRUE)/tau/(u_kk.F+eps))
          # 用 0 初始化时以下不合适，会总被压缩为 0
          # print(paste("case 4: ", sum(diff_v_est < 1e-9), sum(diff_w_est < 1e-9)))
          # diff_v_est[diff_v_est < 1e-9] <- 0
          # diff_w_est[diff_w_est < 1e-9] <- 0
        }
      }
    }
    diff_v_list[[iter]] <- diff_v_est
    diff_w_list[[iter]] <- diff_w_est
    
    # xi,zeta 的更新
    dual_xi_list[[iter]] <- dual_xi_list[[iter-1]] +
      (matrix(H_p%*%as.vector(coef_beta_est),ncol=1) - diff_v_list[[iter]])
    dual_zeta_list[[iter]] <- dual_zeta_list[[iter-1]] +
      (matrix(H_q%*%as.vector(coef_alpha_est),ncol=1) - diff_w_list[[iter]])
    
    # diff_v_list[[iter]] <- matrix(H_p%*%as.vector(coef_beta_est),ncol=1) + dual_xi_est/tau
    # diff_w_list[[iter]] <- matrix(H_q%*%as.vector(coef_alpha_est),ncol=1) + dual_zeta_est/tau
    # dual_xi_list[[iter]] <- dual_xi_list[[iter-1]] +
    #   tau*(matrix(H_p%*%as.vector(coef_beta_est),ncol=1) - diff_v_list[[iter]])
    # dual_zeta_list[[iter]] <- dual_zeta_list[[iter-1]] +
    #   tau*(matrix(H_q%*%as.vector(coef_alpha_est),ncol=1) - diff_w_list[[iter]])
    
    pi_est <- apply(q_c_matrix, 2, sum)/n
    pi_list[[iter]] <- pi_est
    
    # update q_c matrix
    # 用没有重参数化的形式
    q_c_matrix_back <- q_c_matrix
    q_c_matrix <- dnorm(y,
                        mean = data%*%coef_full_ori_list[[iter]],
                        sd = 1/extend_x_to_row(rho_est, n)) *
      extend_x_to_row(pi_list[[iter-1]], n)
    nan_id <- which(apply(q_c_matrix, 1, sum) == 0)
    if(length(nan_id == 0) > 0){ print(paste(iter, nan_id)) }
    q_c_matrix[nan_id,] <- q_c_matrix_back[nan_id,]
    q_c_matrix <- q_c_matrix / apply(q_c_matrix, 1, sum)
    
    
    # coef_diff_std <- (sum(coef_full_list[[iter]]-coef_full_list[[iter-1]])**2)/
    #    (sum(coef_full_list[[iter]])**2 + 1e-4)
    # if(coef_diff_std < 1e-6){
    #   return()
    # }
  }
  cdist <- coef_dist(coef_full_ori_list[[iter]], coef$coef_full)
  # print(coef_full_list[[iter]])
  print(coef_full_ori_list[[iter]])
  print(table(case))
  plot(unlist(lapply(coef_full_ori_list, coef_dist, coef$coef_full)), ylab = "",
       main = paste(q_c_seed, aa, lambda_1, lambda_2, lambda_3, tau))
  cat("q_c_seed", q_c_seed, "a", aa, 
      "l1", lambda_1, "l2", lambda_2, "l3", lambda_3, "tau", tau)
  ci_est <- apply(q_c_matrix, 1, which.max)
  ci_prob <- apply(q_c_matrix, 1, max)
  ci_matrix <- t(apply(q_c_matrix, 1, function(x){as.numeric(x == max(x))}))
  y_hat <- rowSums(ci_matrix * data%*%coef_full_list[[iter]])
  mse <- sum((y-y_hat)^2/n)
  sc_score <- sc(ci_est, ci_sim)
  
  print(cdist)
  print(mse)
  return(list(cdist = cdist, 
              ci_prob_mean = mean(ci_prob),
              # q_c_matrix = q_c_matrix,
              coef_full_ori = coef_full_ori_list[[iter]],
              # coef_full = coef_full_list[[iter]],
              mse = mse,
              sc_score = sc_score))
}