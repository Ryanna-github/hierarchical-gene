
flexmix_init <- function(q_c_seed, minprior_value = 0, tag = "flexmix"){
  set.seed(q_c_seed)
  pi_init <- rep(1/K_up, K_up)
  rho_init <- rep(1, K_up)/sigma_est
  q_c_matrix <- abs(t(kronecker(pi_init, matrix(1, ncol = n))) + rnorm(n*K_up, mean = 0, sd = .1))
  q_c_matrix <- q_c_matrix / apply(q_c_matrix, 1, sum)
  
  m_glm <- tryCatch({
    flexmix(y~cbind(X, Z)-1, cluster = q_c_matrix,
            model = FLXMRglm(),
            control = list(minprior = minprior_value))
  }, error = function(err) {
    cat("Error occurred:", conditionMessage(err), "\n")
    flexmix(y~cbind(X, Z)-1, k = K_up,
            model = FLXMRglm(),
            control = list(minprior = minprior_value))
  })
  K_est <- m_glm@k
  coef_est <- parameters(m_glm)[1:(p+q),]
  row.names(coef_est) <- NULL
  # cdist <- ifelse(ncol(coef$coef_full) == ncol(coef_est), 
  #                 coef_dist(coef_est, coef$coef_full),
  #                 NaN)
  cdist <- tryCatch({
    coef_dist(coef_est, coef$coef_full)
  }, error = function(err) {NaN})
  sc_score <- sc(m_glm@cluster, ci_sim)
  ari_score <- ari(m_glm@cluster, ci_sim)
  print(coef_est)
  return(list(cdist = cdist, 
              coef_full_ori = coef_est,
              est_sub_grn = K_est,
              sc_score = sc_score,
              ari_score = ari_score,
              tag = tag))
}

random_init <- function(q_c_seed, tag = "random"){
  set.seed(q_c_seed)
  pi_init <- rep(1/K_up, K_up)
  rho_init <- rep(1, K_up)/sigma_est
  q_c_matrix <- abs(t(kronecker(pi_init, matrix(1, ncol = n))) + rnorm(n*K_up, mean = 0, sd = .1))
  q_c_matrix <- q_c_matrix / apply(q_c_matrix, 1, sum)
  # coef_est <- coef$coef_full - coef$coef_full + rnorm(prod(dim(coef$coef_full)), 0, 1)
  # 随机初始化,事先不知道 coef$coef_full,应该用 K_up 初始化
  coef_est <- matrix(rnorm(K_up*(p+q), 0, 1), ncol = K_up)
  
  # cdist <- ifelse(ncol(coef$coef_full) == ncol(coef_est), 
  #                 coef_dist(coef_est, coef$coef_full),
  #                 NaN)
  cdist <- tryCatch({
    coef_dist(coef_est, coef$coef_full)
  }, error = function(err) {NaN})
  ci_est <- apply(q_c_matrix, 1, which.max)
  ci_prob_mean <- mean(apply(q_c_matrix, 1, max))
  # sc_score <- sc(ci_est, ci_sim)
  sc_score <- tryCatch({
    sc(ci_est, ci_sim)
  }, error = function(err) {NaN})
  ari_score <- tryCatch({
    ari(ci_est, ci_sim)
  }, error = function(err) {NaN})
  print(coef_est)
  return(list(cdist = cdist,
              ci_prob_mean = ci_prob_mean,
              coef_full_ori = coef_est,
              est_sub_grn = K_up,
              sc_score = sc_score,
              ari_score = ari_score,
              tag = tag))
}

bic_score <- function(q_c_matrix, coef_est, est_main_grn, est_sub_grn, rho_est){
  Cn <- log(n*(p+q))
  pi_est <- apply(q_c_matrix, 2, sum)/n
  fit_matrix <- dnorm(y, mean = data%*%coef_est, sd = 1/extend_x_to_row(rho_est, n)) *
    extend_x_to_row(pi_est, n)
  fit_prob <- apply(fit_matrix, 1, sum)
  fit_sum <- -2 * sum(log(fit_prob + min(fit_prob[fit_prob>0])))
  fit_mean <- fit_sum/n
  penal <- Cn*log(n)/n*(est_main_grn*p+est_sub_grn*q)
  bic_sum <- fit_sum + penal
  bic_mean <- fit_mean + penal
  return(list(fit_sum = fit_sum,
              fit_mean = fit_mean,
              penal = penal,
              bic_sum = bic_sum,
              bic_mean = bic_mean))
}

ADMM_trail <- function(aa, tau, lambda_1, lambda_2, lambda_3, q_c_seed, 
                       coef_full_init, iter_type, iter_max, rho_clip = 2,
                       plot_performance = FALSE,
                       eps = 1e-7, eps_abs = 1e-2, eps_rel = 1e-3){
  # iter_type = "fix" 固定迭代轮次, iter_type = "stop" 使用迭代停止准则判断
  tau <- ifelse(tau == 0, 1e-4, tau) # 防止 tau == 0 导致分母为0情况
  sigma_est = ifelse(sigma_est == 0, 0.01, sigma_est)
  rho_init <- rep(1, K_up)/sigma_est
  
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
  
  resi_dual_list <- list(NULL)
  resi_prim_list <- list(NULL)
  eps_dual_list <- list(NULL)
  eps_prim_list <- list(NULL)
  
  case <- NULL
  cu_tracker <- NULL
  cv_tracker <- NULL
  
  for(iter in 2:iter_max){
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
    
    # rho_k 的更新
    # update_step = ifelse(iter > 10, 0.4, 0)
    update_step = 0
    for(k in 1:K_up){
      W_k <- diag(q_c_matrix[,k])
      onen <- matrix(1, nrow = n, ncol = 1)
      AA <- as.numeric(t(y) %*% W_k %*% y)**2
      BB <- as.numeric(-t(onen) %*% W_k %*% (X%*%coef_beta_est[,k] + Z%*%coef_alpha_est[,k]))
      CC <- as.numeric(-n*t(onen) %*% W_k %*% onen)
      rho_k_est <- (-BB+sqrt(BB**2-4*AA*CC))/(2*AA)
      rho_k_est <- ifelse(rho_k_est <= 1, 1, min(rho_clip, rho_k_est))
      rho_est[k] <- update_step*rho_k_est + (1-update_step)*rho_est[k]
    }
    
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
    
    # stopping criteria
    if(iter_type == "stop"){
      resi_prim <- rbind(matrix(H_p%*%as.vector(coef_beta_est),ncol=1) - diff_v_list[[iter]],
                         matrix(H_q%*%as.vector(coef_alpha_est),ncol=1) - diff_w_list[[iter]]) %>%
        norm(type = '2')
      resi_dual <- rbind(-tau * t(H_p) %*% (diff_v_list[[iter]] - diff_v_list[[iter-1]]),
                         -tau * t(H_q) %*% (diff_w_list[[iter]] - diff_w_list[[iter-1]])) %>%
        norm(type = '2')
      resi_prim_list[[iter]] <- resi_prim
      resi_dual_list[[iter]] <- resi_dual
      
      eps_prim <- sqrt(K_up*(K_up-1)/2*(p+q))*eps_abs + 
        eps_rel*max(norm(rbind(H_p%*%as.numeric(coef_beta_est), H_q%*%as.numeric(coef_alpha_est)), type = '2'), 
                    norm(rbind(diff_v_est, diff_w_est), type = '2'))
      eps_dual <- sqrt(K_up*(p+q))*eps_abs +
        eps_rel*norm(rbind(t(H_p)%*%dual_xi_list[[iter]], t(H_q)%*%dual_zeta_list[[iter]]), type = '2')
      eps_prim_list[[iter]] <- eps_prim
      eps_dual_list[[iter]] <- eps_dual
      
      if((resi_prim < eps_prim) & (resi_dual < eps_dual)){
        break
      }
    }

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
    if(length(nan_id == 0) > 0){ print(paste("[nan]", iter, nan_id)) }
    q_c_matrix[nan_id,] <- q_c_matrix_back[nan_id,]
    q_c_matrix <- q_c_matrix / apply(q_c_matrix, 1, sum)
    
    
    # coef_diff_std <- (sum(coef_full_list[[iter]]-coef_full_list[[iter-1]])**2)/
    #    (sum(coef_full_list[[iter]])**2 + 1e-4)
    # if(coef_diff_std < 1e-6){
    #   return()
    # }
  }
  # 距离真实参数距离，只有类别数相同才能计算该指标
  # cdist <- ifelse(ncol(coef$coef_full) == ncol(coef_full_ori_list[[iter]]), 
  #                 coef_dist(coef_full_ori_list[[iter]], coef$coef_full),
  #                 NaN)
  
  cdist <- tryCatch({
    coef_dist(coef_full_ori_list[[iter]], coef$coef_full)
  }, error = function(err) {NaN})
  # 参数结果（未压缩）
  cat(paste0("****[unsqueezed coef]****:\n"))
  print(coef_full_ori_list[[iter]])
  # case 情况（检查是否落入组别压缩的情况
  case_table_full <- rep(0, 4) # v,w 更新四种情况落入次数记录
  case_table <- table(case)
  for(t in 1:dim(case_table)){
    case_table_full[as.integer(names(case_table)[t])] <- case_table[names(case_table)[t]]
  }
  
  if(plot_performance){
    plot(1:iter, c(resi_prim_list[[2]], unlist(resi_prim_list)), col = "Black", main = "primary")
    points(c(eps_prim_list[[2]], unlist(eps_prim_list)), col = "DeepPink")
    lines(c(resi_prim_list[[2]], unlist(resi_prim_list)), col = "Black")
    lines(c(eps_prim_list[[2]], unlist(eps_prim_list)), col = "DeepPink")
    
    plot(1:iter, c(resi_dual_list[[2]], unlist(resi_dual_list)), col = "Black", main = "dual")
    points(c(eps_dual_list[[2]], unlist(eps_dual_list)), col = "DeepPink")
    lines(c(resi_dual_list[[2]], unlist(resi_dual_list)), col = "Black")
    lines(c(eps_dual_list[[2]], unlist(eps_dual_list)), col = "DeepPink")
    # 估计参数距离真实参数距离的可视化展示
    plot(unlist(lapply(coef_full_ori_list, coef_dist, coef$coef_full)), ylab = "",
         main = paste(q_c_seed, aa, lambda_1, lambda_2, lambda_3, tau))
  }
  
  ci_est <- apply(q_c_matrix, 1, which.max)
  sc_score <- tryCatch({
    sc(ci_est, ci_sim)
  }, error = function(err) {NaN})
  ari_score <- tryCatch({
    ari(ci_est, ci_sim)
  }, error = function(err) {NaN})
  
  ci_prob <- apply(q_c_matrix, 1, max)
  ci_matrix <- t(apply(q_c_matrix, 1, function(x){as.numeric(x == max(x))}))
  y_hat <- rowSums(ci_matrix * data%*%coef_full_list[[iter]])
  mse <- sum((y-y_hat)^2/n)
  
  # group result
  main_group_info <- get_group_num(K_up, 
                                  coef_full_ori_list[[iter]][(1:p),], 
                                  diff_v_list[[iter]], p)
  sub_group_info <- get_group_num(K_up, 
                                  coef_full_ori_list[[iter]][(1+p):(p+q),], 
                                  diff_w_list[[iter]], q)
  # print(main_group_info$capgfl.matrix2)
  # print(sub_group_info$capgfl.matrix2)
  est_main_grn <- main_group_info$gr.num
  est_sub_grn <- sub_group_info$gr.num
  
  cappfl.diff <- main_group_info$capgfl.matrix2 - sub_group_info$capgfl.matrix2
  valid_hier <- ifelse(min(cappfl.diff) >= 0, TRUE, FALSE)
  diag(cappfl.diff) <- 1
  group_detail <- apply(which((cappfl.diff) > 0, arr.ind = TRUE), 1, function(x){paste0("(",x[1],",",x[2],")")})
  group_detail <- paste0(group_detail, collapse = "")
  # BIC.var <- log(mse) + log(n*(p+q))*log(n)*(est_main_grn*(p)+est_sub_grn*(q))/n
  # BIC.o <- log(mse) + log(n)*(est_main_grn*p+est_sub_grn*q)/n
  
  bic_info <- bic_score(q_c_matrix, coef_full_ori_list[[iter]], 
                        est_main_grn, est_sub_grn, rho_est)
  
  cat(paste("**** hypter parameter:", "dt_seed", dt_seed, "q_c_seed", q_c_seed, 
            "l1", lambda_1, "l2", lambda_2, "l3", lambda_3, "tau", tau, "a", aa, "\n"))
  cat(paste("**** result:", "est_main_grn", est_main_grn, "est_sub_grn", est_sub_grn,
            "mse", mse, "sc", sc_score, "bic_sum", bic_info$bic_sum, "bic_mean", bic_info$bic_mean, 
            "valid_hier", valid_hier, "group_detail", group_detail, "\n"))
  print(case_table_full)
  print(coef_full_ori_list[[iter]])
  
  return(list(cdist = cdist, 
              ci_prob_mean = mean(ci_prob),
              q_c_matrix = q_c_matrix,
              coef_full_ori = coef_full_ori_list[[iter]],
              # coef_full = coef_full_list[[iter]],
              mse = mse,
              # BIC.var = BIC.var,
              # BIC.o = BIC.o,
              fit_sum = bic_info$fit_sum,
              fit_mean = bic_info$fit_mean,
              penal = bic_info$penal,
              bic_sum = bic_info$bic_sum,
              bic_mean = bic_info$bic_mean,
              sc_score = sc_score,
              ari_score = ari_score,
              case_table_full = case_table_full,
              est_main_grn = main_group_info$gr.num,
              est_sub_grn = sub_group_info$gr.num,
              rho_init = rho_init[1],
              rho_est = paste0("(", paste(as.character(round(rho_list[[iter]],3)), collapse = ","), ")"),
              case_table_full = case_table_full,
              valid_hier = valid_hier,
              group_detail = group_detail,
              iter_total = iter,
              iter_type = iter_type))
}

# 根据 beta（alpha）v（w）返回组别个数
get_group_num <- function(K_up, coef, diff_v, len, merge.all = F, threshold = 1e-2){
  # diff_v <- diff_w_list[[iter]]
  # coef <- coef_full_ori_list[[iter]][(1+p):(p+q),]
  # len <- q
  # diff_v <- diff_v_list[[iter]]
  # coef <- coef_full_ori_list[[iter]][(1:p),]
  # len <- p
  diff.gfl <- apply(matrix(diff_v, nrow = len, ncol = K_up*(K_up-1)/2), 2, sum)
  squeeze_index <- which(abs(diff.gfl) < threshold)
  
  capgfl.matrix <- matrix(0, nrow = K_up, ncol = K_up)
  if(length(squeeze_index) == 0){
    capgfl.matrix2 <- capgfl.matrix + t(capgfl.matrix)
    diag(capgfl.matrix2) <- 1
    return(list(gr.num = K_up, capgfl.matrix2 = capgfl.matrix2))
  }
  if(length(squeeze_index) == K_up*(K_up-1)/2){
    capgfl.matrix[upper.tri(capgfl.matrix)] <- 1
    capgfl.matrix2 <- capgfl.matrix + t(capgfl.matrix)
    diag(capgfl.matrix2) <- 1
    return(list(gr.num = 1, capgfl.matrix2 = capgfl.matrix2))
  }else{
    sample.index.gfl <- comb_pair[,squeeze_index]
    if(length(which(abs(diff.gfl) < threshold)) == 1){
      capgfl.matrix[sample.index.gfl[1], sample.index.gfl[2]]<-1
    }else{
      for(i in 1:length(squeeze_index)){
        capgfl.matrix[sample.index.gfl[1,i],sample.index.gfl[2,i]]<-1
      }
    }
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = nrow(unique(capgfl.matrix2))
    # 解决如 12 同类 13 同类但 23 不同类的问题
    if(merge.all){
      cap <- capgfl.matrix2
      num_subgroup <- unique(apply(cap, 1, function(a){which(a == 1)}))
      non_inter_list <- list()
      vv <- 1
      non_inter <- c(1:length(num_subgroup))
      repeat{
        a <- num_subgroup[[non_inter[1]]]
        KK_k <- setdiff(non_inter,non_inter[1])
        non_inter <- c()
        i=1
        for (k2 in KK_k) {
          if(length(intersect(a,num_subgroup[[k2]])) > 0){
            a <- union(a,num_subgroup[[k2]])
          } else {
            non_inter[i] <- k2
            i=i+1
          }
        }
        non_inter_list[[vv]] <- a
        vv <- vv+1
        if(length(non_inter) == 0){break}
      }
      
      for (i in 1:dim(cap)[1]) {
        for (k in 1:length(non_inter_list)) {
          if(length(match(cap[i,],non_inter_list[[k]])) > 0){
            cap[i,non_inter_list[[k]]] <- 1
          }
        }
      }
      capgfl.matrix2 <- cap
      group.num.gf = nrow(unique(capgfl.matrix2))
    }
    return(list(gr.num=group.num.gf,capgfl.matrix2=capgfl.matrix2))
  }
}

# tune l3 first then l2
tuning_hyper <- function(l2_seq, l3_seq, fix_para, coef_full_init, grid = TRUE, save_all = FALSE){
  if(save_all){
    result <- NULL
    trail_set <- expand.grid(list(l3 = l3_seq, l2 = l2_seq))
    trail_num <- nrow(trail_set)
    # bic_record <- rep(-Inf, trail_num)
    trail_record <- vector(mode = "list",length = trail_num)
    # 分组实验
    for(trail_idx in 1:trail_num){
      lambda_2 <- trail_set$l2[trail_idx]
      lambda_3 <- trail_set$l3[trail_idx]
      
      trail <- ADMM_trail(aa = fix_para$aa,
                          tau = fix_para$tau,
                          lambda_1 = fix_para$lambda_1,
                          lambda_2 = lambda_2,
                          lambda_3 = lambda_3,
                          q_c_seed = fix_para$q_c_seed,
                          coef_full_init = coef_full_init,
                          iter_type = 'stop',
                          iter_max = 200)
      trail_record[[trail_idx]] <- trail
      # bic_record[trail_idx] <- trail$BIC.var
      
      result <- rbind(result, c(fix_para$dt_seed,
                  fix_para$q_c_seed,
                  fix_para$aa,
                  fix_para$tau,
                  fix_para$lambda_1,
                  lambda_2,
                  lambda_3,
                  trail$cdist,
                  trail$ci_prob_mean,
                  trail$mse,
                  trail$sc_score,
                  trail$ari_score,
                  # trail$BIC.var,
                  trail$fit_sum,
                  trail$fit_mean,
                  trail$penal,
                  trail$bic_sum,
                  trail$bic_mean,
                  trail$est_main_grn,
                  trail$est_sub_grn,
                  trail$rho_init,
                  trail$rho_est,
                  trail$valid_hier,
                  trail$group_detail,
                  trail$case_table_full,
                  trail$iter_total,
                  trail$iter_type,
                  "hier"))
    }
    return(result)
  }
  else if(grid){
    trail_set <- expand.grid(list(l3 = l3_seq, l2 = l2_seq))
    trail_num <- nrow(trail_set)
    bic_record <- rep(-Inf, trail_num)
    trail_record <- vector(mode = "list",length = trail_num)
    for(trail_idx in 1:trail_num){
      lambda_2 <- trail_set$l2[trail_idx]
      lambda_3 <- trail_set$l3[trail_idx]
      
      trail <- ADMM_trail(aa = fix_para$aa,
                          tau = fix_para$tau,
                          lambda_1 = fix_para$lambda_1,
                          lambda_2 = lambda_2,
                          lambda_3 = lambda_3,
                          q_c_seed = fix_para$q_c_seed,
                          coef_full_init = coef_full_init,
                          sigma_est = )
      trail_record[[trail_idx]] <- trail
      bic_record[trail_idx] <- trail$BIC.var
    }
    best_idx <- which(bic_record == min(bic_record))
    trail <- trail_record[[best_idx]]
    lambda_2 <- trail_set$l2[best_idx]
    lambda_3 <- trail_set$l3[best_idx]
  }else{
    # ============================== tune l3 ======================
    trail_set_step1 <- expand.grid(list(l3 = l3_seq, l2 = 0))
    trail_num <- nrow(trail_set_step1)
    bic_record <- rep(-Inf, trail_num)
    trail_record <- vector(mode = "list",length = trail_num)
    for(trail_idx in 1:trail_num){
      lambda_2 <- trail_set_step1$l2[trail_idx]
      lambda_3 <- trail_set_step1$l3[trail_idx]
      
      trail <- ADMM_trail(aa = fix_para$aa,
                          tau = fix_para$tau,
                          lambda_1 = fix_para$lambda_1,
                          lambda_2 = lambda_2,
                          lambda_3 = lambda_3,
                          q_c_seed = fix_para$q_c_seed,
                          coef_full_init = coef_full_init)
      trail_record[[trail_idx]] <- trail
      bic_record[trail_idx] <- trail$BIC.var
    }
    best_idx <- which(bic_record == min(bic_record))
    trail <- trail_record[[best_idx]]
    lambda_3 <- trail_set_step1$l3[best_idx]
    # ============================== tune l2 ======================
    trail_set_step2 <- expand.grid(list(l3 = lambda_3, l2 = l2_seq))
    trail_num <- nrow(trail_set_step2)
    bic_record <- rep(-Inf, trail_num)
    trail_record <- vector(mode = "list",length = trail_num)
    for(trail_idx in 1:trail_num){
      lambda_2 <- trail_set_step2$l2[trail_idx]
      lambda_3 <- trail_set_step2$l3[trail_idx]
      
      trail <- ADMM_trail(aa = fix_para$aa,
                          tau = fix_para$tau,
                          lambda_1 = fix_para$lambda_1,
                          lambda_2 = lambda_2,
                          lambda_3 = lambda_3,
                          q_c_seed = fix_para$q_c_seed,
                          coef_full_init = coef_full_init)
      trail_record[[trail_idx]] <- trail
      bic_record[trail_idx] <- trail$BIC.var
    }
    best_idx <- which(bic_record == min(bic_record))
    trail <- trail_record[[best_idx]]
    lambda_2 <- trail_set_step2$l2[best_idx]
    lambda_3 <- trail_set_step2$l3[best_idx]
  }
  result <- c(fix_para$dt_seed,
              fix_para$q_c_seed,
              fix_para$aa,
              fix_para$tau,
              fix_para$lambda_1,
              lambda_2,
              lambda_3,
              trail$cdist,
              trail$ci_prob_mean,
              trail$mse,
              trail$sc_score,
              # trail$BIC.var,
              trail$fit_sum,
              trail$fit_mean,
              trail$penal,
              trail$bic_sum,
              trail$bic_mean,
              trail$est_main_grn,
              trail$est_sub_grn,
              trail$valid_hier,
              trail$group_detail,
              trail$case_table_full,
              "hier")
  return(result)
}








