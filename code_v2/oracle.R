library(dplyr)
library(Matrix)
library(flexmix)
# source("sim.R")
# source("tools.R")
# source("func.R")
source("hierarchical-gene/code_v2/sim.R")
source("hierarchical-gene/code_v2/tools.R")
source("hierarchical-gene/code_v2/func.R")
library(argparse)

parser <- ArgumentParser()

parser$add_argument("-n", "--num", default=200,  help="sample size")
# parser$add_argument("--dt_seed", default=9,  help="epsilon seed while generating data")
parser$add_argument("-p", default = 8, help = "dim of X")
parser$add_argument("-q", default = 4, help = "dim of Z")
parser$add_argument("-b", default = 1, help = "if data balanced (1: balanced, 2: unbalanced groups with balanced subgroups, 3: balanced groups with unbalanced sbgroups)")
parser$add_argument("-e", "--epsilon_sd", default = 0.5, help = "error")
parser$add_argument("--epsilon_sd_init", default = 0.5, help = "epsilon_sd initiation for estimating rho")
parser$add_argument("--rho_ratio", default = 0.5, help = "update ratio of rho, rho(t+1)=ratio*rho+(1-ratio)rho(t)")
parser$add_argument("-ss", "--signal_size", default = 1, help = "signal size")
parser$add_argument("-bl", "--beta_vlen", default = 3, help = "nonzero length of beta")
parser$add_argument("-al", "--alpha_vlen", default = 2, help = "nonzero length of alpha")
parser$add_argument("--path", default="temp.csv",  help="csv result save path")
parser$add_argument("--K_up", default=4,  help="Upper class number")
parser$add_argument("--cotype", default="En",  help="Covariate matrix of X and Z")



args <- parser$parse_args()
print(str(args))

n <- as.numeric(args$n)
p <- as.numeric(args$p)
q <- as.numeric(args$q)
balance <- as.numeric(args$b)
epsilon_sd <- as.numeric(args$epsilon_sd)
sigma_est <- as.numeric(args$epsilon_sd_init)
rho_ratio <- as.numeric(args$rho_ratio)
signal_size <- as.numeric(args$signal_size)
beta_vlen <- as.numeric(args$beta_vlen)
alpha_vlen <- as.numeric(args$alpha_vlen)
save_path <- args$path
# dt_seed <- as.numeric(args$dt_seed)
K_up <- as.numeric(args$K_up)
cotype <- args$cotype

if(0){
  n <- 500
  p <- 80
  q <- 40
  balance <- 1
  epsilon_sd <- 0.5
  sigma_est <- 0.5
  rho_ratio <- 0
  signal_size <- 1
  beta_vlen <- 3
  alpha_vlen <- 2
  save_path <- "oracle.csv"
  # dt_seed <- as.numeric(args$dt_seed)
  K_up <- 4
  cotype <- "En"
}

oracle_comp <- function(q_c_seed, iter_type = 'stop', iter_max = 200, 
                        rho_ratio = 0, sigma_est = 0.5, K_up_true = 4, 
                        rho_clip = 4, plot_performance = FALSE,
                        eps = 1e-7, eps_abs = 1e-2, eps_rel = 1e-3, tag = "oracle"){
  # iter_type = 'stop'
  # iter_max = 200
  # rho_ratio = 0
  # sigma_est = 0.5
  # K_up_true = 4
  # rho_clip = 4
  # plot_performance = FALSE
  # eps = 1e-7
  # eps_abs = 1e-2
  # eps_rel = 1e-3
  # tag = "oracle"
  tau <- ifelse(tau == 0, 1e-4, tau) # 防止 tau == 0 导致分母为0情况
  sigma_est = ifelse(sigma_est == 0, 0.01, sigma_est)
  rho_init <- rep(1, K_up)/sigma_est
  
  kj <- function(dim) {return((k-1)*dim+j)}
  ks <- function(dim) {return((k-1)*dim+s)}
  rho_list <- list(rho_init)
  pi_init <- rep(1/K_up, K_up)
  pi_list <- list(pi_init)
  set.seed(q_c_seed)
  q_c_matrix <- get_q_c_matrix(n, K_up, ci_sim)
  q_c_matrix <- q_c_matrix / apply(q_c_matrix, 1, sum)
  
  coef_full_init <- matrix(0, nrow = p+q, ncol = K_up)
  coef_full_list <- list(coef_full_init*t(kronecker(rho_init, matrix(1, ncol = p+q))))
  coef_full_ori_list <- list(coef_full_init)
  
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
    
    # rho_k 的更新
    # update_step = ifelse(iter > 10, 0.4, 0)
    rho_ratio = 0
    for(k in 1:K_up){
      W_k <- diag(q_c_matrix[,k])
      onen <- matrix(1, nrow = n, ncol = 1)
      AA <- as.numeric(t(y) %*% W_k %*% y)
      BB <- as.numeric(-t(y) %*% W_k %*% (X%*%coef_beta_est[,k] + Z%*%coef_alpha_est[,k]))
      CC <- as.numeric(-t(onen) %*% W_k %*% onen)
      rho_k_est <- (-BB+sqrt(BB**2-4*AA*CC))/(2*AA)
      rho_k_est <- ifelse(rho_k_est <= 0.5, 1, min(rho_clip, rho_k_est))
      rho_est[k] <- rho_ratio*rho_k_est + (1-rho_ratio)*rho_est[k]
    }
    
    rho_list[[iter]] <- rho_est
    
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
    
    # update q_c matrix
    pi_est <- apply(q_c_matrix, 2, sum)/n
    pi_list[[iter]] <- pi_est
  }
  # ========================= 计算结束整理结果 ==============================
  cdist <- tryCatch({
    coef_dist(coef_full_ori_list[[iter]], coef$coef_full)
  }, error = function(err) {NaN})
  
  # case 情况（检查是否落入组别压缩的情况
  case_table_full <- rep(0, 4) # v,w 更新四种情况落入次数记录
  case_table <- table(case)
  for(t in 1:dim(case_table)){
    case_table_full[as.integer(names(case_table)[t])] <- case_table[names(case_table)[t]]
  }
  
  # 参数距离可视化 & 停止准则和对应临界值
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
    # 类别个数与真实不相同时，以下可视化会出错
    tryCatch(
      plot(unlist(lapply(coef_full_ori_list, coef_dist, coef$coef_full)), ylab = "",
           main = paste(q_c_seed, aa, lambda_1, lambda_2, lambda_3, tau))
      , error = function(e) e)
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
  
  # est_main_grn <- main_group_info$gr.num
  # est_sub_grn <- sub_group_info$gr.num
  
  cappfl.diff <- main_group_info$capgfl.matrix2 - sub_group_info$capgfl.matrix2
  valid_hier <- ifelse(min(cappfl.diff) >= 0, TRUE, FALSE)
  diag(cappfl.diff) <- 1
  
  # ==========================================================================
  # 压缩到相同类后，再计算 sc，ari d等指标得分
  # (1) 找到哪些类别为同一大组，小组
  main_group_info_compact <- find_connected_nodes(main_group_info$capgfl.matrix2)
  sub_group_info_compact <- find_connected_nodes(sub_group_info$capgfl.matrix2)
  est_main_grn <- length(main_group_info_compact)
  est_sub_grn <- length(sub_group_info_compact)
  
  # (2) 计算平均的系数，可以根据样本数进行加权，但因为数值非常接近，可以直接求平均
  coef_beta_ori_comp <- matrix(NaN, p, est_main_grn)
  coef_alpha_ori_comp <- matrix(NaN, q, est_sub_grn)
  coef_full_ori_comp <- matrix(NaN, p+q, est_sub_grn)
  for(k in 1:est_main_grn){
    # 多于两类才需要压缩，否则取原值
    if(length(main_group_info_compact[[k]]) > 1){
      coef_beta_ori_comp[,k] <- rowMeans(coef_full_ori_list[[iter]][1:p,][,main_group_info_compact[[k]]])
    }else{
      coef_beta_ori_comp[,k] <- coef_full_ori_list[[iter]][1:p,][,main_group_info_compact[[k]]]
    }
  }
  for(k in 1:est_sub_grn){
    # 多于两类才需要压缩，否则取原值
    if(length(sub_group_info_compact[[k]]) > 1){
      coef_alpha_ori_comp[,k] <- rowMeans(coef_full_ori_list[[iter]][(p+1):(p+q),][,sub_group_info_compact[[k]]])
    }else{
      coef_alpha_ori_comp[,k] <- coef_full_ori_list[[iter]][(p+1):(p+q),][,sub_group_info_compact[[k]]]
    }
  }
  for(k in 1:K_up){
    coef_full_ori_comp[,which(sapply(sub_group_info_compact, function(x) k %in% x))][(p+1):(p+q)] <- coef_alpha_ori_comp[,which(sapply(sub_group_info_compact, function(x) k %in% x))]
    coef_full_ori_comp[,which(sapply(sub_group_info_compact, function(x) k %in% x))][1:p] <- coef_beta_ori_comp[,which(sapply(main_group_info_compact, function(x) k %in% x))]
  }
  ci_est_main <- sapply(ci_est, function(k){which(sapply(main_group_info_compact, function(x) k %in% x))})
  ci_est_sub <- sapply(ci_est, function(k){which(sapply(sub_group_info_compact, function(x) k %in% x))})
  
  sc_score_main <- sc(ci_est_main, ci_sim_main)
  sc_score_sub <- sc(ci_est_sub, ci_sim_sub)
  ari_score_main <- ari(ci_est_main, ci_sim_main)
  ari_score_sub <- ari(ci_est_sub, ci_sim_sub)
  # per 无法再此处计算，汇总函数中再计算
  cdist_main <- tryCatch({ coef_dist(coef_full_ori_comp[1:p,], coef$coef_beta) }, error = function(err) {NaN})
  if(is.na(cdist_main)){
    cdist_main <- tryCatch({ coef_dist(coef_beta_ori_comp, coef$coef_beta[, !duplicated(t(coef$coef_beta))]) }, error = function(err) {NaN})
  }
  cdist_sub <- tryCatch({ coef_dist(coef_alpha_ori_comp, coef$coef_alpha) }, error = function(err) {NaN})
  # ==========================================================================
  
  group_detail <- paste0(convert_to_parentheses(main_group_info_compact), ";",
                         convert_to_parentheses(sub_group_info_compact))
  # print(coef_full_ori_list[[iter]])
  bic_info <- bic_score(q_c_matrix, coef_full_ori_list[[iter]], 
                        est_main_grn, est_sub_grn, rho_est)
  return(list(dt_seed = dt_seed,
              q_c_seed = q_c_seed,
              tag = tag,
              penal = bic_info$penal,
              bic_sum = bic_info$bic_sum,
              bic_mean = bic_info$bic_mean,
              cdist_main = cdist_main,
              cdist_sub = cdist_sub,
              sc_score_main = sc_score_main,
              sc_score_sub = sc_score_sub,
              ari_score_main = ari_score_main,
              ari_score_sub = ari_score_sub,
              est_main_grn = est_main_grn,
              est_sub_grn = est_sub_grn,
              pi_est = paste0("(", paste(as.character(round(prop.table(table(ci_est_sub)),4)), collapse = ","), ")"),
              mse = mse,
              iter = iter))
}

q_c_seed <- 1
cotype_x <- cotype
cotype_z <- cotype
beta_nonzero <- c(-2, -2, 2, 2)*signal_size # 长度应和真实 group_num_sub 保持一致
alpha_nonzero <- c(-3, -1, 1, 3)*signal_size

q_c_seed_max <- 10
group_num_main <- 2                    
group_num_sub <- 4    
hier_struc <- list(c(1,2),c(3,4))
# prob_sub <- rep(1/group_num_sub, group_num_sub)  
if(balance == 1){
  prob_sub <- rep(1/group_num_sub, group_num_sub)
}else if(balance == 2){
  prob_sub<-c(1/6,1/6,1/3,1/3)
}else{
  prob_sub<-c(1/6,1/3,1/6,1/3)
}
reverse <- FALSE


for(dt_seed in seq(9,99,10)){
  whole.data <- generate_all_data(dt_seed, n, p, q, prob_sub, hier_struc, 
                                  beta_nonzero, alpha_nonzero, beta_vlen, alpha_vlen, 
                                  cotype_x, cotype_z, epsilon_sd, reverse)
  X <- whole.data$data$X
  Z <- whole.data$data$Z
  data <- whole.data$data$data_full
  coef <- whole.data$coef
  # coefv <- lapply(coef, as.vector) # 按照类别拉长
  ci_sim <- whole.data$ci_sim
  ci_sim_main <- whole.data$ci_sim_main
  ci_sim_sub <- whole.data$ci_sim_sub
  y <- matrix(whole.data$y, ncol=1)
  
  # =============================== prepare =================================
  comb_pair <- combn(K_up, 2)
  H_p <- kronecker(t(apply(comb_pair, 2, get_e_mat, K_up)),
                   diag(p)) %>% Matrix(sparse = TRUE)
  H_q <- kronecker(t(apply(comb_pair, 2, get_e_mat, K_up)),
                   diag(q)) %>% Matrix(sparse = TRUE)
  
  iter_max <- 200
  l2_seq <- c(0.5,1,1.2,1.5)
  l3_seq <- c(0.5,1,1.2,1.5)
  # l2_seq <- c(1)
  # l3_seq <- c(1)
  fix_para <- list(dt_seed = dt_seed, q_c_seed = q_c_seed, lambda_1 = 0.3,
                   aa = 1.2, tau = 1)
  trail_set <- expand.grid(list(l3 = l3_seq, l2 = l2_seq))
  
  oraresult <- NULL
  trail_num <- nrow(trail_set)
  for(trail_idx in 1:trail_num){
    aa <- fix_para$aa
    tau <- fix_para$tau
    lambda_1 <- fix_para$lambda_1
    lambda_2 <- trail_set$l2[trail_idx]
    lambda_3 <- trail_set$l3[trail_idx]
    
    trail <- oracle_comp(q_c_seed)
    oraresult <- rbind(oraresult, c(trail$dt_seed, trail$q_c_seed, lambda_2, lambda_3,
                                    trail$tag, trail$penal, trail$bic_sum, trail$bic_mean,
                                    trail$cdist_main, trail$cdist_sub,
                                    trail$sc_score_main, trail$sc_score_sub,
                                    trail$ari_score_main, trail$ari_score_sub,
                                    trail$est_main_grn, trail$est_sub_grn,
                                    trail$pi_est, trail$mse, trail$q_c_seed))
    print(str(trail))
    
    
  }
  oraresult <- data.frame(oraresult)
  colnames(oraresult) <- c("dt_seed", "q_c_seed", "lambda_2", "lambda_3", "tag", 
                           "penal", "bic_sum", "bic_mean", 
                           "cdist_main", "cdist_sub", "sc_score_main",
                           "sc_score_sub", "ari_score_main", "ari_score_sub",
                           "est_main_grn", "est_sub_grn", "pi_est", "mse", "iter")
  oraresult$bic_mean <- as.numeric(oraresult$bic_mean)
  write_csv(oraresult, file=save_path, col_names=!file.exists(save_path), append=TRUE)
}
