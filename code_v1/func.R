# 工具函数
# *********************************************************************************
# 等价相互计算
coef_eta2gamma <- function(coef, p, q){
  coef_beta <- coef[1:p,]
  coef_eta <- coef[(p+q+1):(p+q+p*q),]
  coef_gamma <- coef_eta / kronecker(matrix(1, nrow=q, ncol=1), coef_beta)
  return(rbind(coef[1:(p+q),], coef_gamma))
}
coef_gamma2eta <- function(coef, p, q){
  coef_beta <- coef[1:p,]
  coef_gamma <- coef[(p+q+1):(p+q+p*q),]
  coef_eta <- kronecker(matrix(1, nrow=q, ncol=1), coef_beta) * coef_gamma
  return(rbind(coef[1:(p+q),], coef_eta))
}

# 结果 t(res) %*% beta_k 计算所有 beta_k 相关项
# 返回结果 n*1 矩阵
get_beta_var_fromlast <- function(wi, p, q, gamma_k_old){
  beta_var <- apply(matrix(wi*gamma_k_old, ncol = q), 1, sum)
  return(beta_var)
}
# gamma_k_old <- coef_gamma[,1]
# # res <- apply(whole.data$W, 1, get_beta_var, p, q, gamma_k_old)
# res <- apply(w, 1, get_beta_var_fromlast, p, q, gamma_k_old)
# t(res) %*% beta_k_old

get_gamma_var_fromlast <- function(wi, p, q, beta_k_old){
  gamma_var <- as.numeric(wi * kronecker(matrix(1,nrow=q,ncol=1), beta_k_old))
  return(gamma_var)
}
# res <- apply(w, 1, get_gamma_var, p, q, beta_k_old)
# t(res) %*% gamma_k_old

coef_bind <- function(coef_beta, coef_alpha, coef_last){
  return(rbind(coef_beta, coef_alpha, coef_last))
}
coef_split <- function(coef, p, q){
  coef_beta <- coef[1:p,]
  coef_alpha <- coef[(p+1):(p+q),]
  coef_last <- coef[(p+q+1):(nrow(coef)),]
  return(list(coef_beta, coef_alpha, coef_last))
}
# ***************************************************************************


# 参数初始化
initialize_coef <- function(whole.data, p, q, K_up = 4, int_rm = TRUE) {
  pq <- p * q
  as <- matrix(1, nrow = p + q + pq + 1, ncol = K_up) # activeset parameter in fmrs package
  if (int_rm) { as[1, ] <- 0 }
  # init of init
  coef_rough_init <- rnorm(K_up * (p + q + pq + 1))
  res.mle <- fmrs.mle(y = whole.data$data_y,
                      x = whole.data$X_all,
                      delta = rep(0, length(whole.data$data_y)),
                      nComp = K_up,
                      disFamily = "norm",
                      initCoeff = coef_rough_init,
                      initDispersion = rep(1, K_up),
                      initmixProp = rep(1 / K_up, K_up),
                      nIterNR = 200,
                      activeset = as)
  res.lam <- fmrs.tunsel(y = whole.data$data_y, 
                         x = whole.data$X_all, 
                         delta = rep(0, length(whole.data$data_y)),
                         nComp = K_up, 
                         disFamily = "norm",
                         initCoeff = c(coefficients(res.mle)),
                         initDispersion = dispersion(res.mle),
                         initmixProp = mixProp(res.mle),
                         penFamily = "mcp", 
                         nIterNR = 200,
                         activeset = as)
  # res.lam@lambPen
  res.var <- fmrs.varsel(y = whole.data$data_y, 
                         x = whole.data$X_all, 
                         delta = rep(0, length(whole.data$data_y)),
                         nComp = ncomp(res.mle), 
                         disFamily = "norm",
                         initCoeff = c(coefficients(res.mle)),
                         initDispersion = dispersion(res.mle),
                         initmixProp = mixProp(res.mle),
                         penFamily = "mcp",
                         lambPen = slot(res.lam, "lambPen"), 
                         nIterNR = 200)
  # res.var@dispersion
  # res.var@mixProp
  
  # cluster result
  fmr_weight_ci <- apply(res.var@weights, 1, which.max)
  fmr_res_ci <- apply(abs(res.var@residuals), 1, which.min) 
  fmr_res <- c()
  fmr_weight <- c()
  fmr_weight_y <- c()
  fmr_res_y <- c()
  for (i in 1:n) {
    fmr_res[i] <- res.var@residuals[i, fmr_res_ci[i]]
    fmr_weight[i] <- res.var@weights[i, fmr_weight_ci[i]]
    fmr_res_y[i] <- res.var@fitted[i, fmr_res_ci[i]]
    fmr_weight_y[i] <- res.var@fitted[i, fmr_weight_ci[i]]
  }
  fmr_var <- sqrt(sum(fmr_res^2))
  fmr_bic <- res.var@BIC
  fmr_coe <- res.var@coefficients
  if(int_rm){ fmr_coe <- res.var@coefficients[-1,]}
  coef_init <- as.vector(fmr_coe)
  
  rho_k_weight <- c()
  rho_k_res <- c()
  mu_k_weight <- c()
  mu_k_res <- c()
  for(k in 1:4){
    rho_k_weight[k] <- 1/sd(whole.data$data_y[fmr_weight_ci == k])
    rho_k_res[k] <- 1/sd(whole.data$data_y[fmr_res_ci == k])
    mu_k_weight[k] <- mean(whole.data$data_y[fmr_weight_ci == k])
    mu_k_res[k] <- mean(whole.data$data_y[fmr_res_ci == k])
  }
  
  return(list(ci_weight = fmr_weight_ci,
              ci_res = fmr_res_ci,
              y_weight = fmr_weight_y,
              y_res = fmr_res_y,
              mu_k_weight = mu_k_weight,
              mu_k_res = mu_k_res,
              pi_k = res.var@mixProp,
              rho_k_weight = rho_k_weight,
              rho_k_res = rho_k_res,
              info = list(bic = fmr_bic,
                          coef = fmr_coe,
                          sd = fmr_var),
              model = list(mle = res.mle,
                           lam = res.lam,
                           var = res.var)))
}

# 评价分类效果
sc <- function(ci_est, ci_true){
  if(length(ci_est) != length(ci_true)){
    stop("Wrong len of ci's")
  }
  n <- length(ci_true)
  comb_matrix <- combn(n, 2)
  scs <- NULL
  for(pair_idx in 1:ncol(comb_matrix)){
    scs <- c(scs, (ci_est[comb_matrix[,pair_idx]][1] == ci_est[comb_matrix[,pair_idx]][2]) ==
               (ci_true[comb_matrix[,pair_idx]][1] == ci_true[comb_matrix[,pair_idx]][2]))
  }
  return(mean(scs))
}

# 向量按列填充到矩阵(如果传入矩阵先按列拉伸为向量)
# 将新矩阵按行求和
row_lap <- function(x, period){
  mat <- matrix(x, nrow = period)
  return(apply(mat, 1, sum))
}

# 取正部
positive_part <- function(x){
  x[x < 0] <- 0
  return(x)
}

# beta(kp*1) -> kpq*1，q个p维beta_1,...,q个p维beta_K
beta_reform <- function(bt, p, q){
  K <- length(bt) %/% p
  res <- c()
  for(k in 1:K){
    beta_k <- bt[((k-1)*p+1):(k*p)]
    res <- c(res, kronecker(matrix(1,nrow=q,ncol=1), beta_k))
  }
  return(res)
}

get_q_c_matrix <- function(n, K_up, ci_sim){
  q_c_matrix <- matrix(0, nrow = n, ncol = K_up)
  for(i in 1:n){
    q_c_matrix[i,ci_sim[i]] <- 1
  }
  return(q_c_matrix)
}
# get_q_c_matrix(n, K_up, ci_sim)
get_q_c_matrix_r <- function(q_c_matrix, mean = 0.2, sd = 0.1, seed_num = 999){
  set.seed(seed_num)
  q_c_matrix_r <- q_c_matrix + matrix(rnorm(n*K_up, mean = mean, sd = 0.1), nrow = n)
  q_c_matrix_r[q_c_matrix_r < 0] <- 0
  q_c_matrix_r <- q_c_matrix_r/apply(q_c_matrix_r, 1, sum)
  return(q_c_matrix_r)
}
# get_q_c_matrix_r(get_q_c_matrix(n, K_up, ci_sim))

