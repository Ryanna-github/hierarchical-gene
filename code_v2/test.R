# auto
library(ggplot2)
library(dplyr)
library(Matrix)
library(flexmix)
# library(ncvreg)
# library(mclust)
library(readr)
source("sim.R")
source("tools.R")
source("func.R")
# source("hierarchical-gene/code_v2/sim.R")
# source("hierarchical-gene/code_v2/tools.R")
# source("hierarchical-gene/code_v2/func.R")
# library(argparse)
# #
# # 创建参数解析对象
# parser <- ArgumentParser()
# # Rscript test.R -n 120 --dt_seed 9 -p 8 -q 4 -e 0.5 --path 2023-07-28_eps05_fminit.csv
# 
# parser$add_argument("-n", "--num", default=200,  help="sample size")
# parser$add_argument("--dt_seed", default=9,  help="epsilon seed while generating data")
# parser$add_argument("-p", default = 8, help = "dim of X")
# parser$add_argument("-q", default = 4, help = "dim of Z")
# parser$add_argument("-b", default = 1, help = "if data balanced (1: balanced, 2: unbalanced groups with balanced subgroups, 3: balanced groups with unbalanced sbgroups)")
# parser$add_argument("-e", "--epsilon_sd", default = 0.5, help = "error")
# parser$add_argument("--epsilon_sd_init", default = 0.5, help = "epsilon_sd initiation for estimating rho")
# parser$add_argument("--rho_ratio", default = 0.5, help = "update ratio of rho, rho(t+1)=ratio*rho+(1-ratio)rho(t)")
# parser$add_argument("-ss", "--signal_size", default = 1, help = "signal size")
# parser$add_argument("-bl", "--beta_vlen", default = 3, help = "nonzero length of beta")
# parser$add_argument("-al", "--alpha_vlen", default = 2, help = "nonzero length of alpha")
# parser$add_argument("--path", default="temp.csv",  help="csv result save path")
# parser$add_argument("--K_up", default=4,  help="Upper class number")
# parser$add_argument("--cotype", default="En",  help="Covariate matrix of X and Z")
# # parser$add_argument("-a", "--aa", default = 1.2, help = "penalty para in MCP")
# # parser$add_argument("--lambda_1", default = 0.3, help = "lambda 1")
# # parser$add_argument("--lambda_2", default = 2, help = "lambda 2")
# # parser$add_argument("--lambda_3", default = 5, help = "lambda 3")
# # parser$add_argument("--tau", default = 0.5, help = "ADMM penalty")
# 
# # Rscript --verbose test.R -n 500 --dt_seed 9 -p 8 -q 4 --epsilon_sd 0.5
# # --epsilon_sd_init 0.5 --beta_vlen 3 --alpha_vlen 2 --K_up 4
# 
# args <- parser$parse_args()
# print(str(args))
# 
# n <- as.numeric(args$n)
# p <- as.numeric(args$p)
# q <- as.numeric(args$q)
# balance <- as.numeric(args$b)
# epsilon_sd <- as.numeric(args$epsilon_sd)
# sigma_est <- as.numeric(args$epsilon_sd_init)
# rho_ratio <- as.numeric(args$rho_ratio)
# signal_size <- as.numeric(args$signal_size)
# beta_vlen <- as.numeric(args$beta_vlen)
# alpha_vlen <- as.numeric(args$alpha_vlen)
# save_path <- args$path
# dt_seed <- as.numeric(args$dt_seed)
# K_up <- as.numeric(args$K_up)
# cotype <- args$cotype

if(1){
  n <- 500
  p <- 80
  q <- 40
  balance <- 1
  epsilon_sd <- 0.5
  epsilon_sd_init <- 0.5
  sigma_est <- as.numeric(epsilon_sd_init)
  rho_ratio <- 0.1
  signal_size <- 1
  # beta_vlen <- 3
  beta_vlen <- 30
  # alpha_vlen <- 2
  alpha_vlen <- 20
  save_path <- "temp.csv"
  dt_seed <- 9
  K_up <- 4  # 估计时的最大类别，应该不少于 group_num_sub
  cotype <- "En"
  print("*")
}

# 超参数设定
# n <- 200
# p <- 8
# q <- 4
cotype_x <- cotype
cotype_z <- cotype
# epsilon_sd <- 0.5
beta_nonzero <- c(-2, -2, 2, 2)*signal_size # 长度应和真实 group_num_sub 保持一致
# beta_nonzero <- c(-3, -1, 1, 3) # 长度应和真实 group_num_sub 保持一致
alpha_nonzero <- c(-3, -1, 1, 3)*signal_size
# alpha_nonzero <- c(-2, -2, 2, 2) # 验证最简单情况不需要 alpha
# alpha_nonzero <- c(0, 0, 0, 0) # 验证最简单情况不需要 alpha
# beta_vlen <- 3
# alpha_vlen <- 2

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

# aa <- 1.2
# lambda_1 <- 0.2
# =============================== data =================================
# set.seed(9)
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


# =============================== result =================================
colnames_all <- c("dt_seed", "q_c_seed", "aa", "tau", "l1", "l2", "l3",
                  "cdist", "cdist_main", "cdist_sub", "ci_prob_mean", "mse", 
                  "sc", "sc_main", "sc_sub", 
                  "ari", "ari_main", "ari_sub", "fit_sum", "fit_mean",
                  "penal", "bic_sum", "bic_mean", "main_grn", "sub_grn", 
                  "rho_init", "rho_est", "rho_ratio", "pi_est", "valid_hier", 
                  "group_detail", paste0("case_", 1:4), 
                  "iter_total", "iter_type", "tag")

q_c_seed_max <- 1
for(q_c_seed in 1:10){

  result <- as.data.frame(matrix(NaN, nrow = 2, ncol = length(colnames_all)))
  colnames(result) <- colnames_all
  result$dt_seed <- dt_seed

  # flexmix for initialization
  miniprior <- ifelse(K_up == 4, 0, 0.1)
  flemix_forinit <- tryCatch({
    flexmix_init(q_c_seed, miniprior)
  }, error = function(err) {
    # cat("Random Init\n")
    # random_init(q_c_seed)
    cat(paste0("Fail to init for dt_seed ", dt_seed," q_c_seed ", q_c_seed, ", continue...\n"))
    FALSE
  })
  # 若未初始化成功，直接用下一个随机种子
  if(isFALSE(flemix_forinit)){
    next
  }
  result[1,'q_c_seed'] <- q_c_seed
  result[1,'sub_grn'] <- flemix_forinit$est_sub_grn
  result[1,'sc'] <- flemix_forinit$sc_score
  result[1,'sc_sub'] <- flemix_forinit$sc_score
  result[1,'ari'] <- flemix_forinit$ari_score
  result[1,'ari_sub'] <- flemix_forinit$ari_score
  result[1,'mse'] <- flemix_forinit$mse
  result[1,'cdist'] <- flemix_forinit$cdist
  result[1,'cdist_sub'] <- flemix_forinit$cdist_sub
  result[1,'cdist_main'] <- flemix_forinit$cdist_main
  result[1,'tag'] <- flemix_forinit$tag
  q_c_matrix_init <- flemix_forinit$q_c_matrix
  
  if(K_up != flemix_forinit$est_sub_grn){
    cat(paste("K_up changes from", K_up, "to", flemix_forinit$est_sub_grn))
    K_up <- flemix_forinit$est_sub_grn
    
    comb_pair <- combn(K_up, 2)
    H_p <- kronecker(t(apply(comb_pair, 2, get_e_mat, K_up)),
                     diag(p)) %>% Matrix(sparse = TRUE)
    H_q <- kronecker(t(apply(comb_pair, 2, get_e_mat, K_up)),
                     diag(q)) %>% Matrix(sparse = TRUE)
    # flexmix 大组数目估计过小不需要再继续运算了
    if(K_up < 4){
      next
    }
  }
  
  # flexmix best result (K squeezed)
  flemix_best <- tryCatch({
    flexmix_init(q_c_seed, 0.1)
  }, error = function(err) {
    cat("Random Init\n")
    random_init(q_c_seed)
  })
  result[2,'q_c_seed'] <- q_c_seed
  result[2,'sub_grn'] <- flemix_best$est_sub_grn
  result[2,'sc'] <- flemix_best$sc_score
  result[2,'sc_sub'] <- flemix_best$sc_score
  result[2,'ari'] <- flemix_best$ari_score
  result[2,'ari_sub'] <- flemix_best$ari_score
  result[2,'mse'] <- flemix_best$mse
  result[2,'cdist'] <- flemix_best$cdist
  result[2,'cdist_sub'] <- flemix_best$cdist_sub
  result[2,'cdist_main'] <- flemix_best$cdist_main
  result[2,'tag'] <- flemix_best$tag


  # our method
  l2_seq <- c(0, 0.5, 1, 1.5, 2, 4)
  l3_seq <- c(0, 0.5, 1, 1.5, 2, 4)
  # l2_seq <- c(5,  6,  6.5,  7,  7.5)
  # l3_seq <- c(1.5,2,  2.5,  3)
  # l2_seq <- c(6,7,8,9,10,12,14)
  # l3_seq <- c(2.5,3.5,4.5,5.5,7,8.5)
  lambda_1 = 0.5
  
  if(signal_size == 1 & K_up == 4){
    l2_seq <- c(2.5,3,3.5)
    l3_seq <- c(3.5,4,4.5,5,6)
  }else if(signal_size == 2 & K_up == 4){
    # l2_seq <- c(5,6.5,7,7.5,8)
    # l3_seq <- c(8.5,9,9.5,10)
    l2_seq <- c(6.5,7,7.5,8)
    l3_seq <- c(8.5,9,9.5)
    lambda_1 <- 1
  }else{
    l2_seq <- c(5,6.5,7,7.5,8)
    l3_seq <- c(8,8.5,9,10)
  }
  
  # l2_seq <- c(12, 13, 14)
  # l3_seq <- c(13, 14, 15)
  
  fix_para <- list(dt_seed = dt_seed, q_c_seed = q_c_seed, lambda_1 = lambda_1,
                   aa = 1.2, tau = 1)
  # q_c_matrix 初
  hp <- tuning_hyper(l2_seq, l3_seq, fix_para, flemix_forinit$coef_full_ori,
                     q_c_matrix_init = q_c_matrix_init, save_all = TRUE, add0 = FALSE)
  colnames(hp) <- colnames_all
  result <- rbind(result, hp)
  write_csv(result, file=save_path, col_names=!file.exists(save_path), append=TRUE)
}

print("Done!")

# aa = fix_para$aa
# tau = fix_para$tau
# lambda_1 = fix_para$lambda_1
# lambda_2 = l2_seq[1]
# lambda_3 = l3_seq[1]
# q_c_seed = fix_para$q_c_seed
# coef_full_init = flemix_forinit$coef_full_ori
# rho_ratio = 0.2
# iter_type= "stop"
# iter_max = 200
# rho_clip = 2
# plot_performance = TRUE
# eps = 1e-7
# eps_abs = 1e-2
# eps_rel = 1e-3
result <- data.frame(result)
result$penal2 <- as.numeric(result$penal)/(p+q)
result$bic_mean2 <- as.numeric(result$fit_mean) + result$penal2
