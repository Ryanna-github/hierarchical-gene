
rm(list = ls(all = TRUE))
library(ggplot2)
library(dplyr)
library(Matrix)
library(flexmix)
library(readr)
source("sim.R")
source("tools.R")
source("func.R")


# 超参数设定
n <- 200
p <- 8
q <- 4
cotype_x <- "En"
cotype_z <- "En"
epsilon_sd <- 0.5
beta_nonzero <- c(-2, -2, 2, 2) # 长度应和真实 group_num_sub 保持一致
# beta_nonzero <- c(-3, -1, 1, 3) # 长度应和真实 group_num_sub 保持一致
alpha_nonzero <- c(-3, -1, 1, 3)
# alpha_nonzero <- c(-2, -2, 2, 2) # 验证最简单情况不需要 alpha
# alpha_nonzero <- c(0, 0, 0, 0) # 验证最简单情况不需要 alpha
beta_vlen <- 3
alpha_vlen <- 2
dt_seed <- 9

K_up <- 4  # 估计时的最大类别，应该不少于 group_num_sub
group_num_main <- 2                    
group_num_sub <- 4    
hier_struc <- list(c(1,2),c(3,4))
prob_sub <- rep(1/group_num_sub, group_num_sub)  
reverse <- FALSE

aa <- 1.2
lambda_1 <- 0.2

# set.seed(9)
whole.data <- generate_all_data(dt_seed, n, p, q, prob_sub, hier_struc, 
                                beta_nonzero, alpha_nonzero, beta_vlen, alpha_vlen, 
                                cotype_x, cotype_z, epsilon_sd, reverse)
X <- whole.data$data$X
Z <- whole.data$data$Z
data <- whole.data$data$data_full
coef <- whole.data$coef
coefv <- lapply(coef, as.vector) # 按照类别拉长
ci_sim <- whole.data$ci_sim
y <- matrix(whole.data$y, ncol=1)


comb_pair <- combn(K_up, 2)
H_p <- kronecker(t(apply(comb_pair, 2, get_e_mat, K_up)),
                 diag(p)) %>% Matrix(sparse = TRUE)
H_q <- kronecker(t(apply(comb_pair, 2, get_e_mat, K_up)),
                 diag(q)) %>% Matrix(sparse = TRUE)


para_set <- expand.grid(list(q_c_seed = 1,
                             lambda_1 = 0.3,
                             lambda_2 = c(0.1, 0.2, 0.5),
                             lambda_3 = 2,
                             aa = 1.2,
                             tau = 0.5)) %>% 
  filter(lambda_2 < lambda_3)


result <- NULL
for(i in 1:nrow(para_set[,])){
  para <- para_set[i,]
  # fmres <- flexmix_trail(para$q_c_seed)
  # 改初始化
  # set.seed(99)
  coef_full_init <- fmres[[para$q_c_seed]]$coef_full_ori
  # coef_full_init <- coef$coef_full + rnorm(length(coef$coef_full), 0, 0.5)
  td <- ADMM_trail(aa = para$aa,
                   tau = para$tau,
                   lambda_1 = para$lambda_1,
                   lambda_2 = para$lambda_2,
                   lambda_3 = para$lambda_3,
                   q_c_seed = para$q_c_seed,
                   coef_full_init = coef_full_init)
  result <- rbind(result, c(para$dt_seed,
                            para$q_c_seed,
                            para$aa,
                            para$tau,
                            para$lambda_1,
                            para$lambda_2,
                            para$lambda_3,
                            td$cdist, 
                            td$ci_prob_mean, 
                            td$mse, 
                            td$sc_score))
}


# df <- read.csv("2023-07-24_nonepsilon_around.csv")
# ggplot(df, aes(x = l2, y = l3, fill = cdist, size = 1/cdist)) +
#   geom_point(alpha = 0.3)
# plot(df$l2, df$l3, color = df$cdist)



# 先计算 flexmix 结果
fmres <- list()
for(q_c_seed in 1:5){
  fmres[[q_c_seed]] <- flexmix_trail(q_c_seed)
}



# fmres <- flexmix_trail(1)
# coef_full_init <- fmres$coef_full_ori

save_path <- "2023-07-25_epsilon_fminit.csv"
first_time_write <- TRUE
result <- NULL
for(i in 1:nrow(para_set[,])){
  para <- para_set[i,]
  # fmres <- flexmix_trail(para$q_c_seed)
  # 改初始化
  # set.seed(99)
  coef_full_init <- fmres[[para$q_c_seed]]$coef_full_ori
  # coef_full_init <- coef$coef_full + rnorm(length(coef$coef_full), 0, 0.5)
  td <- ADMM_trail(aa = para$aa,
                   tau = para$tau,
                   lambda_1 = para$lambda_1,
                   lambda_2 = para$lambda_2,
                   lambda_3 = para$lambda_3,
                   q_c_seed = para$q_c_seed,
                   coef_full_init = coef_full_init)
  result <- rbind(result, c(para$dt_seed,
                            para$q_c_seed,
                            para$aa,
                            para$tau,
                            para$lambda_1,
                            para$lambda_2,
                            para$lambda_3,
                            td$cdist, 
                            td$ci_prob_mean, 
                            td$mse, 
                            td$sc_score))
  if(i%%5 == 0){
    print(paste(i, "**********************************"))
    colnames(result) <- c("dt_seed", "q_c_seed", "aa", "tau", "l1", "l2", "l3",
                          "cdist", "ci_prob_mean", "mse", "sc")
    write_csv(data.frame(result), file=save_path, col_names=first_time_write, append=TRUE)
    first_time_write <- FALSE
    result <- NULL
  }
}
if(! is.null(result)){
  print(paste(i, "**********************************"))
  colnames(result) <- c("q_c_seed", "aa", "tau", "l1", "l2", "l3",
                        "cdist", "ci_prob_mean", "mse", "sc")
  write_csv(data.frame(result), file=save_path, col_names=first_time_write, append=TRUE)
  first_time_write <- FALSE
  result <- NULL
}







save_path <- "2023-07-24_nonepsilon_random.csv"
first_time_write <- TRUE
result <- NULL
for(i in 1:nrow(para_set[,])){
  para <- para_set[i,]
  fmres <- flexmix_trail(para$q_c_seed)
  # 改初始化
  set.seed(99)
  coef_full_init <- coef$coef_full - coef$coef_full + 
    rnorm(length(coef$coef_full), 0, 1)
  td <- ADMM_trail(aa = para$aa,
                   tau = para$tau,
                   lambda_1 = para$lambda_1,
                   lambda_2 = para$lambda_2,
                   lambda_3 = para$lambda_3,
                   q_c_seed = para$q_c_seed,
                   coef_full_init = coef_full_init)
  result <- rbind(result, c(para$q_c_seed,
                            para$aa,
                            para$tau,
                            para$lambda_1,
                            para$lambda_2,
                            para$lambda_3,
                            td$cdist, 
                            td$ci_prob_mean, 
                            td$mse, 
                            td$sc_score))
  if(i%%5 == 0){
    print(paste(i, "**********************************"))
    colnames(result) <- c("q_c_seed", "aa", "tau", "l1", "l2", "l3",
                          "cdist", "ci_prob_mean", "mse", "sc")
    write_csv(data.frame(result), file=save_path, col_names=first_time_write, append=TRUE)
    first_time_write <- FALSE
    result <- NULL
  }
}
if(! is.null(result)){
  print(paste(i, "**********************************"))
  colnames(result) <- c("q_c_seed", "aa", "tau", "l1", "l2", "l3",
                        "cdist", "ci_prob_mean", "mse", "sc")
  write_csv(data.frame(result), file=save_path, col_names=first_time_write, append=TRUE)
  first_time_write <- FALSE
  result <- NULL
}


