
rm(list = ls(all = TRUE))
library(ggplot2)
library(dplyr)
library(Matrix)
library(flexmix)
source("sim.R")
source("tools.R")
source("func.R")

# 超参数设定
n <- 200
p <- 8
q <- 4
cotype_x <- "En"
cotype_z <- "En"
epsilon_sd <- 0.0
beta_nonzero <- c(-2, -2, 2, 2) # 长度应和真实 group_num_sub 保持一致
# beta_nonzero <- c(-3, -1, 1, 3) # 长度应和真实 group_num_sub 保持一致
alpha_nonzero <- c(-3, -1, 1, 3)
# alpha_nonzero <- c(-2, -2, 2, 2) # 验证最简单情况不需要 alpha
# alpha_nonzero <- c(0, 0, 0, 0) # 验证最简单情况不需要 alpha
beta_vlen <- 3
alpha_vlen <- 2

K_up <- 4  # 估计时的最大类别，应该不少于 group_num_sub
group_num_main <- 2                    
group_num_sub <- 4    
hier_struc <- list(c(1,2),c(3,4))
prob_sub <- rep(1/group_num_sub, group_num_sub)  
reverse <- FALSE

aa <- 1.2
lambda_1 <- 0.2

set.seed(9)
whole.data <- generate_all_data(n, p, q, prob_sub, hier_struc, 
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


para_set <- expand.grid(list(q_c_seed = 1:5,
                             lambda_1 = seq(0.1, 0.5, 0.1),
                             lambda_2 = seq(0.0, 2, 0.5),
                             lambda_3 = seq(0.0, 2, 0.5),
                             aa = seq(1.2, 2.4, 0.6),
                             tau = seq(0, 2, 0.5)))
result <- NULL
for(i in 1:nrow(para_set[,])){
  if(i%%10 == 0){print(i);print("*************")}
  para <- para_set[i,]
  td <- ADMM_trail(aa = para$aa,
                   tau = para$tau,
                   lambda_1 = para$lambda_1,
                   lambda_2 = para$lambda_2,
                   lambda_3 = para$lambda_3,
                   q_c_seed = para$q_c_seed)
  result <- rbind(result, c(td$cdist, td$ci_prob_mean, td$mse))
}
colnames(result) <- c("cdist", "ci_prob_mean", "mse")
res <- cbind(para_set[1:nrow(result),], result)
write.csv(res, file = "2023-07-15-trail.csv", row.names = F)














