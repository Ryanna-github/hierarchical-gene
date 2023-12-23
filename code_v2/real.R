library(dplyr)
library(Matrix)
library(flexmix)
# source("tools.R")
# source("func.R")
source("hierarchical-gene/code_v2/tools.R")
source("hierarchical-gene/code_v2/func.R")
library(argparse)
parser <- ArgumentParser()

parser$add_argument("-e", "--epsilon_sd", default = 0.5, help = "error")
parser$add_argument("--path", default="temp.csv",  help="csv result save path")
parser$add_argument("--K_up", default=8,  help="Upper class number")

args <- parser$parse_args()
epsilon_sd <- as.numeric(args$e)
sigma_est <- as.numeric(args$e)
save_path <- args$path
K_up <- as.numeric(args$K_up)

# data
# X <- read.csv("../data/realX_image_6.csv") %>% as.matrix()
# Z <- read.csv("../data/realZ_gene_20.csv") %>% as.matrix()
# y <- read.csv("../data/y.csv")$fev1

X <- read.csv("hierarchical-gene/data/realX_image_6.csv") %>% as.matrix()
Z <- read.csv("hierarchical-gene/data/realZ_gene_20.csv") %>% as.matrix()
y <- read.csv("hierarchical-gene/data/y.csv")$fev1
data <- cbind(X, Z)

n <- 160
p <- 6
q <- 20
dt_seed <- NaN
q_c_seed_max <- 10

if(0){
  K_up <- 8
  save_path <- "temp.csv"
  epsilon_sd <- 0.4
  sigma_est <- as.numeric(epsilon_sd)
  q_c_seed <- 9
}


comb_pair <- combn(K_up, 2)
H_p <- kronecker(t(apply(comb_pair, 2, get_e_mat, K_up)),
                 diag(p)) %>% Matrix(sparse = TRUE)
H_q <- kronecker(t(apply(comb_pair, 2, get_e_mat, K_up)),
                 diag(q)) %>% Matrix(sparse = TRUE)


# =============================== result =================================
colnames_all <- c("dt_seed", "q_c_seed", "aa", "tau", "l1", "l2", "l3",
                  "cdist", "ci_prob_mean", "mse", "sc", "fit_sum", "fit_mean",
                  "penal", "bic_sum", "bic_mean", "main_grn", "sub_grn", "valid_hier", 
                  "group_detail", paste0("case_", 1:4), "tag")

for(q_c_seed in 1:q_c_seed_max){
  
  result <- as.data.frame(matrix(NaN, nrow = 2, ncol = length(colnames_all)))
  colnames(result) <- colnames_all
  result$dt_seed <- dt_seed
  
  # flexmix for initialization
  flemix_forinit <- tryCatch({
    flexmix_init(q_c_seed, 0)
  }, error = function(err) {
    cat("Random Init\n")
    random_init(q_c_seed)
  })
  result[1,'q_c_seed'] <- q_c_seed
  result[1,'sub_grn'] <- flemix_forinit$est_sub_grn
  result[1,'sc'] <- flemix_forinit$sc_score
  result[1,'cdist'] <- flemix_forinit$cdist
  result[1,'tag'] <- flemix_forinit$tag
  
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
  result[2,'cdist'] <- flemix_best$cdist
  result[2,'tag'] <- flemix_best$tag
  
  
  # our method
  l2_seq <- c(0, 0.5, 1, 2, 4)
  l3_seq <- c(0, 0.5, 1, 2, 4)
  # l2_seq <- c(0, 1, 3, 5, 7)
  # l3_seq <- c(0, 2, 4, 6, 8, 10)
  # l2_seq <- c(0.5)
  # l3_seq <- c(0.5)
  fix_para <- list(dt_seed = dt_seed, q_c_seed = q_c_seed, lambda_1 = 0.3,
                   aa = 1.2, tau = 1)
  hp <- tuning_hyper(l2_seq, l3_seq, fix_para, flemix_forinit$coef_full_ori,
                     save_all = TRUE)
  colnames(hp) <- colnames_all
  result <- rbind(result, hp)
  write_csv(result, file=save_path, col_names=!file.exists(save_path), append=TRUE)
}





