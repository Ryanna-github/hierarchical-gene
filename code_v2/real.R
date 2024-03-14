library(dplyr)
library(Matrix)
library(readr)
library(glue)
library(flexmix)
source("tools.R")
source("func.R")
# source("hierarchical-gene/code_v2/tools.R")
# source("hierarchical-gene/code_v2/func.R")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-n", default = 160)
parser$add_argument("-p", default = 6)
parser$add_argument("-q", default = 30)
parser$add_argument("-e", "--epsilon_sd", default = 0.5, help = "error")
parser$add_argument("--path", default="temp.csv",  help="csv result save path")
parser$add_argument("--K_up", default = 8,  help="Upper class number")
parser$add_argument("--y_type", default = "fev1", help = "fev1 or dlco")

args <- parser$parse_args()
n <- as.numeric(args$n)
p <- as.numeric(args$p)
q <- as.numeric(args$q)
epsilon_sd <- as.numeric(args$e)
sigma_est <- as.numeric(args$e)
save_path <- args$path
K_up <- as.numeric(args$K_up)
y_type <- args$y_type

print(str(args))

dt_seed <- NaN
q_c_seed_max <- 10
rho_ratio <- 0.2

if(0){
  n <- 160
  p <- 6
  q <- 20
  epsilon_sd <- 0.5
  sigma_est <- 0.5
  save_path <- "temp.csv"
  K_up <- 8
  y_type <- "fev1"
}

# data
# X <- read.csv("../data/realX_image_6.csv") %>% as.matrix()
# Z <- read.csv(glue("../data/realZ_gene_{q}.csv")) %>% as.matrix()
# y <- read.csv("../data/y.csv")$fev1

X <- read.csv("../data/realX_v2_image_6.csv") %>% as.matrix()
Z <- read.csv(glue("../data/realZ_v2_gene_{q}.csv")) %>% as.matrix()
y <- read.csv("../data/realy_v2.csv")

cat("Using", y_type)
col <- rlang::sym(y_type)
y <- y[[col]]
# print(y)


# X <- read.csv("hierarchical-gene/data/realX_image_6.csv") %>% as.matrix()
# Z <- read.csv("hierarchical-gene/data/realZ_gene_20.csv") %>% as.matrix()
# y <- read.csv("hierarchical-gene/data/y.csv")$fev1
data <- cbind(X, Z)

# if(0){
#   K_up <- 8
#   save_path <- "temp.csv"
#   epsilon_sd <- 0.4
#   sigma_est <- as.numeric(epsilon_sd)
#   q_c_seed <- 9
# }

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
  result[1,'sc_sub'] <- flemix_forinit$sc_score
  result[1,'ari'] <- flemix_forinit$ari_score
  result[1,'ari_sub'] <- flemix_forinit$ari_score
  result[1,'mse'] <- flemix_forinit$mse
  result[1,'cdist'] <- flemix_forinit$cdist
  result[1,'cdist_sub'] <- flemix_forinit$cdist_sub
  result[1,'cdist_main'] <- flemix_forinit$cdist_main
  result[1,'tag'] <- flemix_forinit$tag
  q_c_matrix_init <- flemix_forinit$q_c_matrix

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

  l2_seq <- c(0.5, 1, 1.5, 2, 2.5)
  l3_seq <- c(1, 2, 3, 4)
  fix_para <- list(dt_seed = dt_seed, q_c_seed = q_c_seed, lambda_1 = 0.3,
                   aa = 1.2, tau = 1)
  hp <- tuning_hyper(l2_seq, l3_seq, fix_para, flemix_forinit$coef_full_ori,
                     save_all = TRUE)
  colnames(hp) <- colnames_all
  result <- rbind(result, hp)

  write_csv(result, file=save_path, col_names=!file.exists(save_path), append=TRUE)
}





