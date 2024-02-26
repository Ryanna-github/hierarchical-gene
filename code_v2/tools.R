library(gtools) # for permutation in coef_dist

positive_part <- function(x){
  return(ifelse(x>0, x, 0))
}

negative_part <- function(x){
  return(ifelse(x<0, x, 0))
}

soft_threshold <- function(x, t){
  return(sign(x)*(positive_part(abs(x)-t)))
}

mcp_value <- function(x, gamma, lambda, element_wise = TRUE){
  if(!element_wise){x <- norm(matrix(x), type = "2")}
  res <- ifelse(abs(x)<=gamma*lambda, 
                lambda*abs(x) - (x**2)/(2*gamma), 
                (gamma*lambda**2)/(2))
  return(res)
}

# Partial derivatives of MCP penalty
mcp_d <- function(x, gamma, lambda, element_wise = TRUE){
  if(!element_wise){x <- norm(matrix(x), type = "2")}
  res <- ifelse(abs(x)<=gamma*lambda, (lambda - abs(x)/gamma), 0)
  return(res)
}


# standard \theta/2(x-a)^2 + pen(x;lambda,gamma) solution
mcp_solution <- function(theta, a, gamma, lambda){
  key <- norm(matrix(a), type = "2")
  if(key <= gamma*lambda){
    res <- soft_threshold(a, lambda/theta)/(1-1/(gamma*theta))
  }else{
    res <- a
  }
  return(res)
}

# =============================== 评价指标 ================================
## -------------------------------- 分类 ----------------------------------
# adj.rand.index 为 fossil 包中函数，服务器无法下载
adj.rand.index <- function (group1, group2) 
{
  a <- length(table(group1))
  N <- length(group1)
  ctab <- matrix(, a, a)
  for (j in 1:a) {
    for (i in 1:a) {
      ctab[j, i] <- length(which(group2[which(group1 == 
                                                i)] == j))
    }
  }
  sumnij <- sum(choose(ctab, 2))
  sumai <- sum(choose(colSums(ctab), 2))
  sumbj <- sum(choose(rowSums(ctab), 2))
  Ntwo <- choose(N, 2)
  ari <- abs((sumnij - (sumai * sumbj)/Ntwo)/(0.5 * (sumai + 
                                                       sumbj) - (sumai * sumbj)/Ntwo))
  return(ari)
}

# rand.index 为 fossil 包中函数，服务器无法下载
rand.index <- function (group1, group2) 
{
  x <- abs(sapply(group1, function(x) x - group1))
  x[x > 1] <- 1
  y <- abs(sapply(group2, function(x) x - group2))
  y[y > 1] <- 1
  sg <- sum(abs(x - y))/2
  bc <- choose(dim(x)[1], 2)
  ri <- 1 - sg/bc
  return(ri)
}

sc <- function(ci_est, ci_true){
  if(length(ci_est) != length(ci_true)){
    stop("Wrong len of ci's")
  }
  return(rand.index(ci_est, ci_true))
  # n <- length(ci_true)
  # comb_matrix <- combn(n, 2)
  # scs <- NULL
  # for(pair_idx in 1:ncol(comb_matrix)){
  #   scs <- c(scs, (ci_est[comb_matrix[,pair_idx]][1] == ci_est[comb_matrix[,pair_idx]][2]) ==
  #              (ci_true[comb_matrix[,pair_idx]][1] == ci_true[comb_matrix[,pair_idx]][2]))
  # }
  # return(mean(scs))
}

ari <- function(ci_est, ci_true){
  if(length(ci_est) != length(ci_true)){
    stop("Wrong len of ci's")
  }
  # return(adjustedRandIndex(ci_est, ci_true)) # mclust 
  return(adj.rand.index(ci_est, ci_true))
}

## -------------------------------- MSE ----------------------------------
mse <- function(y_true, y_est){
  return(sum((y_true-y_est)**2))
}

mse_2 <- function(q_c_matrix, y, X, coef_beta_est){
  n <- nrow(q_c_matrix)
  K_up <- ncol(q_c_matrix)
  res <- NULL
  for(k in 1:K_up){
    W_k <- diag(q_c_matrix[,k])
    res[k] <- t(y-X%*%coef_beta_est[,k])%*%W_k%*%(y-X%*%coef_beta_est[,k])
  }
  return(sum(res)/n)
}

## ----------------------------- Coef Dist -------------------------------
coef_dist <- function(coef_est, coef_true){
  K_up <- ncol(coef_true)
  if(is.null(K_up)){
    # 单类，数组
    return(sum((coef_true - coef_est)**2)/length(coef_true))
  }
  allposs <- permutations(n = K_up, r = K_up)
  miniest <- 1000 * prod(dim(coef_true))
  for(id_pos in 1:nrow(allposs)){
    mini <- sum((coef_true[,allposs[id_pos,]] - coef_est)**2)
    miniest <- ifelse(mini < miniest, mini, miniest)
  }
  return(miniest/prod(dim(coef_true)))
}


get_e_mat <- function(indx, len){
  e <- matrix(0, len, 1)
  e[indx[1],] <- 1
  e[indx[2],] <- (-1)
  return(e)
}


extend_x_to_row <- function(x, row_to){
  K_up <- length(x)
  return(matrix(kronecker(x, matrix(1, nrow=row_to)), ncol = K_up))
}

# ========================== graph for compact group ==============================

# m1 <- matrix(0,6,6)
# m1[1,3] <- 1
# m1[3,1] <- 1
# m1[1,5] <- 1
# m1[5,1] <- 1
# m1[2,4] <- 1
# m1[4,2] <- 1
# m1[2,6] <- 1
# m1[6,2] <- 1
# m2 <- matrix(0,6,6)
# m2[1,3] <- 1
# m2[3,1] <- 1
# m2[2,4] <- 1
# m2[4,2] <- 1
# main_group_info_compact <- find_connected_nodes(m1)
# sub_group_info_compact <- find_connected_nodes(m2)

dfs <- function(adj_matrix, node, visited) {
  visited[node] <- TRUE
  connected_nodes <- c(node)
  
  for (neighbor in which(adj_matrix[node, ] == 1)) {
    if (!visited[neighbor]) {
      connected_nodes <- c(connected_nodes, dfs(adj_matrix, neighbor, visited))
    }
  }
  return(connected_nodes)
}

find_connected_nodes <- function(adj_matrix) {
  n <- nrow(adj_matrix)
  visited <- rep(FALSE, n)
  connected_components <- list()
  for (i in 1:n) {
    if (!visited[i]) {
      connected_nodes <- dfs(adj_matrix, i, visited)
      connected_components <- c(connected_components, list(sort(unique(connected_nodes))))
    }
  }
  return(unique(connected_components))
}

convert_to_parentheses <- function(lst) {
  result <- sapply(lst, function(x) paste0("(", paste(x, collapse = ","), ")"))
  paste(result, collapse = ",")
}
