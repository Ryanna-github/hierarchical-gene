################################################################################
# Codes for conducting simulation studies
################################################################################
rm(list = ls(all = TRUE))
ls()
###############################

library(Matrix)
library(MASS)
library(fmrs)
source("function.R")
source("sim-func.R")


n <- 120                         # The sample size
q <- 5                           # The dimension of x
p <- 10                          # The dimension of z
pp <- q+p                        # The dimension of each sample
S_E <- 2                         # The dimension of non-zero coefficients of x
S_GE <- 4                        # The dimension of non-zero coefficients of z
E.cat <- 2                       # The dimension of non-zero coefficients of z
beta00 <- 2                      # The signal intensity of non-zero coefficients
balance <- F                     # Balanced/Imbalanced subgroup
gr.main <- 2                     # The number of rough subgroup 
gr.GE <- 4                       # The number of refined subgroup 
num_main <- list(c(1,2),c(3,4))  # Hierarchical subgroup structure
sd.epsi<-0.5                     # The residual error
corZX <- F                       # Correlated/Uncorrelated X and Z
intercorr0<-"AR1"                # The correlation structure of Z


# ------------ Necessary parameters to support algorithm implementation --------
sample.index <- combn(n,2)
a.matrix <- kronecker(Matrix(t(apply(sample.index,2,dMatrixFun)),sparse = T),bdiag(diag(pp)))
a.matrix.main <- kronecker(Matrix(t(apply(sample.index,2,dMatrixFun)),sparse = T),bdiag(diag(q)))
a.matrix.GE <- kronecker(Matrix(t(apply(sample.index,2,dMatrixFun)),sparse = T),bdiag(diag(p)))
one.matrix <- bdiag(rep(list(rep(1,pp)),(n*(n-1)/2)))
one.matrix.main <- bdiag(rep(list(rep(1,q)),(n*(n-1)/2)))
one.matrix.GE <- bdiag(rep(list(rep(1,p)),(n*(n-1)/2)))

# ------------ Generate the regression coefficients -----------
beta_main.vec <- c(-beta00,beta00)
beta_GE.vec <- c(-1.5*beta00,-0.5*beta00,1.5*beta00,0.5*beta00)
if(balance){pr_GE<-rep(1/gr.GE,gr.GE)}else{pr_GE<-c(1/6,1/6,1/3,1/3)}
group.num <- group_GE(n, pr_GE=pr_GE, num_main=num_main)
group.num_main <- group.num$group.num_main
group.num_GE <- group.num$group.num_GE
Beta0 <- Generatebeta(n, q, p, group.num_main, group.num_GE, beta_main.vec, beta_GE.vec, sparse=T, S_GE=S_GE, S_E=S_E)
beta.true.list <- Beta0$beta.true.list


# ------------ Generate data ------------ 
whole.data <- GenerateData(n, p, q, beta.true.list, maincorr="AR_E", intercorr=intercorr0,
                        sd.epsi=sd.epsi, corZX = corZX, E.cat = E.cat)

# ------------ Apply the proposed method ------------ 
beta.init.list <- gen_int_beta(n, p, q, whole.data)
beta.init <- beta.init.list$beta.init
lambda <- genelambda.obo()
result <- tuning.lambda(lambda, whole.data, n, q, p, beta.init)
index.list <- evaluation.tuning(n,q,p, result$admmres, result$abic.n, result$admmres2, result$abic.n2, Beta0, gr.GE=gr.GE, result$bic.var)
index.list$err.s

bt <- result$beta.over
bt <- matrix(bt, ncol = p+q)



data.total <- whole.data
lambda1 <- mean(lambda$lambda1)
lambda2 <- mean(lambda$lambda2)
p_y <- 1
iter.max=1
epsi=0.1;gam.mcp=3;penal.para=1;merge.all=F


pp=q+p
main.n=c(1:q)
GE.n=setdiff(1:pp,main.n)
a.main.n = sort(rep((0:(n-1))*pp,q))+main.n
a.GE.n = sort(rep((0:(n-1))*pp,p))+GE.n
d.main.n =sort(rep((0:(n*(n-1)/2-1))*pp,q)) + main.n
d.GE.n = sort(rep((0:(n*(n-1)/2-1))*pp,p)) + GE.n

cova<-data.total$data_x
data.response<-data.total$data_y

row_vec<-vector(mode="list",length=n)
for(i5 in 1:n){
  row_vec[[i5]]<-matrix(cova[i5,],1,pp)
}
cova_diag<-bdiag(row_vec)

#################################### Initialize beta #########################################
# group
eta.init<-a.matrix%*%beta.init
gamma.init<-Matrix(0,nrow=p_y,ncol=pp*n*(n-1)/2,sparse = T)

# iteration
iter<-1
# coe
beta.est.list<-vector(mode="list",length=iter.max);beta.est.list[[iter]]<-beta.init
# differences of coe
eta.est.list<-vector(mode="list",length=iter.max);eta.est.list[[iter]]<-eta.init
# the dual variables
gamma.est.list<-vector(mode="list",length=iter.max);gamma.est.list[[iter]]<-gamma.init

while(iter<=iter.max){
  iter<-iter+1
  if(iter>iter.max){
    break
  }
  
  # 更新 beta
  beta.est.list[[iter]]<-solve(t(cova_diag)%*%cova_diag+penal.para*t(a.matrix)%*%a.matrix)%*%
    (t(cova_diag)%*%data.response-t(a.matrix)%*%t(gamma.est.list[[iter-1]])+penal.para*(t(a.matrix)%*%eta.est.list[[iter-1]]))
  
  eta.i = a.matrix%*%beta.est.list[[iter]]
  gamma.i = gamma.est.list[[iter-1]]
  eta.tem1<-eta.i+t(gamma.i)/penal.para
  eta.tem1.main<-as.matrix(eta.tem1[d.main.n,],sparse=T)
  eta.tem1.GE<-as.matrix(eta.tem1[d.GE.n,],sparse=T)
  eta.i1<-eta.i
  
  # 2-norm
  eta.tem2norm.main<-sqrt(t(one.matrix.main)%*%(eta.tem1.main^2))
  eta.tem2norm.GE<-sqrt(t(one.matrix.GE)%*%(eta.tem1.GE^2))
  eta.tem2norm<-sqrt(eta.tem2norm.GE^2 + eta.tem2norm.main^2)
  
  # all:non-shrinkage & main:non-shrinkage
  num_n_n<-which(eta.tem2norm > gam.mcp*lambda1 & eta.tem2norm.main > gam.mcp*lambda2)
  
  # all:shrinkage & main:non-shrinkage
  vu.coe<-apply(1-lambda1/penal.para/eta.tem2norm,1,function(x) max(x,0)) / (1-1/(gam.mcp*penal.para))
  vv<- vu.coe * eta.tem2norm.main
  num_s_n<-which(eta.tem2norm <= gam.mcp*lambda1 & vv > gam.mcp*lambda2)
  
  # all:non-shrinkage & main:shrinkage
  vv.coe<-apply(1-lambda2/penal.para/eta.tem2norm.main,1,function(x) max(x,0)) / (1-1/(gam.mcp*penal.para))
  ww<-sqrt(eta.tem2norm.GE^2 + ( vv.coe * eta.tem2norm.main)^2)
  num_n_s<-which(ww > gam.mcp*lambda1 & eta.tem2norm.main <= gam.mcp*lambda2)
  
  # all:shrinkage & main:shrinkage
  num_s_s<-setdiff(1:(n*(n-1)/2),Reduce(union,list(num_n_n,num_s_n,num_n_s)))
  
  if(length(num_n_n)<2){
    num_n_n2<-which(rowSums(as.matrix(one.matrix[,num_n_n]))!=0)
  }else{num_n_n2<-which(rowSums(one.matrix[,num_n_n])!=0)}
  if(length(num_s_n)<2){
    num_s_n2<-which(rowSums(as.matrix(one.matrix[,num_s_n]))!=0)
  }else{num_s_n2<-which(rowSums(one.matrix[,num_s_n])!=0)}
  if(length(num_n_s)<2){
    num_n_s2<-which(rowSums(as.matrix(one.matrix[,num_n_s]))!=0)
  }else{num_n_s2<-which(rowSums(one.matrix[,num_n_s])!=0)}
  if(length(num_s_s)<2){
    num_s_s2<-which(rowSums(as.matrix(one.matrix[,num_s_s]))!=0)
  }else{num_s_s2<-which(rowSums(one.matrix[,num_s_s])!=0)}
  
  # all:non-shrinkage & main:non-shrinkage
  if(length(num_n_n2) > 0){eta.i1[num_n_n2,]<-eta.i[num_n_n2,]}
  length(num_n_n2)
  # all:shrinkage & main:non-shrinkage
  num = num_s_n; num2 = num_s_n2;
  if(length(num) > 0){
    eta.tem3<-as.matrix(vu.coe[num])
    eta.tem4<-as.vector(apply(eta.tem3, 1, function(x) rep(x,pp)))
    eta.i1[num2,]<-eta.tem4*eta.i[num2,]
  }
  
  # all:non-shrinkage & main:shrinkage
  num = num_n_s; num2 = num_n_s2;
  if(length(num) > 0){
    num.main = num2[sort(rep((1:length(num)-1)*pp,q)) + main.n]
    num.GE = num2[sort(rep((1:length(num)-1)*pp,p)) + GE.n]
    eta.tem3<-as.matrix(vv.coe[num])
    eta.tem4<-as.vector(apply(eta.tem3, 1, function(x) rep(x,q)))
    eta.i1[num.main,]<-eta.tem4*eta.i[num.main,]
    eta.i1[num.GE,]<-eta.i[num.GE,]
  }
  
  # all:shrinkage & main:shrinkage
  num = num_s_s; num2 = num_s_s2;
  if(length(num) > 0){
    num.main = num2[sort(rep((1:length(num)-1)*pp,q)) + main.n]
    num.GE = num2[sort(rep((1:length(num)-1)*pp,p)) + GE.n]
    eta.i0norm<-sqrt(t(one.matrix)%*%(eta.est.list[[iter-1]]^2))
    eta.i0norm.main<-sqrt(t(one.matrix.main)%*%(eta.est.list[[iter-1]][d.main.n,]^2))
    mcp.u<-mcp_d(eta.i0norm[num],lambda1,gam.mcp)/penal.para/(eta.i0norm[num]+10^(-7))
    mcp.v<-mcp_d(eta.i0norm.main[num],lambda2,gam.mcp)/penal.para/(eta.i0norm.main[num]+10^(-7))
    eta.tem3.main<-as.matrix(mcp.u+mcp.v)
    eta.tem3.GE<-as.matrix(mcp.u)
    eta.tem4.main<-as.vector(apply(eta.tem3.main, 1, function(x) rep(x,q)))
    eta.tem4.GE<-as.vector(apply(eta.tem3.GE, 1, function(x) rep(x,p)))
    eta.i1[num.main,]<-eta.i[num.main,]/(1+eta.tem4.main)
    eta.i1[num.GE,]<-eta.i[num.GE,]/(1+eta.tem4.GE)
    eta.i1[abs(eta.i1) < 10^(-9)]=0
  }
  
  eta.est.list[[iter]] = eta.i1
  gamma.est.list[[iter]]<-gamma.est.list[[iter-1]]+penal.para*t(a.matrix%*%beta.est.list[[iter]]-eta.est.list[[iter]])
  eps.group = sqrt(sum((a.matrix%*%beta.est.list[[iter]]-eta.est.list[[iter]])^2))
  
  if(eps.group<epsi){
    break
  }
  
  beta.est.list[iter-1]<-list(NULL)
  eta.est.list[iter-1]<-list(NULL)
  gamma.est.list[iter-1]<-list(NULL)
}

if(iter>iter.max){
  beta.gfl<-beta.est.list[[iter-1]]
  eta.gfl<-eta.est.list[[iter-1]]
}else{
  beta.gfl<-beta.est.list[[iter]]
  eta.gfl<-eta.est.list[[iter]]
}

b.main = as.matrix(beta.gfl[a.main.n,])
b.GE = as.matrix(beta.gfl[a.GE.n,])
eta.main = as.matrix(eta.gfl[d.main.n,],sparse=T)
eta.GE = as.matrix(eta.gfl[d.GE.n,],sparse=T)

rm(beta.est.list)
rm(eta.est.list)
rm(gamma.est.list)
gc()
###################
num.main <- get.gr.num(n, b.main, eta.main, q, merge.all=merge.all)
num.GE <- get.gr.num(n, b.GE, eta.GE, p, merge.all=merge.all)
gr.main.est<-num.main$gr.num
gr.GE.est<-num.GE$gr.num
residual <- log(sum((data.response-cova_diag%*%beta.gfl)^2/n))
BIC.var<-residual+log(n*pp)*log(n)*(gr.main.est*(q)+gr.GE.est*(p))/n
BIC.o<-residual+log(n)*(gr.main.est*q+gr.GE.est*p)/n
  
# return(list(beta.gfl=beta.gfl, b.main=b.main, b.GE=b.GE, eta.gfl=eta.gfl, eta.main=eta.main, eta.GE=eta.GE,
              # residual=residual,BIC.var=BIC.var,BIC.o=BIC.o,num.main=num.main,num.GE=num.GE))








size = n; beta.gfl = b.main; eta.gfl = eta.main; q = q; merge.all=F
# 计算每对组合的 eta 之和，如果足够小则说明属于同一类
diff.gfl<-apply(matrix(eta.gfl,nrow=q,ncol=size*(size-1)/2),2,sum)

capgfl.matrix<-matrix(0,nrow=size,ncol=size)
# 如果没有完全为 0 的 eta（单维），则 n 个样本就是 n 类
if(length(which(diff.gfl==0))==0){
  capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
  return(list(gr.num=size, capgfl.matrix2=capgfl.matrix2))
}
# 如果所有 eta 都是 0，则直接压缩为一类
if(length(which(diff.gfl==0))==size*(size-1)/2){
  capgfl.matrix[upper.tri(capgfl.matrix)]=1
  capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
  return(list(gr.num=1, capgfl.matrix2=capgfl.matrix2))
}else{
  # 找出哪一对组合 eta == 0
  sample.index.gfl<-sample.index[,which(diff.gfl==0)]
  # 如果只有一对样本可以归为同一类
  if(length(which(diff.gfl==0))==1){
    capgfl.matrix[sample.index.gfl[1],sample.index.gfl[2]]<-1
  }else{
    # 如果可以压缩的很多
    for(i in 1:length(which(diff.gfl==0))){
      capgfl.matrix[sample.index.gfl[1,i],sample.index.gfl[2,i]]<-1
    }
  }
  # 不会出现大于1的情况，所有压缩的组合都在右上三角阵中
  capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
  # 关键步骤！！注意借鉴这种计算类别个数的方法
  group.num.gf = nrow(unique(capgfl.matrix2))
  
  if(merge.all){
    cap <- capgfl.matrix2
    # 有多少个组，num_subgroup 就是多长的列表
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
    # 对每个样本
    for (i in 1:dim(cap)[1]) {
      # 对每个整合后的类别
      for (k in 1:length(non_inter_list)) {
        if(length(match(cap[i,],non_inter_list[[k]])) > 0){
          cap[i,non_inter_list[[k]]] <- 1
        }
      }
    }
    capgfl.matrix2 <- cap
    group.num.gf = nrow(unique(capgfl.matrix2))
  }
}
