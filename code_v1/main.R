################################################################################
# Codes for conducting simulation studies
################################################################################
rm(list = ls(all = TRUE))
ls()
###############################
setwd("D:/Ryanna/RUC_courses/Y1_2/Group/code")
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

# 原设定下得到的参数是不相同的
# 但是若参数真值的初始化给的不好也不行
# nC <- 4
# b00 = c(rep(-1,(nC*pp+nC)/2),rep(1,nC*pp+nC-floor((nC*pp+nC)/2)))
# b00 = rep(1,(nC*pp+nC))
# b00 = c(rep(3,(nC*pp+nC)/2),rep(-1,nC*pp+nC-floor((nC*pp+nC)/2)))
# res.mle <- fmrs.mle(y = whole.data$data_y, x =whole.data$data_x,delta = rep(0,length(whole.data$data_y)),
#                     nComp =nC, disFamily = "norm",
#                     # initCoeff = rnorm(nC*pp+nC),
#                     initCoeff = b00,
#                     initDispersion = rep(1,nC),
#                     initmixProp = rep(1/nC, nC),nIterNR = 200)
# coefficients(res.mle)
# res.lam <- fmrs.tunsel(y = whole.data$data_y, x =whole.data$data_x,delta = rep(0,length(whole.data$data_y)),
#                        nComp =nC, disFamily = "norm",
#                        initCoeff = c(coefficients(res.mle)),
#                        initDispersion = dispersion(res.mle),
#                        initmixProp = mixProp(res.mle),
#                        penFamily = "adplasso",nIterNR = 200)
# res.var <- fmrs.varsel(y = whole.data$data_y, x =whole.data$data_x, delta = rep(0,length(whole.data$data_y)),
#                        nComp = ncomp(res.mle), disFamily = "norm",
#                        initCoeff=c(coefficients(res.mle)),
#                        initDispersion = dispersion(res.mle),
#                        initmixProp = mixProp(res.mle),
#                        penFamily = "adplasso",
#                        lambPen = slot(res.lam, "lambPen"),nIterNR = 200)
# library(corrplot)
# corrplot(cor(whole.data$data_x))
whole.data
subgroup=c(2,4)
ridge = F
lambda.min=0.0001
pp<-p+q
n<-length(whole.data$data_y)
if(!ridge){
  fmr.bic<-c()
  fmr.res0<-c()
  res.mle0.list<-list()
  nC.vec<-subgroup
  for (ii in 1:length(nC.vec)) {
    nC<-nC.vec[ii]
    b00 = c(rep(-1,(nC*pp+nC)/2),rep(1,nC*pp+nC-floor((nC*pp+nC)/2)))
    res.mle <- fmrs.mle(y = whole.data$data_y, x =whole.data$data_x,delta = rep(0,length(whole.data$data_y)),
                        nComp =nC, disFamily = "norm",
                        # initCoeff = rnorm(nC*pp+nC),
                        initCoeff = b00,
                        initDispersion = rep(1,nC),
                        initmixProp = rep(1/nC, nC),nIterNR = 200)
    res.lam <- fmrs.tunsel(y = whole.data$data_y, x =whole.data$data_x,delta = rep(0,length(whole.data$data_y)),
                           nComp =nC, disFamily = "norm",
                           initCoeff = c(coefficients(res.mle)),
                           initDispersion = dispersion(res.mle),
                           initmixProp = mixProp(res.mle),
                           penFamily = "adplasso",nIterNR = 200)
    res.var <- fmrs.varsel(y = whole.data$data_y, x =whole.data$data_x, delta = rep(0,length(whole.data$data_y)),
                           nComp = ncomp(res.mle), disFamily = "norm",
                           initCoeff=c(coefficients(res.mle)),
                           initDispersion = dispersion(res.mle),
                           initmixProp = mixProp(res.mle),
                           penFamily = "adplasso",
                           lambPen = slot(res.lam, "lambPen"),nIterNR = 200)
    res.mle0<-res.var
    res.mle0.list[[ii]]<-res.mle0
    
    fmr.res=apply(abs(res.mle0@residuals),1,which.min)
    fmr.resj=c()
    for (jj in 1:n) {
      fmr.resj[jj]=res.mle0@residuals[jj,fmr.res[jj]]
    }
    fmr.res0[ii]<-sqrt(sum(fmr.resj^2))
    fmr.bic[ii]<-res.mle0@BIC
  }
  fmr.coe<-res.mle0@coefficients[-1,]
  fmr.w<-apply(res.mle0@weights,1,which.max)
  fmr.r<-apply(abs(res.mle0@residuals),1,which.min)
  
  # weight
  res.mle0 <- res.mle0.list[[which.max(fmr.bic)]]
  nC<-nC.vec[which.max(fmr.bic)]
  beta.int.fmr<-list()
  for (i in 1:n) {
    beta.int.fmr[[i]]<-as.matrix(as.numeric(fmr.coe[,fmr.w[i]]))
  }
  beta.init=do.call(rbind,beta.int.fmr)
}

# ------------ Apply the proposed method ------------ 
beta.init.list <- gen_int_beta(n, p, q, whole.data)
beta.init <- beta.init.list$beta.init
lambda <- genelambda.obo()
result <- tuning.lambda(lambda, whole.data, n, q, p, beta.init)
index.list <- evaluation.tuning(n,q,p, result$admmres, result$abic.n, result$admmres2, result$abic.n2, Beta0, gr.GE=gr.GE, result$bic.var)
index.list$err.s

