###a sample code

###install the R packages WQRADMM and quantreg
install.packages("Rcpp")
library(Rcpp)
install.packages("RcppArmadillo")
library(RcppArmadillo)
install.packages("quantreg")
library(quantreg)
install.packages("devtools")
library(devtools)
install_github("https://github.com/nanlin999/WQRADMM", force = TRUE)
#for macOS users, use the following code to install the WQRADMM package
#install_github("https://github.com/nanlin999/WQRADMMmac", force = TRUE)
library(WQRADMM)

###function for generating the correlation matrix (AR(1) or exchangeable)
gcov = function(p, rho, type){
  if(type == "exchangeable"){
    cov = matrix(rho, p, p)
    diag(cov) = rep(1, p)
  }
  else{
    cov = diag(p)
    for(i in 1:p){
      for(j in 1:p){
        if(i < j) cov[i,j] = rho^{j-i}
        else cov[i,j] = cov[j,i]
      }
    }
  }
  cov
}

###generate synthetic data (Student's t error with d = 0.75*p)
N = 10000
p = 100
n = 10
rep = rep(n, N)
nsum = sum(rep)
d = 0.75*p
rho_X = 0.5
rho_e = 0.5
tau = 0.75
set.seed(999)
X = matrix(rnorm(nsum*p), nsum, p)
cov_X = gcov(p, rho_X, "ar1")
X = X%*%chol(cov_X)
for(i in 1:d){
  X[,i] = pnorm(X[,i])
}
set.seed(999)
e = matrix(rt(N*n, 3), N, n)
cov_e = gcov(n, rho_e, "ar1")
e = as.vector(t(e%*%chol(cov_e)))
sigma = 0.5
e = sigma*e
beta0 = rnorm(1)
beta = rnorm(p)
Y = beta0+X%*%beta+apply(X[,1:d]*e/d, 1, sum)
beta_true = c(beta0, quantile(e/d, tau)+beta[1:d], beta[(d+1):p])

###calculate the WQR estimator by parallel WQR-ADMM
k = 10            ###number of partitions
paraWQR = paraWQRADMM(X, Y, k, rep, tau, TRUE, "WQR")
beta_paraWQR = paraWQR$Estimation_WQR
AE_paraWQR = sum(abs(beta_paraWQR-beta_true))
Time_paraWQR = paraWQR$Time_WQR
Time_paraWQRADMM = paraWQR$Time_total

###output the results
AE_paraWQR        ###estimation error
Time_paraWQR      ###running time for solving WQR problem
Time_paraWQRADMM  ###total computational time for WQR