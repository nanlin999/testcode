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

###calculate the WQR estimator by WQR-ADMM
WQR = WQRADMM(X, Y, rep, tau, TRUE, "WQR")
beta_CQR = WQR$Estimation_CQR
AE_CQR = sum(abs(beta_CQR-beta_true))
Time_CQR = WQR$Time_CQR
beta_WQR = WQR$Estimation_WQR
AE_WQR = sum(abs(beta_WQR-beta_true))
Time_WQR = WQR$Time_WQR
Time_WQRADMM = WQR$Time_total

###calculate the WQR estimator by IP-based method 
Xinter = cbind(1, X)
Time_IP = Time_IPCQR = Time_IPW = Time_IPWQR = 0 
ptm = proc.time()
IP_CQR = rq.fit.fnb(Xinter, Y, tau)
Time_IPCQR = as.numeric((proc.time() - ptm)[3])
beta_IPCQR = IP_CQR$coefficients
AE_IPCQR = sum(abs(beta_IPCQR-beta_true))
IP_W = WS(X, Y, rep, tau, beta_IPCQR, TRUE)
Weight = drop(IP_W$Weight)
Time_IPW = IP_W$Time_W
Xw = Xinter*Weight*N
Yw = Y*Weight*N
ptm = proc.time()
IP_WQR = rq.fit.fnb(Xw, Yw, tau)
Time_IPWQR = as.numeric((proc.time() - ptm)[3])
beta_IPWQR = IP_WQR$coefficients
AE_IPWQR = sum(abs(beta_IPWQR-beta_true))
Time_IP = Time_IPCQR+Time_IPW+Time_IPWQR

###output the results

###comparison of WQR-ADMM and IP for WQR problem 
AE_WQR         ###estimation error of WQR-ADMM for WQR   
AE_IPWQR       ###estimation error of IP for WQR
Time_WQRADMM   ###computational time of WQR-ADMM for WQR
Time_IP        ###computational time of IP for WQR
beta_WQR[1]    ###estimate of beta0 obtained by WQR-ADMM
beta_IPWQR[1]  ###estimate of beta0 obtained by IP
beta0          ###true value of beta0

###comparison of WQR-ADMM and IP for CQR problem 
AE_CQR         ###estimation error of WQR-ADMM for CQR   
AE_IPCQR       ###estimation error of IP for CQR
Time_CQR       ###computational time of WQR-ADMM for CQR
Time_IPCQR     ###computational time of IP for CQR
beta_CQR[1]    ###estimate of beta0 obtained by WQR-ADMM
beta_IPCQR[1]  ###estimate of beta0 obtained by IP
beta0          ###true value of beta0

###comparison of WQR-ADMM and IP for model (1)
beta0 = 100
Y = beta0+X%*%beta+apply(X[,1:d]*e, 1, sum)
beta_true = c(beta0, quantile(e, tau)+beta[1:d], beta[(d+1):p])
###rerun the codes for calculating WQR estimator based on WQR-ADMM and IP, and we can obtain the results for model (1)



