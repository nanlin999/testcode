###install the R package WQRADMM
install.packages("Rcpp")
library(Rcpp)
install.packages("RcppArmadillo")
library(RcppArmadillo)
install.packages("devtools")
library(devtools)
install_github("https://github.com/nanlin999/WQRADMM", force = TRUE)
library(WQRADMM)
#for macOS users, use the following code to install the WQRADMM package
#install_github("https://github.com/nanlin999/WQRADMMmac", force = TRUE)

###read data from .csv file
realdata = read.table("realdata.csv", sep=',', header = TRUE)  ###please put the .csv file directory into the double quotation marks
realdata = as.matrix(realdata)
NSUM = dim(realdata)[1]
p = 26
m = 24

###impute the missing values using the before-and-after data method  
PM2.5 = as.numeric(realdata[,6])
for(i in 1:NSUM){
  j = 0
  k = 0
  if(PM2.5[i]==10000){
    j = i-1
    k = i+1
    while(PM2.5[j]==10000){
      j = j-1
    }
    while(PM2.5[k]==10000){
      k = k+1
    }
    PM2.5[(j+1):(k-1)] = (PM2.5[j]+PM2.5[k])/2
  }
}

TEMP = as.numeric(realdata[,7])
for(i in 1:NSUM){
  j = 0
  k = 0
  if(TEMP[i]==10000){
    j = i-1
    k = i+1
    while(TEMP[j]==10000){
      j = j-1
    }
    while(TEMP[k]==10000){
      k = k+1
    }
    TEMP[(j+1):(k-1)] = (TEMP[j]+TEMP[k])/2
  }
}

PRES = as.numeric(realdata[,8])
for(i in 1:NSUM){
  j = 0
  k = 0
  if(PRES[i]==10000){
    j = i-1
    k = i+1
    while(PRES[j]==10000){
      j = j-1
    }
    while(PRES[k]==10000){
      k = k+1
    }
    PRES[(j+1):(k-1)] = (PRES[j]+PRES[k])/2
  }
}

DEWP = as.numeric(realdata[,9])
for(i in 1:NSUM){
  j = 0
  k = 0
  if(DEWP[i]==10000){
    j = i-1
    k = i+1
    while(DEWP[j]==10000){
      j = j-1
    }
    while(DEWP[k]==10000){
      k = k+1
    }
    DEWP[(j+1):(k-1)] = (DEWP[j]+DEWP[k])/2
  }
}

RAIN = as.numeric(realdata[,10])
for(i in 1:NSUM){
  j = 0
  k = 0
  if(RAIN[i]==10000){
    j = i-1
    k = i+1
    while(RAIN[j]==10000){
      j = j-1
    }
    while(RAIN[k]==10000){
      k = k+1
    }
    RAIN[(j+1):(k-1)] = (RAIN[j]+RAIN[k])/2
  }
}

WSPM = as.numeric(realdata[,11])
for(i in 1:NSUM){
  j = 0
  k = 0
  if(WSPM[i]==10000){
    j = i-1
    k = i+1
    while(WSPM[j]==10000){
      j = j-1
    }
    while(WSPM[k]==10000){
      k = k+1
    }
    WSPM[(j+1):(k-1)] = (WSPM[j]+WSPM[k])/2
  }
}

WD = as.numeric(realdata[,12])
for(i in 1:NSUM){
  j = 0
  k = 0
  if(WD[i]==10000){
    j = i-1
    k = i+1
    while(WD[j]==10000){
      j = j-1
    }
    while(WD[k]==10000){
      k = k+1
    }
    if(WD[j] == 0){
      if(WD[k]>180){
        WD[j] = 360
      }
    }
    if(WD[k] == 0){
      if(WD[j]>180){
        WD[k] = 360
      }
    }
    WD[(j+1):(k-1)] = (WD[j]+WD[k])/2
  }
}


###Standardize the response and continuous variables
y = PM2.5
y = (y-mean(y))/sd(y)
x = matrix(0, NSUM, p)
x[,1] = TEMP
x[,2] = PRES
x[,3] = DEWP
x[,4] = RAIN
x[,5] = WSPM
for(j in 1:5){
  x[,j]=(x[,j]-mean(x[,j]))/sd(x[,j])
}

###Construct the dummy variables for WD, season and monitoring site
for(i in 1:NSUM){
  if((WD[i]>0)&(WD[i]<90)){
    x[i,6] = 1
  }
  if(WD[i]==90){
    x[i,7] = 1
  }
  if((WD[i]>90)&(WD[i]<180)){
    x[i,8] = 1
  }
  if(WD[i]==180){
    x[i,9] = 1
  }
  if((WD[i]>180)&(WD[i]<270)){
    x[i,10] = 1
  }
  if(WD[i]==270){
    x[i,11] = 1
  }
  if((WD[i]>270)&(WD[i]<360)){
    x[i,12] = 1
  }
}

month = as.numeric(realdata[,3])
for(i in 1:NSUM){
  if((month[i]==6)|(month[i]==7)|(month[i]==8)){
    x[i,13] = 1
  }
  if((month[i]==9)|(month[i]==10)|(month[i]==11)){
    x[i,14] = 1
  }
  if((month[i]==12)|(month[i]==1)|(month[i]==2)){
    x[i,15] = 1
  }
}

for(i in 1:NSUM){
  if(realdata[i,13]=="Dongsi"){
    x[i,16] = 1
  }
  if(realdata[i,13]=="Guanyuan"){
    x[i,17] = 1
  }
  if(realdata[i,13]=="Wanshouxigong"){
    x[i,18] = 1
  }
  if(realdata[i,13]=="Tiantan"){
    x[i,19] = 1
  }
  if(realdata[i,13]=="Wanliu"){
    x[i,20] = 1
  }
  if(realdata[i,13]=="Nongzhanguan"){
    x[i,21] = 1
  }
  if(realdata[i,13]=="Gucheng"){
    x[i,22] = 1
  }
  if(realdata[i,13]=="Shunyi"){
    x[i,23] = 1
  }
  if(realdata[i,13]=="Changping"){
    x[i,24] = 1
  }
  if(realdata[i,13]=="Dingling"){
    x[i,25] = 1
  }
  if(realdata[i,13]=="Huairou"){
    x[i,26] = 1
  }
}

###randomly select the training set
xvector = as.vector(t(x))
xnew = matrix(xvector, NSUM/m, p*m, byrow = TRUE)
ynew = matrix(y, NSUM/m, m, byrow = TRUE)
Ntrain = 14400
NSUMtrain = Ntrain*m
rep = rep(m, Ntrain)
K = 12
train = matrix(0, Ntrain/K, K)
set.seed(666)
for(k in 1:K){
  train[,k] = sample((1+(k-1)*(NSUM/(m*K))):((NSUM/(m*K))*k), Ntrain/K)
}
train = as.vector(train)
train = sort(train)
xtrain = xnew[train,]
ytrain = ynew[train,]
X = matrix(as.vector(t(xtrain)), NSUMtrain, p, byrow = TRUE)
Y = as.vector(t(ytrain))

tau = 0.1
WQR = paraWQRADMM(X, Y, K, rep, tau, TRUE, "WQR", FALSE, "ar1")
T_WQR = WQR$Time_WQR                                    ###time cost for WQR estimation 
Est_WQR = WQR$Estimation_WQR[c(1,2,4,9,10,11,14,15)]    ###estimated coefficients of ``Intercept" and seven main variables      

###output the results 
T_WQR        
Est_WQR


