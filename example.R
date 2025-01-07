source("funIHC.R")

library(R.matlab)
library(mclust)
library(fda)
library(cluster) 
library(stats)    
library(proxy) 

testing_data<- readMat('U1505.mat') #Data from Table 2, M=15 sigma = 0.05 
data<-testing_data[["Data"]] #The curves
id<-testing_data[["id"]] #The curves cluster assignment (ground truth)

x = 1:15 #Number of time points (M=15 or 200 in the simulations)
Data = t(data)

funihc = funIHC(x,Data,0.01,type=0) #Run funIHC, type refers to the choice of distance metric: curves (0), first derivative (1) or coefficients (2)
adjustedRandIndex(funihc$label,id) #Compute ARI

