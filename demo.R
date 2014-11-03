library(doMC)
registerDoMC(cores=2)
library(ggplot2)
library(data.table)
library(digest)
library(pROC)

current_path = "/Users/giorgio/Google Drive/Laplacian Mean Map/public_github/"
source(file=paste(current_path,"auc.R", sep=""))
source(file=paste(current_path,"mean.map.R", sep=""))
source(file=paste(current_path,"laplacian.mean.map.R", sep=""))
source(file=paste(current_path,"alternating.mean.map.R", sep=""))
source(file=paste(current_path,"logistic.regression.R", sep=""))


#Hyperparameters
lambda <- 10
gamma <- .01
sigma <- 1

set.seed(123456)

#The dataset is provided in the form of (label, bag, features)
#with the label being still the binary labels.
#This is useful in case of testing the model
load(file=paste(current_path,"demo.heart.R", sep="")) #load heart.data
heart.data <- heart.data[sample(nrow(heart.data)),] #shuffle

test.fold <- 1:(floor(nrow(heart.data)/5))
testset <- heart.data[test.fold,]
trainset <- heart.data[-test.fold,]

N <- length(unique(trainset$bag)) #count the bags into the trainset

#Logistic regression - the labels are the binary ones, not proportions - Oracle
trainset <- trainset <- heart.data[-test.fold,]
w.lr <- logistic.regression(trainset,lambda)
test.X <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.lr))
test.auc <- auc((testset$label+1)/2, test.pred)
print(test.auc)


#Cast label in proportions
for (bag in unique(trainset$bag)) {
  id <- which(trainset$bag==bag)
  trainset$label[id] <- rep(mean((trainset$label[id]+1)/2), length(id))
}

#Mean Map
w.mm <- mean.map(trainset,lambda)
test.X <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.mm))
test.auc <- auc((testset$label+1)/2, test.pred)
print(test.auc)


#Laplacian Mean Map with similarity v^G
#Functions for building the Laplacian matrix are into laplacian.mean.map.R
L <- laplacian(similarity="G,s", trainset, N, sigma=1)
w.lmm <- laplacian.mean.map(trainset, lambda, gamma, L = L)
test.X <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.lmm))
test.auc <- auc((testset$label+1)/2, test.pred)
print(test.auc)


#Laplacian Mean Map with similarity v^{G,s}
L <- laplacian(similarity="G,s", trainset, N, sigma=sigma)
w.lmm <- laplacian.mean.map(trainset, lambda, gamma, L = L)
test.X <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.lmm))
test.auc <- auc((testset$label+1)/2, test.pred)
print(test.auc)


#Laplacian Mean Map with similarity v^NC
L <- laplacian(similarity="NC", trainset, N)
w.lmm <- laplacian.mean.map(trainset, lambda, gamma, L = L)
test.X <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.lmm))
test.auc <- auc((testset$label+1)/2, test.pred)
print(test.auc)


#Alternating Mean Map started with MM
w.amm <- alternating.mean.map(trainset, lambda=lambda, init="MM")
w.amm <- w.amm$theta #the algorithm returns a structure that contains also the number of step until termination
test.X <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.amm))
test.auc <- auc((testset$label+1)/2, test.pred)
print(test.auc)


#Alternating Mean Map started with LMM with similarity v^{G,s}
L <- laplacian(similarity="G,s", trainset, N, sigma=10)
w.amm <- alternating.mean.map(trainset, lambda=lambda, init="LMM", L=L, gamma=gamma)
w.amm <- w.amm$theta #the algorithm returns a structure that contains also the number of step until termination
test.X <- as.matrix(testset[,-c(1,2)])
test.pred <- 1/(1+exp(-2*test.X %*% w.amm))
test.auc <- auc((testset$label+1)/2, test.pred)
print(test.auc)
