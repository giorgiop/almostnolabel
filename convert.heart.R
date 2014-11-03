convert_heart <- function(nbag=10,path){
  load <- read.table(paste(path,"heart.dat.txt",sep=""))
  
  nrow <- nrow(load)
  ncol <- ncol(load)
  
  #Prepare feature set
  data.matrix <- matrix(0,ncol=0,nrow=nrow(load))
  
  #Bagging through clustering
  fit <- kmeans(load[,-ncol],nbag,iter.max=50,nstart=50) 
  bags <- fit$cluster
  
  #Continuous variables
  for (i in 1:(ncol(load)-1)){
    if (!is.factor(load[,i])){
      data.matrix <- cbind(data.matrix, (load[,i] - mean(load[,i]))/sd(load[,i]))
    }
  }
  
  data <- data.frame(label=(2*load[,ncol]-3), bag=bags, x=data.matrix)
  
  #Constant feature
  data <- cbind(data, rep(1,nrow(data)))
  colnames(data)[ncol(data)] <- "constant"
  
  data
}