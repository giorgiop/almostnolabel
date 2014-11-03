
#------------------------------------------------------------------------------
#This implements AMM^{min,max}
#
#The algorithm expects the data organized in data$bag and data$label (the proportions).
#The rest is the features vector.
#
#Value for init: "1"=constance 1s vector, "MM", "LMM", or a number= # random restart and pick the one with lowest loss
#
#
#Parameters gamma and L are relevant only for the initialisation step with LMM
#------------------------------------------------------------------------------


#Precision for the exit condition
ERR = 1e-3



alternating.mean.map <- function(data, L=NULL, minmax=FALSE, lambda=10, gamma=10,
                                 weight=NULL, init="LMM", maxcount=30){
  
  f <- function(w) {
    xw    <- X %*% w
    ai    <- log(exp(xw) + exp(-xw))
    ai    <- ifelse(is.finite(ai), ai, xw) # Handle numerical overflow for log(exp(..))
    lterm <- sum(ai)
    rterm <- w %*% yx
    loss  <- as.numeric(lterm - rterm + (0.5*lambda*(w %*% w)))
  }
  
  g <- function(w) {
    xw    <- as.vector(X %*% w)
    lterm <- colSums(X*tanh(xw))
    lterm - yx + lambda*w
  }
  
  #Eliminate row names from other ordering - this is an issue with sorting in R
  rownames(data) <- NULL
  
  K <- ncol(data)-2
  M <- nrow(data)
  #number of bags
  N <- length(unique(data$bag))
  bags <- 1:N
  
  #To build mapping with original bag numbers. Now select map.bag[j]
  map.bag <- sort(unique(data$bag))
  
  # extract proportions
  proportions <- foreach(bag = map.bag, .combine=rbind) %do% {
    id <- which(data$bag==bag)
    data$label[id[1]]
  }
  
  X <- as.matrix(data[,-(1:2)])
  
  outer_loop = FALSE
  best_loss = +Inf
  n_restart = 1
  
  #Find theta_0
  if (is.numeric(init)){
    outer_loop = TRUE
    n_restart <- nrow(init)
  }else if (init == "LMM"){
    if (is.null(L)){
      print("Matrix L must be specified for using LMM")
      quit(save="no")
    }
    
    theta <- laplacian.mean.map(data, L=L, lambda=lambda, gamma=gamma, weight=weight)
  }else if (init == "MM"){
    theta <- mean.map(data, lambda=lambda, weight=weight)
  }else if (init == "1"){ #Constant vector
    theta <- rep(1,K)
  }
  
  
  for (start in 1:n_restart){ 
    
    if (outer_loop == TRUE){
      theta <- init[start,]
    }
    
    
    loss <- +Inf
    
    theta_best <- theta #Only used for AMM^max
    
    count <- 1
    while (TRUE){
      
      #Find y_t
      prod <- X %*% theta
      prod <- prod / max(abs(prod)) #To solve an issue with sort() that works only up to 100k
      mat <- as.data.frame(cbind(data$label,data$bag,prod)) 
      mat[,1] <- -1 #By default, all labels are negative. Then go to modify the positive ones
      
      
      if (minmax == FALSE){ #default problem max max
        mat <- mat[order(mat[,3], decreasing=TRUE),, drop=FALSE]
      }else{ #Solve minmax problem instead
        mat <- mat[order(mat[,3], decreasing=FALSE),, drop=FALSE]  
      }
      
      for (b in bags){
        id <- which(mat[,2] == map.bag[b])
        #Take the first m_j * pi_j
        n.positive <- proportions[b] * length(id)
        bmat <- mat[id,]
        bmat[1:n.positive,1] <- +1 #Change in positive label only the first \pi_j * |B_j|
        mat[id,] <- bmat
      }
      
      
      #Re-use the original order
      mat <- mat[order(as.numeric(rownames(mat))),]
      
      #Find theta_t
      yx <- colSums(mat[,1] %*% X)
      
      w0 <- rep(0.001,K)
      #print("Time for learning\n")
      
      res <- optim(w0, fn=f, gr=g, method="L-BFGS-B")
      theta <- res$par
      
      #print(f(theta))
      
      #Termination condition - different for minmax since it does not converge in proper sense
      if (minmax == FALSE & abs(f(theta) - loss) < ERR)
        break
      else if (minmax == TRUE){
        
        #Keep track of the best solution
        if (f(theta) < loss)
          theta_best <- theta
        
        if (abs(f(theta) - loss) < ERR | count > maxcount) #After at least maxcount iterations or when it doesn't decrease much anymore
          break
      }
      
      
      count <- count + 1
      loss <- f(theta)
    }
    
    #AMM^max
    if (minmax == TRUE)
      theta <- theta_best
    
    #AMM_random - save the final model with the lowest loss over the n_restart
    if (outer_loop == TRUE && best_loss > loss){
      best_loss <- loss
      theta_random <- theta
    }
    
  } #random restart loop
  
  if (outer_loop == TRUE)
    theta <- theta_random
  
  #print(sprintf("%d iteration to AMM convergence",count))
  data.frame(theta=theta,steps=count)
}
