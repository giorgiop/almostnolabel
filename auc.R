#Compute AUC
auc <- function(test, pred){
  
  N <- length(pred) 
  N_pos <- sum(test == 1) #Number of P
  dfauc <- data.frame(out = test, prob = pred)
  dfauc <- dfauc[order(-dfauc$prob),] #Order by descending score
  
  dfauc$above <- (1:N) - cumsum(dfauc$out)
  return (1- sum( dfauc$above * dfauc$out ) / (N_pos * (N-N_pos)))
}