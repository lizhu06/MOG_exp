##' Generate the data for simulation 1
##' This function generates the data for simulation 1: 
##' 	single level non-overlapping groups
##'
##' @title Generate data for simulation 1
##' @param seed random seed
##' @return a data list with following items
##' 	Y: observation vector of length n;
##'   X: n by p covariate matrix;
##'   U1: p by m1 matrix, indicate feature membership
##'   types: the type of each feature, different type can have different s^2
##'   true_beta: true coefficients
##' @export

geneSimu1 <- function(seed=1){
	set.seed(seed)
	N <- 125   # number of samples
	P <- 200   # number of features         
	m1 <- 10   # number of level-1 groups
	pk <- rep(20, m1)  # number of features in each level-1 group
	k_index <- rep(seq(1,m1),times=pk) # group index for each feature     
	X <- matrix(NA, N, P) # feature matrix
	Y <- rep(NA, N)
	true_group_ind <- c(rep(1, 5), rep(0, 5)) # indicate if group contribute to outcome
	true_beta_nonzero <- rnorm(50, 0, sqrt(5)) # draw all non-zero betas's
	true_beta <- c(c(true_beta_nonzero[1:20]), 
	             c(true_beta_nonzero[21:30], rep(0, 10)),  
	             c(true_beta_nonzero[31:40], rep(0, 10)), 
	             c(true_beta_nonzero[41:45], rep(0, 15)), 
	             c(true_beta_nonzero[46:50], rep(0, 15)), rep(0, 100))
	for(n in 1:N){
	  z <- rnorm(m1, 0, 1)
	  e <- rnorm(P, 0, 1)
	  for(k in 1:m1){
	    X[n, ((k-1)*pk[k]+1):(k*pk[k])] <- 
	    	(z[k]+e[((k-1)*pk[k]+1):(k*pk[k])])/sqrt(1+1) # within group correlation is 0.5
	  }
	}
	Y <- X%*%true_beta + rnorm(N)  # subject to change
	# create feature group (U1) matrix
	U1 <- matrix(0, P, m1)
	for(k in 1:m1){
		if(k==1){
			U1[1:pk[1], 1] <- 1
		}else{
			U1[(sum(pk[1:(k-1)])+1):(sum(pk[1:k])), k] <- 1
		}
	}
	types <- rep(1, P)
	return(list(Y=Y, X=X, U1=U1, types=types, true_beta=true_beta))
}
