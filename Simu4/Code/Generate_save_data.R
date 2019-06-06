rm(list=ls())

setwd("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping")
code_folder <- "/ihome/ctseng/liz86/MOG/Code/MOG_revision1/"

ind <- matrix(NA, 3, 1+100+10)
ind[,1] <- c(-1, -1, 1) # level-0 features
ind[1, 2:101] <- seq(1, 300, by=3)
ind[2, 2:101] <- seq(3, 300, by=3)
ind[3, 2:101] <- sqrt(3)  # level-1 groups
ind[1, 102:111] <- seq(1, 300, by=30) + 300
ind[2, 102:111] <- seq(30, 300, by=30) + 300
ind[3, 102:111] <- sqrt(30) # level-2 groups
write.table(ind, sep=",", 
  quote=FALSE, row.names=FALSE, col.names=FALSE,
  file="Data/ind.csv")

U_seq<-c(0.2,0.5)
U_folder <- c(2, 5)
ns<-100

t0<-proc.time()
for(du in 1:length(U_seq)){
  for(s in 1:ns){
    seed <- 2016+100000*s
    source(paste0(code_folder, "R/GeneSimu3_v2_nonOverlapping.R"))

    data <- geneSimu3_nonOverlapping(U=U_seq[du], seed)
    X <- data$X
    Y <- data$Y
    true_beta <- data$true_beta

    ## split the data to training and testing
    seed <- 2016+100000*s + 1
    train_id <- sample(seq(1, nrow(X)), size = nrow(X)*0.8, replace=FALSE)
    X_train <- X[train_id, ]
    Y_train <- Y[train_id]
    X_test <- X[-train_id, ]
    Y_test <- Y[-train_id]

    inner_fold_id <- sample(rep(seq(1,10), each=nrow(X_train)/10))
    write.table(inner_fold_id, sep=",", file=paste0("Data/U", U_folder[du],
     "/innerFoldID_", s, ".csv"), quote=FALSE, row.names=FALSE, col.names=FALSE)

    write.table(X_train, sep=",", file=paste0("Data/U", U_folder[du],
     "/X_train_", s, ".csv"), quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(Y_train, sep=",", file=paste0("Data/U", U_folder[du], 
      "/Y_train_", s, ".csv"), quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(X_test, sep=",", file=paste0("Data/U", U_folder[du],
     "/X_test_", s, ".csv"), quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(Y_test, sep=",", file=paste0("Data/U", U_folder[du], 
      "/Y_test_", s, ".csv"), quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(true_beta, sep=",", file=paste0("Data/U", U_folder[du], 
     "/true_beta_", s, ".csv"), quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}
t1<-proc.time()-t0



