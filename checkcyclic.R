checkcyclic <- function(amat){
  
####### check_cyclic #######

N<- ncol(amat)

checked_dag <- matrix(0, N, N)
checked_dag[amat==2 & t(amat)==3 ] <- 1
check_cyclic <- pcalg::isValidGraph(checked_dag, type = "dag", verbose = TRUE)

return(check_cyclic)
}
