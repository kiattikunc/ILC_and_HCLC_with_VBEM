construct_aux_latent_matrix <- function(amat,min_latentconfounder,partition_test){
  
  
latentconfounder_matrix_x <- matrix(0, nrow = ncol(amat), ncol = min_latentconfounder)
latentconfounder_matrix_y <- matrix(0, nrow = min_latentconfounder, ncol = ncol(amat)+min_latentconfounder)

latent_name <-c()
for (i in 1:min_latentconfounder)
{
  latent_name<-  c(latent_name, paste0("LV_",i))
  
}  

colnames(latentconfounder_matrix_x ) <- latent_name
rownames (latentconfounder_matrix_y) <- latent_name


amat <- cbind(amat,latentconfounder_matrix_x)
amat <- rbind(amat,latentconfounder_matrix_y)

  
  ######### Add latent confounder ###################
  
  condf_name <- c()
  
  for (partition_test_index in 1:length(partition_test))
  {
    for (child_index in 1:length(partition_test[[partition_test_index]]))
    {
      child <- which(colnames(amat)==partition_test[[partition_test_index]][child_index])

      
      
      amat[ncol(amat)+partition_test_index-min_latentconfounder,child]<-2
      amat[child,ncol(amat)+partition_test_index-min_latentconfounder]<-3
    }
    
  }
  
 
  ########## End adding the auxiliary Latent confounders    #############       
  

  return(amat)
  
}