construct_max_latent_matrix <- function(amat){
  
 
  
  min_latentconfounder <-   length(which(amat ==2 & t(amat)==2))/2
  
  
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
currentbidirectpair <- which(amat ==2 & t(amat)==2, arr.ind=TRUE)
currentbidirectpair <- remove_duplicated_pair(currentbidirectpair)
if (is.null(nrow(currentbidirectpair))) 
{
  currentbidirectpair <- matrix(currentbidirectpair,1,2)
}
  for (partition_test_index in 1:nrow(currentbidirectpair))
  {
   
      
      amat[ncol(amat)+partition_test_index-min_latentconfounder,currentbidirectpair[partition_test_index,1]]<-2
      amat[currentbidirectpair[partition_test_index,1],ncol(amat)+partition_test_index-min_latentconfounder]<-3
    
      amat[ncol(amat)+partition_test_index-min_latentconfounder,currentbidirectpair[partition_test_index,2]]<-2
      amat[currentbidirectpair[partition_test_index,2],ncol(amat)+partition_test_index-min_latentconfounder]<-3
      
    
  }
  
 
  ########## End adding the auxiliary Latent confounders    #############       
  

  return(amat)
  
}