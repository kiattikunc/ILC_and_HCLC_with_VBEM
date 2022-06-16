########## remove the duplicated the locations of latent confounder  #############
remove_duplicated_pair <- function(input){ 
  rows <- c()
for (index_index in 1:nrow(input))
{
  for (index_index_2  in index_index:nrow(input))
  {
    
    if (input[index_index,1]== input[index_index_2,2] & input[index_index,2] ==input[index_index_2,1])
    {
      rows  <- rbind(rows,index_index_2)
    }
    
  }
}

if (length(rows)!=0)
{
  input <- input[-rows, ]
}
  output <- input
  return(output) 
}
########## End to remove the duplicated the locations of latent confounder  #############
