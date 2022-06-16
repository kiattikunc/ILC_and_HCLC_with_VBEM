bidirectededges_to_partition <- function(amat){
 
  library(igraph)
  library(dagitty)
  library(pcalg)
  
  
  mag<-  pcalg2dagitty(amat,  colnames(amat),type="mag")
  
  partition <- list()

  for( v in names( mag ) ){
    spouse <- dagitty::spouses( mag, v )
    if (length(spouse) >1)
    {
      for (index_spouse in spouse)
      {
        cat( v, "<->", index_spouse, "\n" )
        partition[[v]] <- append(partition[[v]] ,v)
      }
      
      comb_index_latentconfounder <- t(combn(length(spouse), 2))
      m_connected_count <- 0
      for ( index in 1:nrow(comb_index_latentconfounder))
      {
        m_connected <- dconnected(mag,  spouse[comb_index_latentconfounder[index,1]],  spouse[comb_index_latentconfounder[index,2]], c() ) 
        cat ("dconnected : ", spouse[comb_index_latentconfounder[index,1]], "and ", spouse[comb_index_latentconfounder[index,2]] ,"= ",m_connected, "\n")
        if (m_connected) m_connected_count = m_connected_count+1
      }
      if (m_connected_count ==nrow(comb_index_latentconfounder))
      {
        for ( index in 1:nrow(comb_index_latentconfounder))
        {
          partition[[v]] <- append(partition[[v]],spouse[comb_index_latentconfounder[index,1]])
          partition[[v]] <- append(partition[[v]],spouse[comb_index_latentconfounder[index,2]])
          partition[[v]] <- partition[[v]][!duplicated(partition[[v]]),drop=F]
        }
        
      }
      else 
      {
        partition[[v]] <-NULL
      }
    }
    if (length(spouse) ==1)
    {
      cat( v, "<->", spouse, "\n" )
      partition[[v]] <- append(partition[[v]] ,v)
     partition[[v]] <- append(partition[[v]] ,spouse)
    }
    
  }
  
  index <- length(partition)
  while(index>=1)
  {
    v1 <-1
    print(index)
    if(index==1)
    {
      while (v1 !=index && v1 <= length(partition) &&  index <= length(partition)){
        if(all(partition[[index]] %in% partition[[v1]]) )
        {
          cat( "Variable : ", names(partition[index]) , "is sub set of ", names(partition[v1])," =", all(partition[[index]] %in% partition[[v1]]), "\n" )  
          
          partition[index] <- NULL
          
        }
        v1<-v1+1
      }
    }
    else
    {
      while ( v1 <= length(partition)){
        if (v1 !=index)
        {
            if(all(partition[[index]] %in% partition[[v1]]) )
            {
              cat( "Variable : ", names(partition[index]) , "is sub set of ", names(partition[v1])," =", all(partition[[index]] %in% partition[[v1]]), "\n" )  
              
              partition[index] <- NULL
              
            }
            if( index > length(partition))
            {
              break()
            }
          }
        v1<-v1+1
      
      }
    }
    
    
    index <-  index-1
    
  }
  return(partition)
}

