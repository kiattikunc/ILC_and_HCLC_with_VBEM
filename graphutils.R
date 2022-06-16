###############################################################################
adjtodot <- function(ugmat, node_list){
  numNodes <- ncol(ugmat)
  varnames <- node_list
  edgelist <- c()

  for (i in 2:numNodes)
  {
    for (j in 1:(i-1))
    {
      if (ugmat[i,j]==3 & ugmat[j,i]==3) 
      {
        edgelist <- c(edgelist,paste(varnames[j],"->",varnames[i],"[arrowtail=none, arrowhead=none];"))
      }
      
      else if (ugmat[i,j]==2 & ugmat[j,i]==3)
      {edgelist <- c(edgelist,
                     paste(varnames[i],
                           "->",
                           varnames[j],"[arrowtail=none, arrowhead=normal];"))
      }
      else if (ugmat[i,j]==3 & ugmat[j,i]==2)
      {edgelist <- c(edgelist,
                     paste(varnames[j],
                           "->",
                           varnames[i],"[arrowtail=none, arrowhead=normal];"))
      }
      
      else if (ugmat[i,j]==1 & ugmat[j,i]==1) 
      {edgelist <- c(edgelist,
                     paste(varnames[i],
                           "->",
                           varnames[j],"[dir=both, arrowtail=odot, arrowhead=odot];"))
      }
      else if (ugmat[i,j]==2 & ugmat[j,i]==1) 
      {edgelist <- c(edgelist,
                     paste(varnames[i],
                           "->",
                           varnames[j]," [dir=both, arrowtail=odot, arrowhead=normal];"))
      }
      else if (ugmat[i,j]==1 & ugmat[j,i]==2) 
      {edgelist <- c(edgelist,
                     paste(varnames[j],
                           "->",
                           varnames[i]," [dir=both, arrowtail=odot, arrowhead=normal];"))
      }
     
      else if (ugmat[i,j]==2 & ugmat[j,i]==2) 
      {edgelist <- c(edgelist,
                     paste(varnames[i],
                           "->",
                           varnames[j]," [dir=both, color=blue, arrowtail=normal, arrowhead=normal ];"))
      }
     
       
    }
    
  }
  

    return( paste0("digraph g {",gsub(",","",toString(edgelist)),"}"))

}
