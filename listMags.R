
listMags <- function(pag, nMags = 500, method=method, Z_i = NULL, opag = pag){
  cat("Entering listMags... \n")
  tag <- pag
  ccomp <- matrix(0,nrow(pag),ncol(pag)) #empty copy of pag
  p <- ncol(pag)
  cat("original pag = ", "\n")
  print(pag)
  
  ind <- which(pag==1)
  lind <- length(ind)
  if(lind==0) return(list(pag)) ### if there are no circles...
  
  ### CREATE THE TAIL AUGEMENTED GRAPH (TAG) ###
  
  for(i in 1:nrow(tag)) {
    for(j in 1:ncol(tag)){
      #double_circle <- 999
      tag[i,j] <- ifelse(tag[i,j]==1 && tag[j,i]!=1,3,tag[i,j]) # circles on o-> and o-- become tails 
      if(tag[i,j]==1 && tag[j,i]==1) ccomp[i,j]<-ccomp[j,i]<-1
    }
  }
  #cat("tag (tail augmented graph) = ", "\n")
  #print(tag)
  #cat("circle component of tag = ", "\n")
  #print(ccomp)
  
  ### END TAG ###
  
  # NEED TO MAKE SURE CCOMP IS A WGTMATRIX NOT JUST DEFAULT MATRIX
  ccomp_g <- as(ccomp,"graphNEL") # make the circle component of the pag into a graph
  ccomp_m <- wgtMatrix(ccomp_g) # get the wgtMatrix of the circle component graph (a pattern/CPDAG)
  cat("Entering allDags method... \n")
  #cat("print circle component: \n")
  #print(ccomp_m)
  listDags <- allDags.fast(ccomp_m,ccomp_m,NULL) # get all possible DAGs from that pattern
  cat("Finished allDags. \n")
  if(is.null(listDags)){
    cat("ERROR! Circle compenent of input PAG cannot be oriented into a DAG! \n")
    return(NULL)
  }
  
  if(FALSE){
    ### this stuff is new as of 12.14.2014 (below)###
    
    if(method=="local"){
      
      ##function which requires igraph package
      connected <- function(ccomp){
        if(!require("igraph")) stop("Package 'igraph' must be installed!")
        ccomp_g2 <- graph.adjacency(ccomp,mode="undirected") # turn ccomp into an igraph object
        return(is.connected(ccomp_g2)) # if the circle component is connected, return TRUE; otherwise, FALSE.
      }
      
      ## if the circle component is NOT connected...
      if(!connected(ccomp)){
        cat("@@@@@@@ not-connected circle component @@@@@@ \n")
        #opag is the original pag
        ccomp_opag <- matrix(0,nrow(opag),ncol(opag))
        for(i in 1:nrow(opag)){
          for(j in 1:ncol(opag)){
            if(opag[i,j]==1 && opag[j,i]==1) ccomp_opag[i,j] <- ccomp_opag[j,i] <- 1
          }
        }
        # now ccomp_opag is the circle component of original pag
        rem <- c()
        for(k in 1:nrow(listDags)){
          ccomp_opag_tmp <- ccomp_opag
          dag <- matrix(listDags[k,],p,p) # current dag in list
          for(i in 1:p){
            for(j in 1:p){
              if(dag[i,j]==1){
                ccomp_opag_tmp[Z_i[[i]],Z_i[[j]]] <- 1
                ccomp_opag_tmp[Z_i[[j]],Z_i[[i]]] <- 0
              }
            }
          }
          # now ccomp_opag_tmp has the dag as a subgraph
          # check if ccomp_opag_tmp is extendable to a dag; if not, then
          # add k to rem (list of dags to remove from allDags)		
          check <- pdag2dag(as(ccomp_opag_tmp,"graphNEL"))	
          if(!check$success) rem <- c(rem,k)
          if(!check$success) cat("@@@@@@@ dag thrown out! @@@@@@ \n")	
        } # for k in listDags
        if(!is.null(rem)) listDags <- listDags[-rem,]
        
      } # if !connected
      
    } # if method=="local"
    
    ### this stuff is new as of 12.14.2014 (above)###
  } # if FALSE temporary as of 11.27.2016
  
  mags <- list() # make a list to fill with MAGs
  for(k in 1:nrow(listDags)){
    mags[[k]] <- tag
    for(i in 1:p){
      for(j in 1:p){ # loop thru the list of DAGs
        if(matrix(listDags[k,],p,p)[i,j]==1) # where there are edges in the DAG...
        {
          mags[[k]][i,j] <-2
          mags[[k]][j,i] <-3 # put those edges in a copy of the original TAG
        }
      }
    }
  }
  
  ###### Now we have a list of MAGs, each of which is a TAG with the circle component oriented as one of the possible DAGs ######
  
  
  ### The following code works roughly like this: 
  ### the outermost loop (with index k) iterates through the mags that have only 
  ### invariant double-headed arrows (these are stored in a list called mags). 
  ### Then for n=1 it looks at all the possible combinations of length n of circle marks
  ### in the original pag. Then it loops through the matrix entries which correspond to 
  ### circles in the original pag. For each of these, test if the current mag under study 
  ### has a tail at the right matrix location. If it does, change it to an arrowhead. 
  ### This is a single mark change from the mag in mags. Then look at all the possible 
  ### combinations of circles, length n=2. Try changing both marks. That would be a two-mark 
  ### different mag. Then try n=3 three mark changes, then 4... Only save those mags that 
  ### haven't been saved before. ###
  
  bigmaglist <- vector(mode="list", nMags) # default value: the list can hold max 500 mags!!!
  for(k in 1:length(mags)){bigmaglist[[k]]<-mags[[k]]}
  i <- k+1 # counter to fill the list called bigmaglist
  mag_t <- tmp <- tmp_old <- matrix(0,nrow(pag),ncol(pag))
  for(k in 1:length(mags)){ # loop over elements in the list of mags
    mag_t <- tmp <- mags[[k]]
    for(n in 1:lind){
      #cat("Entering combn for circle marks... \n")
      comb <- combn(ind,m=n,FUN=NULL,simplify=FALSE)
      #cat("Finished combn. \n")
      ### BUG FIX
      ### if input to combn(x) is an integer (i.e., ind is of length 1), the it returns seq(1:x)
      ### we don't want that! so just let comb <- ind in that case
      if(lind==1) comb <- ind
      ###
      for(y in 1:length(comb)){
        tmp <- mag_t
        for(q in comb[[y]]){
          tmp_old <- tmp
          if(mag_t[q]==3){ # if the mark is tail
            tmp[q] <- 2 # change it to an arrowhead
          } # if mag_t[q]==3
        } # for q
        #n <- n+1		
        #cat("Checking for duplicates... \n")
        if(any(sapply(bigmaglist,identical,tmp))) next
        if(any(sapply(bigmaglist,identical,tmp_old))){
          if(i > nMags){
            print(i)
            cat("MORE MAGS THAN SPECIFIED BY PARAMETER 'nMags' !!!!! ", "\n")
            bigmaglist <- bigmaglist[!sapply(bigmaglist, is.null)]
            cat("Finishing listMags. \n")
            return(bigmaglist)
            #break
          }
          foundpath <- c()
          #cat("Checking conditions 1, 2, and 3 for transformation... \n")
          if(cond1(tmp_old,pag,q) && cond2(tmp_old,pag,q) && cond3(tmp_old,pag,q)) { # conditions for transformational equivalence....
            #cat("Passed all 3 conditions. \n")
            bigmaglist[[i]] <- tmp
            i <- i+1
            #cat("saving graph with k = ",k," and m = ", m, "and q = ", q, "\n")
          } # conditions are TRUE		
        } # if tmp_old is in bigmaglist
      } # for y	
    } # for n
  } # for k
  
  bigmaglist <- bigmaglist[!sapply(bigmaglist, is.null)]
  cat("Finishing listMags. \n")
  return(bigmaglist)
  
} # end function listMags


allDags.fast <- function(gm,a,tmp, verbose=FALSE)
{
  ## Purpose: Find all DAGs for a given PDAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix of initial PDAG; only 0-1 entries
  ##   i -> j iff gm(j,i)=1
  ## - a: copy of gm
  ## - tmp: NULL
  ## ----------------------------------------------------------------------
  ## Value:
  ## - one 0/1 adj.matrix per row
  ## Reversion to graph: as(matrix(res[i,],p,p),"graphNEL")
  ## Reversion to wgtMatrix (i->j iff a[j,i]=1): t(matrix(res[i,],p,p))
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  7 Apr 2008, 14:08
  if (sum(a) == 0) {
    if (verbose) {
      cat("Last Call - Final Graph: \n")
      print(gm)
      cat("#################### \n")
    }
    tmp2 <- rbind(tmp,c(t(gm)))
    if (all(!duplicated(tmp2))) tmp <- tmp2
  } else {
    sinks <- find.sink2(a)
    if (verbose) {
      cat("Main Call: ################## \n")
      print(gm)
      print(a)
      cat("Sinks: ",sinks,"\n")
    }
    for(x in sinks) {
      if (verbose) cat("Try removing", x," in a.\n")
      gm2 <- gm
      a2 <- a
      if (adj.check(a,x)) {
        inc.to.x <- a[, x] == 1 & a[x, ] == 1
        if (any(inc.to.x)) {
          real.inc.to.x <- as.numeric(row.names(a)[inc.to.x])
          real.x <- as.numeric(row.names(a)[x])
          gm2[real.x, real.inc.to.x] <- 1
          gm2[real.inc.to.x, real.x] <- 0
        }
        a2 <- a[-x,-x]
        if (verbose) {
          cat("Removed sink",as.numeric(row.names(a)[x]),
              "in g (", x,"in a).\n")
          cat("New graphs: \n")
          print(gm2)
          print(a)
        }
        tmp <- allDags.fast(gm2,a2,tmp, verbose)
      }
    }
  }
  tmp
}


find.sink2 <- function(gm) {
  ## Purpose: Find sink of an adj matrix; return numeric(0) if there is none;
  ## a sink may have incident undirected edges, but no directed ones
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix (gm_i_j is edge from j to i)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006;  speedup: Martin Maechler, Dec.2013
  ## New speedup: DMalinsky, Feb.2017
  
  uncon <- which(colSums(gm) == 0) # added 2.27.2017 to speed things up
  
  ## treat undirected edges
  gm[gm == t(gm) & gm == 1] <- 0
  ## treat directed edges
  setdiff(which(colSums(gm) == 0),uncon)
}
adj.check <- function(gm,x) {
  ## Purpose:  Return "TRUE", if:
  ## For every vertex y, adj to x, with (x,y) undirected, y is adjacent to
  ## all the other vertices which are adjacent to x
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: adjacency matrix of graph
  ## - x: node number (number)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006;
  ## several smart speedups: Martin Maechler, Dec.2013
  
  gm.1 <- (gm == 1)
  xr <- gm.1[x,]
  xc <- gm.1[,x]
  nx <- which(xr | xc)
  ## undirected neighbors of x
  un <- which(xr & xc)
  for(y in un) {
    adj.x <- setdiff(nx, y)
    adj.y <- setdiff(which(gm.1[y,] | gm.1[,y]), x)
    if(!all(adj.x %in% adj.y))
      return(FALSE)
  }
  TRUE
}



cond1 <- function(tmp_old,pag,q){
  ind3 <- which(pag==1)
  ind3.a <- which(pag==1,arr.ind=TRUE,useNames=FALSE)
  xy3 <- ind3.a[which(ind3==q),]
  a <- xy3[2]
  b <- xy3[1]
  if(is.path(a,b,tmp_old)) return(FALSE)
  else return(TRUE)
} # end function



is.path <- function(a, b, g, internal = FALSE){
  ind1 <- which(g==3, arr.ind=TRUE, useNames=FALSE)
  ind1 <- subset(ind1,ind1[,2]==a) # array of indices for tails out of A
  if(!internal){
    ind1 <- subset(ind1,!(ind1[,1]==b)) # minus A-->B
  }
  if(nrow(ind1)==0) return(FALSE) # addition: if there are no tails at A aside from the A-->B edge
  for(x in 1:nrow(ind1)){ # loop through tails out of A
    if(g[ind1[x,2],ind1[x,1]]==2){ # if there is an arrowhead at the other end of the x-th tail (call this C)
      if(ind1[x,1]==b){
        foundpath <- append(foundpath,TRUE)
        break
      }
      if(any(g[,ind1[x,1]]==3)){ # if there are any tails out of C, i.e., A-->C--*
        a_old <- a
        a2 <- ind1[x,1]
        if(a2==a_old) next
        foundpath <- append(foundpath,is.path(a2,b,g,internal=TRUE))
        if(any(foundpath)==TRUE) break
      }
    } # if there isn't an arrowhead at C - !(A-->C) - don't return anything
  } # for x in 1:nrow(ind1)
  if(any(foundpath)==TRUE) return(TRUE)
  else return(FALSE)
} # end function

cond2 <- function(tmp_old,pag,q){
  ind <- which(pag==1)
  ind.a <- which(pag==1,arr.ind=TRUE,useNames=FALSE)
  xy <- ind.a[which(ind==q),] # a two-element list which is the i,j coordinate of an element in ind
  # xy[1] is the i coordinate
  # xy[2] is the j coordinate
  # can do things like pag[xy[1],xy[2]]
  if (any(tmp_old[,xy[2]]==2)){ # there is at least one C*->A
    cd.ind <- which(tmp_old[,xy[2]]==2) + (nrow(tmp_old)*(xy[2]-1)) # list of positions of variables with arrowheads into A
    cd.ind2 <- which(tmp_old==2, arr.ind=TRUE, useNames=FALSE)
    cd.ind.a <- subset(cd.ind2,cd.ind2[,2]==xy[2])
    #cd.ind.a <- which(tmp_old[,xy[2]]==2, arr.ind=TRUE, useNames=FALSE) # same list as above but in <row,col> form
    for(f in cd.ind){
      cd.xy <- cd.ind.a[which(cd.ind==f),] # cd.xy[1] is i coordinate, cd.xy[2] is j
      if(tmp_old[cd.xy[2],cd.xy[1]]==3){ # if C-->A
        #cat("THERE IS A C-->A *****************", "\n")
        if(!(tmp_old[cd.xy[1],xy[1]]==2 && tmp_old[xy[1],cd.xy[1]]==3)){
          #cat("BUT NO C-->B !!!! *****************", "\n")
          return(FALSE)
        } #return(FALSE) # if not C-->B
        else next ##### addition 4/29
      } # if(tmp_old[cd.xy[2],cd.xy[1]]==3)
      if(tmp_old[cd.xy[2],cd.xy[1]]==2){ # if C<->A
        #cat("THERE IS A C<->A **********************", "\n")
        if(!(tmp_old[cd.xy[1],xy[1]]==2 && (tmp_old[xy[1],cd.xy[1]]==3 || tmp_old[xy[1],cd.xy[1]]==2))){
          #cat("BUT NO C<->B OR C-->B !!!! **********************", "\n")
          return(FALSE)
        } #return(FALSE) # if neither C-->B or C<->B
        else next ##### addition 4/29 
      } # if(tmp_old[cd.xy[2],cd.xy[1]]==3)
      else return(FALSE) # if neither of these if-statements are entered... (something is wrong!)
    }
    return(TRUE) # for(f in cd.ind) ##### also added 4/29
  } # if any(tmp_old[,j]==2)
  else return(TRUE)
} # end of function


cond3 <- function(tmp_old,pag,q){
  #need to specify b, c in col coordinates
  ind4 <- which(pag==1)
  ind4.a <- which(pag==1,arr.ind=TRUE,useNames=FALSE)
  xy4 <- ind4.a[which(ind4==q),]
  b <- xy4[2]
  c <- xy4[1]
  pag <- tmp_old ### !!! since is.discr.path is written with 'pag', replace that with the MAG under condisideration
  indA <- which((pag[b, ] == 2 & pag[, b] != 0) & (pag[c, ] == 3 & pag[, c] == 2))
  if(length(indA)==0) return(TRUE) #if indA is empty return false
  for(a in indA){
    founddpath <- c()
    pathexists <- c()
    pathexists <- append(pathexists, is.discr.path(path=c(c,b,a),pag=pag))
    if(any(pathexists)==TRUE) return(FALSE)
  }
  return(TRUE)
}

founddpath <- c()
is.discr.path <- function (path, pag)
{
  stopifnot((n <- length(path)) >= 3)
  if (n > nrow(pag)) return(FALSE)
  
  pag <- pag
  c <- path[1]
  b <- path[2]
  a <- path[3]
  first.pos <- path[n]
  del.pos <- path[n - 1]
  indD <- which(pag[first.pos, ] != 0 &
                  pag[, first.pos] == 2)
  indD <- setdiff(indD, del.pos)
  for (d in indD)  if(all(d != path)) {
    if (pag[c, d] == 0 && pag[d, c] == 0) {
      ### found discr path There is a discriminating path between d and c for b
      cat("There is a discriminating path between:",
          d, "and", c, "for", b, "\n")
      founddpath <- append(founddpath,TRUE)
      break
    }
    else {
      if (pag[first.pos, d] == 2 && pag[d, c] == 2 && pag[c, d] == 3) {
        founddpath <- append(founddpath, is.discr.path(path = c(path, d), pag = pag))
        if(any(founddpath)==TRUE) break
      } ## else : keep 'pag'
    }
  } ## for( d )
  
  if(any(founddpath)==TRUE) return(TRUE)
  else return(FALSE)
}## {discr.path}


