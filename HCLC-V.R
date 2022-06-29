############################################################################################################################################
#v9 HC the direction of the structure and greedy with states
#9.1 compare all pairs of o-o

options(java.parameters = "-Xmx4096m")
library(pcalg)
library(Rgraphviz)
library(tictoc)
library(igraph)
library(bnlearn)
library(rcausal)
library(rJava)
library(RWeka)
library(DOT)
source("graphutils.R")
source("bidirectededges_to_partition.R")
source("remove_duplicated_pair.R")
source("construct_aux_latent_matrix.R")
source("checkcyclic.R")
############## For preliminary setup #########################
#require(devtools)
#install.packages("stringr")
#library(devtools)
#install_github("bd2kccd/r-causal", INSTALL_opts=c("--no-multiarch"),force = TRUE)
#remotes::install_github("bd2kccd/r-causal")
#install.packages("rJava",type='source',"https://cran.r-project.org/src/contrib/Archive/rJava/rJava_0.9-12.tar.gz") 

###############################################################


############## START #########################
tic.clearlog()
tic("scorePAG time")
start_time <- Sys.time()

######  hyper parameters

######  use PROPERTY dataset #### 
algo <- c('GFCI')
set <- c('PROPERTY')
latent1 <- c("propertyPurchaseValue")
n = 1000  # sample size
max_latentconfounder <- 4
alpha_stat <- 0.05
improve <- TRUE
bestELBO<- -Inf
filepath <- paste0("Input/",set,"/",set,"_",latent1,"_",n)
.jinit(classpath="C:\\Users\\kiattikun\\Documents\\GitHub\\ILC_HCLC_with_VBEM\\VBEM")
training_set <-paste0(filepath,".csv")
foundpath <- c()
time_max<- 20000


################function to convert tetrad format to Matrix format ###############
TetradGetAdjmat <- function(tetradrunner) {
  p <- length(tetradrunner$nodes)
  adjmat <- matrix(0, p, p)
  for (e in edges) {
    edgevec <- unlist(strsplit(e, " "))
    i <- match(edgevec[1], tetradrunner$nodes)
    j <- match(edgevec[3], tetradrunner$nodes)
    if (edgevec[2] == '---') {
      edge <- c(edgevec[1], edgevec[3])
      adjmat[i, j] <- 3
      adjmat[j, i] <- 3
    } else if (edgevec[2] == "-->") {
      adjmat[i, j] <- 2
      adjmat[j, i] <- 3
    }
    else if (edgevec[2] == "<--") {
      adjmat[i, j] <- 3
      adjmat[j, i] <- 2
    }
    else if (edgevec[2] == "o-o") {
      edge <- c(edgevec[1], edgevec[3])
      adjmat[i, j] <- 1
      adjmat[j, i] <- 1
    }
    else if (edgevec[2] == "o->") {
      adjmat[i, j] <- 2
      adjmat[j, i] <- 1
    }
    else if (edgevec[2] == "<-o") {
      adjmat[i, j] <- 1
      adjmat[j, i] <- 2
    }
    else if (edgevec[2] == "<->") {
      edge <- c(edgevec[1], edgevec[3])
      adjmat[i, j] <- 2
      adjmat[j, i] <- 2
    }
    
  }
  return(adjmat)
}

##############################################################

################ Generate synthetic data #####################
#bn_mod <- bn.fit(true_dag, truedag_dataset, method = "bayes")

#set.seed(5*n)
#new_sampling <- rbn(bn_mod, n)
#new_sampling <- new_sampling[ , !(names( new_sampling) %in% latent1)]
#new_sampling <- new_sampling[,order(colnames(new_sampling))]
#write.arff(new_sampling, file = paste0(filepath,".arff") )
#write.csv(new_sampling,paste0(filepath,".csv"),row.names = FALSE)
################ End generating synthetic data #####################


################### load observational dataset ###################
#dataset <- new_sampling
#load observational dataset
dataset<- read.csv(training_set,header = TRUE,na.strings=c(""),check.names = FALSE)

colnames(dataset) <- gsub("[%]", "", colnames(dataset))
colnames(dataset) <- gsub("[+]", "", colnames(dataset))
dataset <- dataset[,order(colnames(dataset))]
dataset_pcalg <- dataset[,order(colnames(dataset))]

cols <- colnames(dataset_pcalg)
for(j in cols){
  dataset_pcalg[[j]]<- as.numeric(dataset_pcalg[[j]])
  dataset_pcalg[[j]]<- dataset_pcalg[[j]]-1
}

suffStat <- list(dm = dataset_pcalg, nlev = sapply(dataset , nlevels),adaptDF = FALSE)

varNames_pcalg <- colnames(dataset_pcalg)

################### identify max of state  ###################
nlev = sapply(dataset , nlevels)
max_state <- max(nlev)

################### run structure learning algorithms ###############

if (algo =='FCI')
{
  
  ################### run FCI ###################
  
  ####### pcalg lib  ###
  varNames <- colnames(dataset)
  fci.est <- fci(suffStat, indepTest = disCItest,skel.method =c("stable"), label=c(varNames_pcalg), alpha=alpha_stat)
  #fci.est <- fci(suffStat, indepTest = gaussCItest, label=c(varNames), alpha=alpha_stat,conservative = FALSE,maj.rule = FALSE  )
  fci.est <- fci.est@amat
  varNames <- colnames(fci.est)
  x<- adjtodot(fci.est, varNames)
  dot(x)
  ############################################
} 
if (algo =='GFCI')
{
  ####### tetrad  ###
  #tetradrunner.getAlgorithmDescription(algoId = 'fci')
  #tetradrunner.listAlgorithms()
  tetradrunner_fci <- tetradrunner(algoId = 'gfci',
                                   df = dataset,
                                   dataType = 'discrete',
                                   alpha=alpha_stat,
                                   faithfulnessAssumed = TRUE,
                                   maxDegree = -1, #argv$maxDegree,
                                   verbose = FALSE)
  edges <- tetradrunner_fci$edges
  fci.est <-  TetradGetAdjmat(tetradrunner_fci)
  
  ################ Plot DoT #######################
  graph_dot <- tetradrunner.tetradGraphToDot(tetradrunner_fci$graph)
  dot(graph_dot)
}
if (algo =='FCI_tetrad')
{
  
  tetradrunner <- tetradrunner(algoId = 'fci',
                               df = dataset,
                               dataType = 'discrete',
                               alpha=alpha_stat,
                               faithfulnessAssumed = TRUE,
                               maxDegree = -1, #argv$maxDegree,
                               verbose = FALSE)
  edges <- tetradrunner$edges
  fci.est <-  TetradGetAdjmat(tetradrunner)
  
  ################ Plot DoT #######################
  graph_dot <- tetradrunner.tetradGraphToDot(tetradrunner$graph)
  dot(graph_dot)
}


colnames(fci.est) <- colnames(dataset)
rownames(fci.est) <- colnames(dataset)

####################################################

################### Incrementally search for MAG  ###################

###### Add the pair #####
fci.est_temp <- fci.est
# latentconfounder for  <-> only
bidirect_pair <- which((fci.est ==2 & t(fci.est)==2)  , arr.ind=TRUE)
if (nrow(bidirect_pair) >0) 
{
  bidirect_pair <- remove_duplicated_pair(bidirect_pair) 
}

# latentconfounder for o-> 
variantdirect_pair <- which((fci.est ==1 & t(fci.est)==2) , arr.ind=TRUE)
if (nrow(variantdirect_pair) >0) 
{
  variantdirect_pair <- remove_duplicated_pair(variantdirect_pair)
}

if (is.null(nrow(variantdirect_pair))) 
{
  variantdirect_pair <- matrix(variantdirect_pair,1,2)
}

# latentconfounder for o-o and their score
variantvariant_pair_temp <- c()
variantvariant_pair <- which((fci.est ==1 & t(fci.est)==1) , arr.ind=TRUE)


variantvariant_pair_temp <- cbind(variantvariant_pair,matrix(-Inf, nrow = nrow(variantvariant_pair), ncol = 1))
variantvariant_pair_temp <- cbind(variantvariant_pair_temp,matrix(0, nrow = nrow(variantvariant_pair), ncol = 1))

if (nrow(variantvariant_pair) >0) 
{
  variantvariant_pair <- remove_duplicated_pair(variantvariant_pair)
}

if (is.null(nrow(variantvariant_pair))) 
{
  variantvariant_pair <- matrix(variantvariant_pair_temp,1,2)
}

# all possilbe latentconfounder for o-o and o-> 
variant_pair <- rbind(variantdirect_pair, variantvariant_pair)

#### scan the <-> for the min latentconfounder
current_bidirectedge <-  length(which(fci.est ==2 & t(fci.est)==2))/2
learnt_bidirectedge <- current_bidirectedge
if (current_bidirectedge > 0) {
  partition_test <-bidirectededges_to_partition(fci.est)
  current_latentconfounder <- length(partition_test)
  min_latentconfounder <- current_latentconfounder
  score_improve = TRUE
} else if (current_bidirectedge ==0 && (nrow(variantdirect_pair)>0 || nrow(variantvariant_pair)>0)) {
  current_bidirectedge<-1
  current_latentconfounder<-1
  min_latentconfounder <- 1
  score_improve = TRUE
} else
{
  print("The learnt graph has no possible latent confounders")
  score_improve = FALSE
}


best_ELBO <- -Inf
best_local_ELBO <- -Inf
best_structure_initial <- fci.est
best_confounders <-c()
rest_confounder <- c()
list_ELBO <-c()
list_confounder<-c()
###### End adding the pair #####
while(score_improve && current_bidirectedge <=max_latentconfounder)
{ 
  print(paste0("Increased the number of bidirected edges =", current_bidirectedge) )
  list_ELBO <- c(list_ELBO,paste0("bidirected= ",current_bidirectedge ))
  list_confounder <- c(list_confounder,paste0("bidirected= ",current_bidirectedge ))
  
  
  ####### use <-> to consider the min latentconfounder and orientate all directions
  if (length(which(fci.est ==2 & t(fci.est)==2))/2 >0  & current_bidirectedge<=length(which(fci.est ==2 & t(fci.est)==2))/2)
  {
    if (is.null(nrow(bidirect_pair)))
    {
      bidirect_pair <- matrix(bidirect_pair,1,2)
    }
    
    for (index_bidirect_pair in 1:nrow(bidirect_pair))
    {
      fci.est_temp[bidirect_pair[index_bidirect_pair,1],bidirect_pair[index_bidirect_pair,2]] <- 2  
      fci.est_temp[bidirect_pair[index_bidirect_pair,2],bidirect_pair[index_bidirect_pair,1]] <- 2  
      best_confounders <- c(best_confounders,index_bidirect_pair)
      
    }
    
    ####change o-> to -> ####
    fci.est_temp[ fci.est_temp ==1 & t(fci.est_temp)==2 ] <- 3
    
    #### remove change o-o  ####
    fci.est_temp[ fci.est_temp ==1 & t(fci.est_temp)==1 ] <- 0
    
    partition_test <-bidirectededges_to_partition(fci.est_temp)
    current_latentconfounder <- length(partition_test)
    print(paste0("Update the number of latentconfounder to ",current_latentconfounder) )
    
    initial_amat <- construct_aux_latent_matrix(fci.est_temp,current_latentconfounder, partition_test)
    
    ###remove <->  between children to be added by A<-LV->B 
    currentbidirectpair <- which(fci.est_temp ==2 & t(fci.est_temp)==2, arr.ind=TRUE)
    currentbidirectpair <- remove_duplicated_pair(currentbidirectpair)
    if (is.null(nrow(currentbidirectpair))) 
    {
      currentbidirectpair <- matrix(currentbidirectpair,1,2)
    }
    for (index_row in 1:nrow(currentbidirectpair))
    {
      initial_amat[currentbidirectpair[index_row,2],currentbidirectpair[index_row,1]]<-0
      initial_amat[currentbidirectpair[index_row,1],currentbidirectpair[index_row,2]]<-0
      
    }
    # consider o-o 
    if (!is.null(nrow(variantvariant_pair_temp))) 
    {
      while (nrow(variantvariant_pair_temp) !=0)
      {  
        tabu_list <- c()
        for(index_variantvariant_pair in 1:nrow(variantvariant_pair_temp))
        {
          if(index_variantvariant_pair==1) 
          {
            best_local_ELBO <- -Inf
          }
          #try A->B
          fci.est_temp_AB <- initial_amat
          fci.est_temp_AB[variantvariant_pair_temp[index_variantvariant_pair,1] ,variantvariant_pair_temp[index_variantvariant_pair,2]]<-2
          fci.est_temp_AB[variantvariant_pair_temp[index_variantvariant_pair,2] ,variantvariant_pair_temp[index_variantvariant_pair,1]]<-3
          
          computeELBOinitial <- .jnew("experiments.real.computeELBOinitial")
          if(checkcyclic(fci.est_temp_AB))
          {
            input <- matrix(as.integer(fci.est_temp_AB),ncol=ncol(fci.est_temp_AB))
            
            ELBO_AB <- .jcall(computeELBOinitial, returnSig="D", method="excute",.jarray(input, dispatch = T), as.integer(current_latentconfounder),as.integer(2),paste0(filepath,".arff"))
            tabu_list <- rbind(tabu_list,c(variantvariant_pair_temp[index_variantvariant_pair,1] ,variantvariant_pair_temp[index_variantvariant_pair,2],ELBO_AB))
            
          }
          else
          {
            print("cyclic")
            ELBO_AB <- -Inf
          }
          variantvariant_pair_temp[index_variantvariant_pair,3]=ELBO_AB
          variantvariant_pair_temp[index_variantvariant_pair,4]=current_latentconfounder
         
        }
        
        
        
        # perform the orientation 
        ######### perform the orientation 
        
        if (max(variantvariant_pair_temp[,3]) !=-Inf)
        {
          
          
          A = variantvariant_pair_temp[ which.max( variantvariant_pair_temp[,3] ),1]
          B = variantvariant_pair_temp[ which.max( variantvariant_pair_temp[,3] ),2]
          
          initial_amat[A ,B]<-2
          initial_amat[B ,A]<-3
          fci.est_temp[A ,B]<-2
          fci.est_temp[B ,A]<-3
          best_structure_initial  <- fci.est_temp
          print(paste0("Best ELBO score: ",max(variantvariant_pair_temp[,3]), " when ", colnames(fci.est_temp_AB)[A],"->",colnames(fci.est_temp_AB)[B]) )
          best_local_ELBO <- max(variantvariant_pair_temp[,3])
          current_ELBO <- best_local_ELBO
          
          list_ELBO <- c(list_ELBO,current_ELBO)
          list_confounder <- c(list_confounder,variantvariant_pair_temp[which.max(variantvariant_pair_temp[,3]),4])
          write.csv(list_ELBO, paste0("output/HCLC_",algo,"_",set,"_listELBO_",latent1,"_",n,".csv"), row.names=F)
          write.csv(list_confounder, paste0("output/HCLC_",algo,"_",set,"_list_confounder_",latent1,"_",n,".csv"), row.names=F)
          
          variantvariant_pair_temp <- variantvariant_pair_temp[-c(which(variantvariant_pair_temp[,1]==B & variantvariant_pair_temp[,2]==A )),]
        
          if (is.null(nrow(variantvariant_pair_temp))) 
          {
            variantvariant_pair_temp <- c()
            break
          }
          else
          {
            variantvariant_pair_temp <- variantvariant_pair_temp[-c(which.max(variantvariant_pair_temp[,3])),]
            
          }
          
        }
        # reverse id
        else
        {
          
          A = variantvariant_pair_temp[ which.max( variantvariant_pair_temp[,3] ),1]
          B = variantvariant_pair_temp[ which.max( variantvariant_pair_temp[,3] ),2]
          
          initial_amat[A ,B]<-3
          initial_amat[B ,A]<-2
          fci.est_temp[A ,B]<-3
          fci.est_temp[B ,A]<-2
          best_structure_initial  <- fci.est_temp
        }
        
        
        
        
        if (!is.null(nrow(variantvariant_pair_temp))) 
        {
          if (nrow(variantvariant_pair_temp)>0) 
          { 
            
            for ( row_index in 1:nrow(variantvariant_pair_temp))
            {
              
              if (variantvariant_pair_temp[row_index,1]==B && variantvariant_pair_temp[row_index,2]==A )
              {
                variantvariant_pair_temp <- variantvariant_pair_temp[-c(row_index),]
                break
              }
              
            }
            
          }
        }
        
        ######### Ending perform the orientation 
        
        
      }
    }
    if (nrow(variantvariant_pair) ==0)
    {
      computeELBOinitial <- .jnew("experiments.real.computeELBOinitial")
      input <- matrix(as.integer(initial_amat),ncol=ncol(initial_amat))
      current_ELBO  <- .jcall(computeELBOinitial, returnSig="D", method="excute",.jarray(input, dispatch = T), as.integer(current_latentconfounder),as.integer(2),paste0(filepath,".arff"))
      fci.est_temp<- initial_amat
      if (current_ELBO>best_ELBO)
      {
        best_ELBO <- current_ELBO
        best_structure <- initial_amat
        write.csv(best_ELBO, paste0("output/HCLC_",algo,"_",set,"_bestELBO_",latent1,"_",n,".csv"), row.names=F)
        write.csv(best_structure, paste0("output/HCLC_",algo,"_",set,"_bestdag_",latent1,"_",n,".csv"), row.names=F)
        
      }
      
    }
    
    
    #revert back to learnt <-> for the next 
    if(nrow(bidirect_pair) !=0)
    {
      for (index_row in 1:nrow(bidirect_pair))
      {
        fci.est_temp[bidirect_pair[index_row,1],bidirect_pair[index_row,2]] <- 2  
        fci.est_temp[bidirect_pair[index_row,2],bidirect_pair[index_row,1]] <- 2
      }
    }
    
    #revert back to initial dimention for the next 
    fci.est_temp <- fci.est_temp[-c(ncol(dataset)+1:ncol(fci.est_temp)),-c(ncol(dataset)+1:ncol(fci.est_temp))] 
    
    rest_confounder <- best_confounders
    best_local_confounder <- c()
    best_structure_initial  <- fci.est_temp
    best_structure <- initial_amat
    best_ELBO <- best_local_ELBO
    
  }
  ####### if no more <->, incrementally add a latentconfounder from o-> or o-o 
  else if (length(which(fci.est ==2 & t(fci.est)==2))/2 < current_bidirectedge)
  {
    
    # consider o-o 
    ####change o-> to -> ####
    fci.est_temp[ fci.est_temp ==1 & t(fci.est_temp)==2 ] <- 3
    
    
    if(nrow(bidirect_pair)==0)
    {
      initial_amat <-fci.est_temp
      if (!is.null(nrow(variantvariant_pair_temp))) 
      {
        while (nrow(variantvariant_pair_temp) !=0)
        {  
          tabu_list <- c()
          for(index_variantvariant_pair in 1:nrow(variantvariant_pair_temp))
          {
            if(index_variantvariant_pair==1) 
            {
              best_local_ELBO <- -Inf
            }
            #try A->B
            fci.est_temp_AB <- initial_amat
            fci.est_temp_AB[variantvariant_pair_temp[index_variantvariant_pair,1] ,variantvariant_pair_temp[index_variantvariant_pair,2]]<-2
            fci.est_temp_AB[variantvariant_pair_temp[index_variantvariant_pair,2] ,variantvariant_pair_temp[index_variantvariant_pair,1]]<-3
            
            computeELBOinitial <- .jnew("experiments.real.computeELBOinitial")
            if(checkcyclic(fci.est_temp_AB))
            {
              input <- matrix(as.integer(fci.est_temp_AB),ncol=ncol(fci.est_temp_AB))
              
              ELBO_AB <- .jcall(computeELBOinitial, returnSig="D", method="excute",.jarray(input, dispatch = T), as.integer(current_latentconfounder),as.integer(2),paste0(filepath,".arff"))
              tabu_list <- rbind(tabu_list,c(variantvariant_pair_temp[index_variantvariant_pair,1] ,variantvariant_pair_temp[index_variantvariant_pair,2],ELBO_AB))
              
            }
            else
            {
              ELBO_AB <- -Inf
            }
            variantvariant_pair_temp[index_variantvariant_pair,3]=ELBO_AB
            variantvariant_pair_temp[index_variantvariant_pair,4]=current_latentconfounder
          }
          
          
          
          # perform the orientation 
          ######### perform the orientation 
          
          if (max(variantvariant_pair_temp[,3]) !=-Inf)
          {
            
            
            A = variantvariant_pair_temp[ which.max( variantvariant_pair_temp[,3] ),1]
            B = variantvariant_pair_temp[ which.max( variantvariant_pair_temp[,3] ),2]
            
            initial_amat[A ,B]<-2
            initial_amat[B ,A]<-3
            fci.est_temp[A ,B]<-2
            fci.est_temp[B ,A]<-3
            best_structure_initial  <- fci.est_temp
            print(paste0("Best ELBO score: ",max(variantvariant_pair_temp[,3]), " when ", colnames(fci.est_temp_AB)[A],"->",colnames(fci.est_temp_AB)[B]) )
            best_local_ELBO <- max(variantvariant_pair_temp[,3])
            current_ELBO <- best_local_ELBO
            
            list_ELBO <- c(list_ELBO,current_ELBO)
            list_confounder <- c(list_confounder,variantvariant_pair_temp[which.max(variantvariant_pair_temp[,3]),4])
            
            write.csv(list_ELBO, paste0("output/HCLC_",algo,"_",set,"_listELBO_",latent1,"_",n,".csv"), row.names=F)
            
            write.csv(list_confounder, paste0("output/HCLC_",algo,"_",set,"_list_confounder_",latent1,"_",n,".csv"), row.names=F)
            
            
            
            variantvariant_pair_temp <- variantvariant_pair_temp[-c(which(variantvariant_pair_temp[,1]==B & variantvariant_pair_temp[,2]==A )),]
            
            if (is.null(nrow(variantvariant_pair_temp))) 
            {
              variantvariant_pair_temp <- c()
              break
            }
            else
            {
              variantvariant_pair_temp <- variantvariant_pair_temp[-c(which.max(variantvariant_pair_temp[,3])),]
              
            }
            
          }
          # reverse id
          else
          {
            
            A = variantvariant_pair_temp[ which.max( variantvariant_pair_temp[,3] ),1]
            B = variantvariant_pair_temp[ which.max( variantvariant_pair_temp[,3] ),2]
            
            initial_amat[A ,B]<-3
            initial_amat[B ,A]<-2
            fci.est_temp[A ,B]<-3
            fci.est_temp[B ,A]<-2
            best_structure_initial  <- fci.est_temp
          }
          
          
          
          
          if (!is.null(nrow(variantvariant_pair_temp))) 
          {
            if (nrow(variantvariant_pair_temp)>0) 
            { 
              
              for ( row_index in 1:nrow(variantvariant_pair_temp))
              {
                
                if (variantvariant_pair_temp[row_index,1]==B && variantvariant_pair_temp[row_index,2]==A )
                {
                  variantvariant_pair_temp <- variantvariant_pair_temp[-c(row_index),]
                  break
                }
                
              }
              
            }
          }
          
          ######### Ending perform the orientation 
          
          
        }
        
      }
    }
    # incrementally add <-> 
    # loop over all o-> and o-o to be a <-> 
    
    for(index_variantdirect in 1:nrow(variant_pair))
    {
      
      
      fci.est_temp <- best_structure_initial 
      # reset the memory for the learnt bidirected 
      if ( current_bidirectedge - nrow(bidirect_pair) ==1)
      {
        for (bidirect_pair_index in 1:nrow(bidirect_pair) )
        {
          rest_confounder <- rest_confounder[-1]
          print("update the rest confounders")
        }
        fci.est_temp <- best_structure_initial  
      }
      
      #revert back to learnt <-> for the next 
      if(nrow(bidirect_pair) !=0)
      {
        for (index_row in 1:nrow(bidirect_pair))
        {
          fci.est_temp[bidirect_pair[index_row,1],bidirect_pair[index_row,2]] <- 2  
          fci.est_temp[bidirect_pair[index_row,2],bidirect_pair[index_row,1]] <- 2
        }
      }
      
      #revert back to initial dimention for the next 
      fci.est_temp <- fci.est_temp[-c(ncol(dataset)+1:ncol(fci.est_temp)),-c(ncol(dataset)+1:ncol(fci.est_temp))] 
      
      
      
      if(index_variantdirect==1) 
      {
        best_local_ELBO <- -Inf
        best_local_confounder <- c()
      }
      
      #generate <-> from the best confounder
      for (index_rest_confounder in 1:length(rest_confounder))
      {
        fci.est_temp[variant_pair[rest_confounder[index_rest_confounder],1],variant_pair[rest_confounder[index_rest_confounder],2]] <- 2  
        fci.est_temp[variant_pair[rest_confounder[index_rest_confounder],2],variant_pair[rest_confounder[index_rest_confounder],1]] <- 2  
        
      }
      
      
      
      if( !all(index_variantdirect %in% rest_confounder))
      {  
        
        
        if (current_bidirectedge>learnt_bidirectedge && nrow(bidirect_pair)==0)
        {
          
          for (index_rest_confounder in 1:length(rest_confounder))
          {
            fci.est_temp[variant_pair[rest_confounder[index_rest_confounder],1],variant_pair[rest_confounder[index_rest_confounder],2]] <- 2  
            fci.est_temp[variant_pair[rest_confounder[index_rest_confounder],2],variant_pair[rest_confounder[index_rest_confounder],1]] <- 2  
            
          }
        }
        
        fci.est_temp[variant_pair[index_variantdirect,1],variant_pair[index_variantdirect,2]] <- 2  
        fci.est_temp[variant_pair[index_variantdirect,2],variant_pair[index_variantdirect,1]] <- 2
        
        ####change o-> to -> ####
        fci.est_temp[ fci.est_temp ==1 & t(fci.est_temp)==2 ] <- 3
        
        partition_test <-bidirectededges_to_partition(fci.est_temp)
        current_latentconfounder <- length(partition_test)
        print(paste0("Update the number of latentconfounder to ",current_latentconfounder) )
        
        initial_amat <- construct_aux_latent_matrix(fci.est_temp,current_latentconfounder, partition_test)
        
        ###remove <->  between children to be added by A<-LV->B 
        currentbidirectpair <- which(fci.est_temp ==2 & t(fci.est_temp)==2, arr.ind=TRUE)
        currentbidirectpair <- remove_duplicated_pair(currentbidirectpair)
        if (is.null(nrow(currentbidirectpair))) 
        {
          currentbidirectpair <- matrix(currentbidirectpair,1,2)
        }
        for (index_row in 1:nrow(currentbidirectpair))
        {
          initial_amat[currentbidirectpair[index_row,2],currentbidirectpair[index_row,1]]<-0
          initial_amat[currentbidirectpair[index_row,1],currentbidirectpair[index_row,2]]<-0
          
        }
        
        
        
        
        
      }
      
      
      if(checkcyclic(initial_amat))
      {
        
        
        ######## compute the local score for the <-> index
        computeELBOinitial <- .jnew("experiments.real.computeELBOinitial")
        input <- matrix(as.integer(initial_amat),ncol=ncol(initial_amat))
        current_ELBO  <- .jcall(computeELBOinitial, returnSig="D", method="excute",.jarray(input, dispatch = T), as.integer(current_latentconfounder),as.integer(2),paste0(filepath,".arff"))
        list_ELBO <- c(list_ELBO,current_ELBO)
        write.csv(list_ELBO, paste0("output/HCLC_",algo,"_",set,"_listELBO_",latent1,"_",n,".csv"), row.names=F)
        list_confounder <- c(list_confounder,current_latentconfounder)
        
        write.csv(list_confounder, paste0("output/HCLC_",algo,"_",set,"_list_confounder_",latent1,"_",n,".csv"), row.names=F)
        
        
        
      }
      
      
      if (current_ELBO>best_local_ELBO)
      {
        print("Score improve locally")
        best_local_ELBO <- current_ELBO
        best_local_confounder <- c(index_variantdirect)
      }
      if (current_ELBO>best_ELBO)
      {
        best_structure <- initial_amat
        print("Improved score globlally")
        write.csv(best_ELBO, paste0("output/HCLC_",algo,"_",set,"_bestELBO_",latent1,"_",n,".csv"), row.names=F)
        write.csv(best_structure, paste0("output/HCLC_",algo,"_",set,"_bestdag_",latent1,"_",n,".csv"), row.names=F)
        
      }
      
      print(paste0("the best confounder =",rest_confounder))
      print(paste0("the best ELBO =",best_ELBO))
      
      #keep the best confounder at the first iteration
      if(index_variantdirect==1 && current_bidirectedge==1)
      {
        best_structure <- initial_amat
        best_local_ELBO <- current_ELBO
      }
      
      
      
      ######## ending compute the local score for the <-> index
      ####### break for the time limite ########
      end_time <- Sys.time()
      time_lab = difftime(end_time ,start_time, units = "secs")[[1]]
      print(paste0("Elap time is :", time_lab))
      if  (time_lab > time_max)
      {
        break
      }
      ####### end break for the time limite ########
      
    }
    
  }
  
  
  
  
  
  else print("Error !!")
  
  
  ####### break for the time limite ########
  end_time <- Sys.time()
  time_lab = difftime(end_time ,start_time, units = "secs")[[1]]
  print(paste0("Elap time is :", time_lab))
  if  (time_lab > time_max)
  {
    break
  }
  ####### end break for the time limite ########
  
  print (paste0("best_local_ELBO=",best_local_ELBO , " best_ELBO=",best_ELBO))
  if(best_local_ELBO >= best_ELBO )
  {
    best_ELBO <- best_local_ELBO
    rest_confounder <- c(rest_confounder,best_local_confounder)
  }
  else
  {
    score_improve = FALSE
  }  
  
  
  current_bidirectedge=current_bidirectedge+1
  
}

################### END Incrementally search for MAG  ##############





########## Compute ELBO BIC and LL for the DAG with the additional states of auxiliary Latent confounders #############

ELBO <- list()
BIC <- list()
LL <- list()
D <- list()

current_latentconfounder =  ncol(best_structure) - ncol(fci.est)
computeELBO <- .jnew("experiments.real.computeELBOgreedylocalState")

input <- matrix(as.integer(best_structure),ncol=ncol(best_structure))

print(paste0("Calculating number of states the final DAG"))

bestELBO <- .jcall(computeELBO, returnSig="D", method="excute",.jarray(input, dispatch = T), as.integer(current_latentconfounder),as.integer(max_state),paste0(filepath,".arff"))

#### Stop computing ELBO BIC and LL for each DAG with the auxiliary Latent confounders ##


########## write outputs ##########
print(paste0("Writing the model"))

.jcall(computeELBO, method="exportGenieModel", returnSig="V",paste0("output/","HCLC_",algo,"_",set,"_",latent1,"_",n))

best_BIC <- .jcall(computeELBO, returnSig="D", method="getBIC")
best_LL <- .jcall(computeELBO, returnSig="D", method="getLL")
best_D <- .jcall(computeELBO, returnSig="I", method="getdim")

write.csv(bestELBO, paste0("output/HCLC_",algo,"_",set,"_bestELBO_",latent1,"_",n,".csv"), row.names=F)
write.csv(best_BIC, paste0("output/HCLC_",algo,"_",set,"_bestBIC_",latent1,"_",n,".csv"), row.names=F)
write.csv(best_LL, paste0("output/HCLC_",algo,"_",set,"_bestLL_",latent1,"_",n,".csv"), row.names=F)  
write.csv(best_D, paste0("output/HCLC_",algo,"_",set,"_dim_",latent1,"_",n,".csv"), row.names=F)  
write.csv(best_structure, paste0("output/HCLC_",algo,"_",set,"_bestdag_",latent1,"_",n,".csv"), row.names=F)
########## End writing outputs ##########

############## STOP #########################
toc(log = TRUE, quiet = TRUE)
log.txt <- tic.log(format = TRUE)
#############################################
write.csv(log.txt, paste0("output/HCLC_",algo,"_",set,"_runtime_",latent1,"_",n,".csv"), row.names=F)

############################################################################################################################################





