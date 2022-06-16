#v7 HC the structure and greedy with states
#v5 greedy with the structure and state
#v4 use m-sep to indentify the min of latent confounders
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
source("listMags.R")
############## For preliminary setup #########################
#require(devtools)
#install_version("rJava", version = "0.9-12", repos = "http://cran.us.r-project.org")
#install.packages("stringr")
#install.packages("rJava", type="source")
#library(devtools)
#install_github("bd2kccd/r-causal", INSTALL_opts=c("--no-multiarch"),force = TRUE)
#remotes::install_github("bd2kccd/r-causal")
#install.packages("rJava",type='source',"https://cran.r-project.org/src/contrib/Archive/rJava/rJava_0.9-12.tar.gz") 
#Sys.setenv(JAVA_HOME="C:\\Program Files\\Java\\jdk-13.0.1")
#.jinit()
#.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
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

#latent1 <- c("HYPOVOLEMIA", "LVFAILURE","ERRCAUTER","PULMEMBOLUS","INTUBATION","KINKEDTUBE")
#latent1 <- c("HYPOVOLEMIA", "LVFAILURE","ERRCAUTER")
#latent1 <- c("PULMEMBOLUS","INTUBATION","KINKEDTUBE")
#latent1 <- c("INTUBATION") #example of 3 latent conf.
#latent1 <- c("ERRCAUTER") #example of1 latent conf.



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
#write.arff(new_sampling, file = paste0(filepath,"_",n,".arff") )
#write.csv(new_sampling,paste0(filepath,"_",n,".csv"),row.names = FALSE)

################ End generating synthetic data #####################


################### load observational dataset ###################
#dataset <- new_sampling
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

##############################################################

if (algo =='FCI')
{
  
  ################### run FCI ###################
  
  ####### pcalg lib  ###
  varNames <- colnames(dataset)
  fci.est <- fci(suffStat, indepTest = disCItest,skel.method =c("stable"), label=c(varNames_pcalg), alpha=alpha_stat)
  #fci.est <- fci(suffStat, indepTest = gaussCItest, label=c(varNames), alpha=alpha_stat,conservative = FALSE,maj.rule = FALSE  )
  fci.est <- fci.est@amat
  
  
  varNames <- colnames(fci.est)
  #varNames <- c(varNames,"LV1")
  x<- adjtodot(fci.est, varNames)
  dot(x)
  ############################################
} 
if (algo =='GFCI')
{
  ####### tetrad  ###
  tetradrunner.getAlgorithmDescription(algoId = 'fci')
  #tetradrunner.getAlgorithmParameters(algoId = 'fges',scoreId = 'sem-bic')
  #tetradrunner.getAlgorithmParameters(algoId = 'gfci',testId = 'chi-square-test',scoreId = 'sem-bic')
  #tetradrunner.listIndTests()
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
                               faithfulnessAssumed = FALSE,
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

# latentconfounder for o-o 
variantvariant_pair <- which((fci.est ==1 & t(fci.est)==1) , arr.ind=TRUE)
if (nrow(variantdirect_pair) >0) 
{
  variantvariant_pair <- remove_duplicated_pair(variantvariant_pair)
}


if (is.null(nrow(variantvariant_pair))) 
{
  variantvariant_pair <- matrix(variantvariant_pair,1,2)
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
best_structure <- c()
best_confounders <-c()
rest_confounder <- c()
list_ELBO <-c()
###### End adding the pair #####


#### list all MAG and remove MAG that contains <->
fci.est_temp <-fci.est
listMag <-list()
listMag <- listMags(fci.est_temp, nMags = 1000, method="local")
computeELBOinitial <- .jnew("experiments.real.computeELBOinitial")

while(score_improve && current_bidirectedge <=max_latentconfounder)
{ 
  print(paste0("Increased the number of bidirected edges =", current_bidirectedge) )
  
  ####### use <-> to consider the min latentconfounder
  if (length(which(fci.est ==2 & t(fci.est)==2))/2 >0  & current_bidirectedge<=length(which(fci.est ==2 & t(fci.est)==2))/2)
  {
        
        index_listMag_toremove <-c()
        listMag_temp <- listMag
        for(index_listMag in 1:length(listMag_temp))
        {
          
          no_bidirectedge <- length(which(listMag_temp[[index_listMag]] ==2 & t(listMag_temp[[index_listMag]])==2))/2
          
           if(current_bidirectedge!=no_bidirectedge)
          {
            index_listMag_toremove <- append(index_listMag_toremove, index_listMag)
          }
          
        }
        if ( is.null(index_listMag_toremove) == FALSE)
        {
          for(index_toremove in length(index_listMag_toremove):1)
          {
          
            listMag_temp[[index_listMag_toremove[index_toremove]]] <-NULL
          }
        }
        list_ELBO <- c(list_ELBO,paste0("bidirected= ",current_bidirectedge ," total=",length(listMag_temp)))
        for(index_listMag in 1:length(listMag_temp))
        {
          write.csv(length(listMag_temp), paste0("output/ILC_",algo,"_",set,"_nummag_",latent1,"_",n,".csv"), row.names=F)
          
          partition_test <-bidirectededges_to_partition(listMag_temp[[index_listMag]])
          current_latentconfounder <- length(partition_test)
          print(paste0("Update the number of latentconfounder to ",current_latentconfounder) )
          initial_amat <- construct_aux_latent_matrix(listMag_temp[[index_listMag]],current_latentconfounder, partition_test)
          input <- matrix(as.integer(initial_amat),ncol=ncol(initial_amat))
          current_ELBO  <- .jcall(computeELBOinitial, returnSig="D", method="excute",.jarray(input, dispatch = T), as.integer(current_latentconfounder),as.integer(2),paste0(filepath,".arff"))
          list_ELBO <- c(list_ELBO,current_ELBO)
       
          write.csv(list_ELBO, paste0("output/ILC_",algo,"_",set,"_listELBO_",latent1,"_",n,".csv"), row.names=F)
          if (current_ELBO>best_ELBO)
          {
            best_ELBO <- current_ELBO
            best_structure <-initial_amat
            write.csv(best_ELBO, paste0("output/ILC_",algo,"_",set,"_bestELBO_",latent1,"_",n,".csv"), row.names=F)
            write.csv(best_structure, paste0("output/ILC_",algo,"_",set,"_bestdag_",latent1,"_",n,".csv"), row.names=F)
            
            
          }
          
          
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
        best_local_ELBO <- best_ELBO

      
  }
  ####### if no more <->, incrementally add a latentconfounder from o-> or o-o 
  else if (length(which(fci.est ==2 & t(fci.est)==2))/2 < current_bidirectedge)
  {
    best_local_ELBO <- -Inf
    index_listMag_toremove <-c()
    listMag_temp <- listMag
  
      
    for(index_listMag in 1:length(listMag_temp))
    {
      
      no_bidirectedge <- length(which(listMag_temp[[index_listMag]] ==2 & t(listMag_temp[[index_listMag]])==2))/2
      print(no_bidirectedge)
      if(current_bidirectedge!=no_bidirectedge)
      {
        index_listMag_toremove <- append(index_listMag_toremove, index_listMag)
      }
      
    }
    if ( is.null(index_listMag_toremove) == FALSE)
    {
      for(index_toremove in length(index_listMag_toremove):1)
      {
       
        listMag_temp[[index_listMag_toremove[index_toremove]]] <-NULL
      }
    }
    list_ELBO <- c(list_ELBO,paste0("bidirected= ",current_bidirectedge ," total=",length(listMag_temp)))
    
    if(length(listMag_temp)!=0)
    {
    for(index_listMag in 1:length(listMag_temp))
    {
      write.csv(length(listMag_temp), paste0("output/ILC_",algo,"_",set,"_nummag_",latent1,"_",n,".csv"), row.names=F)
      
      partition_test <-bidirectededges_to_partition(listMag_temp[[index_listMag]])
      current_latentconfounder <- length(partition_test)
      print(paste0("Update the number of latentconfounder to ",current_latentconfounder) )
      initial_amat <- construct_aux_latent_matrix(listMag_temp[[index_listMag]],current_latentconfounder, partition_test)
      input <- matrix(as.integer(initial_amat),ncol=ncol(initial_amat))
      current_ELBO  <- .jcall(computeELBOinitial, returnSig="D", method="excute",.jarray(input, dispatch = T), as.integer(current_latentconfounder),as.integer(2),paste0(filepath,".arff"))
      list_ELBO <- c(list_ELBO,current_ELBO)
      write.csv(list_ELBO, paste0("output/ILC_",algo,"_",set,"_listELBO_",latent1,"_",n,".csv"), row.names=F)
      
      if (current_ELBO>best_ELBO)
      {
        best_ELBO <- current_ELBO
        best_structure <-initial_amat
        write.csv(best_ELBO, paste0("output/ILC_",algo,"_",set,"_bestELBO_",latent1,"_",n,".csv"), row.names=F)
        write.csv(best_structure, paste0("output/ILC_",algo,"_",set,"_bestdag_",latent1,"_",n,".csv"), row.names=F)
        
      }
      if(current_ELBO > best_local_ELBO)
      {
        best_local_ELBO <- current_ELBO
      }  
      
      
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
    else break
   
  }

  ####### break for the time limite ########
  end_time <- Sys.time()
  time_lab = difftime(end_time ,start_time, units = "secs")[[1]]
  print(paste0("Elap time is :", time_lab))
  if  (time_lab > time_max)
  {
    break
  }
  ####### end break for the time limite ########
  
  if(best_ELBO > best_local_ELBO)
  {
    score_improve = FALSE
  } 
  
  current_bidirectedge=current_bidirectedge+1
  
  
  
  
}



##################### 


########## Add the auxiliary Latent confounders    #############

ELBO <- list()
BIC <- list()
LL <- list()
D <- list()


#best_structure<-  read.matrix("output/ILC_GFCI_PROPERTY_bestdag_propertyPurchaseValue.csv",header = TRUE,na.strings=c(""),check.name=FALSE)
#best_structure <- matrix(unlist(best_structure), ncol =ncol(best_structure), nrow =ncol(best_structure))
########## Compute ELBO BIC and LL for the DAG with the additional states of auxiliary Latent confounders #############

current_latentconfounder =  ncol(best_structure) - ncol(fci.est)
computeELBO <- .jnew("experiments.real.computeELBOgreedylocalState")

input <- matrix(as.integer(best_structure),ncol=ncol(best_structure))

print(paste0("Calculating number of states the final DAG"))

bestELBO <- .jcall(computeELBO, returnSig="D", method="excute",.jarray(input, dispatch = T), as.integer(current_latentconfounder),as.integer(max_state),paste0(filepath,".arff"))


########## write ##########
print(paste0("Writing the model"))

.jcall(computeELBO, method="exportGenieModel", returnSig="V",paste0(filepath,"_hc_",algo,"_",n))

best_BIC <- .jcall(computeELBO, returnSig="D", method="getBIC")
best_LL <- .jcall(computeELBO, returnSig="D", method="getLL")
best_D <- .jcall(computeELBO, returnSig="I", method="getdim")


write.csv(bestELBO, paste0("output/ILC_",algo,"_",set,"_bestELBO_",latent1,"_",n,".csv"), row.names=F)
write.csv(best_BIC, paste0("output/ILC_",algo,"_",set,"_bestBIC_",latent1,"_",n,".csv"), row.names=F)
write.csv(best_LL, paste0("output/ILC_",algo,"_",set,"_bestLL_",latent1,"_",n,".csv"), row.names=F)  
write.csv(best_D, paste0("output/ILC_",algo,"_",set,"_dim_",latent1,"_",n,".csv"), row.names=F)  
write.csv(best_structure, paste0("output/ILC_",algo,"_",set,"_bestdag_",latent1,"_",n,".csv"), row.names=F)

############## STOP #########################
toc(log = TRUE, quiet = TRUE)
log.txt <- tic.log(format = TRUE)
#############################################
write.csv(log.txt, paste0("output/ILC_",algo,"_",set,"_runtime_",latent1,"_",n,".csv"), row.names=F)

#  ############################




#### Stop computing ELBO BIC and LL for each DAG with the auxiliary Latent confounders ##


