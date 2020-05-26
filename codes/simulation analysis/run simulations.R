# Author: Giovanni Nattino 
# Simulation analysis comparing NN matching and 
#  3-way conditional matching in 3 groups. 
# Reference: Triplet Matching for Estimating Causal Effects with Three Treatment Arms: 
#  A Comparative Study of Mortality by Trauma Center Level
######################################################################################

library(optmatch)
library(doSNOW)
library(parallel)
library(doRNG)

source("lib/functionsMatching_7.R")
source("lib/functionsDataGeneration_5.R")

##############
# Parameters #
##############

#Parameters for matching algorithms
#----------------
#3-way optimal matching - change of starting edge
startingEdges <- c("1-2","2-3","1-3")

#Variable(s) to match on
variablesMatch <- "varMatching"
variableTreatment <- "treatment"

################
# SIMULATION 1 #
################

#Number of iterations
N <- 1000

set.seed(13112016)
listResults_3wayCond <- list()

start <- Sys.time()
for(iter in 1:N) {
  
  #H0 data
  #--------------
  n1 <- 100
  n2 <- 200
  n3 <- 300
  mu <- 0; sd <- 1
  data <- generateData_H0(n1, n2, n3, mu, sd)
  
  temp <- diff(sort(data$varMatching))
  data$varMatching <- data$varMatching/min(temp[temp>0])
  
  #3-way constrained optimal matching
  #----------------------
  result3wayConstr <- applyModifiedOptMatching(data, startingEdges, variablesMatch = variablesMatch,
                                               variableTreatment = variableTreatment)
  
  listResults_3wayCond[[iter]] <- result3wayConstr$results
}

end <- Sys.time()
end-start

save(listResults_3wayCond, file = paste("result simulations/resultSimulations_H0_n1.100_n2.200_n3.300.Rdata",sep=""))

################
# SIMULATION 2 #
################

#Total sample sizes: 1,500-2,000-2,500
#-------------------------------------
#Number of iterations
N <- 1000   

listParmaters <- list(c(alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 3, alpha2 = 3, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 4, alpha2 = 4, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 5, alpha2 = 5, beta2 = 2, alpha3 = 2, beta3 = 2))
listSampleSizes <- list(c(n1 = 500, n2 = 500, n3 = 500),
                        c(n1 = 1000, n2 = 1000, n3 = 500),
                        c(n1 = 500, n2 = 500, n3 = 1000),
                        c(n1 = 1000, n2 = 500, n3 = 500))

set.seed(123456)

cl <- makeCluster(6)
registerDoSNOW(cl)


startOVERALL <- Sys.time()

for(parameters in listParmaters) {
  for(sampleSizes in listSampleSizes) {
    
    startParam <- Sys.time()
    
    print(parameters)
    print(sampleSizes)
    
    pb <- txtProgressBar(max = N, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    resultTemp <- foreach(iter=1:N, 
                          .packages = c("optmatch"), .options.snow = opts) %dorng% {
                            
                            
                            data <- generateData_Betas(sampleSizes[["n1"]], sampleSizes[["n2"]], sampleSizes[["n3"]], 
                                                       parameters[["alpha1"]], parameters[["beta1"]], 
                                                       parameters[["alpha2"]], parameters[["beta2"]], 
                                                       parameters[["alpha3"]], parameters[["beta3"]]) 
                            
                            #3-way constrained optimal matching
                            #----------------------
                            result3wayConstr <- applyModifiedOptMatching(data, startingEdges, variablesMatch = variablesMatch,
                                                                         variableTreatment = variableTreatment)
                            
                            #Nearest Neighbors
                            #----------------------
                            resultsNN <- applyNearNeigh(data, variablesMatch = variablesMatch,
                                                        variableTreatment = variableTreatment)
                            
                            
                            #If 3-way constrained optimal matching worse than NN: apply 3 way constr on result of NN
                            if(min(result3wayConstr$results$overallDistance) > resultsNN$results$matchingSchemesDistances) {
                              
                              resultsNN$matchedData$distance <- NULL
                              result3wayConstr_postNN <- applyModifiedOptMatching_startMatch(data = resultsNN$matchedData, 
                                                                                             startingEdges = startingEdges, 
                                                                                             variablesMatch = variablesMatch,
                                                                                             variableTreatment = variableTreatment, 
                                                                                             varInitialMatch = "indexMatch")
                              
                            } else {
                              
                              result3wayConstr_postNN <- list(results = NULL,
                                                              matchedData=NULL)
                              
                            }
                            
                            list(resultsNN = resultsNN$results, 
                                 result3wayConstr = result3wayConstr$results,
                                 result3wayConstr_postNN = result3wayConstr_postNN$results)
                          }
    
    listResults_NN <- lapply(resultTemp, FUN = "[[", "resultsNN")
    listResults_3wayCond <- lapply(resultTemp, FUN = "[[", "result3wayConstr")
    listResults_3wayCond_postNN <- lapply(resultTemp, FUN = "[[", "result3wayConstr_postNN")
    
    save(listResults_3wayCond, file = paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],
                                            "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                            "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                            "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                            "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    save(listResults_NN, file = paste("result simulations/result_NN_alpha1.",parameters[["alpha1"]],
                                      "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                      "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                      "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                      "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    save(listResults_3wayCond_postNN, file = paste("result simulations/result_3wayCond_postNN_alpha1.",parameters[["alpha1"]],
                                                   "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                                   "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                                   "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                                   "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    endParam <- Sys.time()
    print("")
    print(endParam-startParam)
    print("********************************************************************************")
    print("")
  }
}

endOVERALL <- Sys.time()
endOVERALL-startOVERALL
stopCluster(cl)


#Total sample sizes: 300-400-500
#-------------------------------

#Number of iterations
N <- 1000   

listParmaters <- list(c(alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 3, alpha2 = 3, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 4, alpha2 = 4, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 5, alpha2 = 5, beta2 = 2, alpha3 = 2, beta3 = 2))
listSampleSizes <- list(c(n1 = 100, n2 = 100, n3 = 100),
                        c(n1 = 200, n2 = 200, n3 = 100),
                        c(n1 = 100, n2 = 100, n3 = 200),
                        c(n1 = 200, n2 = 100, n3 = 100))

listResults_2waySeq <- list()
listResults_3wayCond <- list()
listResults_NN <- list()

set.seed(13112016)

start <- Sys.time()

for(parameters in listParmaters) {
for(sampleSizes in listSampleSizes) {
for(iter in 1:N) {
  
  data <- generateData_Betas(sampleSizes[["n1"]], sampleSizes[["n2"]], sampleSizes[["n3"]], 
                             parameters[["alpha1"]], parameters[["beta1"]], 
                             parameters[["alpha2"]], parameters[["beta2"]], 
                             parameters[["alpha3"]], parameters[["beta3"]]) 
  
  #3-way constrained optimal matching
  #----------------------
  result3wayConstr <- applyModifiedOptMatching(data, startingEdges, variablesMatch = variablesMatch,
                                               variableTreatment = variableTreatment)

  listResults_3wayCond[[iter]] <- result3wayConstr$results
  
  #Nearest Neighbors
  #----------------------
  resultsNN <- applyNearNeigh(data, variablesMatch = variablesMatch,
                              variableTreatment = variableTreatment)
  listResults_NN[[iter]] <- resultsNN$results
  
  if(iter %% 10 == 1) { print(iter); }
}

save(listResults_3wayCond, file = paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],
                                        "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                        "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                        "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                        "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
save(listResults_NN, file = paste("result simulations/result_NN_alpha1.",parameters[["alpha1"]],
                                       "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                       "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                       "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                       "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))

}
}

end <- Sys.time()
end-start


#Total sample sizes: 2 x 300-400-500
#-------------------------------
#Number of iterations
N <- 100

listParmaters <- list(c(alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 3, alpha2 = 3, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 4, alpha2 = 4, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 5, alpha2 = 5, beta2 = 2, alpha3 = 2, beta3 = 2))
listSampleSizes <- list(c(n1 = 200, n2 = 200, n3 = 200),
                        c(n1 = 400, n2 = 400, n3 = 200),
                        c(n1 = 200, n2 = 200, n3 = 400),
                        c(n1 = 400, n2 = 200, n3 = 200))

listResults_2waySeq <- list()
listResults_3wayCond <- list()
listResults_NN <- list()

set.seed(13112016)

start <- Sys.time()

for(parameters in listParmaters) {
  for(sampleSizes in listSampleSizes) {
    for(iter in 1:N) {
      
      data <- generateData_Betas(sampleSizes[["n1"]], sampleSizes[["n2"]], sampleSizes[["n3"]], 
                                 parameters[["alpha1"]], parameters[["beta1"]], 
                                 parameters[["alpha2"]], parameters[["beta2"]], 
                                 parameters[["alpha3"]], parameters[["beta3"]]) 
      
      #3-way constrained optimal matching
      #----------------------
      result3wayConstr <- applyModifiedOptMatching(data, startingEdges, variablesMatch = variablesMatch,
                                                   variableTreatment = variableTreatment)
      
      listResults_3wayCond[[iter]] <- result3wayConstr$results
      
      #Nearest Neighbors
      #----------------------
      resultsNN <- applyNearNeigh(data, variablesMatch = variablesMatch,
                                  variableTreatment = variableTreatment)
      listResults_NN[[iter]] <- resultsNN$results
      
      if(iter %% 10 == 1) { print(iter); }
    }
    
    save(listResults_3wayCond, file = paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],
                                            "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                            "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                            "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                            "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    save(listResults_NN, file = paste("result simulations/result_NN_alpha1.",parameters[["alpha1"]],
                                      "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                      "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                      "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                      "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    
  }
}

end <- Sys.time()
end-start

#Total sample sizes: 3 x 300-400-500
#-------------------------------
#Number of iterations
N <- 100

listParmaters <- list(c(alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 3, alpha2 = 3, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 4, alpha2 = 4, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 5, alpha2 = 5, beta2 = 2, alpha3 = 2, beta3 = 2))
listSampleSizes <- list(c(n1 = 300, n2 = 300, n3 = 300),
                        c(n1 = 600, n2 = 600, n3 = 300),
                        c(n1 = 300, n2 = 300, n3 = 600),
                        c(n1 = 600, n2 = 300, n3 = 300))

listResults_2waySeq <- list()
listResults_3wayCond <- list()
listResults_NN <- list()

set.seed(13112016)

start <- Sys.time()

for(parameters in listParmaters) {
  for(sampleSizes in listSampleSizes) {
    for(iter in 1:N) {
      
      data <- generateData_Betas(sampleSizes[["n1"]], sampleSizes[["n2"]], sampleSizes[["n3"]], 
                                 parameters[["alpha1"]], parameters[["beta1"]], 
                                 parameters[["alpha2"]], parameters[["beta2"]], 
                                 parameters[["alpha3"]], parameters[["beta3"]]) 
      
      #3-way constrained optimal matching
      #----------------------
      result3wayConstr <- applyModifiedOptMatching(data, startingEdges, variablesMatch = variablesMatch,
                                                   variableTreatment = variableTreatment)
      
      listResults_3wayCond[[iter]] <- result3wayConstr$results
      
      #Nearest Neighbors
      #----------------------
      resultsNN <- applyNearNeigh(data, variablesMatch = variablesMatch,
                                  variableTreatment = variableTreatment)
      listResults_NN[[iter]] <- resultsNN$results
      
      if(iter %% 10 == 1) { print(iter); }
    }
    
    save(listResults_3wayCond, file = paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],
                                            "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                            "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                            "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                            "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    save(listResults_NN, file = paste("result simulations/result_NN_alpha1.",parameters[["alpha1"]],
                                      "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                      "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                      "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                      "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    
  }
}

end <- Sys.time()
end-start

#Total sample sizes: 4 x 300-400-500
#-------------------------------
#Number of iterations
N <- 100

listParmaters <- list(c(alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 3, alpha2 = 3, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 4, alpha2 = 4, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 5, alpha2 = 5, beta2 = 2, alpha3 = 2, beta3 = 2))
listSampleSizes <- list(c(n1 = 400, n2 = 400, n3 = 400),
                        c(n1 = 800, n2 = 800, n3 = 400),
                        c(n1 = 400, n2 = 400, n3 = 800),
                        c(n1 = 800, n2 = 400, n3 = 400))

listResults_2waySeq <- list()
listResults_3wayCond <- list()
listResults_NN <- list()

set.seed(13112016)

start <- Sys.time()

for(parameters in listParmaters) {
  for(sampleSizes in listSampleSizes) {
    for(iter in 1:N) {
      
      data <- generateData_Betas(sampleSizes[["n1"]], sampleSizes[["n2"]], sampleSizes[["n3"]], 
                                 parameters[["alpha1"]], parameters[["beta1"]], 
                                 parameters[["alpha2"]], parameters[["beta2"]], 
                                 parameters[["alpha3"]], parameters[["beta3"]]) 
      
      #3-way constrained optimal matching
      #----------------------
      result3wayConstr <- applyModifiedOptMatching(data, startingEdges, variablesMatch = variablesMatch,
                                                   variableTreatment = variableTreatment)
      
      listResults_3wayCond[[iter]] <- result3wayConstr$results
      
      #Nearest Neighbors
      #----------------------
      resultsNN <- applyNearNeigh(data, variablesMatch = variablesMatch,
                                  variableTreatment = variableTreatment)
      listResults_NN[[iter]] <- resultsNN$results
      
      if(iter %% 10 == 1) { print(iter); }
    }
    
    save(listResults_3wayCond, file = paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],
                                            "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                            "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                            "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                            "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    save(listResults_NN, file = paste("result simulations/result_NN_alpha1.",parameters[["alpha1"]],
                                      "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                      "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                      "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                      "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    
  }
}

end <- Sys.time()
end-start


#Total sample sizes: 10 x 300-400-500
#-------------------------------

#Number of iterations
N <- 100

listParmaters <- list(c(alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 3, alpha2 = 3, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 4, alpha2 = 4, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 5, alpha2 = 5, beta2 = 2, alpha3 = 2, beta3 = 2))
listSampleSizes <- list(c(n1 = 1000, n2 = 1000, n3 = 1000),
                        c(n1 = 2000, n2 = 2000, n3 = 1000),
                        c(n1 = 1000, n2 = 1000, n3 = 2000),
                        c(n1 = 2000, n2 = 1000, n3 = 1000))

set.seed(123456)

cl <- makeCluster(1)
registerDoSNOW(cl)


startOVERALL <- Sys.time()

for(parameters in listParmaters) {
  for(sampleSizes in listSampleSizes) {
    
    startParam <- Sys.time()
    
    print(parameters)
    print(sampleSizes)
    
    pb <- txtProgressBar(max = N, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    resultTemp <- foreach(iter=1:N, 
                          .packages = c("optmatch"), .options.snow = opts) %dorng% {
                            
                            
                            data <- generateData_Betas(sampleSizes[["n1"]], sampleSizes[["n2"]], sampleSizes[["n3"]], 
                                                       parameters[["alpha1"]], parameters[["beta1"]], 
                                                       parameters[["alpha2"]], parameters[["beta2"]], 
                                                       parameters[["alpha3"]], parameters[["beta3"]]) 
                            
                            #3-way constrained optimal matching
                            #----------------------
                            result3wayConstr <- applyModifiedOptMatching(data, startingEdges, variablesMatch = variablesMatch,
                                                                         variableTreatment = variableTreatment)
                            
                            #Nearest Neighbors
                            #----------------------
                            resultsNN <- applyNearNeigh(data, variablesMatch = variablesMatch,
                                                        variableTreatment = variableTreatment)
                            
                            
                            #If 3-way constrained optimal matching worse than NN: apply 3 way constr on result of NN
                            if(min(result3wayConstr$results$overallDistance) > resultsNN$results$matchingSchemesDistances) {
                            
                            resultsNN$matchedData$distance <- NULL
                            result3wayConstr_postNN <- applyModifiedOptMatching_startMatch(data = resultsNN$matchedData, 
                                                                                            startingEdges = startingEdges, 
                                                                                            variablesMatch = variablesMatch,
                                                                                            variableTreatment = variableTreatment, 
                                                                                            varInitialMatch = "indexMatch")
                            
                            } else {
                              
                              result3wayConstr_postNN <- list(results = NULL,
                                                              matchedData=NULL)
                            
                            }
                            
                            list(resultsNN = resultsNN$results, 
                                 result3wayConstr = result3wayConstr$results,
                                 result3wayConstr_postNN = result3wayConstr_postNN$results)
                          }
    
    listResults_NN <- lapply(resultTemp, FUN = "[[", "resultsNN")
    listResults_3wayCond <- lapply(resultTemp, FUN = "[[", "result3wayConstr")
    listResults_3wayCond_postNN <- lapply(resultTemp, FUN = "[[", "result3wayConstr_postNN")
    
    save(listResults_3wayCond, file = paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],
                                            "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                            "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                            "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                            "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    save(listResults_NN, file = paste("result simulations/result_NN_alpha1.",parameters[["alpha1"]],
                                      "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                      "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                      "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                      "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    save(listResults_3wayCond_postNN, file = paste("result simulations/result_3wayCond_postNN_alpha1.",parameters[["alpha1"]],
                                      "_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],
                                      "_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],
                                      "_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],
                                      "_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    endParam <- Sys.time()
    print("")
    print(endParam-startParam)
    print("********************************************************************************")
    print("")
  }
}

endOVERALL <- Sys.time()
endOVERALL-startOVERALL
stopCluster(cl)

############################
# EVALUATION OF COMPLEXITY #
############################

#Results not shown in the paper, only mentioned in the discussion.
set.seed(13112016)

#Number of iterations
N <- 20

#Only one starting edge
startingEdges <- c("1-2")

#Parameters simulations: symmetric scenario with all peaks equal
alpha1 <- 2  
beta1 <- 2
alpha2 <- 2
beta2 <- 2
alpha3 <- 2
beta3 <- 2

#Original size of treatment groups
n1_orig <- 100
n2_orig <- 100
n3_orig <- 100

#Factors for sample sizes: from 100 to 2,000 per group
vectorK <- c(.5,1:10,12,14,16,18,20,25)

for(indexK in 1:length(vectorK)) {
  
  k <- vectorK[indexK]
  
  listResultK <- list(k=vectorK[[indexK]])
  listResultK_iter <- list()

  n1 <- n1_orig*k
  n2 <- n2_orig*k
  n3 <- n3_orig*k
  
  print(n1)
  
  for(iter in 1:N) {
   
    data <- generateData_Betas(n1, n2, n3, 
                               alpha1, beta1, 
                               alpha2, beta2, 
                               alpha3, beta3) 
    
    start <- Sys.time()
  
    #3-way constrained optimal matching
    #----------------------
    result3wayConstr <- applyModifiedOptMatching(data, startingEdges, variablesMatch = variablesMatch,
                                                 variableTreatment = variableTreatment)
    
    end <- Sys.time()
    time3way <- as.numeric(difftime(end,start, units = "secs"))
  
    time3way
    
    #Pair optimal matching on largest treatment groups
    #-------------------------------------------------
    dataRed <- data[data$treatment %in% c(1,2),]
    dataRed$treatment <- as.factor(dataRed$treatment)
    
    start <- Sys.time()
    resultPairMatch <- pairmatch(treatment~varMatching, data = dataRed)
    end <- Sys.time()
    
    timePair <- as.numeric(difftime(end,start, units = "secs"))
    
    listResultK_iter[[iter]] <- list(time3way = time3way,
                                     timePair = timePair,
                                     numberIterations = result3wayConstr$results$iterations)
    
    print(iter)
  } 
  
  listResultK$results <- listResultK_iter
  
  save(listResultK, file = paste0("resultSimulations_",k,".Rdata",sep=""))
  
  print("---------------------------------------------------------")
}


