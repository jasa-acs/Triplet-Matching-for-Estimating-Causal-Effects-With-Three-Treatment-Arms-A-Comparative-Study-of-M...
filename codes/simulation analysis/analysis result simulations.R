# Author: Giovanni Nattino 
# File to analyze the results of the simulation analysis comparing NN matching and 
#  3-way conditional matching in 3 groups. 
# Reference: Triplet Matching for Estimating Causal Effects with Three Treatment Arms: 
#  A Comparative Study of Mortality by Trauma Center Level
######################################################################################

library(optmatch)

source("lib/functionsMatching_7.R")
source("lib/functionsDataGeneration_5.R")

################
# SIMULATION 1 #
################

# Existence of "perfect" matched triplets in the sample
#-------------------------------------------------------------------------------
load(file="result simulations/resultSimulations_H0_n1.100_n2.200_n3.300.Rdata")

minTotalDistances_3wayCond <- as.numeric(sapply(sapply(listResults_3wayCond, FUN="[", "overallDistance"), FUN=min))
maxTotalDistances_3wayCond <- as.numeric(sapply(sapply(listResults_3wayCond, FUN="[", "overallDistance"), FUN=max))
summary(minTotalDistances_3wayCond)
summary(maxTotalDistances_3wayCond)
#The total distances of the matched samples are 0 across all of the simulations, 
#hence the perfect-matchable triplet are identified all of the times

################
# SIMULATION 2 #
################

#Plots of distributions used to sample matching variable
#--------------------------------------------------------
par(mfrow = c(1,4))
curve(dbeta(x,2,2), from = 0, to = 1, ylim = c(0,2.5),
      xlab = "", ylab ="", main = "Group 1: beta(2,2), Group 2: beta(2,2),\n Group 3: beta(2,2)")
curve(dbeta(x,2,2), from = 0, to = 1, add=T, lty = 3, lwd = 2)

curve(dbeta(x,3,2), from = 0, to = 1, ylim = c(0,2.5), lty = 3,lwd=2,
      xlab = "", ylab ="", main = "Group 1: beta(2,3), Group 2: beta(3,2),\n Group 3: beta(2,2)")
curve(dbeta(x,2,2), from = 0, to = 1, add=T, lty = 5)
curve(dbeta(x,2,3), from = 0, to = 1, add=T)

curve(dbeta(x,4,2), from = 0, to = 1, ylim = c(0,2.5),lty = 3,lwd=2,
      xlab = "", ylab ="", main = "Group 1: beta(2,4), Group 2: beta(4,2),\n Group 3: beta(2,2)")
curve(dbeta(x,2,2), from = 0, to = 1, add=T, lty = 5)
curve(dbeta(x,2,4), from = 0, to = 1, add=T)

curve(dbeta(x,5,2), from = 0, to = 1, ylim = c(0,2.5),lty = 3,lwd=2,
      xlab = "", ylab ="", main = "Group 1: beta(2,5), Group 2: beta(5,2),\n Group 3: beta(2,2)")
curve(dbeta(x,2,2), from = 0, to = 1, add=T, lty = 5)
curve(dbeta(x,2,5), from = 0, to = 1, add=T)


#Summary simulations for Table 1 paper (total sample size between 1,500 and 2,500)
#---------------------------------------------------------------------------------

#List of parameters considered in simulations
listParmaters <- list(c(alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 3, alpha2 = 3, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 4, alpha2 = 4, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 5, alpha2 = 5, beta2 = 2, alpha3 = 2, beta3 = 2))

#List of sample sizes in main simulation of the paper
listSampleSizes <- list(c(n1 = 500, n2 = 500, n3 = 500),
                        c(n1 = 1000, n2 = 1000, n3 = 500),
                        c(n1 = 500, n2 = 500, n3 = 1000),
                        c(n1 = 1000, n2 = 500, n3 = 500))

for(parameters in listParmaters) {
  for(sampleSizes in listSampleSizes) {
    
    load(file=paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],"_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],"_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],"_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],"_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    load(file=paste("result simulations/result_NN_alpha1.",parameters[["alpha1"]],"_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],"_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],"_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],"_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
    
    totalDistances_3wayCond <- as.numeric(sapply(sapply(listResults_3wayCond, FUN="[", "overallDistance"), FUN=min))
    totalDistances_NN <- as.numeric(sapply(listResults_NN, FUN="[[", "matchingSchemesDistances"))
    
    print(parameters)
    print(sampleSizes)
    summaryDiff<- sprintf("%.2f",summary((totalDistances_NN - totalDistances_3wayCond)/totalDistances_NN)*100)
    sdDiff<- sprintf("%.2f",sd((totalDistances_NN - totalDistances_3wayCond)/totalDistances_NN)*100)
    
    minMax <- paste(summaryDiff[c(1,6)], collapse ="-")
    meanSd <- paste(summaryDiff[4]," (",sdDiff,") ", sep="")
    medianIQR <- paste(summaryDiff[3]," (",summaryDiff[2],"-",summaryDiff[5],") ", sep="")
    
    print(paste(minMax," & ", medianIQR," & ", meanSd))
    
    print("")
    print("---------------------------------------------------------------------")
  }
}



#Summary simulations for Tables in Supplementary Material
#--------------------------------------------------------

listParmaters <- list(c(alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 3, alpha2 = 3, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 4, alpha2 = 4, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 5, alpha2 = 5, beta2 = 2, alpha3 = 2, beta3 = 2))

listSampleSizes <- list(c(n1 = 100, n2 = 100, n3 = 100),
                        c(n1 = 200, n2 = 200, n3 = 100),
                        c(n1 = 100, n2 = 100, n3 = 200),
                        c(n1 = 200, n2 = 100, n3 = 100))

#List of multipliers for sample sizes to obtain the results of simulations presented in the Appendix
vectorMultipliers <- c(1,2,3,4,10)

for(multiplier in vectorMultipliers) {
  print(paste0("multiplier: ",  multiplier))
  
  for(parameters in listParmaters) {
    for(sampleSizes in listSampleSizes) {
      
      sampleSizes <- sampleSizes*multiplier
      load(file=paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],"_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],"_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],"_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],"_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
      load(file=paste("result simulations/result_NN_alpha1.",parameters[["alpha1"]],"_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],"_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],"_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],"_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
      
      totalDistances_3wayCond <- as.numeric(sapply(sapply(listResults_3wayCond, FUN="[", "overallDistance"), FUN=min))
      totalDistances_NN <- as.numeric(sapply(listResults_NN, FUN="[[", "matchingSchemesDistances"))
      
      print(parameters)
      print(sampleSizes)
      summaryDiff<- sprintf("%.2f",summary((totalDistances_NN - totalDistances_3wayCond)/totalDistances_NN)*100)
      sdDiff<- sprintf("%.2f",sd((totalDistances_NN - totalDistances_3wayCond)/totalDistances_NN)*100)
      
      minMax <- paste(summaryDiff[c(1,6)], collapse ="-")
      meanSd <- paste(summaryDiff[4]," (",sdDiff,") ", sep="")
      medianIQR <- paste(summaryDiff[3]," (",summaryDiff[2],"-",summaryDiff[5],") ", sep="")
      
      print(paste(minMax," & ", medianIQR," & ", meanSd))
      
      print("")
      print("----    ----    ----   ----    ----    ----    ----    ----")
    }
  }
  print("--------------------------------------------------------")
}


#Check: in one scenario, NN attained a smaller total distance than our matching algorithm
# in a few of the 1,000 simulations. If our algorithm is applied on the result of NN, 
# can it improve the total distance attained by the NN?
#-----------------------------------------------------------------------

#Scenario where NN attained smaller total distance:
# - parameters 1: alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2
# - sample size 2: n1 = 2000, n2 = 2000, n3 = 1000
parameters <- listParmaters[[1]]
sampleSizes <- listSampleSizes[[2]]*10

#Load result of simulations.
load(file=paste("result simulations/result_3wayCond_postNN_alpha1.",parameters[["alpha1"]],"_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],"_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],"_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],"_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
load(file=paste("result simulations/result_NN_alpha1.",parameters[["alpha1"]],"_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],"_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],"_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],"_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
load(file=paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],"_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],"_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],"_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],"_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))

#NN did better than 3way cond in: 41 simulations (out of 100)
indexNNbetter <- (!sapply(listResults_3wayCond_postNN, FUN = is.null))
sum(indexNNbetter) #41

#Verify this:
distNN <- sapply(listResults_NN, FUN = "[[", "matchingSchemesDistances")
dist3way <- sapply(listResults_3wayCond, FUN = function(x) {min(x[,"overallDistance"])})
sum(distNN<dist3way) #OK, 41

#Check results of conditional 3-way after NN
dist3way_postNN_sel <- sapply(listResults_3wayCond_postNN[indexNNbetter], FUN = function(x) {min(x[,"overallDistance"])})
dist3way_sel <- sapply(listResults_3wayCond[indexNNbetter], FUN = function(x) {min(x[,"overallDistance"])})
distNN_sel <- sapply(listResults_NN[indexNNbetter], FUN = "[[", "matchingSchemesDistances")

#3way post NN improved the distance 100% of the times:
sum(dist3way_postNN_sel < distNN_sel) #Improvement in all 41

#How much is the improvement over NN?
summaryDiff<- sprintf("%.2f",summary((distNN_sel - dist3way_postNN_sel)/distNN_sel)*100)
sdDiff<- sprintf("%.2f",sd((distNN_sel - dist3way_postNN_sel)/distNN_sel)*100)

#compute and print summary statistics of improvement:
minMax <- paste(summaryDiff[c(1,6)], collapse ="-")
meanSd <- paste(summaryDiff[4]," (",sdDiff,") ", sep="")
medianIQR <- paste(summaryDiff[3]," (",summaryDiff[2],"-",summaryDiff[5],") ", sep="")
print(paste(minMax," & ", medianIQR," & ", meanSd)) 


#Verify number of iterations
#----------------------------

listParmaters <- list(c(alpha1 = 2, beta1 = 2, alpha2 = 2, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 3, alpha2 = 3, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 4, alpha2 = 4, beta2 = 2, alpha3 = 2, beta3 = 2),
                      c(alpha1 = 2, beta1 = 5, alpha2 = 5, beta2 = 2, alpha3 = 2, beta3 = 2))

listSampleSizes <- list(c(n1 = 100, n2 = 100, n3 = 100),
                        c(n1 = 200, n2 = 200, n3 = 100),
                        c(n1 = 100, n2 = 100, n3 = 200),
                        c(n1 = 200, n2 = 100, n3 = 100))

#List of multipliers for sample sizes to obtain the results of simulations presented in the Appendix
vectorMultipliers <- c(1,2,3,4,5,10)

summaryIterations <- data.frame(smallestGroup=100*vectorMultipliers,
                                meanIter = NA,
                                sdIter = NA,
                                medianIter = NA,
                                q1Iter = NA,
                                q3Iter = NA,
                                minIter = NA,
                                maxIter = NA)


for(multiplier in vectorMultipliers) {
  
  vectorIter_12 <- NULL
  vectorIter_23 <- NULL
  vectorIter_13 <- NULL
  
  for(parameters in listParmaters) {
    
    for(sampleSizes in listSampleSizes) {
      
        sampleSizes <- sampleSizes*multiplier
        
        load(file=paste("result simulations/result_3wayCond_alpha1.",parameters[["alpha1"]],"_beta1.",parameters[["beta1"]],"_alpha2.",parameters[["alpha2"]],"_beta2.",parameters[["beta2"]],"_alpha3.",parameters[["alpha3"]],"_beta3.",parameters[["beta3"]],"_n1.",sampleSizes[["n1"]],"_n2.",sampleSizes[["n2"]],"_n3.",sampleSizes[["n3"]],".Rdata",sep=""))
        
        vectorIter_12_temp <- sapply(listResults_3wayCond, function(x) {x$iterations[x$scheme=="1-2"]})
        vectorIter_23_temp <- sapply(listResults_3wayCond, function(x) {x$iterations[x$scheme=="2-3"]})
        vectorIter_13_temp <- sapply(listResults_3wayCond, function(x) {x$iterations[x$scheme=="1-3"]})
        
        vectorIter_12 <- c(vectorIter_12, vectorIter_12_temp)
        vectorIter_23 <- c(vectorIter_23, vectorIter_23_temp)
        vectorIter_13 <- c(vectorIter_13, vectorIter_13_temp)
        
      }
      
  }
  
  vectorIterAll <- c(vectorIter_12,vectorIter_23,vectorIter_13)
  smallestGroup <- min(sampleSizes)
  
  summaryIterations$meanIter[summaryIterations$smallestGroup==smallestGroup] <- mean(vectorIterAll)
  summaryIterations$medianIter[summaryIterations$smallestGroup==smallestGroup] <- median(vectorIterAll)
  summaryIterations$sdIter[summaryIterations$smallestGroup==smallestGroup] <- sd(vectorIterAll)
  summaryIterations$q1Iter[summaryIterations$smallestGroup==smallestGroup] <- quantile(vectorIterAll, probs = .25)
  summaryIterations$q3Iter[summaryIterations$smallestGroup==smallestGroup] <- quantile(vectorIterAll, probs = .75)
  summaryIterations$minIter[summaryIterations$smallestGroup==smallestGroup] <- min(vectorIterAll)
  summaryIterations$maxIter[summaryIterations$smallestGroup==smallestGroup] <- max(vectorIterAll)
}      


summaryIterations