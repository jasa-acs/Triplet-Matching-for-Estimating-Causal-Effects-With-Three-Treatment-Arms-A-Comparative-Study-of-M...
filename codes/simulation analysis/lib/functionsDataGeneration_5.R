# Author: Giovanni Nattino 
# Functions for data generation in simulation analysis
# Reference: Triplet Matching for Estimating Causal Effects with Three Treatment Arms: 
#  A Comparative Study of Mortality by Trauma Center Level
#####################################################################################

generateData_logitPS <- function(n, mu1, mu2, sd1, sd2) {
  
  for(i in 1:2) {
    logitTemp <- rnorm(n,evaluate("mu",i),evaluate("sd",i))
    assign(paste("logit",i,sep=""), logitTemp)
    rm(logitTemp)
  }
  
  prob3 <- 1/(1 + exp(logit1) + exp(logit2))
  prob1 <- exp(logit1)/(1 + exp(logit1) + exp(logit2))
  prob2 <- exp(logit2)/(1 + exp(logit1) + exp(logit2))
  
  treatment <- rep(NA, n)
  for(j in 1:n) {
    treatment[j] <- which(rmultinom(1, size = 1, 
                                    prob = c(prob1[j], prob2[j], prob3[j]))==1)
  }
  
  data <- data.frame(prob1 = prob1, 
                     prob2 = prob2, 
                     prob3 = prob3, 
                     logit1 = logit1,
                     logit2 = logit2,
                     treatment = treatment)
  
  data$logitProb1 <- logit(data$prob1)
  data$logitProb2 <- logit(data$prob2)
  data$logitProb3 <- logit(data$prob3)
  
  return(data)
  
}

generateData_shift <- function(n, shift) {
  
  varMatching1 <- 1:n
  varMatching2 <- varMatching1 + shift
  varMatching3 <- varMatching1 - shift
  
  indexTrueMatches3 <- indexTrueMatches2 <- indexTrueMatches1 <- seq(from = 1, to = n)
  
  treatment1 <- rep(1, n)
  treatment2 <- rep(2, n)
  treatment3 <- rep(3, n)
  
  data <- data.frame(treatment = c(treatment1, treatment2,  treatment3),
                     indexTrueMatches = c(indexTrueMatches1, indexTrueMatches2, indexTrueMatches3),
                     varMatching = c(varMatching1, varMatching2,  varMatching3))
  
  data <- data[sample(nrow(data), nrow(data), replace = F),]
  row.names(data) <- seq(from = 1, to = nrow(data))
  
  
  return(data)
  
}

generateData_H0 <- function(n1, n2, n3, mu, sd) {
  
  varMatching3 <- varMatching2 <- varMatching1 <- rnorm(n1, mu, sd)
  
  indexTrueMatches3 <- indexTrueMatches2 <- indexTrueMatches1 <- seq(from = 1, to = n1)
  
  treatment1 <- rep(1, n1)
  treatment2 <- rep(2, n1)
  treatment3 <- rep(3, n1)
  
  data <- data.frame(treatment = c(treatment1, treatment2,  treatment3),
                     indexTrueMatches = c(indexTrueMatches1, indexTrueMatches2, indexTrueMatches3),
                     varMatching = c(varMatching1, varMatching2,  varMatching3))
  
  if(n2 != n1 | n3 != n1) {
    largerGroups <- which(c(n1,n2,n3) > min(c(n1,n2,n3)))
    nSmallerGroups <- min(c(n1,n2,n3))
    for(groupTemp in largerGroups){
      sizeGroupTemp <- eval(parse(text = paste("n",groupTemp,sep="")))
      dataTemp <- data.frame(treatment = rep(groupTemp, sizeGroupTemp - nSmallerGroups),
                             indexTrueMatches = rep(NA, sizeGroupTemp - nSmallerGroups),
                             varMatching = rnorm(sizeGroupTemp - nSmallerGroups, mu, sd),
                             stringsAsFactors = F)
      data <- rbind(data,dataTemp)
    }
    
  }
  
  data <- data[sample(nrow(data), nrow(data), replace = F),]
  row.names(data) <- seq(from = 1, to = nrow(data))
  
  return(data)
  
}

generateData_H0andShift <- function(n1, n2, n3, mu, sd, nShift, shift) {
  
  varMatching3 <- varMatching2 <- varMatching1 <- rnorm(n1, mu, sd)
  
  indexTrueMatches3 <- indexTrueMatches2 <- indexTrueMatches1 <- seq(from = 1, to = n1)
  
  treatment1 <- rep(1, n1)
  treatment2 <- rep(2, n1)
  treatment3 <- rep(3, n1)
  
  data <- data.frame(treatment = c(treatment1, treatment2,  treatment3),
                     indexTrueMatches = c(indexTrueMatches1, indexTrueMatches2, indexTrueMatches3),
                     varMatching = c(varMatching1, varMatching2,  varMatching3))
  
  
  if(n2 != n1 | n3 != n1) {
    largerGroups <- which(c(n1,n2,n3) > min(c(n1,n2,n3)))
    nSmallerGroups <- min(c(n1,n2,n3))
    for(groupTemp in largerGroups){
      sizeGroupTemp <- eval(parse(text = paste("n",groupTemp,sep="")))
      dataTemp <- data.frame(treatment = rep(groupTemp, sizeGroupTemp - nSmallerGroups),
                             indexTrueMatches = rep(NA, sizeGroupTemp - nSmallerGroups),
                             varMatching = rnorm(sizeGroupTemp - nSmallerGroups, mu, sd),
                             stringsAsFactors = F)
      data <- rbind(data,dataTemp)
    }
    
  }
  
  varMatching1shift <- seq(from = -2, to = 2, length = nShift)
  varMatching2shift <- varMatching1shift + shift*4/nShift
  varMatching3shift <- varMatching1shift - shift*4/nShift
  
  indexTrueMatches1shift <- seq(from = (nrow(data)+1), to = (nrow(data)+nShift))
  indexTrueMatches3shift <- indexTrueMatches2shift <- indexTrueMatches1shift
  
  treatment1shift <- rep(1, nShift)
  treatment2shift <- rep(2, nShift)
  treatment3shift <- rep(3, nShift)
  
  dataShift <- data.frame(treatment = c(treatment1shift,treatment2shift,treatment3shift),
                          indexTrueMatches = c(indexTrueMatches1shift,indexTrueMatches2shift,indexTrueMatches3shift),
                          varMatching = c(varMatching1shift,varMatching2shift,varMatching3shift),
                          stringsAsFactors = F)
  data <- rbind(data,dataShift)
  
  data <- data[sample(nrow(data), nrow(data), replace = F),]
  row.names(data) <- seq(from = 1, to = nrow(data))
  
  return(data)
  
}

generateData_Peaks <- function(nMid, nPeaks) {
  
  if(nMid < nPeaks) {
    varMatching2 <- rbeta(nPeaks, 5, 2)
    varMatching3 <- rbeta(nPeaks, 2, 5)
    pooledObs <-c(varMatching2,varMatching3)
    varMatching1 <- sample(pooledObs, nMid, prob = (1 - (pooledObs-0.5)^2))
  } 
  
  if(nMid >= nPeaks) {
    varMatching1 <- rbeta(nMid, 2, 2)
    varMatching2 <- sample(varMatching1, nPeaks, prob = varMatching1^2)
    varMatching3 <- sample(varMatching1, nPeaks, prob = (1 - varMatching1^2))
  }
  
  indexTrueMatches3 <- rep(NA,nPeaks)
  indexTrueMatches2 <- rep(NA,nPeaks)
  indexTrueMatches1 <- rep(NA,nMid)
  
  treatment1 <- rep(1, nMid)
  treatment2 <- rep(2, nPeaks)
  treatment3 <- rep(3, nPeaks)
  
  data <- data.frame(treatment = c(treatment1, treatment2,  treatment3),
                     indexTrueMatches = c(indexTrueMatches1, indexTrueMatches2, indexTrueMatches3),
                     varMatching = c(varMatching1, varMatching2,  varMatching3))
  
  data <- data[sample(nrow(data), nrow(data), replace = F),]
  row.names(data) <- seq(from = 1, to = nrow(data))
  
  return(data)
  
}

generateData_neighbor <- function(nGrid, epsilon, n1byInt, n2byInt, n3byInt) {
  
  grid <- 1:nGrid
  
  noise1 <- runif(n1byInt * nGrid, min = -epsilon, max = epsilon)
  noise2 <- runif(n2byInt * nGrid, min = -epsilon, max = epsilon)
  noise3 <- runif(n3byInt * nGrid, min = -epsilon, max = epsilon)
  
  noise1Mat <- matrix(noise1, ncol = n1byInt)
  noise2Mat <- matrix(noise2, ncol = n2byInt)
  noise3Mat <- matrix(noise3, ncol = n3byInt)
  
  varMatching1 <- as.vector(sapply(1:n1byInt, function(x){return(grid)}) + noise1Mat)
  varMatching2 <- as.vector(sapply(1:n2byInt, function(x){return(grid)}) + noise2Mat)
  varMatching3 <- as.vector(sapply(1:n3byInt, function(x){return(grid)}) + noise3Mat)
  
  indexTrueMatches1 <- rep(NA,length(varMatching1))
  indexTrueMatches2 <- rep(NA,length(varMatching2))
  indexTrueMatches3 <- rep(NA,length(varMatching3))
  
  treatment1 <- rep(1, length(varMatching1))
  treatment2 <- rep(2, length(varMatching2))
  treatment3 <- rep(3, length(varMatching3))
  
  data <- data.frame(treatment = c(treatment1, treatment2,  treatment3),
                     indexTrueMatches = c(indexTrueMatches1, indexTrueMatches2, indexTrueMatches3),
                     varMatching = c(varMatching1, varMatching2,  varMatching3))
  
  data <- data[sample(nrow(data), nrow(data), replace = F),]
  row.names(data) <- seq(from = 1, to = nrow(data))
  
  return(data)
  
}

generateDataNCH <- function(n, alfa1, alfa2, OR1_13, OR1_23, 
                            OR21_13, OR21_23, OR22_13, OR22_23, 
                            OR3_13, OR3_23) {
  
  #X1 Normal(0,1)
  x1 <- rnorm(n,mean = 0, sd = 1)
  #X2 Multinomial(0.2,0.55,0.25)
  x2 <- t(rmultinom(n, size = 1, prob = c(0.20, 0.55, 0.25)))
  x2_1 <-x2[,1]
  x2_2 <-x2[,2]
  #X3 Bernoulli(0.55)
  x3 <- rbinom(n, size = 1, prob = 0.55)
  
  logit_13 <- alfa1+log(OR1_13)*x1+log(OR21_13)*x2_1+log(OR22_13)*x2_2+log(OR3_13)*x3;
  logit_23 <- alfa2+log(OR1_23)*x1+log(OR21_23)*x2_1+log(OR22_23)*x2_2+log(OR3_23)*x3;
  
  #calculate p1/p3 and p2/p3 by taking exponential of the logits;
  p1_p3 <- exp(logit_13);
  p2_p3 <- exp(logit_23);
  
  #calculate p1,p2,and p3 by deviding p1/p3 and p2/3 with (p1/p3)+(p2/p3)+1 and assuming p1+p2+p3=1;
  prob1 <- p1_p3/(p1_p3+p2_p3+1);
  prob2 <- p2_p3/(p1_p3+p2_p3+1);
  prob3 <- 1-prob1-prob2;
  
  treatment <- rep(NA, n)
  for(j in 1:n) {
    treatment[j] <- which(rmultinom(1, size = 1, 
                                    prob = c(prob1[j], prob2[j], prob3[j]))==1)
  }
  
  data <- data.frame(prob1 = prob1, 
                     prob2 = prob2, 
                     prob3 = prob3, 
                     treatment = treatment)
  return(data)
}

generateData_Betas <- function(n1, n2, n3, alpha1, beta1, alpha2, beta2, alpha3, beta3) {
  
  varMatching1 <- rbeta(n1, alpha1, beta1)
  varMatching2 <- rbeta(n2, alpha2, beta2)
  varMatching3 <- rbeta(n3, alpha3, beta3)
  
  indexTrueMatches1 <- rep(NA,n1)
  indexTrueMatches2 <- rep(NA,n2)
  indexTrueMatches3 <- rep(NA,n3)
  
  treatment1 <- rep(1, n1)
  treatment2 <- rep(2, n2)
  treatment3 <- rep(3, n3)
  
  data <- data.frame(treatment = c(treatment1, treatment2,  treatment3),
                     indexTrueMatches = c(indexTrueMatches1, indexTrueMatches2, indexTrueMatches3),
                     varMatching = c(varMatching1, varMatching2,  varMatching3))
  
  data <- data[sample(nrow(data), nrow(data), replace = F),]
  row.names(data) <- seq(from = 1, to = nrow(data))
  
  return(data)
  
}
