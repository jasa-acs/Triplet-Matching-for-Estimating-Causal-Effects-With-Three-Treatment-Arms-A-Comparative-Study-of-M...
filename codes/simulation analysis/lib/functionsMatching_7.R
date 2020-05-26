# Author: Giovanni Nattino 
# Functions implementing 3-way conditional matching 
# Reference: Triplet Matching for Estimating Causal Effects with Three Treatment Arms: 
#  A Comparative Study of Mortality by Trauma Center Level
#####################################################################################

evaluate <- function(text,j) {
  return(eval(parse(text = paste(text,j,sep=""))))
}

logistic<- function(x) {1/(1+exp(-x))}
logit<- function(x) {log(x/(1-x))}


euclideanDistance <- function(A,B) {
  return(sqrt(sum((A-B)^2)))
}

perimeterTriangle <- function(dataTriangle) {
  dataTriangle <- as.matrix(dataTriangle)
  output <- (euclideanDistance(dataTriangle[1,],dataTriangle[2,]) +
             euclideanDistance(dataTriangle[1,],dataTriangle[3,]) +
             euclideanDistance(dataTriangle[2,],dataTriangle[3,]))
  return(output)
}

evaluateMatching <- function(data, stringMatchedIndex, variablesMatch) {
  
  withinTripletDistances <- by(data[,variablesMatch], 
                               INDICES = data[,stringMatchedIndex], 
                               FUN = perimeterTriangle)

  dfTempToExport <- data.frame(indexMatch = names(withinTripletDistances),
                               distance = as.vector(withinTripletDistances),
                               stringsAsFactors = F)
  
  names(dfTempToExport)[names(dfTempToExport) == "indexMatch"] <- stringMatchedIndex
  
  return(list(overallMatchedDistance = sum(withinTripletDistances),
              distanceByTriplet = dfTempToExport))
}

smallestValueNon0 <- function(x){
  temp <- diff(sort(x))
  return(min(temp[temp>0]))
}

constrained3wayMatching <- function(data, stringMatchesStep1, treatStep1, otherTreat,
                                    variablesMatch, variableTreatment) {
  
  #Functions for constrained 3way:
  #-------------------------------
  applyPersonalDistance <- function(index, data, z) { 
    indexDf <- as.data.frame(index, stringsAsFactors = F)
    groupTreated <- unique(data[z,variableTreatment])
    groupControls <- unique(data[!z,variableTreatment])
    names(indexDf) <- paste("group",c(groupTreated,groupControls), sep ="")
    
    longIndexStep1 <- dataAll[!is.na(dataAll[,stringMatchesStep1]),c(variableTreatment,stringMatchesStep1)]
    longIndexStep1$value <- row.names(longIndexStep1) 
    wideIndexStep1 <- reshape(longIndexStep1, direction = "wide", idvar = c(stringMatchesStep1), 
                              timevar = variableTreatment)
    wideIndexStep1[,stringMatchesStep1]<- NULL
    names(wideIndexStep1) <- gsub("value.","group",names(wideIndexStep1))
    
    matchingVector <- match(indexDf[,paste("group",treatStep1[1],sep="")],wideIndexStep1[,paste("group",treatStep1[1],sep="")])
    indexDf[,paste("group",treatStep1[2],sep="")] <- wideIndexStep1[matchingVector,paste("group",treatStep1[2],sep="")]
    
    #Merge with logit prob for each group
    varsByGroup <- paste("group",outer(c(treatStep1,otherTreat), variablesMatch, FUN = paste, sep="_"),sep="")
    
    indexDf[,varsByGroup] <- NA
    indexDf[,paste("group",treatStep1[1],"_",variablesMatch,sep="")] <- dataAll[indexDf[,paste("group",treatStep1[1],sep="")],variablesMatch]
    indexDf[,paste("group",treatStep1[2],"_",variablesMatch,sep="")] <- dataAll[indexDf[,paste("group",treatStep1[2],sep="")],variablesMatch]
    indexDf[,paste("group",otherTreat,"_",variablesMatch,sep="")] <- dataAll[indexDf[,paste("group",otherTreat,sep="")],variablesMatch]
     
    distances <- (sqrt(rowSums(as.matrix((indexDf[,paste("group",treatStep1[1],"_",variablesMatch,sep="")] - indexDf[,paste("group",treatStep1[2],"_",variablesMatch,sep="")])^2))) +
                  sqrt(rowSums(as.matrix((indexDf[,paste("group",treatStep1[1],"_",variablesMatch,sep="")] - indexDf[,paste("group",otherTreat,"_",variablesMatch,sep="")])^2))) +
                  sqrt(rowSums(as.matrix((indexDf[,paste("group",treatStep1[2],"_",variablesMatch,sep="")] - indexDf[,paste("group",otherTreat,"_",variablesMatch,sep="")])^2))))
    
    return(distances)
  }
  
  #---------------------------------------------------------------------------------------------
  
  #In the selection, selecting the treatStep1[1] or treatStep1[2] is the same.
  #However, the same choice made here must be reported also into the function
  #'applyPersonalDistance'
  selectionSecondStep <- (data[,variableTreatment] %in% otherTreat |
                            (data[,variableTreatment] %in% treatStep1[1] & 
                               !is.na(data[,stringMatchesStep1])))
  
  stringTreatStep1 <- paste("treatment",paste(sort(treatStep1), collapse = ""),sep ="")
  stringTreatStep2 <- paste(stringTreatStep1,"_", otherTreat, sep = "")
  
  data[,stringTreatStep2] <- NA
  data[data[,variableTreatment] %in% treatStep1,stringTreatStep2] <- paste(sort(treatStep1), collapse = "")
  data[data[,variableTreatment] %in% otherTreat,stringTreatStep2] <- otherTreat
  data[,stringTreatStep2] <- factor(data[,stringTreatStep2], 
                                   levels = c(paste(sort(treatStep1), collapse = ""),
                                              otherTreat))
  data[,stringTreatStep2] <- relevel(data[,stringTreatStep2], 
                                    names(table(data[selectionSecondStep,stringTreatStep2]))[which.max(table(data[selectionSecondStep,stringTreatStep2]))])
  #Make treatment variable binary. Package optmatch deprecated factor type for treatment variables.
  data[,stringTreatStep2] <- (data[,stringTreatStep2] == (levels(data[,stringTreatStep2])[2]))*1
  
  #Global variables to be used in the function 'applyPersonalDistance'
  dataAll <-  data
  distanceStep2 <- match_on(applyPersonalDistance,
                            z = data[selectionSecondStep,stringTreatStep2],
                            data = data[selectionSecondStep,])
  
  matchSecondStep <- pairmatch(distanceStep2,
                               controls = 1,
                               data = data[selectionSecondStep,]) 
  
  data$matchesStep2 <- NA
  data$matchesStep2[selectionSecondStep] <- matchSecondStep
  
  matchBothSides <- !is.na(data$matchesStep2) & !is.na(data[,stringMatchesStep1])
  
  # Match groups 1-2-3
  #-----------------
  
  data$matches123 <- NA
  matchedCounter <- 0
  for (i in which(matchBothSides)) {
    matchedCounter <- matchedCounter + 1
    data$matches123[data[,stringMatchesStep1] %in% data[i,stringMatchesStep1] | 
                          data$matchesStep2 %in% data$matchesStep2[i]] <- matchedCounter    
  }
  
  restultEvaluation <- evaluateMatching(data, "matches123", variablesMatch)
  return(list(data=data, distances = restultEvaluation$distanceByTriplet$distance))
}



applyModifiedOptMatching <- function(data, startingEdges, variablesMatch, variableTreatment) {

  #Initialize dataframe to store results for starting edges
  results <- data.frame(scheme = startingEdges,
                        iterations = NA,
                        sizeMatchedSample = NA,
                        overallDistance = NA)

  #Standardize the variable to match on
  #-------------------------------------
  # If the values were too small, optimal matching
  # was not distinguishng among real optimal matching.
  variablesMatchStndzd <- variablesMatch

  #ID to be used to merge the results to original dataset
  data$idTempForMatching <- 1:nrow(data)
  
  #Keep track of the best scheme, store matched data only if improvement
  #in overall distance
  minOverallDistanceEdges <- Inf

  for (j in 1:length(startingEdges)){
    
    startingEdge <- startingEdges[j]
    treatStep1 <- unlist(strsplit(startingEdge,split="-"))
    
    dataStep <- data[,c(variablesMatch,variablesMatchStndzd,
                        variableTreatment, "idTempForMatching")]
    
    
    # Step 1
    #-----------------
    
    #Definition of selection
    selectionFirstStep <- dataStep[,variableTreatment] %in% treatStep1
    
    #Definition of matching criteria
    stringTreatStep1 <- paste("treatment",paste(sort(treatStep1), collapse = ""),sep ="")
    stringStep1 <- paste(stringTreatStep1, " ~ ", paste(variablesMatchStndzd, collapse = " + "), sep ="")
    
    #Definition of treatment variable (the first level must be the one with more observations)
    dataStep[,stringTreatStep1] <- factor(dataStep[,variableTreatment], levels = treatStep1)
    dataStep[,stringTreatStep1] <- relevel(dataStep[,stringTreatStep1], 
                                           names(table(dataStep[selectionFirstStep,stringTreatStep1]))[which.max(table(dataStep[selectionFirstStep,stringTreatStep1]))])
    #Make treatment variable binary. Package optmatch deprecated factor type for treatment variables.
    dataStep[,stringTreatStep1] <- (dataStep[,stringTreatStep1] == (levels(dataStep[,stringTreatStep1])[2]))*1
      
    matchFirstStep <- pairmatch(formula(stringStep1), #
                                controls = 1,
                                method = "euclidean",
                                data = dataStep[selectionFirstStep,]) 
    
    dataStep$matchesStep1 <- NA
    dataStep$matchesStep1[selectionFirstStep] <- matchFirstStep

    # Step 2
    #-----------------

    otherTreat <- setdiff(names(table(dataStep[,variableTreatment])),treatStep1)
    outputStep2_1 <- constrained3wayMatching(data = dataStep, stringMatchesStep1 = "matchesStep1", 
                                      treatStep1 = treatStep1, otherTreat = otherTreat,
                                      variablesMatch = variablesMatchStndzd, 
                                      variableTreatment = variableTreatment) 
    distStep2_1 <- sum(outputStep2_1$distances)
    dataStep2_1 <- outputStep2_1$data
    
    for(iter in 1:50) {
    
      #Scheme 2
      treatStep1_2 <- sort(c(treatStep1[1],otherTreat))
      otherTreat_2 <- setdiff(names(table(dataStep[,variableTreatment])),treatStep1_2)
      
      dataStep2_1$matchesStep1_2 <- dataStep2_1$matches123
      dataStep2_1$matchesStep1_2[dataStep2_1[,variableTreatment] %in% otherTreat_2] <- NA
      outputStep2_2 <- constrained3wayMatching(data = dataStep2_1, stringMatchesStep1 = "matchesStep1_2", 
                                        treatStep1 = treatStep1_2, otherTreat_2,
                                        variablesMatch = variablesMatchStndzd, 
                                        variableTreatment = variableTreatment) 
      distStep2_2 <- sum(outputStep2_2$distances)
      
      #Scheme 3
      treatStep1_3 <- sort(c(treatStep1[2],otherTreat))
      otherTreat_3 <- setdiff(names(table(dataStep[,variableTreatment])),treatStep1_3)
      
      dataStep2_1$matchesStep1_3 <- dataStep2_1$matches123
      dataStep2_1$matchesStep1_3[dataStep2_1[,variableTreatment] %in% otherTreat_3] <- NA
      
      outputStep2_3 <- constrained3wayMatching(data = dataStep2_1, stringMatchesStep1 = "matchesStep1_3", 
                                               treatStep1 = treatStep1_3, otherTreat = otherTreat_3,
                                               variablesMatch = variablesMatchStndzd, 
                                               variableTreatment = variableTreatment) 
      distStep2_3 <- sum(outputStep2_3$distances)
      
      if( min(distStep2_2,distStep2_3) < distStep2_1) {
        
        labelBest <- c(2,3)[which.min(c(distStep2_2,distStep2_3))]
        output <- eval(parse(text=paste("outputStep2_",labelBest,sep="")))
        dataStep2_1 <- output$data
        distStep2_1 <- sum(output$distances)

        treatStep1 <- eval(parse(text=paste("treatStep1_",labelBest,sep="")))
        otherTreat <- eval(parse(text=paste("otherTreat_",labelBest,sep="")))
        
      } else {
        
        break
        
      }
    
    }

    # Evaluate distance of matched triplets
    #--------------------------------------
    matchingEvaluation <- evaluateMatching(data = dataStep2_1,
                                           stringMatchedIndex = "matches123", 
                                           variablesMatch = variablesMatch)
    
    results$sizeMatchedSample[j] <- nrow(matchingEvaluation[["distanceByTriplet"]])
    results$overallDistance[j] <- matchingEvaluation$overallMatchedDistance
    results$iterations[j] <- iter
    
    if(matchingEvaluation$overallMatchedDistance < minOverallDistanceEdges){
      
      dataToMerge <- matchingEvaluation[["distanceByTriplet"]]
      dataToMerge[,variableTreatment] <-  levels(factor(data[,variableTreatment]))[1]
      resultFirstMerge <- merge(dataStep2_1, dataToMerge, 
                                all.x = T, all.y = T, 
                                by = c("matches123",variableTreatment))
      resultFirstMerge$indexMatch <- resultFirstMerge$matches123
      
      resultSecondMerge <- merge(data, 
                                 resultFirstMerge[,c("idTempForMatching","indexMatch","distance")],
                                 by = "idTempForMatching") 
      resultSecondMerge <- resultSecondMerge[,! (names(resultSecondMerge) %in% c("idTempForMatching", variablesMatchStndzd))]
      
      #New best matching scheme:
      minOverallDistanceEdges <- matchingEvaluation$overallMatchedDistance
      
    }
    
  }
  
  return(list(results= results,
              matchedData = resultSecondMerge))
}


#####################
# Nearest Neighbors #
#####################


applyNearNeigh <- function(data, variablesMatch, variableTreatment) {
  permutation <- sample(nrow(data), replace = F)
  data <- data[permutation,]
  
  #Standardize the variable to match on. If the values were too small, optimal matching
  #was not distinguishng among real optimal matching.
  variablesMatchStndzd <- paste(variablesMatch,"Stndzd",sep="")
  
  for(k in 1:length(variablesMatch)) {
    meanK <- mean(data[,variablesMatch[k]])
    sdK <- sd(data[,variablesMatch[k]])
    data[,variablesMatchStndzd[k]] <- (data[,variablesMatch[k]] - meanK)/sdK*100
  }
  
  orderedTable <- table(data[,variableTreatment])[order(table(data[,variableTreatment]))]
  smallestGroup <- names(orderedTable)[1]
  largestGroup <- names(orderedTable)[3]
  
  data$picked <- 0
  data$index <- 1:nrow(data)
  data$matches123 <- NA
  
  dataS <- data[data[,variableTreatment] %in% smallestGroup, ]
  dataL <- data[data[,variableTreatment] %in% largestGroup, ]
  dataM <- data[! (data[,variableTreatment] %in% c(largestGroup,smallestGroup)), ]
  
  matchedCounter <- 0
  overallMatchedDistance <- 0
  
  for(iS in dataS$index[dataS$picked %in% 0]) {
    
    matchedCounter <- matchedCounter + 1
    
    availableL <- dataL$picked %in% 0
    availableM <- dataM$picked %in% 0
    
    numAvailableL <- sum(availableL)
    numAvailableM <- sum(availableM)

    cumulativeDistTot <- rep(0,numAvailableL*numAvailableM)

    indexesL <- rep(dataL$index[availableL],rep(numAvailableM,numAvailableL))
    indexesM <- rep(dataM$index[availableM],numAvailableL)
    
    A <- NULL
    B <- NULL
    C <- NULL
    for(variable in variablesMatchStndzd) {
      a <- rep(dataL[availableL,variable],rep(numAvailableM,numAvailableL))
      b <- rep(dataM[availableM,variable],numAvailableL)
      c <- rep(dataS[dataS$index %in% iS, variable], length(a))
      
      A <- cbind(A,a)
      B <- cbind(B,b)
      C <- cbind(C,c)
    }
    
    
    cumulativeDistTot <- sqrt(rowSums((A-B)^2)) + sqrt(rowSums((A-C)^2)) + sqrt(rowSums((B-C)^2))
    
    indexBest <- which.min(cumulativeDistTot)
    
    indexLbest <- indexesL[indexBest]
    indexMbest <- indexesM[indexBest]

    mdBest <- cumulativeDistTot[indexBest]
    
    dataS$picked[dataS$index %in% iS] <- 1
    dataL$picked[dataL$index %in% indexLbest] <- 1
    dataM$picked[dataM$index %in% indexMbest] <- 1
    
    data$matches123[data$index %in% iS] <- matchedCounter
    data$matches123[data$index %in% indexLbest] <- matchedCounter
    data$matches123[data$index %in% indexMbest] <- matchedCounter
    
    
    overallMatchedDistance <- overallMatchedDistance + mdBest 
    
  }
  
  matchingEvaluation <- evaluateMatching(data = data, 
                                         stringMatchedIndex = "matches123",
                                         variablesMatch = variablesMatch)
  
  
  results <- data.frame(sizeMatchedSample = nrow(matchingEvaluation[["distanceByTriplet"]]),
                        matchingSchemesDistances = matchingEvaluation$overallMatchedDistance)
  
  dataToMerge <- matchingEvaluation[["distanceByTriplet"]]
  dataToMerge[,variableTreatment] <-  levels(factor(data[,variableTreatment]))[1]
  resultMerge <- merge(data, dataToMerge, 
                       all.x = T, all.y = T, 
                       by = c("matches123",variableTreatment))
  
  resultMerge$indexMatch <- resultMerge$matches123
  resultMerge <- resultMerge[,! (names(resultMerge) %in% c("index", variablesMatchStndzd, "picked", "matches123"))]
  
  return(list(results=results, permutation=permutation,
              matchedData = resultMerge))
  
}


############################################################
# 3-way optimal matching starting from any set of triplets #
############################################################

applyModifiedOptMatching_startMatch <- function(data, startingEdges, variablesMatch, 
                                                variableTreatment, varInitialMatch) {
  
  #Initialize dataframe to store results for starting edges
  results <- data.frame(scheme = startingEdges,
                        iterations = NA,
                        sizeMatchedSample = NA,
                        overallDistance = NA)
  
  
  #ID to be used to merge the results to original dataset
  data$idTempForMatching <- 1:nrow(data)
  
  #Standardize the variable to match on. If the values were too small, optimal matching
  #was not distinguishng among real optimal matching.
  variablesMatchStndzd <- variablesMatch
  
  matchingEvaluation <- evaluateMatching(data = data, 
                                         stringMatchedIndex = varInitialMatch,
                                         variablesMatch = variablesMatch)
  
  initialOverallMatchedDistance <- matchingEvaluation$overallMatchedDistance
  
  minOverallDistanceEdges <- Inf
  
  for (j in 1:length(startingEdges)){
    
    dataStep <- data
    startingEdge <- startingEdges[j]
    treatStep1 <- unlist(strsplit(startingEdge,split="-"))
    otherTreat <- setdiff(names(table(dataStep[,variableTreatment])),treatStep1)
    
    distStep2_1 <- initialOverallMatchedDistance
    dataStep2_1 <- dataStep
    dataStep2_1$matches123 <- dataStep2_1[,varInitialMatch]
    
    for(iter in 1:50) {
      
      #Scheme 2
      treatStep1_2 <- sort(c(treatStep1[1],otherTreat))
      otherTreat_2 <- setdiff(names(table(dataStep[,variableTreatment])),treatStep1_2)
      
      dataStep2_1$matchesStep1_2 <- dataStep2_1$matches123
      dataStep2_1$matchesStep1_2[dataStep2_1[,variableTreatment] %in% otherTreat_2] <- NA
      outputStep2_2 <- constrained3wayMatching(data = dataStep2_1, stringMatchesStep1 = "matchesStep1_2", 
                                               treatStep1 = treatStep1_2, otherTreat_2,
                                               variablesMatch = variablesMatchStndzd,
                                               variableTreatment = variableTreatment) 
      distStep2_2 <- sum(outputStep2_2$distances)
      
      #Scheme 3
      treatStep1_3 <- sort(c(treatStep1[2],otherTreat))
      otherTreat_3 <- setdiff(names(table(dataStep[,variableTreatment])),treatStep1_3)
      
      dataStep2_1$matchesStep1_3 <- dataStep2_1$matches123
      dataStep2_1$matchesStep1_3[dataStep2_1[,variableTreatment] %in% otherTreat_3] <- NA
      
      outputStep2_3 <- constrained3wayMatching(data = dataStep2_1, stringMatchesStep1 = "matchesStep1_3", 
                                               treatStep1 = treatStep1_3, otherTreat = otherTreat_3,
                                               variablesMatch = variablesMatchStndzd,
                                               variableTreatment = variableTreatment) 
      distStep2_3 <- sum(outputStep2_3$distances)
      
      if( min(distStep2_2,distStep2_3) < distStep2_1) {
        
        labelBest <- c(2,3)[which.min(c(distStep2_2,distStep2_3))]
        output <- eval(parse(text=paste("outputStep2_",labelBest,sep="")))
        dataStep2_1 <- output$data
        distStep2_1 <- sum(output$distances)
        
        treatStep1 <- eval(parse(text=paste("treatStep1_",labelBest,sep="")))
        otherTreat <- eval(parse(text=paste("otherTreat_",labelBest,sep="")))
        
      } else {
        
        break
        
      }
      
    }
    
    # Evaluate distance of matched triplets
    #--------------------------------------
    matchingEvaluation <- evaluateMatching(data = dataStep2_1,
                                           stringMatchedIndex = "matches123", 
                                           variablesMatch = variablesMatch)
    
    results$sizeMatchedSample[j] <- nrow(matchingEvaluation[["distanceByTriplet"]])
    results$overallDistance[j] <- matchingEvaluation$overallMatchedDistance
    results$iterations[j] <- iter
    
    if(matchingEvaluation$overallMatchedDistance < minOverallDistanceEdges){
      
      dataToMerge <- matchingEvaluation[["distanceByTriplet"]]
      dataToMerge[,variableTreatment] <-  levels(factor(data[,variableTreatment]))[1]
      resultFirstMerge <- merge(dataStep2_1, dataToMerge, 
                                all.x = T, all.y = T, 
                                by = c("matches123",variableTreatment))
      resultFirstMerge$indexMatch <- resultFirstMerge$matches123
      
      resultSecondMerge <- merge(data, 
                                 resultFirstMerge[,c("idTempForMatching","indexMatch","distance")],
                                 by = "idTempForMatching") 
      resultSecondMerge <- resultSecondMerge[,! (names(resultSecondMerge) %in% c("idTempForMatching", variablesMatchStndzd))]
      
      #New best matching scheme:
      minOverallDistanceEdges <- matchingEvaluation$overallMatchedDistance
      
    }
    
  }
  
  return(list(results= results,
              matchedData = resultSecondMerge))
}
