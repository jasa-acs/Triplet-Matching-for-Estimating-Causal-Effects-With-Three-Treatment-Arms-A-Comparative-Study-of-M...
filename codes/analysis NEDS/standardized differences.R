library(sas7bdat)
library(nnet)
library(optmatch)

setwd("U:/")

stdzDifferenceCont <- function(data, selectionMatched, contVar, factorVar, levelsFactor){
  
  #For difference in means, matched data
  dataMeans <- data[data[,factorVar] %in% levelsFactor & selectionMatched, ]
  means <- tapply(dataMeans[,contVar], INDEX = dataMeans[,factorVar], FUN = mean)
  
  #For variances, unmatched data
  dataSds <- data[data[,factorVar] %in% levelsFactor, ]
  sds <- tapply(dataSds[,contVar], INDEX = dataSds[,factorVar], FUN = sd )
  vars <- sds^2
  
  stdzDiff <- (means[names(means)==levelsFactor[1]] - means[names(means)==levelsFactor[2]])/(
                sqrt((vars[names(vars)==levelsFactor[1]] + vars[names(vars)==levelsFactor[2]])/2))
  return(stdzDiff)
  
}

stdzDifferenceCat <- function(data, selectionMatched, catVar, factorVar, levelsFactor){
  
  #For difference in means, matched data
  dataMeans <- data[data[,factorVar] %in% levelsFactor & selectionMatched, ]
  means <- tapply(dataMeans[,catVar], INDEX = dataMeans[,factorVar], FUN = mean)
  
  #For variances, unmatched data
  dataMeansUnm <- data[data[,factorVar] %in% levelsFactor, ]
  meansUnm <- tapply(dataMeansUnm[,catVar], INDEX = dataMeansUnm[,factorVar], FUN = mean)
  vars <- meansUnm*(1-meansUnm)
  
  stdzDiff <- (means[names(means)==levelsFactor[1]] - means[names(means)==levelsFactor[2]])/(
              sqrt((vars[names(vars)==levelsFactor[1]] + vars[names(vars)==levelsFactor[2]])/2))
  return(stdzDiff)
  
}

#Load data
load("NCH data/bestResult_3wayconstr_allData_1-2.Rdata")

############################
# Standardized differences #
############################
results <- matrix(NA, ncol = 3, nrow = 21)

results[1,1] <- stdzDifferenceCont(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch), contVar = "AGE", factorVar = "HOSP_TRAUMA",
                                   levelsFactor = c(0,1))
results[1,2] <- stdzDifferenceCont(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch), contVar = "AGE", factorVar = "HOSP_TRAUMA",
                                   levelsFactor = c(1,2))
results[1,3] <- stdzDifferenceCont(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch), contVar = "AGE", factorVar = "HOSP_TRAUMA",
                                   levelsFactor = c(0,2))

results[2,1] <- stdzDifferenceCont(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch), contVar = "iss", factorVar = "HOSP_TRAUMA",
                                   levelsFactor = c(0,1))
results[2,2] <- stdzDifferenceCont(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch), contVar = "iss", factorVar = "HOSP_TRAUMA",
                                   levelsFactor = c(1,2))
results[2,3] <- stdzDifferenceCont(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch), contVar = "iss", factorVar = "HOSP_TRAUMA",
                                   levelsFactor = c(0,2))


results[3,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"chronic","HOSP_TRAUMA", c(0,1))
results[3,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"chronic","HOSP_TRAUMA", c(1,2))
results[3,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"chronic","HOSP_TRAUMA", c(0,2))

results[4,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_2","HOSP_TRAUMA", c(0,1))
results[4,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_2","HOSP_TRAUMA", c(1,2))
results[4,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_2","HOSP_TRAUMA", c(0,2))

results[5,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_3","HOSP_TRAUMA", c(0,1))
results[5,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_3","HOSP_TRAUMA", c(1,2))
results[5,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_3","HOSP_TRAUMA", c(0,2))

results[6,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_4","HOSP_TRAUMA", c(0,1))
results[6,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_4","HOSP_TRAUMA", c(1,2))
results[6,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_4","HOSP_TRAUMA", c(0,2))

results[7,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_1","HOSP_TRAUMA", c(0,1))
results[7,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_1","HOSP_TRAUMA", c(1,2))
results[7,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_1","HOSP_TRAUMA", c(0,2))

results[8,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_2","HOSP_TRAUMA", c(0,1))
results[8,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_2","HOSP_TRAUMA", c(1,2))
results[8,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_2","HOSP_TRAUMA", c(0,2))

results[9,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_4","HOSP_TRAUMA", c(0,1))
results[9,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_4","HOSP_TRAUMA", c(1,2))
results[9,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_4","HOSP_TRAUMA", c(0,2))

results[10,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_5","HOSP_TRAUMA", c(0,1))
results[10,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_5","HOSP_TRAUMA", c(1,2))
results[10,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_5","HOSP_TRAUMA", c(0,2))

results[11,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_6","HOSP_TRAUMA", c(0,1))
results[11,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_6","HOSP_TRAUMA", c(1,2))
results[11,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_6","HOSP_TRAUMA", c(0,2))


results[12,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_2","HOSP_TRAUMA", c(0,1))
results[12,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_2","HOSP_TRAUMA", c(1,2))
results[12,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_2","HOSP_TRAUMA", c(0,2))

results[13,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_3","HOSP_TRAUMA", c(0,1))
results[13,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_3","HOSP_TRAUMA", c(1,2))
results[13,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_3","HOSP_TRAUMA", c(0,2))

results[14,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_4","HOSP_TRAUMA", c(0,1))
results[14,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_4","HOSP_TRAUMA", c(1,2))
results[14,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_4","HOSP_TRAUMA", c(0,2))

results[15,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_5","HOSP_TRAUMA", c(0,1))
results[15,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_5","HOSP_TRAUMA", c(1,2))
results[15,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_5","HOSP_TRAUMA", c(0,2))

results[16,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_6","HOSP_TRAUMA", c(0,1))
results[16,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_6","HOSP_TRAUMA", c(1,2))
results[16,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_6","HOSP_TRAUMA", c(0,2))

results[17,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"multiple_injury","HOSP_TRAUMA", c(0,1))
results[17,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"multiple_injury","HOSP_TRAUMA", c(1,2))
results[17,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"multiple_injury","HOSP_TRAUMA", c(0,2))

results[18,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"FEMALE","HOSP_TRAUMA", c(0,1))
results[18,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"FEMALE","HOSP_TRAUMA", c(1,2))
results[18,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"FEMALE","HOSP_TRAUMA", c(0,2))

results[19,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_1","HOSP_TRAUMA", c(0,1))
results[19,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_1","HOSP_TRAUMA", c(1,2))
results[19,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"income_1","HOSP_TRAUMA", c(0,2))

results[20,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_1","HOSP_TRAUMA", c(0,1))
results[20,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_1","HOSP_TRAUMA", c(1,2))
results[20,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"nchs_1","HOSP_TRAUMA", c(0,2))

results[21,1] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_3","HOSP_TRAUMA", c(0,1))
results[21,2] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_3","HOSP_TRAUMA", c(1,2))
results[21,3] <- stdzDifferenceCat(result3wayConstr$matchedData, selectionMatched = !is.na(result3wayConstr$matchedData$indexMatch),"pay1_3","HOSP_TRAUMA", c(0,2))
