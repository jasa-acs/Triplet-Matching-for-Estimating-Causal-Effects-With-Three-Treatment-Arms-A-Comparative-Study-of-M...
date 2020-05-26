library(sas7bdat)
library(nnet)
library(optmatch)

setwd("U:/")

source("simulation analysis/functionsMatching_6.R")

########
# DATA #
########
#Load data and construction analytic file
datH1 <- read.sas7bdat(file = "trauma_25_75_h1.sas7bdat")
datH2 <- read.sas7bdat(file = "trauma_25_75_h2.sas7bdat")

table(datH1$HOSP_TRAUMA,exclude = NULL)
table(datH2$HOSP_TRAUMA,exclude = NULL)

dat <- rbind(datH1,datH2[datH2$HOSP_TRAUMA %in% 0,])

#Treatment variable
table(dat$HOSP_TRAUMA,exclude = NULL)
dat$HOSP_TRAUMA <- factor(dat$HOSP_TRAUMA, levels = c(0, 1, 2))

# summary(dat$iss)
# summary(dat$AGE)
# summary(dat$chronic)
# summary(dat$income_2)
# summary(dat$income_3)
# summary(dat$income_4)
# summary(dat$pay1_1)
# summary(dat$pay1_2)
# summary(dat$pay1_4)
# summary(dat$pay1_5)
# summary(dat$pay1_6)
# 
# summary(dat$nchs_2)
# summary(dat$nchs_3)
# summary(dat$nchs_4)
# summary(dat$nchs_5)
# summary(dat$nchs_6)
# summary(dat$multiple_injury)
# #summary(dat$MULTINJURY)
# summary(dat$FEMALE)

####################
# Propensity score #
####################
formulaPropScore <- formula(HOSP_TRAUMA ~ iss + AGE + 
                              income_2 + income_3 + chronic + multiple_injury + income_4 + 
                              pay1_1 + pay1_2 + pay1_4 + pay1_5 + pay1_6 + 
                              nchs_2 + nchs_3 + nchs_4 + nchs_5 + nchs_6+
                              FEMALE) 
propScoreModel <- multinom(formulaPropScore, family=binomial(link=logit),data=dat)
probabilitiesPS <- predict(propScoreModel, type = "probs")

dat$logit1vs0 <- log(probabilitiesPS[,"1"]/probabilitiesPS[,"0"])
dat$logit2vs0 <- log(probabilitiesPS[,"2"]/probabilitiesPS[,"0"])

dat$prob0 <- probabilitiesPS[,"0"]
dat$prob1 <- probabilitiesPS[,"1"]
dat$prob2 <- probabilitiesPS[,"2"]

rm(propScoreModel, probabilitiesPS)
rm(datH1, datH2)
gc(reset=T)

#Info about propensity score
#----------------------------
# plot(density(dat$prob0))
# plot(density(dat$prob1))
# plot(density(dat$prob2))

############
# Matching #
############

variableTreatment <- "HOSP_TRAUMA"
variablesMatch <- c("logit1vs0","logit2vs0")

#3-way constrained optimal matching
#----------------------

#3-way optimal matching - change of starting edge
start <- Sys.time()
startingEdges <- c("0-2","1-2","0-1")
options("optmatch_max_problem_size" = Inf)
dat$HOSP_TRAUMA <- as.numeric(as.character(dat$HOSP_TRAUMA))

for(iter in 1:3) {
  
  result3wayConstr <- applyModifiedOptMatching(data = dat, startingEdges[iter], 
                                                 variablesMatch = variablesMatch, 
                                                 variableTreatment = variableTreatment)
  endTemp <- Sys.time()
  print(endTemp - start)
  save(result3wayConstr, file = paste("bestResult_3wayconstr_allData_",startingEdges[iter],".Rdata",sep=""))
  
}
end <- Sys.time()
end - start
