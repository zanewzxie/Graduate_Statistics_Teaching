# +=========================================+
# | NSS Lifespan dataset: Propensity score  |
# +=========================================+
# | LOG                                     |
# | Updated - 07/23/14 - Weizhen Xie        |
# | University of California, Riverside     |
# +=========================================+

# Check installation for necessary toolbox. 
checkInstall <- function(packageName)
{
  packages <- installed.packages();
  if (packageName %in% packages[,1] == FALSE)
  {
    install.packages(packageName);
  };  
}
checkInstall("Matching");
checkInstall("rbound");
checkInstall("psych");
checkInstall("foreign");
library(Matching);
library(psych);
library(rbounds);
library(foreign);

#Read in the first data set... 
#Note:This is what it looks like on my MacBook 
root <-  "/Users/Zane/Dropbox/RUC_CAS_UCR_Research/Chan_CAS_NSS/Data_Analysis/";
filename = "Match_1_3.csv"
#filename = "Match_1_2.csv"
#filename = "Match_1_4.csv"
#filename = "Match_1_5.csv"

#filename = "Match_3_2.csv"
#filename = "Match_3_4.csv"
#filename = "Match_3_5.csv"

#filename = "Match_2_4.csv"
#filename = "Match_5_2.csv"
#filename = "Match_5_4.csv"

#filename = "Match_Chronicity_First.csv"
#filename = "Match_Chronicity_Healthy.csv"
#filename = "Match_First_Healthy.csv"


filepath <- paste(root,filename,sep=''); # treatment group always in the second digit
data <- read.csv(filepath,sep=",",header=TRUE);

# Check the data. 
str(data)

# Estimate the propensity model, based on Age, Gendr, IQ, ReEdu
glm1  <- glm(Trt ~ Age + Gender + IQEdu, family = binomial, data=data)

#save data objects
X  <- glm1$fitted
Y  <- data$NSS
Tr  <- data$Trt  # this is a vector to determine which group is the treatment group
                # and which is the control. Usually, the one with smaller sample size is the treatment group

#=============== 1. Exact match with replacement ================
# nearest neighbour exact matching without replacement.
# Note: Estimand should change into ATC in interested in control effect. 
rr1  <- Match(Y=Y, Tr=Tr, X=X, exact = TRUE, replace = TRUE,  Weight =2,tolerance=sqrt(.Machine$double.eps))
summary(rr1,full=TRUE)

# We need these following informaiton for each match
# 1. balance check for X covariates, mean value, standard deviation, standardized mean difference

# Before matching 
describeBy(data$Age,Tr)
describeBy(data$Gender,Tr)
describeBy(data$IQEdu,Tr)
describeBy(data$NSS,Tr)

# After matching 
#mb1  <- MatchBalance(Trt~ Age + Gender + IQ + ReEdu, data=data, match.out=rr1, nboots=500)
mb1  <- MatchBalance(Trt~ Age + Gender + IQEdu, data=data, match.out=rr1, nboots=500)

mb1$BeforeMatching[[1]]$tt$statistic
mb1$AfterMatching[[1]]$tt$statistic

mb1$BeforeMatching[[3]]$tt$statistic
mb1$AfterMatching[[3]]$tt$statistic


# The t-test for before and after matching is avaliable in this function. 
# sort the results as below to put into a table. 
uniqueSubject <- unique(rr1$index.treated)
trtgender <- data$Gender[rr1$index.treated]*rr1$weights
congender <- data$Gender[rr1$index.control]*rr1$weights
trtsex <- NA
consex <- NA
for (i in 1:length(uniqueSubject)) 
{trtsex[i] <- sum(trtgender[rr1$index.treated== uniqueSubject[i]])
 consex[i] <- sum(congender[rr1$index.treated== uniqueSubject[i]])}
trmale <- sum(trtsex<1.5)
conmale <- sum(consex<1.5)
x <- matrix(c(trmale,length(uniqueSubject)- trmale, conmale, length(uniqueSubject)- conmale), ncol = 2)
(chigendertest <- chisq.test(x))
trtmaleper <- trmale/length(uniqueSubject)
conmaleper <- conmale/length(uniqueSubject)

Age  <- balanceUV(data$Age[rr1$index.treated],data$Age[rr1$index.control], weights= rr1$weights,match = TRUE,paired = TRUE)
outputname <- c('Trmean','Trstd', 'ConMean','Contrt')
outputdata <-  c(Age$mean.Tr,sqrt(Age$var.Tr),Age$mean.Co,sqrt(Age$var.Co))
output <- data.frame(outputname,outputdata)
output
Age$tt

IQEdu <- balanceUV(data$IQEdu[rr1$index.treated],data$IQEdu[rr1$index.control], weights= rr1$weights,match = TRUE)
outputname <- c('Trmean','Trstd', 'ConMean','Contrt')
outputdata <-  c(IQEdu$mean.Tr,sqrt(IQEdu$var.Tr),IQEdu$mean.Co,sqrt(IQEdu$var.Co))
output <- data.frame(outputname,outputdata)
output
IQEdu$tt

ReEdu <- balanceUV(data$ReEdu[rr1$index.treated],data$ReEdu[rr1$index.control], weights= rr1$weights,match = TRUE)
outputname <- c('Trmean','Trstd', 'ConMean','Contrt')
outputdata <-  c(ReEdu$mean.Tr,sqrt(ReEdu$var.Tr),ReEdu$mean.Co,sqrt(ReEdu$var.Co))
output <- data.frame(outputname,outputdata)
output

Outcome <- balanceUV(data$NSS[rr1$index.treated],data$NSS[rr1$index.control], weights= rr1$weights,match = TRUE,paired=TRUE)
outputname <- c('Trmean','Trstd', 'ConMean','Contrt')
outputdata <-  c(Outcome$mean.Tr,sqrt(Outcome$var.Tr),Outcome$mean.Co,sqrt(Outcome$var.Co))
output <- data.frame(outputname,outputdata)
output
Outcome$tt

savedataname = paste(root,"Exact_",filename,sep='')
savedata <- cbind(rr1$index.control,rr1$index.treated, rr1$weights)
write.csv(savedata, file = savedataname,row.names = TRUE)

# Look at the outcomes 

#=============== Caliper, oversampling, with replacement ================
# Caliper matching, with replacement and allow oversampling (find the best match 2 data points)
rr2  <- Match(Y=Y, Tr=Tr, X=X, M = 1, caliper = 0.1, estimand="ATT", replace = TRUE, Weight =2,tolerance=sqrt(.Machine$double.eps))
summary(rr2,full=TRUE)

# We need these following informaiton for each match
# 1. balance check for X covariates, mean value, standard deviation, standardized mean difference

# Before matching 
describeBy(data$Age,Tr)
describeBy(data$Gender,Tr)
describeBy(data$IQ,Tr)

# After matching 
#mb2  <- MatchBalance(Trt~ Age + Gender + IQ + ReEdu, data=data, match.out=rr2, nboots=500)
mb2  <- MatchBalance(Trt~ Age + Gender + IQEdu, data=data, match.out=rr2, nboots=500)

mb2$BeforeMatching[[1]]$tt$statistic
mb2$AfterMatching[[1]]$tt$statistic

mb2$BeforeMatching[[3]]$tt$statistic
mb2$AfterMatching[[3]]$tt$statistic

# The t-test for before and after matching is avaliable in this function. 
# The t-test for before and after matching is avaliable in this function. 
# sort the results as below to put into a table. 
uniqueSubject <- unique(rr2$index.treated)
trtgender <- data$Gender[rr2$index.treated]*rr2$weights
congender <- data$Gender[rr2$index.control]*rr2$weights
trtsex <- NA
consex <- NA
for (i in 1:length(uniqueSubject)) 
{trtsex[i] <- sum(trtgender[rr2$index.treated== uniqueSubject[i]])
 consex[i] <- sum(congender[rr2$index.treated== uniqueSubject[i]])}
trmale <- sum(trtsex<=1.5)
conmale <- sum(consex<=1.5)
x <- matrix(c(trmale,length(uniqueSubject)- trmale, conmale, length(uniqueSubject)- conmale), ncol = 2)
(chigendertest <- chisq.test(x))
trtmaleper <- trmale/length(uniqueSubject)
conmaleper <- conmale/length(uniqueSubject)

Age  <- balanceUV(data$Age[rr2$index.treated],data$Age[rr2$index.control], weights= rr2$weights,match = TRUE,paired = TRUE)
outputname <- c('Trmean','Trstd', 'ConMean','Contrt')
outputdata <-  c(Age$mean.Tr,sqrt(Age$var.Tr),Age$mean.Co,sqrt(Age$var.Co))
output <- data.frame(outputname,outputdata)
output
Age$tt

IQEdu <- balanceUV(data$IQEdu[rr2$index.treated],data$IQ[rr2$index.control], weights= rr2$weights,match = TRUE)
outputname <- c('Trmean','Trstd', 'ConMean','Contrt')
outputdata <-  c(IQEdu$mean.Tr,sqrt(IQEdu$var.Tr),IQEdu$mean.Co,sqrt(IQEdu$var.Co))
output <- data.frame(outputname,outputdata)
output
IQEdu$tt

Outcome <- balanceUV(data$NSS[rr2$index.treated],data$NSS[rr2$index.control], weights= rr2$weights,match = TRUE,paired=TRUE)
outputname <- c('Trmean','Trstd', 'ConMean','Contrt')
outputdata <-  c(Outcome$mean.Tr,sqrt(Outcome$var.Tr),Outcome$mean.Co,sqrt(Outcome$var.Co))
output <- data.frame(outputname,outputdata)
output
Outcome$tt

savedataname = paste(root,"Caliper_",filename,sep='')
savedata <- cbind(rr2$index.control, rr2$index.treated, rr2$weights)
write.csv(savedata, file = savedataname,row.names = TRUE)

# Second matching procedure: Capliper matching procedure, and 





# Let's check the covariate balance
# 'nboots' is set to small values in the interest of speed.
# Please increase to at least 500 each for publication quality p-values.  

ks  <- balanceUV(Y[rr$index.treated], Y[rr$index.control], nboots=500)
summary(ks,full=TRUE)


qqout  <- qqstats(Y[rr$index.treated], Y[rr$index.control])
print(qqout)


#Extract Matched Outome Data
ctrl <- Y[rr$index.control]
trt <- Y[rr$index.treated]

#Extract age Doses
age.trt <- data$Age[rr$index.treated]
age.ctrl <- data$Age[rr$index.control]

#Run Sensitivity Analsyis
iv_sens(trt, ctrl, age.trt, age.ctrl, Gamma=2.5, GammaInc=.1)


# Sensitivity test
tmp <- data.prep(rr, group.size=3)
mcontrol(tmp$Y, tmp$id, tmp$treat, group.size = 3)


