# +=========================================+
# | NSS Lifespan dataset: Regression of Age |
# +=========================================+
# | LOG                                     |
# | Updated - 07/26/14 - Weizhen Xie        |
# | University of California, Riverside     |
# +=========================================+


# # Check installation for necessary toolbox. 
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
checkInstall("lavaan");
checkInstall("semPlot");
checkInstall("MASS");
checkInstall("DAAG");
checkInstall("bootstrap");

library(bootstrap)
library(Matching);
library(psych);
library(rbounds);
library(foreign);
library(lavaan);
library(semPlot);
library(MASS);
library(DAAG);
library(ggplot2);

# ====================================================================================

# #Read in the first data set... 
#Note:This is what it looks like on my MacBook 
root <-  "/Users/Zane/Dropbox/RUC_CAS_UCR_Research/Chan_CAS_NSS/Data_Analysis/";
#filename = "Healthy_Data.csv"
#filename = "Healthy_14_62.csv"

#filename = "Schiz_Data.csv"
filename = "Health_Schiz_combine.csv"

filepath <- paste(root,filename,sep=''); # treatment group always in the second digit
data <- read.csv(filepath,sep=",",header=TRUE);
# Note: We will be working with Standardized Regression Models (We want Betas for Path Analysis)
# data = as.data.frame(scale(data));
# In this data set, we actually have reason not to rescale the data

# Check the data. 
str(data)

# See the corrlation matrix
cor(data);


# center the data
data$Age_center <- as.vector(scale(data$Age,scale=F));
data$IQEdu_center <- as.vector(scale(data$IQEdu,scale=F));
data$Age_centersq  <- data$Age_center^2;
data$Age_centercubic  <- data$Age_center^3;

#zscore everthing, and then get the square of age. 
#data = as.data.frame(scale(data));


#============= Model  age only linear =========================# 
model <- lm(NSS ~ Age_center, data);
(model.sum <- summary(model));
sqrt(model.sum$r.square);

library(bootstrap)
# define functions 
theta.fit <- function(x,y){lsfit(x,y)}
theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef} 

# matrix of predictorsuj
X <- as.matrix(data[c("Age_center")])
# vector of predicted values
y <- as.matrix(data[c("NSS")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=2)
cor(y, model$fitted.values) # raw R
cor(y,results$cv.fit) # cross-validated R
cor(results$cv.fit, model$fitted.values) # intercorrelation of R
length(results$groups[[1]])
summary(results)


#============= Model1  age only quadradic  =========================# 
# Model 1 age  quadradic
model1 <- lm(NSS ~ Age_center + Age_centersq, data);
(model1.sum <- summary(model1));
anova(model1)
anova(model, model1)

# matrix of predictorsuj
X <- as.matrix(data[c("Age_center","Age_centersq")])
# vector of predicted values
y <- as.matrix(data[c("NSS")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=2)
cor(y, model1$fitted.values) # raw R
cor(y,results$cv.fit) # cross-validated R
cor(results$cv.fit, model1$fitted.values) # intercorrelation of R
length(results$groups[[1]])
summary(results)


# ================ Model 2 age gender ============= # 
model2 <- lm(NSS ~ Age_center + Age_centersq + Gender, data);
(model2.sum <- summary(model2));
sqrt(model2.sum$r.square);
anova(model1, model2)

# matrix of predictorsuj
X <- as.matrix(data[c("Age_center","Age_centersq","Gender")])
# vector of predicted values
y <- as.matrix(data[c("NSS")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=2)
cor(y, model2$fitted.values) # raw R
cor(y,results$cv.fit) # cross-validated R
cor(results$cv.fit, model2$fitted.values) # intercorrelation of R
length(results$groups[[1]])
summary(results)


# ==============  Model 3 age IQEdu ===============# 
model3 <- lm(NSS ~ Age_center + Age_centersq +Gender + IQEdu_center, data);
(model3.sum <- summary(model3));
sqrt(model3.sum$r.square);
anova(model1, model3)
anova(model3)
# matrix of predictorsuj
X <- as.matrix(data[c("Age_center","Age_centersq","Gender","IQEdu_center")])
# vector of predicted values
y <- as.matrix(data[c("NSS")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=2)
cor(y, model3$fitted.values) # raw R
cor(y,results$cv.fit) # cross-validated R
cor(results$cv.fit, model3$fitted.values) # intercorrelation of R
length(results$groups[[1]])
summary(results)


data = as.data.frame(scale(data));
