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



# =============  Model 4 age IQEdu, and grou ================= # 
model4 <- lm(NSS ~ Age_center + Age_centersq + Gender + IQEdu_center + Group, data);
(model4.sum <- summary(model4));
sqrt(model4.sum$r.square);
anova(model3, model4)

# matrix of predictorsuj
X <- as.matrix(data[c("Age_center","Age_centersq","Gender","IQEdu_center","Group")])
# vector of predicted values
y <- as.matrix(data[c("NSS")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=2)
cor(y, model4$fitted.values) # raw R
cor(y,results$cv.fit) # cross-validated R
cor(results$cv.fit, model4$fitted.values) # intercorrelation of R
length(results$groups[[1]])
summary(results)


# =================  Model 5 one interaction ============= 
data$Age_centerGroup <- as.vector(data$Age_center*data$Group);
data$Age_centersqGroup <- as.vector(data$Age_centersq*data$Group);

model5 <- lm(NSS ~ Age_center + Age_centersq + Gender + IQEdu_center + Group + 
               I(Group*Age_center) + I(Group*Age_centersq), data);
(model5.sum <- summary(model5));
sqrt(model5.sum$r.square);
anova(model4, model5)
# matrix of predictorsuj
X <- as.matrix(data[c("Age_center","Age_centersq","Gender","IQEdu_center","Group","Age_centerGroup","Age_centersqGroup")])
# vector of predicted values
y <- as.matrix(data[c("NSS")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=2)
cor(y, model5$fitted.values) # raw R
cor(y,results$cv.fit) # cross-validated R
cor(results$cv.fit, model5$fitted.values) # intercorrelation of R
length(results$groups[[1]])
summary(results)


# =================  Model6 one interaction ============= 
data$groupgender = as.vector(data$Group*data$Gender);

model6 <- lm(NSS ~ Age_center + Age_centersq + Gender + IQEdu_center + Group + 
               I(Group*Age_center) + I(Group*Age_centersq) + I(Group*Gender), data);
(model6.sum <- summary(model6));
sqrt(model6.sum$r.square);
anova(model5, model6)

# matrix of predictorsuj
X <- as.matrix(data[c("Age_center","Age_centersq","Gender","IQEdu_center","Group","Age_centerGroup","Age_centersqGroup","groupgender")])
# vector of predicted values
y <- as.matrix(data[c("NSS")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=2)
cor(y, model6$fitted.values) # raw R
cor(y,results$cv.fit) # cross-validated R
cor(results$cv.fit, model6$fitted.values) # intercorrelation of R
length(results$groups[[1]])
summary(results)

# =================  Model7 one interaction ============= 
data$GroupIQEdu = as.vector(data$Group*data$IQEdu_center);

model7 <- lm(NSS ~ Age_center + Age_centersq + Gender + IQEdu_center + Group + 
               I(Group*Age_center) + I(Group*Age_centersq) + I(Group*Gender) + I(Group*IQEdu_center), data);
(model7.sum <- summary(model7));
sqrt(model7.sum$r.square);
anova(model6, model7)


# matrix of predictorsuj
X <- as.matrix(data[c("Age_center","Age_centersq","Gender","IQEdu_center","Group","Age_centerGroup","Age_centersqGroup","groupgender","GroupIQEdu")])
# vector of predicted values
y <- as.matrix(data[c("NSS")]) 

results <- crossval(X,y,theta.fit,theta.predict,ngroup=2)
cor(y, model7$fitted.values) # raw R
cor(y,results$cv.fit) # cross-validated R
cor(results$cv.fit, model7$fitted.values) # intercorrelation of R
length(results$groups[[1]])
summary(results)



#=============Interaction Effect of Gender * Group ================ # 
modelsexgroup <- lm(NSS ~ Gender + Group + I(Group*Gender), data);
(modelsexgroup.sum <- summary(modelsexgroup));
sqrt(modelsexgroup.sum$r.square);

modelIQEdu <- lm(NSS ~ IQEdu + Group + I(Group*IQEdu), data);
(modelIQEdu.sum <- summary(modelIQEdu));
sqrt(modelIQEdu.sum$r.square);







# ========================================================================================================================================= # 





# Note: Adding I(Group*Age_centersq) did not help improve the overall fit, and it is not significant as well


# Model 5 age IQEdu, and two interaction
model5 <- lm(NSS ~ Age_center + Age_centersq + IQEdu_center + I(IQEdu_center*Age_center) + I(IQEdu_center* Age_centersq) , data);
(model5.sum <- summary(model5));
sqrt(model5.sum$r.square);
anova(model4, model5)
data$yhat <- predict(model5);
ggplot(data, aes(x=Age_center, y=yhat)) +
  geom_point(shape=1,size=4)+
  geom_smooth();

CVlm(df=nihills, form.lm=formula(log(time)~log(climb)+log(dist)),
     plotit="Observed")


model <- lm( IQEdu ~ NSS+Age_center + Age_centersq, data);
(model.sum <- summary(model));
sqrt(model.sum$r.square);

