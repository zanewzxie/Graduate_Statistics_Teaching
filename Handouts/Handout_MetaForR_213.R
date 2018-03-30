
#Set the directory to your data folder
setwd("....")
#This is what it looks like on my MacBook 

# OPTIONAL: Clear R's memory and graphics:
rm(list=ls())  # Careful! This clears all of R's memory!
graphics.off() # This closes all of R's graphics windows.

# load the metafor toolbox/function
library(metafor)# Get the functions loaded into R's working memory:

# read your data
dat <- read.csv("xxx.csv",sep=",",header=TRUE);

# analyzie fixed effect
fix_res <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="FE") ### use r-to-z transformed correlations

# analyzie random effect
ran_res <- rma(measure="ZCOR", ri=dat$ri, ni=dat$ni, method="REML") ### use r-to-z transformed correlations

# plot the forest plot
forest(fix_res, slab=paste(dat$author, dat$year, sep=", "),
       xlim=c(-3,2),alim=c(-1,1), ylim=c(-2, 22) )  # forest plot based on fixed-effect results
addpoly(ran_res) # add the random-effect results. 
#text(-3,                21, "Study ID. Author(s) and Year",  pos=4)
#text(1.8,                21, "Zr[95% CI]", pos=2)


# plot the funnal plot
funnel(fix_res, yaxis="sei")
funnel(res, main="Standard Error")









