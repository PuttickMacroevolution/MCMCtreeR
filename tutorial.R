##################################################################################
#########################	 SET UP MCMCTREE DATING ANALYSIS  ####################
##################################################################################
######    written by Mark Puttick 2016 (mark.puttick@bristol.ac.uk)         ######
##################################################################################


# This tutorial is designed to show how to estimate node distribution for input into the software mcmcTree. By default the code here estimates the parameters that are needed for a distribution the spans from the a minimum bound (lower age) to an maximum bound (upper age). The distributions are treated so there is a 'hard' minimum and 'soft maximum', and 97.5% of the distribution falls between the minimum and maximum ages. 

# The code can estimate paramaters for the caucy, skew-t, skew-normal, and gamma distirbutions shown in the PAML manual on page 49. 

# Once these parameters are estimated the code can plot them and write them as node labels on a tree that is output in the correct format for mcmcTree. Additionally, the the code can write out node constraintst that are fixed, and have a lower and upper bound respectively.

# The code can be used to simultanaeouly estimate the parameter values and write files. This can be done across the whole tree with a specific type of distribution placed on each node, or by placing different constraints on different nodes.


##################################################################################
#########################	 set working directory and obtain data    ############
##################################################################################

rm(list=ls(all=T))
library(ape)
setwd("")

#install.packages("devtools")
library(devtools)
install_github("PuttickMacroevolution/MCMCTreeR")
library(MCMCTreeR)

# with the library installed we can check the help for any function: i.e, 

?estimateCauchy
		
# we'll get some data for apes: a phylogeny, a list containing members of each clade of a node of interest, and vectors of maximum and minimum times

data(apeData)
attach(apeData)
names(apeData)

# [1] "minimumTimes"       "maximumTimes"      
# [3] "monophyleticGroups" "apeTree"  

# minimumTimes - vector of minimum bounds
# maximumTimes - vector of maximum bounds
# monophyleticGroups - list containing species in each clade
# apeTree - species phhylogeny of apes in APE format


##################################################################################
#########################	 Estimate parameters  ################################
##################################################################################

# We have our phylogeny and node ages. We can now estimate the appropriate parameters for a skew t distribution. In this estimateSkewT function the default assumes we want to estimate scale with a given shape value (the default shape value is 50). 

##################################################################################
# details of the function
#########################	 Skew t distribution  ################################

#estimateSkewT(minAge, maxAge, monoGroups, phy, shape=50, scale=1.5, df=1, addMode=0, maxProb=0.975, minProb=0.003, estimateScale=TRUE, estimateShape=FALSE, estimateMode=F, plot=FALSE, pdfOutput="skewTPlot.pdf", writeMCMCTree=FALSE, mcmcTreeName="skewTInput.tre")

# minAge - vector of minimum age bounds for nodes matching order in monoGroups
# maxAge - vector of maximum age bounds for nodes matching order in monoGroups
# monoGroups - list  with each element containing species that define a node of interest
# phy - fully resolved phylogeny in ape format
# shape - shape value for skew-t distribution (default = 50)
# scale- scale value for skew-t distribution (default = 1.5)
# df - degrees of freedom for skew-t distribution (default = 1)
# addMode - addition to the minimum age to give the location
# maxProb - probability of right tail (maximum bound. default = 0.975) 
# estimateScale - logical specifying whether to estimate scale with a given shape value (default = TRUE)
# estimateShape - logical specifying whether to estimate shape with a given scale value (default = TRUE)
# estimateMode - logical speciftying whether to estimate the scale that produces probabilities of each tail that corresponds roughly to the values given by minProb (lower tail) and maxProb (upper tail)
# plot - logical specifying whether to plot to PDF
# pdfOutput - pdf file name
# writeMCMCTree - logical whether to write tree in format that is compatible with mcmcTree to file
# mcmcTreeName - mcmcTree input file name
##################################################################################

#########################	 estimate scale with a given shape

# The function will take input minimum input times, and estimate the value of the scale that will produce a distribution with the 97.5% cumulative probability of the distribution at the maximum age.

skewT_results <- estimateSkewT(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, phy=apeTree, plot=F)

skewT_results$mcmctree

#we can look at the results. the parameters show the input for the mcmctree format - the location (lower node age), scale, shape, and degrees of freedom

par(mfrow=c(2,2), family="Times")
for(i in 1:4) plotMCMCTree(skewT_results$parameters[i,], method="skewT", title=paste0("node ", i), upperTime=max(maximumTimes))

#if we're happy with this we can write the file to the input format of mcmctree

write.table(skewT_results$mcmctree, "mcmcInputTree.tre", quote=F, row.names=F, col.names=F)

#alternatively we could have set this as default in the estimateSkewt function using the "writeMCMCtree" argument, and set the name in the mcmcTreeName argument. Additionally, we could have requested a pdf output of the estimate distributions by using 'plot=TRUE' and specifying the file of this name using the pdfOutput

skewT_results <- estimateSkewT(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, phy=apeTree, plot=T, pdfOutput="skewTPlot.pdf", writeMCMCTree=T, mcmcTreeName="skewTInput.tre")

# we can change the scale and shape (and in skew t the df arguments from the default values)

# we can change the value of shape to be different for each distribution

skewT_results <- estimateSkewT(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, shape=c(9, 10, 8, 10),  phy=apeTree, plot=T, pdfOutput="skewTPlot.pdf", writeMCMCTree=T, mcmcTreeName="skewTInput.tre")
skewT_results$parameters


#########################	 estimate shape with a given scale

skewT_results <- estimateSkewT(minAge=minimumTimes[2], maxAge=maximumTimes[2], monoGroups=monophyleticGroups, scale=0.05, estimateShape=T, estimateScale=F,  phy=apeTree, plot=T, pdfOutput="skewTPlot.pdf", writeMCMCTree=T, mcmcTreeName="skewTInput.tre")
skewT_results$parameters


#########################	 estimate the location to produce tails of 3% lower and 97.5% upper probability at the minimum and maximum ages 

skewT_results <- estimateSkewT(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, shape=20, estimateShape=F, estimateScale=F, estimateMode=T, maxProb=0.975, minProb=0.03,  phy=apeTree, plot=T, pdfOutput="skewTPlot.pdf", writeMCMCTree=T, mcmcTreeName="skewTInput.tre")
skewT_results$parameters
par(mfrow=c(2,2), family="Times")
for(i in 1:4) plotMCMCTree(skewT_results$parameters[i,], method="skewT", title=paste0("node ", i), upperTime=max(maximumTimes))


#########################	 produce tree with given shape and given scale

skewT_results <- estimateSkewT(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, shape=10,  scale=0.5, estimateShape=F, estimateScale=F, phy=apeTree, plot=T, pdfOutput="skewTPlot.pdf", writeMCMCTree=T, mcmcTreeName="skewTInput.tre")
skewT_results$parameters

##################################################################################
# details of the function
#########################	 Skew normal distribution  ###########################

#estimateSkewNormal(minAge, maxAge, monoGroups, phy, shape=50, scale=1.5, addMode=0,  maxProb=0.975, minProb=0.003, estimateScale=T, estimateShape=F, estimateMode=F, plot=FALSE,  pdfOutput="skewNormalPlot.pdf", writeMCMCTree=FALSE, mcmcTreeName="skewNormalInput.tre")	

# minAge - vector of minimum age bounds for nodes matching order in monoGroups
# maxAge - vector of maximum age bounds for nodes matching order in monoGroups
# monoGroups - list  with each element containing species that define a node of interest
# phy - fully resolved phylogeny in ape format
# shape - shape value for skew normal distribution (default = 50)
# scale- scale value for skew normal distribution (default = 1.5)
# addMode - addition to the minimum age to give the location
# maxProb - probability of right tail (maximum bound. default = 0.975) 
# estimateScale - logical specifying whether to estimate scale with a given shape value (default = TRUE)
# estimateShape - logical specifying whether to estimate shape with a given scale value (default = TRUE)
# estimateMode - logical speciftying whether to estimate the scale that produces probabilities of each tail that corresponds roughly to the values given by minProb (lower tail) and maxProb (upper tail)
# plot - logical specifying whether to plot to PDF
# pdfOutput - pdf file name
# writeMCMCTree - logical whether to write tree in format that is compatible with mcmcTree to file
# mcmcTreeName - mcmcTree input file name
##################################################################################

#########################	 estimate scale with a given shape
# The function will take input minimum input times, and estimate the value of the scale that will produce a skew normal distribution with the 97.5% cumulative probability of the distribution at the maximum age.

skewNormal_results <- estimateSkewNormal(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, addMode=0.05, phy=apeTree, plot=F)
skewNormal_results$parameters

par(mfrow=c(2,2), family="Times")
for(i in 1:4) plotMCMCTree(skewNormal_results$parameters[i,], method="skewNormal", title=paste0("node ", i), upperTime=max(maximumTimes))


##################################################################################
# details of the function
#########################	Cauchy distribution  #################################
# This is the parameter specified by 'L(tL, p, c, pL)' in the PAML manual page 49 


#estimateCauchy <- function(minAge, maxAge, phy, monoGroups, scaleInput=0.2, offsetInput=0.1, estimateScale=T, minProb=0, maxProb=0.975, plot=FALSE, pdfOutput="cauchyPlot.pdf", writeMCMCTree=FALSE, mcmcTreeName="cauchyInput.tre")

# minAge - vector of minimum age bounds for nodes matching order in monoGroups
# maxAge - vector of maximum age bounds for nodes matching order in monoGroups. By default the 97.5% right tail probability will fall here
# monoGroups - list  with each element containing species that define a node of interest
# phy - fully resolved phylogeny in ape format
# offset - offset value for cauchy distribution (default = 50) (p in PAML manual page 49)
# scale - scale value for cauchy distribution (default = 1.5) (c in PAML manual page 49)
# estimateScale - logical specifying whether to estimate scale with a given shape value (default = TRUE)
# minProb - probability of left tail (minimum bound) - default to hard minimum (minProb=0)
# maxProb - probability of right tail (maximum bound. default = 0.975) 
# plot - logical specifying whether to plot to PDF
# pdfOutput - pdf file name
# writeMCMCTree - logical whether to write tree in format that is compatible with mcmcTree to file
# mcmcTreeName - mcmcTree input file name
##################################################################################

### plot example on page 50 of PAML manual

example_page_50 <- estimateCauchy(minAge=1, maxAge=4.32, monoGroups=monophyleticGroups[[1]],   phy=apeTree,  offset=0.5, minProb=0.025, plot=F)[[1]]
dev.new()
plotMCMCTree(example_page_50, method="cauchy", title=paste0("node ", i), upperTime=max(maximumTimes))

#########################	 estimate scale with a given shape
# The function will take input minimum input times, and estimate the value of the scale that will produce a cauchy distribution with the 97.5% cumulative probability of the distribution at the maximum age.

cauchy_results <- estimateCauchy(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups,  offset=0.5, phy=apeTree, plot=F)
cauchy_results$parameters

par(mfrow=c(2,2), family="Times")
for(i in 1:4) plotMCMCTree(cauchy_results$parameters[i,], method="cauchy", title=paste0("node ", i), upperTime=max(maximumTimes))

# we may have constrained our distribution too much for the 2nd, 3rd, and 4th distribution so we can modify that to allow for smaller offset

cauchy_results <- estimateCauchy(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups,  offset=c(0.5, 0.1, 0.1, 0.05), phy=apeTree, plot=F)
cauchy_results$parameters

par(mfrow=c(2,2), family="Times")
for(i in 1:4) plotMCMCTree(cauchy_results$parameters[i,], method="cauchy", title=paste0("node ", i), upperTime=maximumTimes[i])



#########################	 Uniform distribution  ###################################

#estimateBound(minAge, maxAge, minProb=0.025, rightTail=0.025, phy, monoGroups,  writeMCMCTree=FALSE, mcmcTreeName="bound.tre", pdfOutput="uniformPlot.pdf)

# minAge - vector of minimum age bounds for nodes matching order in monoGroups
# maxAge - vector of maximum age bounds for nodes matching order in monoGroups. By default the 97.5% right tail probability will fall here
# monoGroups - list  with each element containing species that define a node of interest
# phy - fully resolved phylogeny in ape format
# minProb - lower tail probabiliity
# right Tail - right tail probability
# plot - logical specifying whether to plot to PDF
# pdfOutput - pdf file name
# writeMCMCTree - logical whether to write tree in format that is compatible with mcmcTree to file
# mcmcTreeName - mcmcTree input file name
##################################################################################

uniform_results <- estimateBound(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, phy=apeTree)
uniform_results$parameters

par(mfrow=c(2,2), family="Times")
for(i in 1:4) plotMCMCTree(uniform_results$parameters[i,], method="bound", title=paste0("node ", i), upperTime=maximumTimes[i]+1)



#########################	 Gamma distribution  ###################################
# This is very rough!!!

#function(minAge, maxAge, phy, monoGroups, alphaInput=188, betaInput=2690, offset=0.1, estimateAlpha=TRUE, estimateBeta=F,  plot=FALSE, pdfOutput="gammaPlot.pdf", writeMCMCTree=FALSE, mcmcTreeName="gammaInput.tre") 

# minAge - vector of minimum age bounds for nodes matching order in monoGroups
# maxAge - vector of maximum age bounds for nodes matching order in monoGroups. By default the 97.5% right tail probability will fall here
# monoGroups - list  with each element containing species that define a node of interest
# phy - fully resolved phylogeny in ape format
# alpha - alpha value for gamma distribution (default = 188) (p in PAML manual page 49)
# beta - betaInput value for gamma distribution (default = 2690) (c in PAML manual page 49)
# offset - distance of mean value from minimum bound
# estimateAlpha - logical specifying whether to estimate alpha with a given beta value (default = TRUE)
# estimateShape - logical specifying whether to estimate beta with a given alpha value (default = FALSE)
# plot - logical specifying whether to plot to PDF
# pdfOutput - pdf file name
# writeMCMCTree - logical whether to write tree in format that is compatible with mcmcTree to file
# mcmcTreeName - mcmcTree input file name
##################################################################################

gamma_results <- estimateGamma(minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, alpha=188, beta=2690, offset=0.1, phy=apeTree, plot=F)
gamma_results$parameters

par(mfrow=c(2,2), family="Times")
for(i in 1:4) plotMCMCTree(gamma_results$parameters[i,], method="gamma", title=paste0("node ", i), upperTime=maximumTimes[i])


#########################	 Upper Age and fixed ages  ###################################

# we can also specify an upper age bound

upper_results <- estimateUpper(maxAge=maximumTimes, monoGroups=monophyleticGroups, rightTail=0.025, phy=apeTree)
upper_results$parameters

# we can also specify a fixed age for a baseML analysis

fixed_results <- estimateFixed(minAge=minimumTimes[1], monoGroups=monophyleticGroups[[1]], phy=apeTree)


##################################################################################
#########################	 Different parameters on different nodes  ############
##################################################################################

# if we want different parameters on different nodes we can specify this by using the function 'mcmcTreePhy'

# we will set a fixed root (node 1), skew normal (node 2), gamma (node 3), and upper distribution (node 4) to our tree

# for each input we give the associated parameter values in a vector in the order of nodes. i.e - for the 'min prob' on four nodes we want nodes 1, 2, 4 to be 1e-8 and  node 3 to be 0.025
# c(1e-8, 1e-8, 0.025, 1e-8)


#the mcmcTreePhy function

# mcmcTreePhy <- function(phy, minAges, maxAges, monoGroups, method=c("cauchy", "upper", "lowerUpper", "gamma", "skewNormal", "skewT", "fixed"), offset=0.1, df=1, shape=50, scale=1.5, minProb=1e-8, maxProb=0.975, rightTail=0.025, rate=100, plotDist=TRUE, estimateScale=T, estimateShape=F, estimateAlpha=TRUE, estimateBeta=FALSE)	

methodApe <- c("skewT", "cauchy", "gamma", "upper")

outputFull <- mcmcTreePhy(phy=apeTree, minAge=minimumTimes, maxAge=maximumTimes, monoGroups=monophyleticGroups, method=methodApe)

# we can fine-tune this. If for example we want the method to estimate beta not gamma for the 3rd node

estimateAlphaApe=c(FALSE, F, F, F)
estimateBetaApe=c(T, T, T, T)

outputFull <- mcmcTreePhy(phy=apeTree, minAges=minimumTimes, maxAges=maximumTimes, monoGroups=monophyleticGroups, method=methodApe, estimateAlpha=estimateAlphaApe, estimateBeta=estimateBetaApe, alphaInput=188, betaInput=100)


#### we can also use the output from individual methods to input to subsequent node estimation. This allows for easier fine-tuning


skewNormal_results_nodeOne <- estimateSkewNormal(minAge=minimumTimes[1], maxAge=maximumTimes[1], monoGroups=monophyleticGroups[[1]], addMode=0.05, phy=apeTree, plot=F)
skewNormal_results_nodeOne$apePhy

cauchy_results_nodeTwo <- estimateCauchy(minAge=minimumTimes[2], maxAge=maximumTimes[2], monoGroups=monophyleticGroups[[2]],  offset=0.5, phy=skewNormal_results_nodeOne$apePhy, plot=F)
cauchy_results_nodeTwo$apePhy

write.table(cauchy_results_nodeTwo$mcmctree, quote=F, row.names=F, col.names=F)
# or to file
#write.table(cauchy_results_nodeTwo$mcmctree, "myMCMCTreeInput.tre", quote=F, row.names=F, col.names=F)




