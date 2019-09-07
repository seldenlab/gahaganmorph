# Robert Z. Selden, Jr.
# Comparing Gahagan Bifaces from Gahagan, Mounds Plantation, and George C Davis
# March 2018
library(devtools)
devtools::install_github("geomorphR/geomorph", ref = "Stable")
library(geomorph)

# set working directory
setwd(getwd())
source('readmulti.csv.r')

# Read in .csv files from data folder
setwd("./data")  # move to subdirectory to read landmark data
filelist <- list.files(pattern = ".csv")
coords<-readmulti.csv(filelist)
setwd("../")  # back out of the subdirectory and into main dirctory

# read-in qualitative data and linear measures
qdata<-read.csv("qdata.csv",header=TRUE, row.names=1)
qdata<-qdata[match(dimnames(coords)[[3]],rownames(qdata)),]  # re-order rows of qdata to match order in landmark dataset

# GPA
Y.gpa<-gpagen(coords, PrinAxes = TRUE)  #GPA of entire shape: NOTE: not a symmetry analysis
plot(Y.gpa)  

# note: this section not needed in 3d - only to pull out X,Y and rotate PCA plots
Y.gpa$coords<-simplify2array(lapply(1:dim(coords)[3], function(j) Y.gpa$coords[,1:2,j] )) #note: ignore this line if multi-curve data is used
  rot<-matrix(c(1,0,0,-1),ncol=2, byrow=T)
Y.gpa$coords<-simplify2array(lapply(1:dim(coords)[3], function(j) Y.gpa$coords[,,j] %*% rot )) 
plotAllSpecimens(Y.gpa$coords)

# statistical analysis
gdf<-geomorph.data.frame(shape=Y.gpa$coords, size=Y.gpa$Csize, site=qdata$SiteName)
# note: add variables to a geomorph data frame column by column - do not 'nest' data frames within each other

  # exploratory statistics
 plotTangentSpace(Y.gpa$coords,label = TRUE,warpgrids = TRUE)
 
#Plot PCA without warpgrids, and generate warpgrids by hand
     PCA<-plotTangentSpace(Y.gpa$coords, warpgrids=FALSE, groups = qdata$SiteName)
     M<-mshape(Y.gpa$coords)
     #shapes along PC1 using Procrustes values -0.1 and 0.1, corresponding to the x-axis values in the plot window
     preds<-shape.predictor(Y.gpa$coords, x= NULL, Intercept = FALSE, pred1 = -0.1, pred2 = 0.1)
     plotRefToTarget(M, preds$pred)
     #shapes at the minima and maxima of PC1 and PC2 using min and max, corresponding to the x/y values in the plot
     PC<-PCA$pc.scores[,1] #x-axis
     preds<-shape.predictor(Y.gpa$coords, x=PC, Intercept = FALSE, pred1 = min(PC), pred2 = max(PC))
     plotRefToTarget(M, preds$pred1) #minimum x-axis
     plotRefToTarget(M, preds$pred2) #maximum x-axis
     PC<-PCA$pc.scores[,2] #y-axis
     preds<-shape.predictor(Y.gpa$coords, x=PC, Intercept = FALSE, pred1 = min(PC), pred2 = max(PC))
     plotRefToTarget(M, preds$pred1) #minimum y-axis
     plotRefToTarget(M, preds$pred2) #maximum y-axis

  # inferential statistics
# allometry: does shape change with size?  
results.allom<-procD.lm(shape~size,data=gdf,iter = 9999) #allometry = regression
summary(results.allom)
plot(results.allom,type="regression" ,predictor=(gdf$size),reg.type = "RegScore")  #could also use reg.type = "CRC"

# PredLine allometry--compare allometry between types
botAllometry<-procD.allometry(shape~size, ~site, data=gdf, iter=9999, RRPP=TRUE)
summary(botAllometry)
plot(botAllometry, method = "PredLine", warpgrids = TRUE)

# anova: does shape differ by group?
res.anova<-procD.lm(shape~site,data=gdf,iter = 9999)  
summary(res.anova)

# anova pt 2: which groups differ from one another?
res.pair<-advanced.procD.lm(shape~site,f2=~1,groups=~site,data=gdf,iter = 9999) #pairwise comparisons of shape between groups in 'Type'
summary(res.pair)
# gives pairise differences in group means and p-values

#Disparity: do any groups display greater shape variation among individuals relative to other groups?
morphol.disparity(shape~size,groups=~site,data=gdf,iter = 9999) #table gives amount of disparity(variation) per group, pairwise differences, pairwise p-values

# graphical depiction of results (stats and GMM plots)

# plot group means for 'SiteName'
PCA<-plotTangentSpace(Y.gpa$coords,label = TRUE,warpgrids = TRUE,groups = qdata$SiteName,legend=TRUE)  #PCA plot groups colored
# PCA barplot
pvar <- (PCA$sdev^2)/(sum(PCA$sdev^2))
names(pvar) <- seq(1:length(pvar))
barplot(pvar, xlab= "Principal Components", ylab = "% Variance")

mean<-mshape(Y.gpa$coords)
lsMeans<-arrayspecs(res.pair$LS.means,p=46,k=2)
#plot mean shapes by site/context
plotRefToTarget(mean,lsMeans[,,1]) #"GCD - F119"
plotRefToTarget(mean,lsMeans[,,2]) #"GCD - F134"
plotRefToTarget(mean,lsMeans[,,3]) #"GM - burial pit 2"
plotRefToTarget(mean,lsMeans[,,4]) #"GM - burial pit 3"
plotRefToTarget(mean,lsMeans[,,5]) #"MP - burial pit 2"
plotRefToTarget(mean,lsMeans[,,6]) #"MP - burial pit 5"
#compare mean shapes by site/context::
plotRefToTarget(lsMeans[,,1],lsMeans[,,2], method="vector",mag=1) #GCD-F119 v GCD-F134
plotRefToTarget(lsMeans[,,1],lsMeans[,,3], method="vector",mag=1) #GCD-F119 v GM BP2
plotRefToTarget(lsMeans[,,1],lsMeans[,,4], method="vector",mag=1) #GCD-F119 v GM BP3
plotRefToTarget(lsMeans[,,1],lsMeans[,,5], method="vector",mag=1) #GCD-F119 v MP BP2
plotRefToTarget(lsMeans[,,1],lsMeans[,,6], method="vector",mag=1) #GCD-F119 v MP BP5
plotRefToTarget(lsMeans[,,2],lsMeans[,,3], method="vector",mag=1) #GCD-F134 v GM BP2
plotRefToTarget(lsMeans[,,2],lsMeans[,,4], method="vector",mag=1) #GCD-F134 v GM BP3
plotRefToTarget(lsMeans[,,2],lsMeans[,,5], method="vector",mag=1) #GCD-F134 v MP BP2
plotRefToTarget(lsMeans[,,2],lsMeans[,,6], method="vector",mag=1) #GCD-F134 v MP BP5
plotRefToTarget(lsMeans[,,3],lsMeans[,,4], method="vector",mag=1) #GM BP2 v GM BP3
plotRefToTarget(lsMeans[,,3],lsMeans[,,5], method="vector",mag=1) #GM BP2 v MP BP2
plotRefToTarget(lsMeans[,,3],lsMeans[,,6], method="vector",mag=1) #GM BP2 v MP BP5
plotRefToTarget(lsMeans[,,4],lsMeans[,,5], method="vector",mag=1) #GM BP3 v MP BP2
plotRefToTarget(lsMeans[,,4],lsMeans[,,6], method="vector",mag=1) #GM BP3 v MP BP5
plotRefToTarget(lsMeans[,,5],lsMeans[,,6], method="vector",mag=1) #MP BP2 v MP BP5
# morphological integration::
land.gps<-c("A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A",
            "B","B","B","B","B","B","B","B","B","A","A","A","A","A","A","A","A","A",
            "A","A","A","A","A","A","A","A","A")
IT<-integration.test(Y.gpa$coords, partition.gp = land.gps, iter = 9999, print.progress = FALSE)
summary(IT)
plot(IT)

#end of script