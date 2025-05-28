#setting working directory
setwd("~/Desktop/big data/big data")
#loading temp data
temperature <- read.csv("temperatureTimeSeries.csv", header=FALSE)
#plotting temp data
plot(temperature)
#making temp plot neater, flipping x axis, changing dots to line, renaming axis, smoothing it out
finaltemp <- smooth(smooth(temperature$V2))
plot(finaltemp)
plot(temperature$V1, finaltemp, xlab='Time')
plot(temperature$V1,finaltemp,  xlab='Time', xlim=c(220,0))
plot(temperature$V1,finaltemp, xlab='Time', xlim=c(220,0), type = 'l')
plot(temperature$V1,finaltemp, xlab='Time (mya)', xlim=c(220,0), type = 'l', ylab = "Temperature (measured by proxy)")
#loading sea level data
seaLevel <- read.csv("seaLevelTimeSeries.csv",header=TRUE)
#plotting sea level data, flipping axis, changing dots to line
plot(seaLevel)
plot(seaLevel, xlim=c(250,0))
plot(seaLevel, xlim=c(250,0), type = 'l')
plot(seaLevel, xlim=c(250,0), type = 'l', ylab = "Sea Level (measured by proxy)", xlab = "Time (mya)")
#loading libraries: phytools and strap
library(phytools)
library(strap)
#loading phylogeny data
tree <- read.tree("fossilCrocPhylogeny.tre")
#plotting phylogenic tree 
plot(tree)
#checking tree data
head(tree$edge)
head(tree$edge.length)
head(tree$Nnode)
head(tree$tip.label)
head(tree$root.edge)
#loading habitat data
habitatdata <- read.csv("HabitatData.csv", header=T, stringsAsFactors = FALSE)
#checking habitat data
head(habitatdata)
#subsetting terrestrial taxa
TerrestrialTaxa <- subset(habitatdata, habitatdata$Habitat=='Terrestrial')$Taxon
#subsetting marine taxa 
MarineTaxa <- subset(habitatdata, habitatdata$Habitat=='Marine')$Taxon
#checking that subsetting data worked correctly 
head(MarineTaxa)
head(TerrestrialTax)
#setting for plotting terrestrial data phylo 
treeT <- keep.tip(tree, TerrestrialTaxa)
#checking data for phylo terrestrial 
treeT
#plotting terrestrial data and making it neat 
plot(treeT)
plot(treeT, cex=0.2)
#saving terestrial data pylo in tre file 
write.tree(treeT, file = 'TerrestrialTaxa.tre')
#Terrestrial data: find root for plotting first by checking ages of nodes (nodeHeights)
lengths <- nodeHeights(treeT)
#checking lengths
head(lengths)
#setting largest number aka root node
root.time <- max(lengths)
#setting root for plotting
treeT$root.time <- root.time
#grabbing OTUs (Operational Taxonomic Units = taxa)
all_otus <- treeT$tip.label
#Creating an empty matrix containing the terrestrial taxa as required for strap
all_otudates <- matrix(0, nrow = length(all_otus), ncol=2)
#Turning the matrix into a data frame
all_otudates <- data.frame(all_otudates)
#set the row names to the terrestrial taxa (OTUs)
row.names(all_otudates) <- all_otus
# set column names to FAD (First Appearance Datum) and LAD (Last Appearance Datum)
colnames(all_otudates) <- c('FAD','LAD')
#graphing terrestial taxa against geological timescale: final graph
geoscalePhylo(treeT,ages=all_otudates, cex.tip=0.1, lwd=1, quat.rm=T, units=c("Period", "Epoch"), boxes="Epoch")
#repeating for marine data,since the values set ie. root.time, all_otus, etc aren't used again, i just overwrote  them for the marine data
lengths <- nodeHeights(treeM)
root.time <- max(lengths)
treeM$root.time <- root.time
all_otus <- treeM$tip.label
all_otudates <- matrix(0, nrow = length(all_otus), ncol=2)
all_otudates <- data.frame(all_otudates)
row.names(all_otudates) <- all_otus
colnames(all_otudates) <- c('FAD','LAD')
#graphing marine taxa phylogeny against geological timescale: final graph 
geoscalePhylo(treeM,ages=all_otudates, cex.tip=0.1, lwd=1, quat.rm=T, units=c("Period", "Epoch"), boxes="Epoch")
#loading libraries
library(BAMMtools)
library(phytools)
#extracting rates through time
edata <- getEventData(tree, eventdata ='fossilCrocDiversificationData.txt', burnin=0.1)
#extracting terrestrial subtree
streeTerrestrial <- subtreeBAMM(edata, tips=TerrestrialTaxa)
#extracting rates through time for terrestrial taxa
rtt_T <- getRateThroughTimeMatrix(streeTerrestrial)
#checking rtt_T and its elements 
summary(rtt_T)
rtt_T$times
rtt_T$lambda
rtt_T$mu
rtt_T$type
#plotting for terrestrial speciaton
plotRateThroughTime(streeTerrestrial, ratetype='speciation', avgCol="red", ylim=c(0,0.5), cex.axis=2, intervalCol='red', intervals=c(0.05, 0.95), opacity=0.1)
#plotting for terrestrial extinction 
plotRateThroughTime(streeTerrestrial, ratetype='extinction', avgCol="red", ylim=c(0,0.5), cex.axis=2, intervalCol='red', intervals=c(0.05, 0.95), opacity=0.1)
#extracting marine subtree
rtt_M <- getRateThroughTimeMatrix(streeTerrestrial)
#checking rtt_M and its elements 
summary(rtt_M)
rtt_M$times
rtt_M$lambda
rtt_M$mu
rtt_M$type
#plotting for marine speciaton
plotRateThroughTime(streeMarine, ratetype='speciation', avgCol="blue", ylim=c(0,0.5), cex.axis=2, intervalCol='blue', intervals=c(0.05, 0.95), opacity=0.1)
#plotting for marine extinction
plotRateThroughTime(rtt_M, ratetype='extinction', avgCol="blue", ylim=c(0,0.5), cex.axis=2, intervalCol='blue', intervals=c(0.05, 0.95), opacity=0.1)
#got an error for marine extinction so i did this 
which(is.na(rtt_M$mu))
rtt_M$mu[is.na(rtt_M$mu)] <- 0
plotRateThroughTime(rtt_M, ratetype='extinction', avgCol="blue", ylim=c(0,0.5), cex.axis=2, intervalCol='blue', intervals=c(0.05, 0.95), opacity=0.1)
#loading libraries for correlation
library(ggplot2)
#loading DCCA function for analysis
DCCA <- function(x,y,s){
  xx<-cumsum(x)
  yy<-cumsum(y)
  t<-1:length(xx)
  F2sj_xy<-runif(floor(length(xx)/s))
  F2sj_xx<-F2sj_xy
  F2sj_yy<-F2sj_xy
  for(ss in seq(1,(floor(length(xx)/s)*s),by=s)){
    F2sj_xy[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_xx[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_yy[(ss-1)/s+1]<-sum((summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
  }
  rho<-mean(F2sj_xy)/sqrt(mean(F2sj_xx)*mean(F2sj_yy))
  return(c(rho,1/sqrt(length(xx)),1-pnorm(abs(rho),mean=0,sd=1/sqrt(length(xx)))))
}
#flipping time for terrestrial 
times = abs(rtt_T$times-max(rtt_T$times))
# Calculate number of simulations in our diversification data
numberOfSims = length(rtt_T$lambda)/length(rtt_T$times)
#setting number of samples 
numberOfSamples = 100
#store the correlation coefficients
cors_temp_T <- rep(NA, numberOfSamples)
cors_seaLevel_T <- rep(NA, numberOfSamples)
#setting random seed
set.seed(1)
# generating the random sample
samples = sample(1:numberOfSims, numberOfSamples, replace = FALSE)
#set count
count = 1
#running correlation
 for (i in 1:numberOfSims ) {
  +     if (i %in% samples){
  + interpdiv = approx(times, rtt_T$lambda[i,], temperature$V1, method='linear', rule=1)
  + end = which(is.na(interpdiv$y))
  + if (length(end) == 0) {
  +     div_rates = interpdiv$y
  +     ft = finaltemp
  + } else {
  +     div_rates = interpdiv$y[-end]
  +     ft = finaltemp[-end]
  + }
  + 
  + c = DCCA(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/10)
  + cors_temp_T[count] = c[1]
  + count = count+1

#plotting cor. temp speciation histogram for terrestrial taxa
plot=qplot(cors_temp_T, geom="histogram", bins=30)
plot
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot2
ggsave(plot2,file="SpeciationTerrestrialTemperature.pdf")
plot3 = plot2 + xlab("Correlation Coefficient (Temperature and Speciation)")
plot3
ggsave(plot2,file="SpeciationTerrestrialTemperature.pdf")
#stat test
sink(file="SpeciationTerrestrialTemperatureStats.txt")
print(summary(cors_temp_T))
print(quantile(cors_temp_T, c(0.025, 0.975)))
stats <- (wilcox.test(cors_temp_T, mu=0.0, paired = FALSE))
stats$p.value
sink()
#running correlation for sea level and speciation for terrestial taxa
for (i in 1:numberOfSims ) {
+     if (i %in% samples){
+         interpdiv = approx(times, rtt_T$lambda[i,], seaLevel$Age, method='linear', rule=1)
+         end = which(is.na(interpdiv$y))
+ if (length(end) == 0) {
+    div_rates = interpdiv$y
+    ft = seaLevel$SL
+} else {
+    div_rates = interpdiv$y[-end]
+    ft = seaLevel$SL[-end]
+ }
+ c = DCCA(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/10)
+cors_seaLevel_T[count] = c[1]
+count = count+1
+ }
+ }
#plotting histogram cor. speciation sea level terrestrial histogram
plotS1=qplot(cors_seaLevel_T, geom="histogram", bins=30)
plotS1
plotS2 = plotS1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
plotS2
plotS3 = plotS2 + xlab("Correlation Cofficient (Sea Level and Speciation)")
plotS3
ggsave(plotS3,file="SpeciationSeaLevelT.pdf")
#setting up for stat test
sink(file="SpeciationTerrestrialSeaLevelStats.txt")
print(summary(cors_seaLevel_T))
print(quantile(cors_seaLevel_T, c(0.025, 0.975))) 
stats <- (wilcox.test(cors_temp_T, mu=0.0, paired = FALSE))
stats$p.value
sink()

#setting up for marine data 
DCCA <- function(x,y,s){
  xx<-cumsum(x)
  yy<-cumsum(y)
  t<-1:length(xx)
  F2sj_xy<-runif(floor(length(xx)/s))
  F2sj_xx<-F2sj_xy
  F2sj_yy<-F2sj_xy
  for(ss in seq(1,(floor(length(xx)/s)*s),by=s)){
    F2sj_xy[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_xx[(ss-1)/s+1]<-sum((summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(xx[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
    F2sj_yy[(ss-1)/s+1]<-sum((summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals)*(summary(lm(yy[ss:(ss+s-1)]~t[ss:(ss+s-1)]))$residuals))/(s-1)
  }
  rho<-mean(F2sj_xy)/sqrt(mean(F2sj_xx)*mean(F2sj_yy))
  return(c(rho,1/sqrt(length(xx)),1-pnorm(abs(rho),mean=0,sd=1/sqrt(length(xx)))))
}
#flipping time, i reset the values cuz they wont be used again 
times = abs(rtt_M$times-max(rtt_M$times))
numberOfSims = length(rtt_M$lambda)/length(rtt_M$times)
#correlation analysis marine temp and speciation 
numberOfSamples = 100
cors_temp_M <- rep(NA, numberOfSamples)
cors_seaLevel_M <- rep(NA, numberOfSamples)
#SETTING RANDOM SEED
set.seed(1)
#GENERATING RANDOM DAMPLE
samples = sample(1:numberOfSims, numberOfSamples, replace = FALSE)
#SETTING COUNT TO 1
count = 1
#RUNNING ANAYLSIS FOR MARINE SPECIATION TEMP
> for (i in 1:numberOfSims ) {
     if (i %in% samples){
 interpdiv = approx(times, rtt_M$lambda[i,], temperature$V1, method='linear', rule=1)
 end = which(is.na(interpdiv$y))
 if (length(end) == 0) {
     div_rates = interpdiv$y
     ft = finaltemp
 } else {
     div_rates = interpdiv$y[-end]
     ft = finaltemp[-end]
 }
 c = DCCA(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/10)
 cors_temp_M[count] = c[1]
 count = count+1
#PLOTTING HISTOGRAM FOR MATINE SPECIATION AND TEMP
plot=qplot(cors_temp_T, geom="histogram", bins=30)
plot
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot2
plot3 = plot2 + xlab("Correlation Coefficient (Temperature and Speciation)")
plot3
ggsave(plot3,file="SpeciationMarineTemperature.pdf")
#STAT ANALYSIS Wilcox
sink(file="SpeciationMarineTemperatureStats.txt")
print(summary(cors_temp_M))
print(quantile(cors_temp_M, c(0.025, 0.975)))
stats <- (wilcox.test(cors_temp_M, mu=0.0, paired = FALSE))
stats$p.value
sink()
#SETTING UP FOR MARINE SPECIATION AND SEA LEVEL ANALYSIS
count = 1
for (i in 1:numberOfSims ) {
  +     if (i %in% samples){
    +         interpdiv = approx(times, rtt_M$lambda[i,], seaLevel$Age, method='linear', rule=1)
    +         end = which(is.na(interpdiv$y))
    + if (length(end) == 0) {
      +    div_rates = interpdiv$y
      +    ft = seaLevel$SL
      +} else {
        +    div_rates = interpdiv$y[-end]
        +    ft = seaLevel$SL[-end]
        + }
    + c = DCCA(as.numeric(unlist(div_rates)),as.numeric(unlist(ft)),length(ft)/10)
    +cors_seaLevel_M[count] = c[1]
    +count = count+1
    + }
  + }
#plotting graph 
plot=qplot(cors_seaLevel_M, geom="histogram", bins=30)
plot
plot2 = plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot2
plot3 = plot2 + xlab("Correlation Coefficient (Sea Level and Speciation)")
plot3
ggsave(plot3,file="SpeciationMarineSeaLevel.pdf")
#stat test wilcox 
sink(file="SpeciationMarineSeaLevelStats.txt")
print(summary(cors_seaLevel_M))
print(quantile(cors_seaLevel_M, c(0.025, 0.975))) 
stats <- (wilcox.test(cors_temp_M, mu=0.0, paired = FALSE))
stats$p.value
sink()










