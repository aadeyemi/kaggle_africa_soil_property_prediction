library(stats)
library(caret)
library(reshape2)
library(signal)
library(tictoc)
library(earth)
library(mda)
library(Metrics)
suppressWarnings(suppressMessages(library(doParallel)))
## Data Loading

data0 <- read.csv(unz("train.zip","training.csv"), na.strings=c("NA",""))
#testing <- read.csv(unz("test.zip","sorted_test.csv"), na.strings=c("NA",""))

a <- grep("m2379.76",names(data0))
b <- grep("m2352.76",names(data0))

data <- data0[,-seq(a,b)]

## extract infrared columns for exploration
infraredCols <- grep("^m[0-9].*?\\.[0-9].*",names(data))
infraredData <- data.frame(data[,infraredCols])

## create plot data for exploration
plotData <- data.frame(data[,infraredCols],
                       sample=seq(1,nrow(infraredData)))

## compare data from multiple sites
pdata0 = data.frame(Spectral.Index=1:ncol(plotData),
                    sample1=as.vector(t(plotData[101,])),
                    sample2=as.vector(t(plotData[876,])),
                    sample3=as.vector(t(plotData[1074,])))
pdata1 <- melt(pdata0,id.vars = "Spectral.Index")

gg1 <- ggplot(pdata1,aes(x=Spectral.Index,y=value))
gg1 <- gg1 + geom_line(aes(color=variable))
gg1 <- gg1 + coord_cartesian(ylim = c(0,2))
gg1 <- gg1 + labs(x="Spectral Index", y="Value")
plot(gg1)

## compare results from decimating site data at different rates
site <- 1
u0 = decimate(as.vector(t(plotData[site,])),1)
x0 = seq(1,length(as.vector(t(plotData[site,]))),by=1)
m0 = rep(1,length(u0))
u1 = decimate(as.vector(t(plotData[site,])),4)
x1 = seq(1,length(as.vector(t(plotData[site,]))),by=4)
m1 = rep(4,length(u1))
u2 = decimate(as.vector(t(plotData[site,])),10)
x2 = seq(1,length(as.vector(t(plotData[site,]))),by=10)
m2 = rep(10,length(u2))
u3 = decimate(as.vector(t(plotData[site,])),20)
x3 = seq(1,length(as.vector(t(plotData[site,]))),by=20)
m3 = rep(20,length(u3))

dx = c(x0,x1,x2,x3)
du = c(u0,u1,u2,u3)
dm = c(m0,m1,m2,m3)

pdata2 = data.frame(index = dx, value = du, Decimation = factor(dm))
gg2 <- ggplot(pdata2,aes(x=index,y=value))
gg2 <- gg2 + geom_line(aes(color=Decimation))
gg2 <- gg2 + coord_cartesian(ylim = c(0,2))
gg2 <- gg2 + labs(x="Spectral Index", y="Value")
plot(gg2)

## plot of signal power
## power is the sum of absolute squares of sinal samples divided by the signal length
signalPower <- apply(infraredData^2,1,sum)/ncol(infraredData)
signalPowerCut <- factor(cut(signalPower, quantile(signalPower, probs=0:4/4), include.lowest=TRUE))

outcomeCols <- grep("Sand|Ca|ph|^P$|SOC",names(data),ignore.case = TRUE)
pat <- "BSAN|BSAS|BSAV|CTI|ELEV|EVI|LSTD|LSTN|Ref1|Ref2|Ref3|Ref7|Reli|TMAP|TMFI|Depth"
varCols <- grep(pat,names(data),ignore.case = TRUE)

properties <- data.frame(data[,outcomeCols],data[,varCols],
                         Power=signalPower,Power.Quantile=signalPowerCut)

gg3 <- ggplot(properties,aes(x=Sand,y=pH))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Sand,y=P))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + coord_cartesian(xlim=c(-1.5,2.5),ylim = c(-0.5,5))
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Sand,y=Ca))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Sand,y=SOC))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Ca,y=pH))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Ca,y=P))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
gg3 <- gg3 + coord_cartesian(xlim=c(-0.5,2.5),ylim = c(-0.5,2.5))
plot(gg3)

gg3 <- ggplot(properties,aes(x=Ca,y=SOC))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
plot(gg3)

gg3 <- ggplot(properties,aes(x=pH,y=P))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + coord_cartesian(xlim=c(-2,2.5),ylim = c(-0.5,5))
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
plot(gg3)

gg3 <- ggplot(properties,aes(x=pH,y=SOC))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + coord_cartesian(xlim=c(-2,2.5),ylim = c(-1,4))
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
plot(gg3)

gg3 <- ggplot(properties,aes(x=P,y=SOC))
gg3 <- gg3 + geom_point(aes(colour=Power.Quantile),alpha = 0.5)
gg3 <- gg3 + coord_cartesian(xlim=c(-0.5,2),ylim = c(-1,4))
gg3 <- gg3 + facet_grid(Depth ~ .)
gg3 <- gg3 + scale_colour_discrete(name ="Spectral Power")
plot(gg3)


## violin plots of distributions in spectral power quantiles
gg3 <- ggplot(properties,aes(x=Power.Quantile,y=Sand))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Power.Quantile,y=log10(Ca)))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Power.Quantile,y=log10(P)))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Power.Quantile,y=pH))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Power.Quantile,y=SOC))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

## explore existence of relationships with 
## additional provided properties
outcomeCols <- grep("Sand|Ca|ph|^P$|SOC",names(data),ignore.case = TRUE)
pat <- "BSAN|BSAS|BSAV|CTI|ELEV|EVI|LSTD|LSTN|Ref1|Ref2|Ref3|Ref7|Reli|TMAP|TMFI"
varCols <- grep(pat,names(data),ignore.case = TRUE)

n           <- length(outcomeCols)*length(varCols)
outcomeName <- character(n)
varName     <- character(n)
corVal      <- numeric(n)

k <- 1
for (i in outcomeCols) {
    for (j in varCols) {
        outcomeName[k] <- names(data)[i]
        varName[k]     <- names(data)[j]
        corVal[k]      <- cor(data[,i],data[,j])
        k <- k + 1
    }
}
dmat <- cbind(outcomeName,varName,corVal)
row.names(dmat) <- c(1:nrow(dmat))
df <- data.frame(dmat)
names(df) <- c("Outcome","Variable","Correlation")
df$Correlation <- as.numeric(as.character(df$Correlation))
print(df)

gg3 <- ggplot(df,aes(x=Variable,y=Correlation))
gg3 <- gg3 + geom_point(aes(color=Outcome,size=abs(Correlation)))
plot(gg3)

## individual variables
# BSAN
gg4 <- ggplot(properties,aes(x=BSAN,y=SOC))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# BSAS
gg4 <- ggplot(properties,aes(x=BSAS,y=SOC))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# BSAV
gg4 <- ggplot(properties,aes(x=BSAV,y=Ca))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# # CTI
gg4 <- ggplot(properties,aes(x=CTI,y=Sand))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# ELEV
gg4 <- ggplot(properties,aes(x=ELEV,y=Sand))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# EVI
gg4 <- ggplot(properties,aes(x=EVI,y=pH))
gg4 <- gg4 + geom_point(alpha=0.5)
# gg4 <- gg4 +
#   coord_cartesian(xlim=c(-1,2),ylim = c(-1,3))
gg4 <- gg4 + geom_smooth()
plot(gg4)
# LSTD
gg4 <- ggplot(properties,aes(x=LSTD,y=pH))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# LSTN
gg4 <- ggplot(properties,aes(x=LSTN,y=pH))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# REF1
gg4 <- ggplot(properties,aes(x=REF1,y=SOC))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# REF2
gg4 <- ggplot(properties,aes(x=REF2,y=SOC))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# REF3
gg4 <- ggplot(properties,aes(x=REF3,y=SOC))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# REF7
gg4 <- ggplot(properties,aes(x=REF7,y=SOC))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# RELI
gg4 <- ggplot(properties,aes(x=RELI,y=Sand))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# TMAP
gg4 <- ggplot(properties,aes(x=TMAP,y=Sand))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# TMFI
gg4 <- ggplot(properties,aes(x=TMFI,y=Sand))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)

## data preparation
infraredColsData <- data.frame(data[,infraredCols])
decimatedInfrared <- data.frame(t(apply(tmp,1,decimate,4)))
outcomeCols <- grep("Sand|Ca|ph|^P$|SOC|pidn",names(data),ignore.case = TRUE)
depthCol    <- grep("Depth",names(data),ignore.case = TRUE)
finalData   <- cbind(data[,outcomeCols],Depth=data[,depthCol],decimatedInfrared)

getFinalData <- function (infraredColsData,decimationFactor) {
    decimatedInfrared <- data.frame(t(apply(infraredColsData,1,decimate,decimationFactor)))
    outcomeCols <- grep("Sand|Ca|ph|^P$|SOC|pidn",names(data),ignore.case = TRUE)
    depthCol    <- grep("Depth",names(data),ignore.case = TRUE)
    finalData   <- cbind(data[,outcomeCols],Depth=data[,depthCol],decimatedInfrared)
    finalData
}

createModelingData <- function(outcomeName,infraredColsData,decimationFactor) {
    preparedData <- getFinalData(infraredColsData,decimationFactor)
    
    if (grepl("sand",outcomeName,ignore.case = TRUE)) 
        removeCols <- grep("Ca|ph|^P$|SOC|pidn",names(preparedData),ignore.case = TRUE)
    
    else if (grepl("Ca",outcomeName,ignore.case = TRUE))
        removeCols <- grep("sand|ph|^P$|SOC|pidn", names(preparedData),ignore.case = TRUE)
    
    else if (grepl("ph",outcomeName,ignore.case = TRUE))
        removeCols <- grep("sand|Ca|^P$|SOC|pidn", names(preparedData),ignore.case = TRUE)
    
    else if (grepl("^P$",outcomeName,ignore.case = TRUE))
        removeCols <- grep("sand|Ca|ph|SOC|pidn",  names(preparedData),ignore.case = TRUE)
    
    else if (grepl("SOC",outcomeName,ignore.case = TRUE))
        removeCols <- grep("sand|Ca|ph|^P$|pidn",  names(preparedData),ignore.case = TRUE)
    
    procData <- preparedData[,-removeCols]
}

## select pertinent columns for each outcome
preparedData   <- getFinalData(infraredColsData,10)
removeColsSand <- grep("Ca|ph|^P$|SOC|pidn",   names(preparedData),ignore.case = TRUE)
removeColsCa   <- grep("sand|ph|^P$|SOC|pidn", names(preparedData),ignore.case = TRUE)
removeColsPh   <- grep("sand|Ca|^P$|SOC|pidn", names(preparedData),ignore.case = TRUE)
removeColsP    <- grep("sand|Ca|ph|SOC|pidn",  names(preparedData),ignore.case = TRUE)
removeColsSoc  <- grep("sand|Ca|ph|^P$|pidn",  names(preparedData),ignore.case = TRUE)

## create datasets containing pertinent data for each outcome variable
procDataSand <- preparedData[,-removeColsSand]
procDataCa   <- preparedData[,-removeColsCa]
procDataPh   <- preparedData[,-removeColsPh]
procDataP    <- preparedData[,-removeColsP]
procDataSoc  <- preparedData[,-removeColsSoc]

## create data partition vectors
inTrainSand <- createDataPartition(y=procDataSand$Sand, p=0.75, list=FALSE)
inTrainCa   <- createDataPartition(y=procDataCa$Ca,     p=0.75, list=FALSE)
inTrainPh   <- createDataPartition(y=procDataPh$pH,     p=0.75, list=FALSE)
inTrainP    <- createDataPartition(y=procDataP$P,       p=0.75, list=FALSE)
inTrainSoc  <- createDataPartition(y=procDataSoc$SOC,   p=0.75, list=FALSE)

## partition data
trainingSand <- procDataSand[inTrainSand,]
testingSand  <- procDataSand[-inTrainSand,]

## train using random forests
## 1) Sand

## with decimation of 10 times

tic()
suppressWarnings(suppressMessages(library(doParallel)))
#cl <- makeCluster(detectCores())
cl <- makeCluster(2)
registerDoParallel(cl)
# modelFitSand.bagEarth.1     <- bagEarth(Sand ~ ., data=trainingSand) # better result
# modelFitSand.bagEarth.10    <- bagEarth(Sand ~ ., data=trainingSand)
# modelFitSand.lm.10          <- train(Sand ~ ., method="lm",          data=trainingSand)
# modelFitSand.bagEarth.10    <- train(Sand ~ ., method="bagEarth",    data=trainingSand)
# modelFitSand.bagEarthGCV.10 <- train(Sand ~ ., method="bagEarthGCV", data=trainingSand)
# modelFitSand.earth.10       <- train(Sand ~ ., method="earth",       data=trainingSand)
# modelFitSand.gcvEarth.10    <- train(Sand ~ ., method="gcvEarth",    data=trainingSand)
# modelFitSand.rf.10 <- train(Sand ~ ., method="rf", data=trainingSand)
stopCluster(cl)
timed <- toc()


## with decimation of 4 times
tic()
suppressWarnings(suppressMessages(library(doParallel)))
cl <- makeCluster(2)
registerDoParallel(cl)
modelFitSand.lm.4          <- train(Sand ~ ., method="lm",          data=trainingSand)
timed <- toc()
modelFitSand.bagEarth.4    <- train(Sand ~ ., method="bagEarth",    data=trainingSand)
timed <- toc()
modelFitSand.bagEarthGCV.4 <- train(Sand ~ ., method="bagEarthGCV", data=trainingSand)
timed <- toc()
modelFitSand.earth.4       <- train(Sand ~ ., method="earth",       data=trainingSand)
timed <- toc()
modelFitSand.gcvEarth.4    <- train(Sand ~ ., method="gcvEarth",    data=trainingSand)
timed <- toc()
modelFitSand.rf.4          <- train(Sand ~ ., method="rf", data=trainingSand)
stopCluster(cl)
timed <- toc()


## plots with decimation of 10
pred.lm.10 <- predict(modelFitSand.lm.10 ,newdata = testingSand)
predT.lm.10 <- predict(modelFitSand.lm.10 ,newdata = trainingSand)
plot(pred.lm.10 ,testingSand$Sand)

pred.bagEarth.10 <- predict(modelFitSand.bagEarth.10 ,newdata = testingSand)
predT.bagEarth.10 <- predict(modelFitSand.bagEarth.10 ,newdata = trainingSand)
plot(pred.bagEarth.10 ,testingSand$Sand)

pred.bagEarthGCV.10 <- predict(modelFitSand.bagEarthGCV.10 ,newdata = testingSand)
predT.bagEarthGCV.10 <- predict(modelFitSand.bagEarthGCV.10 ,newdata = trainingSand)
plot(pred.bagEarthGCV.10 ,testingSand$Sand)

pred.earth.10 <- predict(modelFitSand.earth.10 ,newdata = testingSand)
predT.earth.10 <- predict(modelFitSand.earth.10 ,newdata = trainingSand)
plot(pred.earth.10 ,testingSand$Sand)

pred.gcvEarth.10 <- predict(modelFitSand.gcvEarth.10 ,newdata = testingSand)
predT.gcvEarth.10 <- predict(modelFitSand.gcvEarth.10 ,newdata = trainingSand)
plot(pred.gcvEarth.10 ,testingSand$Sand)

pred.rf.10 <- predict(modelFitSand.rf.10 ,newdata = testingSand)
predT.rf.10 <- predict(modelFitSand.rf.10 ,newdata = trainingSand)
plot(pred.rf.10 ,testingSand$Sand)

# plot errors with respect to regression method used
n      <- 6
method <- c("lm","bagEarth","bagEarthGCV","earth","gcvEarth","rf")
is_error  <- c(rmse(trainingSand$Sand,predT.lm.10),
               rmse(trainingSand$Sand,predT.bagEarth.10),
               rmse(trainingSand$Sand,predT.bagEarthGCV.10),
               rmse(trainingSand$Sand,predT.earth.10),
               rmse(trainingSand$Sand,predT.gcvEarth.10),
               rmse(trainingSand$Sand,predT.rf.10))
os_error  <- c(rmse(testingSand$Sand,pred.lm.10),
               rmse(testingSand$Sand,pred.bagEarth.10),
               rmse(testingSand$Sand,pred.bagEarthGCV.10),
               rmse(testingSand$Sand,pred.earth.10),
               rmse(testingSand$Sand,pred.gcvEarth.10),
               rmse(testingSand$Sand,pred.rf.10))

emat            <- cbind(method,os_error,is_error)
row.names(emat) <- c(1:nrow(emat))
edf             <- data.frame(emat)
names(edf)      <- c("Method","OutOfSampleRMSE","InSampleRMSE")
edf$OutOfSampleRMSE <- as.numeric(as.character(edf$OutOfSampleRMSE))
edf$InSampleRMSE    <- as.numeric(as.character(edf$InSampleRMSE))

print(edf)

gg5 <- ggplot(edf,aes(x=InSampleRMSE,y=OutOfSampleRMSE))
gg5 <- gg5 + geom_point(aes(fill=Method,stroke=1,shape=Method),size=3)
gg5 <- gg5 + geom_abline(intercept=0, slope=1, linetype=2)
gg5 <- gg5 + coord_cartesian(xlim=c(0,1),ylim = c(0,1.03))
gg5 <- gg5 + labs(x="In-sample Error (RMSE)", y="Out-of-sample Error (RMSE)")
gg5 <- gg5 + scale_colour_brewer(palette="Set1")
plot(gg5)

## plots with desimation of 4
pred.lm.4 <- predict(modelFitSand.lm.4 ,newdata = testingSand)
predT.lm.4 <- predict(modelFitSand.lm.4 ,newdata = trainingSand)
plot(pred.lm.4 ,testingSand$Sand)

pred.bagEarth.4 <- predict(modelFitSand.bagEarth.4 ,newdata = testingSand)
predT.bagEarth.4 <- predict(modelFitSand.bagEarth.4 ,newdata = trainingSand)
plot(pred.bagEarth.4 ,testingSand$Sand)

pred.bagEarthGCV.4 <- predict(modelFitSand.bagEarthGCV.4 ,newdata = testingSand)
predT.bagEarthGCV.4 <- predict(modelFitSand.bagEarthGCV.4 ,newdata = trainingSand)
plot(pred.bagEarthGCV.4 ,testingSand$Sand)

pred.earth.4 <- predict(modelFitSand.earth.4 ,newdata = testingSand)
predT.earth.4 <- predict(modelFitSand.earth.4 ,newdata = trainingSand)
plot(pred.earth.4 ,testingSand$Sand)

pred.gcvEarth.4 <- predict(modelFitSand.gcvEarth.4 ,newdata = testingSand)
predT.gcvEarth.4 <- predict(modelFitSand.gcvEarth.4 ,newdata = trainingSand)
plot(pred.gcvEarth.4 ,testingSand$Sand)

pred.rf.4 <- predict(modelFitSand.rf.4 ,newdata = testingSand)
predT.rf.4 <- predict(modelFitSand.rf.4 ,newdata = trainingSand)
plot(pred.rf.4 ,testingSand$Sand)

# plot errors with respect to regression method used
n      <- 6
method <- c("lm","bagEarth","bagEarthGCV","earth","gcvEarth","rf")
is_error  <- c(rmse(trainingSand$Sand,predT.lm.4),
               rmse(trainingSand$Sand,predT.bagEarth.4),
               rmse(trainingSand$Sand,predT.bagEarthGCV.4),
               rmse(trainingSand$Sand,predT.earth.4),
               rmse(trainingSand$Sand,predT.gcvEarth.4),
               rmse(trainingSand$Sand,predT.rf.4))
os_error  <- c(rmse(testingSand$Sand,pred.lm.4),
               rmse(testingSand$Sand,pred.bagEarth.4),
               rmse(testingSand$Sand,pred.bagEarthGCV.4),
               rmse(testingSand$Sand,pred.earth.4),
               rmse(testingSand$Sand,pred.gcvEarth.4),
               rmse(testingSand$Sand,pred.rf.4))

emat            <- cbind(method,os_error,is_error)
row.names(emat) <- c(1:nrow(emat))
edf             <- data.frame(emat)
names(edf)      <- c("Method","OutOfSampleRMSE","InSampleRMSE")
edf$OutOfSampleRMSE <- as.numeric(as.character(edf$OutOfSampleRMSE))
edf$InSampleRMSE    <- as.numeric(as.character(edf$InSampleRMSE))

print(edf)

gg6 <- ggplot(edf,aes(x=InSampleRMSE,y=OutOfSampleRMSE))
gg6 <- gg6 + geom_point(aes(fill=Method,stroke=1,shape=Method),size=3)
gg6 <- gg6 + geom_abline(intercept=0, slope=1, linetype=2)
gg6 <- gg6 + coord_cartesian(xlim=c(0,1),ylim = c(0,1.03))
gg6 <- gg6 + labs(x="In-sample Error (RMSE)", y="Out-of-sample Error (RMSE)")
gg6 <- gg6 + scale_colour_brewer(palette="Set1")
plot(gg6)

## principal components with decimation of 4 times
preProc        <- preProcess(decimatedInfrared,method="pca",pcaComp=ncol(decimatedInfrared))
trainPC        <- predict(preProc,decimatedInfrared)

finalData      <- cbind(data[,outcomeCols],Depth=data[,depthCol],trainPC)

preparedData   <- finalData
removeColsSand <- grep("Ca|ph|^P$|SOC|pidn",names(preparedData),ignore.case = TRUE)
procDataSand   <- preparedData[,-removeColsSand]
inTrainSand    <- createDataPartition(y=procDataSand$Sand, p=0.75, list=FALSE)
trainingSand   <- procDataSand[inTrainSand,]
testingSand    <- procDataSand[-inTrainSand,]

tic()
cl <- makeCluster(3)
registerDoParallel(cl)
modelFitSand.rf.4.pc <- train(Sand ~ ., method="rf", data=trainingSand)
stopCluster(cl)
timed <- toc()

pred.rf.4.pc <- predict(modelFitSand.rf.4.pc ,newdata = testingSand)
predT.rf.4.pc <- predict(modelFitSand.rf.4.pc ,newdata = trainingSand)
plot(pred.rf.4.pc ,testingSand$Sand)
plot(predT.rf.4.pc ,trainingSand$Sand)
rmse(testingSand$Sand,pred.rf.4.pc)