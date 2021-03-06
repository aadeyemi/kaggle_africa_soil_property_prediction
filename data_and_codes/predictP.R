library(stats)
library(reshape2)
library(signal)
library(tictoc)
library(earth)
library(mda)
library(Metrics)
suppressWarnings(suppressMessages(library(caret)))
suppressWarnings(suppressMessages(library(doParallel)))

## Load Data
data0 <- read.csv(unz("train.zip","training.csv"), na.strings=c("NA",""))

# Remove spectra CO2 bands 
a <- grep("m2379.76",names(data0))
b <- grep("m2352.76",names(data0))

data <- data0[,-seq(a,b)]

## Create spectra-only data subset
infraredCols <- grep("^m[0-9].*?\\.[0-9].*",names(data))
infraredData <- data.frame(data[,infraredCols])


## data preparation
# create final data frame based on provided decimation factor
getFinalData <- function (infraredColsData,decimationFactor) {
    decimatedInfrared <- data.frame(t(apply(infraredColsData,1,decimate,decimationFactor)))
    outcomeCols <- grep("Sand|Ca|ph|^P$|SOC|pidn",names(data),ignore.case = TRUE)
    depthCol    <- grep("Depth",names(data),ignore.case = TRUE)
    finalData   <- cbind(data[,outcomeCols],Depth=data[,depthCol],decimatedInfrared)
    finalData
}

# create modeling data for specified outcome variable
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

##------------------------------------------------------------------
## load model for P prediction given data decimation factor of 50
##------------------------------------------------------------------

load(file = "saved_models/fitP.lm.50.RData")
load(file = "saved_models/fitP.bagEarth.50.RData")
load(file = "saved_models/fitP.bagEarthGCV.50.RData")
load(file = "saved_models/fitP.earth.50.RData")
load(file = "saved_models/fitP.gcvEarth.50.RData")
load(file = "saved_models/fitP.rf.50.RData")

##------------------------------------------------------------------
## load model for P prediction given data decimation factor of 10
##------------------------------------------------------------------

load(file = "saved_models/fitP.lm.10.RData")
load(file = "saved_models/fitP.bagEarth.10.RData")
load(file = "saved_models/fitP.bagEarthGCV.10.RData")
load(file = "saved_models/fitP.earth.10.RData")
load(file = "saved_models/fitP.gcvEarth.10.RData")
load(file = "saved_models/fitP.rf.10.RData")

##------------------------------------------------------------------
## load model for P prediction given data decimation factor of 5
##------------------------------------------------------------------

load(file = "saved_models/fitP.lm.5.RData")
load(file = "saved_models/fitP.bagEarth.5.RData")
load(file = "saved_models/fitP.bagEarthGCV.5.RData")
load(file = "saved_models/fitP.earth.5.RData")
load(file = "saved_models/fitP.gcvEarth.5.RData")
load(file = "saved_models/fitP.rf.5.RData")

##------------------------------------------------------------------
## PLOTS
##------------------------------------------------------------------
getRmsVector <- function(fit.lm,fit.bagEarth,
                         fit.bagEarthGCV,fit.earth,
                         fit.gcvEarth,fit.rf,
                         outcome,infrared.data,
                         decim.factor,in.train) {
    
    proc.data <- createModelingData(outcome,infrared.data,decim.factor)
    
    training <- proc.data[in.train,]
    testing  <- proc.data[-in.train,]
    
    in.sample.data      <- training[,1]
    out.sample.data     <- testing[,1]
    
    in.sample.lm.pred   <- predict(fit.lm,newdata=training)
    out.sample.lm.pred  <- predict(fit.lm,newdata=testing)
    
    in.sample.bagEarth.pred   <- predict(fit.bagEarth,newdata=training)
    out.sample.bagEarth.pred  <- predict(fit.bagEarth,newdata=testing)
    
    in.sample.bagEarthGCV.pred   <- predict(fit.bagEarthGCV,newdata=training)
    out.sample.bagEarthGCV.pred  <- predict(fit.bagEarthGCV,newdata=testing)
    
    in.sample.earth.pred   <- predict(fit.earth,newdata=training)
    out.sample.earth.pred  <- predict(fit.earth,newdata=testing)
    
    in.sample.gcvEarth.pred   <- predict(fit.gcvEarth,newdata=training)
    out.sample.gcvEarth.pred  <- predict(fit.gcvEarth,newdata=testing)
    
    in.sample.rf.pred   <- predict(fit.rf,newdata=training)
    out.sample.rf.pred  <- predict(fit.rf,newdata=testing)
    
    in.lm.rms          <- rmse(in.sample.data,in.sample.lm.pred)
    in.bagEarth.rms    <- rmse(in.sample.data,in.sample.bagEarth.pred)
    in.bagEarthGCV.rms <- rmse(in.sample.data,in.sample.bagEarthGCV.pred)
    in.earth.rms       <- rmse(in.sample.data,in.sample.earth.pred)
    in.gcvEarth.rms    <- rmse(in.sample.data,in.sample.gcvEarth.pred)
    in.rf.rms          <- rmse(in.sample.data,in.sample.rf.pred)
    
    out.lm.rms          <- rmse(out.sample.data,out.sample.lm.pred)
    out.bagEarth.rms    <- rmse(out.sample.data,out.sample.bagEarth.pred)
    out.bagEarthGCV.rms <- rmse(out.sample.data,out.sample.bagEarthGCV.pred)
    out.earth.rms       <- rmse(out.sample.data,out.sample.earth.pred)
    out.gcvEarth.rms    <- rmse(out.sample.data,out.sample.gcvEarth.pred)
    out.rf.rms          <- rmse(out.sample.data,out.sample.rf.pred)
    
    vec.in <- c(in.lm.rms,in.bagEarth.rms,in.bagEarthGCV.rms,
                in.earth.rms,in.gcvEarth.rms,in.rf.rms)
    vec.out <- c(out.lm.rms,out.bagEarth.rms,out.bagEarthGCV.rms,
                 out.earth.rms,out.gcvEarth.rms,out.rf.rms)
    
    type <- c("lm","bagEarth","bagEarthGCV","earth","gcvEarth","rf")
    
    errors <- data.frame(in.sample.error=vec.in,
                         out.sample.error=vec.out,
                         decimation=factor(decim.factor),
                         method=type)
}
set.seed(1234)
procData <- createModelingData("P",infraredData,50)
inTrain  <- createDataPartition(y=procData$P, p=0.75, list=FALSE)
errors50 <- getRmsVector(fitP.lm.50,fitP.bagEarth.50,
                         fitP.bagEarthGCV.50,fitP.earth.50,
                         fitP.gcvEarth.50,fitP.rf.50,
                         "P",infraredData,50,inTrain)
set.seed(1234)
procData <- createModelingData("P",infraredData,10)
inTrain  <- createDataPartition(y=procData$P, p=0.75, list=FALSE)
errors10 <- getRmsVector(fitP.lm.10,fitP.bagEarth.10,
                         fitP.bagEarthGCV.10,fitP.earth.10,
                         fitP.gcvEarth.10,fitP.rf.10,
                         "P",infraredData,10,inTrain)
set.seed(1234)
procData <- createModelingData("P",infraredData,5)
inTrain  <- createDataPartition(y=procData$P, p=0.75, list=FALSE)
errors5 <- getRmsVector(fitP.lm.5,fitP.bagEarth.5,
                         fitP.bagEarthGCV.5,fitP.earth.5,
                         fitP.gcvEarth.5,fitP.rf.5,
                         "P",infraredData,5,inTrain)

errors <- rbind(errors50,errors10,errors5)

gg1 <- ggplot(errors,aes(x=in.sample.error,y=out.sample.error)) 
gg1 <- gg1 + geom_point(aes(shape=method,color=decimation),size=3)
gg1 <- gg1 + geom_abline(intercept=0,slope=1,linetype=2)
gg1 <- gg1 + xlim(c(0,1.1)) + ylim(c(0,1.1))
plot(gg1)

# set.seed(1234)
# proc.data <- createModelingData("P",infraredData,10)
# inTrain  <- createDataPartition(y=proc.data$P, p=0.75, list=FALSE)
# training <- proc.data[inTrain,]
# testing  <- proc.data[-inTrain,]
# in.sample.rf.pred   <- predict(fitP.lm.10,newdata=training)
# out.sample.rf.pred  <- predict(fitP.lm.10,newdata=testing)
# 
# plot(training$P,in.sample.rf.pred)
# plot(testing$P,out.sample.rf.pred)
