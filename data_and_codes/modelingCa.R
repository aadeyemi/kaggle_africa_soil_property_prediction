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
## train model for Ca prediction given data decimation factor of 50
##------------------------------------------------------------------
set.seed(1234)
procData <- createModelingData("Ca",infraredData,50)
inTrain  <- createDataPartition(y=procData$Ca, p=0.75, list=FALSE)

trainingCa <- procData[inTrain,]
testingCa  <- procData[-inTrain,]

cl <- makeCluster(3)
registerDoParallel(cl)

tic(); fitCa.lm.50          <- train(Ca ~ ., method="lm",          data=trainingCa); toc()
tic(); fitCa.bagEarth.50    <- train(Ca ~ ., method="bagEarth",    data=trainingCa); toc()
tic(); fitCa.bagEarthGCV.50 <- train(Ca ~ ., method="bagEarthGCV", data=trainingCa); toc()
tic(); fitCa.earth.50       <- train(Ca ~ ., method="earth",       data=trainingCa); toc()
tic(); fitCa.gcvEarth.50    <- train(Ca ~ ., method="gcvEarth",    data=trainingCa); toc()
tic(); fitCa.rf.50          <- train(Ca ~ ., method="rf",          data=trainingCa); toc()

stopCluster(cl)

# save models to disk
save(fitCa.lm.50,          file = "saved_models/fitCa.lm.50.RData")
save(fitCa.bagEarth.50,    file = "saved_models/fitCa.bagEarth.50.RData")
save(fitCa.bagEarthGCV.50, file = "saved_models/fitCa.bagEarthGCV.50.RData")
save(fitCa.earth.50,       file = "saved_models/fitCa.earth.50.RData")
save(fitCa.gcvEarth.50,    file = "saved_models/fitCa.gcvEarth.50.RData")
save(fitCa.rf.50,          file = "saved_models/fitCa.rf.50.RData")

##------------------------------------------------------------------
## train model for Ca prediction given data decimation factor of 10
##------------------------------------------------------------------
set.seed(1234)
procData <- createModelingData("Ca",infraredData,10)
inTrain  <- createDataPartition(y=procData$Ca, p=0.75, list=FALSE)

trainingCa <- procData[inTrain,]
testingCa  <- procData[-inTrain,]

cl <- makeCluster(2)
registerDoParallel(cl)

tic(); fitCa.lm.10          <- train(Ca ~ ., method="lm",          data=trainingCa); toc()
tic(); fitCa.bagEarth.10    <- train(Ca ~ ., method="bagEarth",    data=trainingCa); toc()
tic(); fitCa.bagEarthGCV.10 <- train(Ca ~ ., method="bagEarthGCV", data=trainingCa); toc()
tic(); fitCa.earth.10       <- train(Ca ~ ., method="earth",       data=trainingCa); toc()
tic(); fitCa.gcvEarth.10    <- train(Ca ~ ., method="gcvEarth",    data=trainingCa); toc()
tic(); fitCa.rf.10          <- train(Ca ~ ., method="rf",          data=trainingCa); toc()

stopCluster(cl)

# save models to disk
save(fitCa.lm.10,          file = "saved_models/fitCa.lm.10.RData")
save(fitCa.bagEarth.10,    file = "saved_models/fitCa.bagEarth.10.RData")
save(fitCa.bagEarthGCV.10, file = "saved_models/fitCa.bagEarthGCV.10.RData")
save(fitCa.earth.10,       file = "saved_models/fitCa.earth.10.RData")
save(fitCa.gcvEarth.10,    file = "saved_models/fitCa.gcvEarth.10.RData")
save(fitCa.rf.10,          file = "saved_models/fitCa.rf.10.RData")


##------------------------------------------------------------------
## train model for Ca prediction given data decimation factor of 5
##------------------------------------------------------------------
set.seed(1234)
procData <- createModelingData("Ca",infraredData,5)
inTrain  <- createDataPartition(y=procData$Ca, p=0.75, list=FALSE)

trainingCa <- procData[inTrain,]
testingCa  <- procData[-inTrain,]

cl <- makeCluster(2)
registerDoParallel(cl)

tic(); fitCa.lm.5          <- train(Ca ~ ., method="lm",          data=trainingCa); toc()
tic(); fitCa.bagEarth.5    <- train(Ca ~ ., method="bagEarth",    data=trainingCa); toc()
tic(); fitCa.bagEarthGCV.5 <- train(Ca ~ ., method="bagEarthGCV", data=trainingCa); toc()
tic(); fitCa.earth.5       <- train(Ca ~ ., method="earth",       data=trainingCa); toc()
tic(); fitCa.gcvEarth.5    <- train(Ca ~ ., method="gcvEarth",    data=trainingCa); toc()
tic(); fitCa.rf.5          <- train(Ca ~ ., method="rf",          data=trainingCa); toc()

stopCluster(cl)

# save models to disk
save(fitCa.lm.5,          file = "saved_models/fitCa.lm.5.RData")
save(fitCa.bagEarth.5,    file = "saved_models/fitCa.bagEarth.5.RData")
save(fitCa.bagEarthGCV.5, file = "saved_models/fitCa.bagEarthGCV.5.RData")
save(fitCa.earth.5,       file = "saved_models/fitCa.earth.5.RData")
save(fitCa.gcvEarth.5,    file = "saved_models/fitCa.gcvEarth.5.RData")
save(fitCa.rf.5,          file = "saved_models/fitCa.rf.5.RData")


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

load(file = "saved_models/fitCa.lm.50.RData")
load(file = "saved_models/fitCa.bagEarth.50.RData")
load(file = "saved_models/fitCa.bagEarthGCV.50.RData")
load(file = "saved_models/fitCa.earth.50.RData")
load(file = "saved_models/fitCa.gcvEarth.50.RData")
load(file = "saved_models/fitCa.rf.50.RData")

errors50 <- getRmsVector(fitCa.lm.50,fitCa.bagEarth.50,
                         fitCa.bagEarthGCV.50,fitCa.earth.50,
                         fitCa.gcvEarth.50,fitCa.rf.50,
                         "Ca",infraredData,50,inTrain)

load(file = "saved_models/fitCa.lm.10.RData")
load(file = "saved_models/fitCa.bagEarth.10.RData")
load(file = "saved_models/fitCa.bagEarthGCV.10.RData")
load(file = "saved_models/fitCa.earth.10.RData")
load(file = "saved_models/fitCa.gcvEarth.10.RData")
load(file = "saved_models/fitCa.rf.10.RData")

errors10 <- getRmsVector(fitCa.lm.10,fitCa.bagEarth.10,
                         fitCa.bagEarthGCV.10,fitCa.earth.10,
                         fitCa.gcvEarth.10,fitCa.rf.10,
                         "Ca",infraredData,10,inTrain)

load(file = "saved_models/fitCa.lm.5.RData")
load(file = "saved_models/fitCa.bagEarth.5.RData")
load(file = "saved_models/fitCa.bagEarthGCV.5.RData")
load(file = "saved_models/fitCa.earth.5.RData")
load(file = "saved_models/fitCa.gcvEarth.5.RData")
load(file = "saved_models/fitCa.rf.5.RData")

errors5 <- getRmsVector(fitCa.lm.5,fitCa.bagEarth.5,
                         fitCa.bagEarthGCV.5,fitCa.earth.5,
                         fitCa.gcvEarth.5,fitCa.rf.5,
                         "Ca",infraredData,5,inTrain)

errors <- rbind(errors50,errors10,errors5)

gg1 <- ggplot(errors,aes(x=in.sample.error,y=out.sample.error)) 
gg1 <- gg1 + geom_point(aes(shape=method,color=decimation),size=3)
gg1 <- gg1 + geom_abline(intercept=0,slope=1,linetype=2)
gg1 <- gg1 + xlim(c(0,1)) + ylim(c(0,1))
plot(gg1)

