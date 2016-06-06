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

## train model for sand prediction given data decimation factor of 50
set.seed(1234)
procData <- createModelingData("sand",infraredData,50)
inTrain  <- createDataPartition(y=procData$Sand, p=0.75, list=FALSE)

trainingSand <- procData[inTrain,]
testingSand  <- procData[-inTrain,]

cl <- makeCluster(3)
registerDoParallel(cl)

tic(); fitSand.lm.50          <- train(Sand ~ ., method="lm",          data=trainingSand); toc()
tic(); fitSand.bagEarth.50    <- train(Sand ~ ., method="bagEarth",    data=trainingSand); toc()
tic(); fitSand.bagEarthGCV.50 <- train(Sand ~ ., method="bagEarthGCV", data=trainingSand); toc()
tic(); fitSand.earth.50       <- train(Sand ~ ., method="earth",       data=trainingSand); toc()
tic(); fitSand.gcvEarth.50    <- train(Sand ~ ., method="gcvEarth",    data=trainingSand); toc()
tic(); fitSand.rf.50          <- train(Sand ~ ., method="rf",          data=trainingSand); toc()

stopCluster(cl)

# save models to disk
save(fitSand.lm.50,          file = "saved_models/fitSand.lm.50.RData")
save(fitSand.bagEarth.50,    file = "saved_models/fitSand.bagEarth.50.RData")
save(fitSand.bagEarthGCV.50, file = "saved_models/fitSand.bagEarthGCV.50.RData")
save(fitSand.earth.50,       file = "saved_models/fitSand.earth.50.RData")
save(fitSand.gcvEarth.50,    file = "saved_models/fitSand.gcvEarth.50.RData")
save(fitSand.rf.50,          file = "saved_models/fitSand.rf.50.RData")



## train model for sand prediction given data decimation factor of 10
set.seed(1234)
procData <- createModelingData("sand",infraredData,10)
inTrain  <- createDataPartition(y=procData$Sand, p=0.75, list=FALSE)

trainingSand <- procData[inTrain,]
testingSand  <- procData[-inTrain,]

cl <- makeCluster(3)
registerDoParallel(cl)

tic(); fitSand.lm.10          <- train(Sand ~ ., method="lm",          data=trainingSand); toc()
tic(); fitSand.bagEarth.10    <- train(Sand ~ ., method="bagEarth",    data=trainingSand); toc()
tic(); fitSand.bagEarthGCV.10 <- train(Sand ~ ., method="bagEarthGCV", data=trainingSand); toc()
tic(); fitSand.earth.10       <- train(Sand ~ ., method="earth",       data=trainingSand); toc()
tic(); fitSand.gcvEarth.10    <- train(Sand ~ ., method="gcvEarth",    data=trainingSand); toc()
tic(); fitSand.rf.10          <- train(Sand ~ ., method="rf",          data=trainingSand); toc()

stopCluster(cl)

# save models to disk
save(fitSand.lm.10,          file = "saved_models/fitSand.lm.10.RData")
save(fitSand.bagEarth.10,    file = "saved_models/fitSand.bagEarth.10.RData")
save(fitSand.bagEarthGCV.10, file = "saved_models/fitSand.bagEarthGCV.10.RData")
save(fitSand.earth.10,       file = "saved_models/fitSand.earth.10.RData")
save(fitSand.gcvEarth.10,    file = "saved_models/fitSand.gcvEarth.10.RData")
save(fitSand.rf.10,          file = "saved_models/fitSand.rf.10.RData")


## train model for sand prediction given data decimation factor of 5
set.seed(1234)
procData <- createModelingData("sand",infraredData,5)
inTrain  <- createDataPartition(y=procData$Sand, p=0.75, list=FALSE)

trainingSand <- procData[inTrain,]
testingSand  <- procData[-inTrain,]

cl <- makeCluster(3)
registerDoParallel(cl)

tic(); fitSand.lm.5          <- train(Sand ~ ., method="lm",          data=trainingSand); toc()
tic(); fitSand.bagEarth.5    <- train(Sand ~ ., method="bagEarth",    data=trainingSand); toc()
tic(); fitSand.bagEarthGCV.5 <- train(Sand ~ ., method="bagEarthGCV", data=trainingSand); toc()
tic(); fitSand.earth.5       <- train(Sand ~ ., method="earth",       data=trainingSand); toc()
tic(); fitSand.gcvEarth.5    <- train(Sand ~ ., method="gcvEarth",    data=trainingSand); toc()
tic(); fitSand.rf.5          <- train(Sand ~ ., method="rf",          data=trainingSand); toc()

stopCluster(cl)

# save models to disk
save(fitSand.lm.5,          file = "saved_models/fitSand.lm.5.RData")
save(fitSand.bagEarth.5,    file = "saved_models/fitSand.bagEarth.5.RData")
save(fitSand.bagEarthGCV.5, file = "saved_models/fitSand.bagEarthGCV.5.RData")
save(fitSand.earth.5,       file = "saved_models/fitSand.earth.5.RData")
save(fitSand.gcvEarth.5,    file = "saved_models/fitSand.gcvEarth.5.RData")
save(fitSand.rf.5,          file = "saved_models/fitSand.rf.5.RData")


## train model for sand prediction given data decimation factor of 3
set.seed(1234)
procData <- createModelingData("sand",infraredData,3)
inTrain  <- createDataPartition(y=procData$Sand, p=0.75, list=FALSE)

trainingSand <- procData[inTrain,]
testingSand  <- procData[-inTrain,]

cl <- makeCluster(3)
registerDoParallel(cl)

tic(); fitSand.lm.3          <- train(Sand ~ ., method="lm",          data=trainingSand); toc()
tic(); fitSand.bagEarth.3    <- train(Sand ~ ., method="bagEarth",    data=trainingSand); toc()
tic(); fitSand.bagEarthGCV.3 <- train(Sand ~ ., method="bagEarthGCV", data=trainingSand); toc()
tic(); fitSand.earth.3       <- train(Sand ~ ., method="earth",       data=trainingSand); toc()
tic(); fitSand.gcvEarth.3    <- train(Sand ~ ., method="gcvEarth",    data=trainingSand); toc()
tic(); fitSand.rf.3          <- train(Sand ~ ., method="rf",          data=trainingSand); toc()

stopCluster(cl)

# save models to disk
save(fitSand.lm.3,          file = "saved_models/fitSand.lm.3.RData")
save(fitSand.bagEarth.3,    file = "saved_models/fitSand.bagEarth.3.RData")
save(fitSand.bagEarthGCV.3, file = "saved_models/fitSand.bagEarthGCV.3.RData")
save(fitSand.earth.3,       file = "saved_models/fitSand.earth.3.RData")
save(fitSand.gcvEarth.3,    file = "saved_models/fitSand.gcvEarth.3.RData")
save(fitSand.rf.3,          file = "saved_models/fitSand.rf.3.RData")


## train model for sand prediction given data decimation factor of 1
set.seed(1234)
procData <- createModelingData("sand",infraredData,50)
inTrain  <- createDataPartition(y=procData$Sand, p=0.75, list=FALSE)

trainingSand <- procData[inTrain,]
testingSand  <- procData[-inTrain,]

cl <- makeCluster(3)
registerDoParallel(cl)

# tic(); fitSand.lm.1          <- train(Sand ~ ., method="lm",          data=trainingSand); toc()
# tic(); fitSand.bagEarth.1    <- train(Sand ~ ., method="bagEarth",    data=trainingSand); toc()
# tic(); fitSand.bagEarthGCV.1 <- train(Sand ~ ., method="bagEarthGCV", data=trainingSand); toc()
# tic(); fitSand.earth.1       <- train(Sand ~ ., method="earth",       data=trainingSand); toc()
# tic(); fitSand.gcvEarth.1    <- train(Sand ~ ., method="gcvEarth",    data=trainingSand); toc()
# tic(); fitSand.rf.1          <- train(Sand ~ ., method="rf",          data=trainingSand); toc()

stopCluster(cl)

# save models to disk
# save(fitSand.lm.1,          file = "saved_models/fitSand.lm.1.RData")
# save(fitSand.bagEarth.1,    file = "saved_models/fitSand.bagEarth.1.RData")
# save(fitSand.bagEarthGCV.1, file = "saved_models/fitSand.bagEarthGCV.1.RData")
# save(fitSand.earth.1,       file = "saved_models/fitSand.earth.1.RData")
# save(fitSand.gcvEarth.1,    file = "saved_models/fitSand.gcvEarth.1.RData")
# save(fitSand.rf.1,          file = "saved_models/fitSand.rf.1.RData")