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
