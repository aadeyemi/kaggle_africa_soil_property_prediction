Introduction
------------

I have explored the data set [Africa Soil Property Prediction Challenge
data set](!https://www.kaggle.com/c/afsis-soil-properties/data) (*see
the data exploration report*) and see evidence to support that the
provided data set can be used to build prediction models for the desired
outcome variables (**Sand**, **pH**, **Ca**, **P**, and **SOC**). The
goal of this report is to detail my thought process as I explore
different models and select the best option based on prediction error,
which in this case is the root-mean-squared-error (RMSE). The final
pre-processed data set consists of 3579 predictor variables and they
are:

-   **m7497.96 - m599.76:** 3563 mid-infrared absorbance measurements
-   **Depth:** Depth of the soil sample (2 categories:
    "Topsoil", "Subsoil")
-   BSA: average long-term Black Sky Albedo measurements from MODIS
    satellite images (**BSAN** = near-infrared, **BSAS** = shortwave,
    **BSAV** = visible)
-   **CTI:** compound topographic index calculated from Shuttle Radar
    Topography Mission elevation data
-   **ELEV:** Shuttle Radar Topography Mission elevation data
-   **EVI:** average long-term Enhanced Vegetation Index from MODIS
    satellite images.
-   LST: average long-term Land Surface Temperatures from MODIS
    satellite images (**LSTD** = day time temperature, **LSTN** = night
    time temperature)
-   REF: average long-term Reflectance measurements from MODIS satellite
    images (**REF1** = blue, **REF2** = red, **REF3** = near-infrared,
    **REF7** = mid-infrared)
-   **RELI:** topographic Relief calculated from Shuttle Radar
    Topography mission elevation data
-   **TMAP & TMFI:** average long-term Tropical Rainfall Monitoring
    Mission data (TMAP = mean annual precipitation, TMFI = modified
    Fournier index)

After data exploration, I found that the infrared absorbance
measurements are the most correlated with our outcome variables while
the auxiliary predictor variables (BSA, CTI, ELEV, EVI, LST, REF, RELI,
TMAP and TMFI) are weakly or not correlated to the outcome variables
except for pH where moderate correlations exists with the predictor
variables EVI, TMAP, LSTD, and Ref7.

Data Loading
------------

Load and prepare data for pre-processing.

    ## Load Data
    data0 <- read.csv(unz("data_and_codes/train.zip","training.csv"), na.strings=c("NA",""))

    # Remove spectra CO2 bands 
    a <- grep("m2379.76",names(data0))
    b <- grep("m2352.76",names(data0))

    data <- data0[,-seq(a,b)]

    ## Create spectra-only data subset
    infraredCols <- grep("^m[0-9].*?\\.[0-9].*",names(data))
    infraredData <- data.frame(data[,infraredCols])

Modeling
--------

Since we will be separating data into predictors variables and outcomes
multiple times, define helper functions to accomplish this task. The
first function *getFinalData* builds the final data frame from its
components: (1) the outcome variable, (2) the Depth column, (3) the
decimated infrared measurements column. The second function
*createModelingData* performs data creation for the specified outcome
variable and decimation parameter using *getFinalData*. In creating the
Sand prediction model, I have excluded the auxiliary predictor
variables.

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

### 1. Sand

I tested modeling methods: lm (Linear), earth (Multivariate Adaptive
Regression Spline), gcvEarth (Multivariate Adaptive Regression Splines),
bagEarth (bagged MARS), bagEarthGCV (bagged MARS using gCV pruning), rf
(Random Forests). Below is the code sequence I used for modeling Sand
outcome with a spectral data decimation factor of 10. I show only the
code sequence for this run for brevity because all other similar
decimation and modeling method tests follow the same sequence.
Subsequently, I will only show summary results.

    ## train model for sand prediction given data decimation factor of 10
    set.seed(1234)
    procData <- createModelingData("sand",infraredData,10)
    inTrain  <- createDataPartition(y=procData$Sand, p=0.75, list=FALSE)

    trainingSand <- procData[inTrain,]
    testingSand  <- procData[-inTrain,]

    fitSand.lm.10          <- train(Sand ~ ., method="lm",          data=trainingSand)
    fitSand.bagEarth.10    <- train(Sand ~ ., method="bagEarth",    data=trainingSand)
    fitSand.bagEarthGCV.10 <- train(Sand ~ ., method="bagEarthGCV", data=trainingSand)
    fitSand.earth.10       <- train(Sand ~ ., method="earth",       data=trainingSand)
    fitSand.gcvEarth.10    <- train(Sand ~ ., method="gcvEarth",    data=trainingSand)
    fitSand.rf.10          <- train(Sand ~ ., method="rf",          data=trainingSand)

![](prediction_models_files/figure-markdown_strict/plot%20models-1.png)<!-- -->

We can see from Figure 1 that the errors vary more with training method
than decimation factors. In fact, for rf and bagEarthGCV method, there
isn't much difference in prediction rms errors. Random forest and linear
regression methods produced results with in-sample errors not reflective
of out-of-sample errors.

Example prediction plots are shown below.

![](prediction_models_files/figure-markdown_strict/plot%20sand%20outcome%20example%201a-1.png)<!-- -->

![](prediction_models_files/figure-markdown_strict/plot%20sand%20outcome%20example%201b-1.png)<!-- -->

![](prediction_models_files/figure-markdown_strict/plot%20sand%20outcome%20example%202a-1.png)<!-- -->

![](prediction_models_files/figure-markdown_strict/plot%20sand%20outcome%20example%202b-1.png)<!-- -->

### 2. SOC

Just as with the Sand content, I tested several decimation factors and
modeling methods for the SOC measure. Results are shown below.

![](prediction_models_files/figure-markdown_strict/plot%20soc%20models-1.png)<!-- -->

Figure 6 shows in-sample errors less than 0.4 and out-of-sample errors
less than 0.6 with the exception of the linear regression model for the
SOC variable. Example prediction plots are shown below.

Example prediction plots are shown below.

![](prediction_models_files/figure-markdown_strict/plot%20soc%20outcome%20example%201a-1.png)<!-- -->

![](prediction_models_files/figure-markdown_strict/plot%20soc%20outcome%20example%201b-1.png)<!-- -->

### 3. Ca

Prediction error plots showing outcome of decimation tests for Ca are
shown below.

![](prediction_models_files/figure-markdown_strict/plot%20ca%20models-1.png)<!-- -->

Figure 9 shows in-sample errors less than 0.4 and out-of-sample errors
less than 0.6 with the exception of the linear regression model for the
Ca variable. Example prediction plots are shown below.

Example prediction plots are shown below.

![](prediction_models_files/figure-markdown_strict/plot%20ca%20outcome%20example%201a-1.png)<!-- -->

![](prediction_models_files/figure-markdown_strict/plot%20ca%20outcome%20example%201b-1.png)<!-- -->

### 4. P

Prediction error plots of decimation tests for P are shown below.

![](prediction_models_files/figure-markdown_strict/plot%20p%20models-1.png)<!-- -->

It appears from Figure 12 that the predictor variables are not good
predictors of P. Out-of-sample errors are lower than in-sample errors in
some cases which does not make sense given the data set. Let's take a
look at some prediction plots.

![](prediction_models_files/figure-markdown_strict/plot%20p%20outcome%20example%201a-1.png)<!-- -->

![](prediction_models_files/figure-markdown_strict/plot%20p%20outcome%20example%201b-1.png)<!-- -->

As expected, our prediction model does not explain the data to an
acceptable degree. One way to solve this problem is to fine tune our
modeling parameters or seek additional data to better predict P from the
spectral data provided.

### 5. pH

In modeling pH, we will consider the case where auxiliary predictor
variables (EVI, TMAP, LSTD, REF7) are included in training and the case
where they are not. During exploratory analysis, moderate correlation
between pH and the auxiliary variables was seen.

#### 5a. Without auxilliary variables

![](prediction_models_files/figure-markdown_strict/plot%20ph%20models-1.png)<!-- -->

We can see in Figure 15 that with the exception of random forest and
linear regression methods, all methods produce results with
out-of-sample errors that mimic the in-sample errors. That's a good
generalization. In addition, points cluster around the same location
regardless of decimation factor for earth, bagEarth, bagEarthGCV, and
gcvEarth methods.

![](prediction_models_files/figure-markdown_strict/plot%20ph%20outcome%20example%201a-1.png)<!-- -->

![](prediction_models_files/figure-markdown_strict/plot%20ph%20outcome%20example%201b-1.png)<!-- -->

#### 5b. With auxilliary variables

![](prediction_models_files/figure-markdown_strict/plot%20ph.aux%20models-1.png)<!-- -->

Figure 18 shows improvements in in-sample and out-of-sample error due to
the inclusion of the auxiliary variables in predicting pH. This
observation was expected as we saw the moderate correlation of pH with
the auxiliary variables during exploratory data analysis.

Below are examples of predicted data using the updated model.

![](prediction_models_files/figure-markdown_strict/plot%20ph.aux%20outcome%20example%201a-1.png)<!-- -->

![](prediction_models_files/figure-markdown_strict/plot%20ph.aux%20outcome%20example%201b-1.png)<!-- -->

Conclusions
-----------

In conclusion, I have shown through numerous tests that diffuse
reflectance infrared spectroscopy can be used to predict sand content,
Mehlich-3 extractable calcium, pH values and organic carbon in soils.
However, I was unable to model relationships between the spectral data
and Mehlich-3 extractable phosphorus. From initial exploratory data
analysis, it was evident that any relationship between them will be
tough to model.

Secondly, by decimating the data in such a way that aliasing is avoided,
it was possible to reduce the number of predictor variables by up to a
factor of 50 and still obtain reasonable results. This would be
impossible if aliasing effects were not considered. Reducing the number
of predictor variables from over 3500 to less than 100 definitely saves
computation time.

Thirdly, for the pH outcome variable, including auxiliary predictor
variables EVI (average long-term Enhanced Vegetation Index from MODIS
satellite images), TMAP (mean annual precipitation), LSTD (day time
temperature), and REF7 (average long-term Reflectance measurements from
MODIS satellite images - mid-infrared), helped improve prediction
errors. The possibility of this improvement was first observed during
exploratory data analysis through computed cross-correlation values.
