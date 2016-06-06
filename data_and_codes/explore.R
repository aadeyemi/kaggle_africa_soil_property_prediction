library(stats)
library(reshape2)
library(signal)

## Load Data
data0 <- read.csv(unz("train.zip","training.csv"), na.strings=c("NA",""))

# Remove spectra CO2 bands 
a <- grep("m2379.76",names(data0))
b <- grep("m2352.76",names(data0))

data <- data0[,-seq(a,b)]

## Create spectra-only data subset
infraredCols <- grep("^m[0-9].*?\\.[0-9].*",names(data))
infraredData <- data.frame(data[,infraredCols])

## Create new data frame for plotting
vec <- seq(1,nrow(infraredData))
plotData <- data.frame(infraredData,Index=vec)

## compare spectral data from multiple sites
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

## compare spectral characteristics obtained after decimating original data
siteNum <- 1
newVec  <- newIndex <- newVecId <- vector()
vec     <- as.vector(t(plotData[siteNum,]))

for (fac in c(1,4,10,20)) {
    newVec   <- c(newVec,decimate(vec,fac))
    newIndex <- c(newIndex,seq(1,length(vec),by=fac))
    newVecId <- c(newVecId,rep(fac,length(seq(1,length(vec),by=fac))))
}

pdata2 <- data.frame(index=newIndex,value=newVec,Decimation=factor(newVecId))
gg2 <- ggplot(pdata2,aes(x=index,y=value))
gg2 <- gg2 + geom_line(aes(color=Decimation))
gg2 <- gg2 + coord_cartesian(ylim = c(0,2))
gg2 <- gg2 + labs(x="Spectral Index", y="Value")
plot(gg2)

## plots of spectral signal power 
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
# pH vs LSTD
gg4 <- ggplot(properties,aes(x=LSTD,y=pH))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# pH vs EVI
gg4 <- ggplot(properties,aes(x=EVI,y=pH))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# pH vs TMAP
gg4 <- ggplot(properties,aes(x=TMAP,y=pH))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# Sand vs RELI
gg4 <- ggplot(properties,aes(x=RELI,y=Sand))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)
# P vs CTI
gg4 <- ggplot(properties,aes(x=CTI,y=P))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + coord_cartesian(xlim=c(-1,3.5),ylim = c(-0.5,2.5))
gg4 <- gg4 + geom_smooth()
plot(gg4)
# SOC vs REF2
gg4 <- ggplot(properties,aes(x=REF2,y=SOC))
gg4 <- gg4 + geom_point(alpha=0.5)
gg4 <- gg4 + geom_smooth()
plot(gg4)





gg3 <- ggplot(properties,aes(x=Power,y=SOC))
gg3 <- gg3 + geom_point(alpha = 0.5)
gg3 <- gg3 + geom_smooth(method="lm")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Power,y=pH))
gg3 <- gg3 + geom_point(alpha = 0.5)
gg3 <- gg3 + geom_smooth(method="lm")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Power,y=P))
gg3 <- gg3 + geom_point(alpha = 0.5)
gg3 <- gg3 + geom_smooth(method="lm")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Power,y=Sand))
gg3 <- gg3 + geom_point(alpha = 0.5)
gg3 <- gg3 + geom_smooth(method="lm")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Power,y=Ca))
gg3 <- gg3 + geom_point(alpha = 0.5)
gg3 <- gg3 + geom_smooth(method="lm")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Depth,y=Sand))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Depth,y=pH))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Depth,y=log10(P)))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Depth,y=log10(Ca)))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

gg3 <- ggplot(properties,aes(x=Depth,y=SOC))
gg3 <- gg3 + geom_violin(fill="deeppink4")
plot(gg3)

