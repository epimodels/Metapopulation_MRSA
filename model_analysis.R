library(readr)
library(sm)
library(vioplot)

data <- read_csv("model_data.csv")
data$SST_norm <- (data$SST/6570)*1000
data$Nurse_MD_norm <- (data$Nurse_MD/6570)*1000
data$Meta_norm <- (data$MetaPop/6570)*1000

par(mfrow=c(1,1))
vioplot(data$Meta_norm,data$Nurse_MD_norm,data$SST_norm,names=c("Metapopulation","Nurse-MD","Single Staff Type"),
        col="Grey85",drawRect=FALSE,ylim=c(0,15))

title(main="")
title(ylab="MRSA Acquisitions per 1,000 Patient-Days",cex.lab=1.35)
title(xlab="Model",cex.lab=1.35)
segments(0.75,mean(data$Meta_norm),1.25,mean(data$Meta_norm),col="black",lwd=2)
segments(1.75,mean(data$Nurse_MD_norm),2.25,mean(data$Nurse_MD_norm),col="black",lwd=2)
segments(2.75,mean(data$SST_norm),3.25,mean(data$SST_norm),col="black",lwd=2)

legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)

kruskal.test(x=list(data$SST,data$Nurse_MD,data$MetaPop))
