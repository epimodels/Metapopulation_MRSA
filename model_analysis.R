library(readr)
library(sm)
library(vioplot)

data <- read_csv("model_data.csv")

par(mfrow=c(1,1))
vioplot(data$MetaPop,data$Nurse_MD,data$Random,names=c("Metapopulation","Nurse-MD","Single Staff Type"),
        col="Grey85",drawRect=FALSE,ylim=c(0,100))

title(main="")
title(ylab="MRSA Acquisitions",cex.lab=1.35)
title(xlab="Model",cex.lab=1.35)
segments(0.75,mean(data$MetaPop),1.25,mean(data$MetaPop),col="black",lwd=2)
segments(1.75,mean(data$Nurse_MD),2.25,mean(data$Nurse_MD),col="black",lwd=2)
segments(2.75,mean(data$Random),3.25,mean(data$Random),col="black",lwd=2)

legend("topleft", c("Model Mean"), lwd=2,col=c("black"),lty=1, bty='n', cex=1.2)

kruskal.test(x=list(data$Random,data$Nurse_MD,data$MetaPop))
