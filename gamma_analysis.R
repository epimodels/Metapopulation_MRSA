#### Gamma Sensitivity Analysis ####

library(segmented)

data_gamma <- read.csv("gamma_data.csv", header = T, sep = ",")
attach(data_gamma)
data_gamma$denom <- 18*365

str(data_gamma)
summary(data_gamma)


gammasweep <- glm(cases ~ gamma + I(gamma^2)+offset(log(data_gamma$denom)),family="poisson")
summary(gammasweep)

segmod <- segmented(gammasweep,seg.Z = ~ gamma, psi=0.50)
summary(segmod)
fitted <- fitted(segmod)
changepoint <- segmod$psi[2]
ucl <- changepoint + (1.96*segmod$psi[3])
lcl <- changepoint - (1.96*segmod$psi[3])

modeldata <- data.frame(gamma=data_gamma$gamma,fit=fitted)
modeldata <- modeldata[order(gamma),]

par(mar=c(5,4.75,4,2)+0.1)
plot(gamma,(cases/6570)*1000,pch=21,col="grey90",bg="grey90",xlab=expression(gamma ~ ": Proportion of Time Nurses Spend with Assigned Patients"),ylab="MRSA Acquisitions per 1,000 Patient-Days",cex.lab=1.45)
rect(xleft=lcl,xright=ucl,ybottom=0,ytop=200,col="lightcoral",density=90)
lines(modeldata$gamma,(modeldata$fit/6570)*1000, type="l",lwd=5)
abline(v=changepoint,lty=2,lwd=2)
legend("topright", c("Segmented Poisson Fit",expression(gamma*"*"),expression(gamma*"*"~ ~"95% Confidence Interval")), 
       lwd=c(3,2,NA),pch=c(NA,NA,15),pt.cex=c(1,1,3),col=c("black","black","lightcoral"),lty=c(1,2), bty='n', cex=1.1)
