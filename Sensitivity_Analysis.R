######### Global Sensitivitity Analysis #########

# Base Parameter Model Analysis
base <- read.csv("sens_base_psi.csv", header = T, sep = ",")

#Create a new variable of normalized cases
base$norm_cases <- base$cases/(median(base$cases))

# Separated Parameter Model Analysis
separate <- read.csv("sens_sep_psi.csv", header = T, sep = ",")

#Create a new variable of normalized cases
separate$md_norm <- separate$md_cases/(median(separate$md_cases))
separate$meta_norm <- separate$meta_cases/(median(separate$meta_cases)) 

########## Linear Regression Models ##########

#### Base Model ####

par(mfrow=c(1,1))

par(mfrow=c(2,2))

base_model <- lm(norm_cases ~ iota + mu + tau + theta + nu + psi + rho, data = base)

summary(base_model)

## Plot Sensitivity ##

coefficients <- base_model$coefficients[2:8]*0.1

cols <- c("grey90", "grey50")[(coefficients < 0) + 1]  

greek_names <- c(expression(iota),expression(mu),expression(tau),expression(theta),expression(nu),
                 expression(psi),expression(rho))

barplot(coefficients,horiz=TRUE,names.arg=greek_names,main="(A) Single Staff Type Model",ylab="Parameter",
        xlab="Change in Cumulative Acquisitions",cex.main=1.25,cex.lab=1.25,cex.axis=1.25,cex.names=1.25,
        xlim = c(-0.5,0.5),col=cols)


#### Nurse-MD Model ####

md_model <- lm(md_norm ~ iota_N + iota_D + mu + tau_N + tau_D + theta + nu + psi
               + rho_N + rho_D, data = separate)

summary(md_model)

# Plot Sensitivity

coefficients1 <- md_model$coefficients[2:11]*0.1

cols1 <- c("grey90", "grey50")[(coefficients1 < 0) + 1]  

greek_names1 <- c(expression(iota[N]),expression(iota[D]),expression(mu),
                  expression(tau[N]),expression(tau[D]),expression(theta),
                  expression(nu),expression(psi),expression(rho[N]),expression(rho[D]))

barplot(coefficients1,horiz=TRUE,names.arg=greek_names1,main="(B) Nurse-MD Model",ylab="Parameter",
        xlab="Change in Cumulative Acquisitions",cex.lab=1.25,cex.axis=1.25,cex.names=1.25,
        xlim = c(-0.5,0.5),col=cols1)

#### Meta-Pop Model ####

metapop_model <- lm(meta_norm ~ iota_N + iota_D + mu + tau_N + tau_D + theta + nu + psi
                    + rho_N + rho_D, data = separate)

summary(metapop_model)

# Plot Sensitivity

coefficients2 <- metapop_model$coefficients[2:11]*0.1

cols2 <- c("grey90", "grey50")[(coefficients2 < 0) + 1]  

barplot(coefficients2,horiz=TRUE,names.arg=greek_names1,main="(C) Metapopulation Model",ylab="Parameter",
        xlab="Change in Cumulative Acquisitions",cex.lab=1.25,cex.axis=1.25,cex.names=1.25,
        xlim = c(-0.5,0.5),col=cols2)

#### MD vs Meta-Pop Difference Chart ####

coefficients3 <- coefficients1-coefficients2

cols3 <- c("grey90", "grey50")[(coefficients3 < 0) + 1] 

barplot(coefficients3,horiz=TRUE,names.arg=greek_names1,main="(D) Difference between Model B & C",ylab="Parameter",
        xlab="Change in Cumulative Acquisitions",cex.lab=1.25,cex.axis=1.25,cex.names=1.25,
        xlim = c(-0.5,0.5),col=cols3)

