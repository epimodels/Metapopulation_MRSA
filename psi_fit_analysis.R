
sst <- read.csv("psi_fit_SST.csv", sep="")
nurse_md <- read.csv("psi_fit_nurse_md.csv", sep="")
meta <- read.csv("psi_fit_meta.csv", sep="")


quantile(sst$data, probs = c(0.025,0.50,0.975))
quantile(nurse_md$data, probs = c(0.025,0.50,0.975))
quantile(meta$data, probs = c(0.025,0.50,0.975))