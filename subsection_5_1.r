#_______________________________________________________________________________________________________________________#
#Online accompaniement of the paper "Fraud risk assessment within block chain transactions"
#File: subsection_5_1.r
#Author: Pierre-O Goffard
#_______________________________________________________________________________________________________________________#

#Libraries
library(fitdistrplus)
library(moments)
library(xtable)

#File directory
my_path <- "C:/Users/goffard/Dropbox/goffard/BitCoins/Mathematica_Rcode/"
#Name of the file containing the data
data_file <- "times_and_targets_interp2.txt"

#Importation of the data
#data frame with the original time stamps and the preprocessed time stamp
#The file is available on Rhys Bowden github page
time_stamps <- read.table(paste(my_path, data_file, sep = ""), quote = "\"", comment.char = "")
#Vector of inter-arrival times of blocks
inter_block_times <- time_stamps[[2]][-1] - time_stamps[[2]][-length(time_stamps[[2]])]

#Definition of a new time stamps data frame containing the inter-block times
time_stamps_new <- data.frame(
  block_time_stamp = time_stamps[[1]][2:length(time_stamps[[2]])],
  block_time = time_stamps[[2]][2:length(time_stamps[[2]])],
  inter_block_time = c(0, inter_block_times[2:length(inter_block_times)]) / 60)#Inter-block time in minutes

#Empirical mean and standard deviation of the inter-block times
mean(time_stamps_new$inter_block_time[-1])
sd(time_stamps_new$inter_block_time[-1])
n=length(time_stamps_new$inter_block_time)#Number of blocks

#Time series of the interblock times

#File directory for the outputs
my_path_output <- "C:/Users/goffard/Dropbox/goffard/BitCoins/"

#Time series of the inter-block times
file_name <- "InterBlockTimeChronoSerie.png"
png(paste(my_path_output, file_name, sep = ""))
plot(time_stamps_new$inter_block_time, 
     type = "l", main = "", xlab = "Block Number", ylab = "Inter-Block times (in minutes)", cex.lab = 1.5, cex.axis = 1.5)
dev.off()

#Time series of the inter-block times starting from the 400000th onward
file_name <- "InterBlockTimeChronoSerieOver40000.png"
png(paste(my_path_output, file_name, sep = ""))
plot(400000:length(time_stamps_new$inter_block_time),
     time_stamps_new$inter_block_time[400000:length(time_stamps_new$inter_block_time)],
     type = "l", xlab = "Block Number", ylab = "Inter-Block times (in minutes)"
     ,cex.lab = 1.5,cex.axis = 1.5)
dev.off()

#Inter block time selected (i.e. the most recent ones) starting from the 400,000th block onward
inter_block_time_recent <- time_stamps_new$inter_block_time[400000:length(time_stamps_new$inter_block_time)]
#Number of data points, empirical mean and standard deviation
n <- length(inter_block_time_recent)
mean(inter_block_time_recent)
sd(inter_block_time_recent)

#Histogram of the inter-block time with benchmarked by the exponential p.d.f.
lam_hat <- 1 / mean(inter_block_time_recent)#MME of the rate parameter of the exponential distribution
file_name <- "InterBlockHist.png"
png(paste(my_path_output, file_name, sep = ""))
h <- hist(inter_block_time_recent, plot = FALSE) #generate hist
plot(h, col = "grey", main = "", xlab = "Time (in minutes)",cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5) #plot hist
xlines <- seq(min(h$breaks), max(h$breaks), length.out = 100) #seq of x for pdf
lines(x = xlines, y = dexp(xlines, rate = lam_hat) * length(inter_block_time_recent) * diff(h$breaks)[1], lwd = 2)
dev.off()

#QQ-Plot against exponential distribution quantile
file_name <- "InterBlockQQPlot.png"
png(paste(my_path_output, file_name, sep = ""))
plot(sapply(seq(0.01, 0.99, 0.01), function(x) quantile(inter_block_time_recent, x)),
  sapply(seq(0.01, 0.99, 0.01),function(x) qexp(x, rate = lam_hat)),
  xlab = "Theoretical quantiles", ylab = "Empirical quantiles", type = "p",
  cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
lines(sapply(seq(0.01, 0.99, 0.01),function(x) qexp(x, rate = lam_hat)), 
      sapply(seq(0.01, 0.99, 0.01),function(x) qexp(x, rate = lam_hat)), type='l')
dev.off()

#Inference and goodness of fit test for the inter block times
## Exponential distribution
fit_exp <- fitdist(inter_block_time_recent, "exp", method = 'mme')
## Gamma Distribution
fit_gamma <-fitdist(inter_block_time_recent, "gamma", method = 'mme')#This is an exponential distribution

##Function to compute the MME of the parameters of a Weibull distribution
mme_weibull <- function(U){
  f <- function(x){
    gamma(1 + 2 / x) / gamma(1 + 1 / x)^(2) - moment(U, order = 2, central = FALSE) / mean(U)^(2)
    }
  roots <- uniroot(f, c(0.1, 5))
  shape <- roots$root[1]
  scale <- mean(U) / gamma(1 + 1 / shape)
  return(c(shape, scale))
}

##Weibull distribution
fit_weibull <- mme_weibull(inter_block_time_recent)#Which is equivalent to an exponential distribution

##Adequacy of the exponential distribution
gof_exp <- gofstat(fit_exp,fitnames = c("exp"))
#Exponential is definitely the distribution of those inter-block times

#Parametric bootstrap routine to assess the critical point of the K-S test statistic
ks_bootstrap <- {}
i=1
for(i in 1:1000){
  bootstrap_inter_block_time <- rexp(n, rate = fit_exp$estimate[1])
  fit_exp_bootstrap <- fitdist(bootstrap_inter_block_time, "exp", method = 'mme')
  gof_bootstrap <- gofstat(fit_exp_bootstrap, fitnames = c("exp"))
  ks_bootstrap[i]=gof_bootstrap$ks[1]
}

#Rejection
quantile(ks_bootstrap,0.95) > gof_exp$ks

#Distribution of the inter-block time within the malicious chain
#Transformation of the distribution in the interblock time of the public chain as DeltaT = exp(lam)
#Assume that the block discovery in the malicious chain is r times slower than in the honest chain

#Scaling transformation: DeltaS = r * DeltaT = Exp(lam / r = mu)
##PDF of the double spending time for exponential DeltaS
pdf_tau_exp <- function(K, z, lam, mu, t){
  sum(
    sapply(
      0:K, function(n) z / (z + n) * gamma(2 * n + z) / gamma(n + 1) / gamma( n + z)*
               (mu / (mu + lam))^(n + z) * (lam / (mu + lam))^(n)*dgamma(t, shape = 2 * n + z, rate = lam + mu)
    )
  )
}
##CDF of the double spending time for exponential DeltaS
cdf_tau_exp <- function(K, z, lam, mu, t){
  sum(
    sapply(
      0:K, function(n) z / (z + n) * gamma(2 * n + z) / gamma(n + 1) / gamma(n + z)*
               (mu / (mu + lam))^(n + z)*(lam / (mu + lam))^(n)*pgamma(t, shape = 2 * n + z, rate = lam + mu)
    )
  )
}
##Probability of a successful double spending for exponential DeltaS
p_tau_exp <- function(z, lam, mu){
  (mu / lam)^(z)
}

#Compounding transformation: DeltaS = DeltaT_1 + ... + DeltaT_r = gamma(r, lam)
##PDF of the double spending time for gamma DeltaS
pdf_tau_gamma <- function(K, z, lam, r, t){
  sum(
    sapply(
      0:K, function(n) (z / (z + n)) * gamma((n + z) * r + n) / 2^((n + z) * r + n)/gamma(n + 1)
             /gamma((n + z) * r) * dgamma(t, shape=(n + z) * r + n,rate = 2 * lam)
    )
  )
}
##CDF of the double spending time for gamma DeltaS
cdf_tau_gamma <- function(K, z, lam, r, t){
  sum(
    sapply(
      0:K, function(n) (z / (z + n))*gamma((n + z) * r + n) / 2^((n + z) * r + n) / gamma(n + 1) / 
               gamma((n + z) * r) * pgamma(t, shape=(n + z) * r + n,rate = 2 * lam)
    )
  )
}
##Probability of a successful double spending for gamma DeltaS
p_tau_gamma <- function(z, lam, r, gamInit){
  gam <- uniroot(function(x) log(lam^(r + 1) / (lam - x) / (lam + x)^(r)),c(gamInit, lam))$root
  return((lam / (lam + gam))^(r * (z - 1))*(lam - gam) / lam)
}

##Table of probability of a successful double spending
z_vec <- c(1, 2, 3, 4, 5) #Various values of z
r_vec <- c(2, 3, 4, 5) #Various values for the repartitrion of the computing power
gamInit <- 0.01

#NUmber of digits in the table
options(digits = 5)
p_exp <- sapply(r_vec, function(r) sapply(z_vec, function(z) p_tau_exp(z, lam_hat, lam_hat/r)))
p_gamma <- sapply(r_vec, function(r) sapply(z_vec, function(z) p_tau_gamma(z, lam_hat, r, gamInit)))
xtable(cbind(p_exp, p_gamma), include.rownames = FALSE, include.colnames = TRUE, digits = 4)

##Plot of the PDF and the CDF of the double spending time
#Truncation order
K <- 50
#The malicious chain is z block behind
z <- 1
#The honest miners mine twice as fast as the dishonest miners
r <- 2
#Limit of the X-axis which cporesponds to the time axis here
t_max <- 200

###PDF
file_name <- "PDFTau.png"
png(paste(my_path_output, file_name, sep = ""))
plot(seq(0, t_max, t_max / 1000),
     sapply(seq(0, t_max, t_max / 1000),function(t) pdf_tau_gamma(K, z, lam_hat, r, t))
     ,type = 'l',xlab = "Time (in minutes)", ylab = "", main= "", lwd = 2, ylim = c(0, 0.06),lty = 'twodash', cex.axis = 1.5,
     cex.lab = 1.5)
lines(seq(0, t_max, t_max / 1000),
      sapply(seq(0, t_max, t_max / 1000), function(t) pdf_tau_exp(K, z, lam_hat, lam_hat / r, t)),
      lty = 'solid', lwd = 2)
dev.off()

###CDF
file_name <- "CDFTau.png"
png(paste(my_path_output, file_name, sep = ""))
plot(seq(0, t_max, t_max / 1000), 
     sapply(seq(0, t_max, t_max / 1000), function(t) cdf_tau_exp(K, z, lam_hat, lam_hat / r, t)),
     type = 'l', xlab = "Time (in minutes)", ylab = "",main = "",lwd = 2, cex.axis = 1.5, cex.lab = 1.5,cex.main = 1.5)
abline(h = p_tau_exp(z, lam_hat, lam_hat / r))
lines(seq(0, t_max, t_max / 1000),
      sapply(seq(0, t_max,t_max / 1000),function(t) cdf_tau_gamma(K, z, lam_hat, r, t)),
      lty = 'twodash', lwd = 2)
abline(h = p_tau_gamma(z, lam_hat, r, gamInit),lty = 'twodash')
dev.off()

  
  





