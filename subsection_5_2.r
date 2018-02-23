#_______________________________________________________________________________________________________________________#
#Online accompaniement of the paper "Fraud risk assessment within block chain transactions"
#File: subsection_5_2.r
#Author: Pierre-O Goffard
#_______________________________________________________________________________________________________________________#

#Libraries
library(fitdistrplus)
library(moments)
library(xtable)

#The length of the chains are modeled by non-homogeneous Poisson process

##Integrated intensity function associated to the non-homogeneous Poisson processes
###The global hashrate has a parametric form with exp(at+b)
###L corresponds to the difficulty of the cryptopuzzle
Lambda=function(a, b, L, t){
  return( exp(b) * L/ 2^(256) / a * (exp(a * t)- 1 ) )
}

#The evaluation of the survival function of the double spending time is achieved via Monte-Carlo simulations.
#The estimator is named Appell Polynomial Monte-Carlo APMC estimator
#It requires the recursive computation of Appell Polynomials at x=1 for any order

##Computation of the Appell polynomial at x=1, denoted by Am(1|U)
### m is the degree
### U is a sequence of real numbers which uniquely define a sequence of Appell polynomials
A1=function(m,U){
  A0 <- c(1)
  A1 <- c(1)
  for(k in 1:m){
    A0[k+1]=-sum(sapply(1:k, function(j) choose(k,j) * (U[k])^(j) * A0[k-j+1])
    )
    A1[k+1]=sum(sapply(0:k, function(j) choose(k,j) * A0[k-j+1]))
  }
  return(A1[m+1])      
}


##Computation of the APMC estimator of the probability P(\tau_z\in(t1,t2))
### a,b,L relate to the intensity function of the block arrival
### We compute the probability that the double spending attack occurs during (t1,t2)
### 0 < p < 1 reflects the repartition of the computing powers among the miners
### The malicious chain is assumed to be z>0 blocks behind
### R is the number of trajectories generated for the Monte-Carlo evaluation of the expectation

APMC=function(a, b, L, t1, t2, p, z, R){
  #Length of the public chain
  Nt <- rpois( R, Lambda(a, b, p*L, t2) - Lambda(a, b, p*L, t1) )
  #Length of the dishonest chain
  Mt <- rpois( R, Lambda(a, b, (1-p)*L, t2) - Lambda(a, b, (1-p)*L, t1) )
  
  return(mean(sapply(1:R, function(r) if( Mt[r] >= (Nt[r] + z) ){#The dishonest chain upcrossed z+Nt
     0
    }else if(Mt[r]<z){#The dishonest chain is lower than z
       1
      }else{
          if(z == 1){
            A1(Mt[r], sort(runif(Nt[r], 0, 1))[1:Mt[r]])
          }else{#When z>1 some 0's must appended to the sequence U
             A1(Mt[r], c(rep(0,z - 1), sort(runif(Nt[r], 0, 1))[1:(Mt[r] + 1 - z)]))
            }
          }
  )))
}

#Computation of the APMC estimator of the probability P(\tau_z\in(t1,t2)) in the special case where z==1

APMC_z1=function(a, b, L, t1, t2, p, R){
  #Length of the public chain
  Nt <- rpois(R, Lambda(a, b, p*L, t2) - Lambda(a, b, p*L, t1))
  #Length of the private chain
  Mt <- rpois(R, Lambda(a, b, (1-p)*L, t2) - Lambda(a, b, (1-p)*L, t1))
  
  return(
    mean(sapply(1:R, function(r) if(Mt[r] >= (Nt[r] + 1)){
      0
    }else if(Mt[r] < 1){
        1
    }else{
        (Nt[r] - Mt[r] + 1) / (1 + Nt[r])
      }
    ))
  )
}
#____________________________________________________________________________________________________________________#

#Number of digits
options(digits = 15)

#Path for the output files 
my_path <- "C:/Users/goffard/Dropbox/goffard/BitCoins/"

#Parameters of the global hashrate function, see the first row of Table 1 in Bowden et al.
a <- -9.44 * 10^(-9)
b <- 27.1

#Horizon of projection in seconds
t <- 1209600 #which corresponds to 2 weeks
NbBlocks <- 2016 #Targeted number of blocks mined within two weeks

#Difficulty of the cryptopuzzles
L <- 2^(256) * a / exp(b) / (exp(a * 1209600) - 1) * NbBlocks #Maintain 2016 blocks mined in two weeks

#Plot of the integrated intensity function over time for varying computing power distribution (parameter p)
file_name <-"IntensityFunction.png"

png(paste(my_path, file_name, sep = ""))
par(mfrow = c(2, 2))

p <- 0.6
plot(seq(1, t, 100), sapply(seq(1, t, 100), function(t) Lambda(a, b, p * L, t)),
     type="l", ylim = c(0,2016),xlab = "Time (in seconds)", ylab="Intensity function", main=paste("p = ", p, sep = ""))
lines(seq(1, t, 100), sapply(seq(1,t,100), function(t) Lambda(a, b, (1-p) * L, t)),
      type = "l", lty = 'twodash')

p <- 0.7
plot(seq(1, t, 100), sapply(seq(1, t, 100), function(t) Lambda(a, b, p * L, t)),
     type = "l", ylim = c(0, 2016), xlab = "Time (in seconds)", ylab = "Intensity function", main = paste("p=", p, sep = ""))
lines(seq(1, t, 100), sapply(seq(1, t, 100), function(t) Lambda(a, b, (1-p) * L, t)),
      type = "l",lty = 'twodash')

p <- 0.8
plot(seq(1 , t, 100), sapply(seq(1, t, 100), function(t) Lambda(a, b, p * L , t)),
     type = "l",ylim = c(0, 2016),xlab = "Time (in seconds)",ylab = "Intensity function",main = paste("p = ", p, sep = ""))
lines(seq(1, t, 100), sapply(seq(1, t, 100), function(t) Lambda(a, b, (1-p) * L, t)),
      type = "l", lty = 'twodash')

p <- 0.9
plot(seq(1, t, 100), sapply(seq(1, t, 100), function(t) Lambda(a, b, p * L, t)),
     type = "l", ylim = c(0, 2016), xlab = "Time (in seconds)", ylab = "Intensity function", main = paste("p = ", p, sep=""))
lines(seq(1, t, 100), sapply(seq(1, t, 100), function(t) Lambda(a, b, (1-p) * L, t)),
      type = "l", lty = 'twodash')

dev.off()

#Probability of a successful double spending attack over time spans of 3 hours, for z = 1
##Vector of time steps 3 hours = 10800 over the course of 2 weeks = 1209600 seconds
time_steps <- seq(0,1209600,10800)
n_steps <- length(time_steps)

par(mfrow = c(1, 1))

#Number of simulated trajectory 10^5
R <- 10000
file_name <- "SuccesfulDoubleSpendingAttackNonHomogeneous.png"

png(paste(my_path, file_name, sep = ""))

p <- 0.6
plot(time_steps[2:n_steps], 
     sapply(2:n_steps, function(t)  1 - APMC_z1(a, b, L, time_steps[t-1], time_steps[t], p, R)),
     type = 'l', ylim = c(0, 0.8), xlab = 'Time (in seconds)', ylab = '',main = '',
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,lwd = 2)

p <- 0.7
lines(time_steps[2:n_steps],
      sapply(2:n_steps, function(t) 1 - APMC_z1(a, b, L, time_steps[t-1], time_steps[t], p, R)),
      type = 'l',lty = 'twodash', lwd = 2)

p <- 0.8
lines(time_steps[2:n_steps],
      sapply(2:n_steps, function(t) 1 - APMC_z1(a, b, L, time_steps[t-1], time_steps[t], p, R)),
      type = 'l',lty = 'dotted',lwd = 2)

p <- 0.9
lines(time_steps[2:length(time_steps)],
      sapply(2:length(time_steps), function(t) 1 - APMC_z1(a, b, L, time_steps[t-1], time_steps[t], p, R)),
      type = 'l',lty = 'dotdash',lwd = 2)

dev.off()




     