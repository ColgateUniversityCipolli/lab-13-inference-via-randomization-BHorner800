###########################################
# Ben Horner
# Lab 13
# Math 240
###########################################
###########################################
# Libraries
###########################################
library(tidyverse)
library(patchwork)
library(e1071)

###########################################
# Question 1
###########################################
finch.dat = read.csv('./zebrafinches.csv')

#potential error in p-value
n = length(finch.dat$further)
t.further = t.test(finch.dat$further, alternative = 'less',
                   mu = 0)$statistic[1]
gauss.pdf  = dnorm(finch.dat$further)
finch.skew = skewness(finch.dat$further)
  
error = gauss.pdf*finch.skew*(2*t.further^2 + 1)/(sqrt(n)*6)

#t-statistic error from -10 to 10
t.errors = c()
for (i in -10:10){
  t = i
  error = gauss.pdf*finch.skew*(2*t^2 + 1)/(sqrt(n)*6)
  t.errors = append(t.errors, mean(error))
}
plot(-10:10, t.errors, 
     xlab = 't statistic',
     ylab = 'error in t',
     main = 'A changing error in t statistic')


