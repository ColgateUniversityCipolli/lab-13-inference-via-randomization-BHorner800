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
gauss.pdf  = dnorm(t.further)
finch.skew = skewness(finch.dat$further)
  
error = gauss.pdf*finch.skew*(2*t.further^2 + 1)/(sqrt(n)*6)

#t-statistic error from -10 to 10
t.errors = c()
for (i in seq(-10, 10, by = 0.01)){
  t = i
  gauss.pdf  = dnorm(t)
  error = gauss.pdf*finch.skew*(2*t^2 + 1)/(sqrt(n)*6)
  t.errors = append(t.errors, error)
}
plot(seq(-10, 10, by = 0.01), t.errors, 
     xlab = 't statistic',
     ylab = 'error in t',
     main = 'A changing error in t statistic')

#Sample size needed
a= 0.05
gauss.pdf  = dnorm(t.further)
n = ((finch.skew/(6*0.10*a))*(2*t.further^2 + 1)*gauss.pdf)^2
n
