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
library(boot)

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
par(mar = c(4, 4, 2, 1))
plot(seq(-10, 10, by = 0.01), t.errors, 
     xlab = 't statistic',
     ylab = 'error in t',
     main = 'A changing error in t statistic')

#Sample size needed
a= 0.05
gauss.pdf  = dnorm(t.further)
n = ((finch.skew/(6*0.10*a))*(2*t.further^2 + 1)*gauss.pdf)^2
n

t.further
gauss.pdf
finch.dat$further


###########################################
# Question 2
###########################################
R <- 10000
resamples <- tibble(xbars = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finch.dat$`closer`,
                          size = nrow(finch.dat),
                          replace = T)
  
  resamples$xbars[i] <- mean(curr.resample)
}

# Confidence Interval
quantile(resamples$xbars, c(0.025, 0.975))

library(boot)
boot.mean <- function(d, i){
  mean(d[i])
}
boots <- boot(data = finch.dat$`closer`,
              statistic = boot.mean,
              R = R)
boot.ci(boots, type="bca")

#Shifted Resamples
delta = mean(resamples$xbars)
resamples.null.closer = data.frame(xbars.shifted = resamples$xbars - delta)

bootstrap.plot <- ggplot() +
  # Plot Resampling Distribution
  stat_density(data=resamples, aes(x=xbars), geom="line", color="lightgrey")+
  stat_density(data=resamples.null.closer, aes(x=xbars.shifted), geom="line", color="black")+
  geom_hline(yintercept=0)+

  scale_color_manual("",
                     values=c("black", "red"))+
  # Plot Observation
  # Clean up
  theme_bw()+

  ylab("Density")+
  ggtitle("Bootstrap Hypothesis Test")


