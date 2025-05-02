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
library(boot.pval)

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
  curr.resample <- sample(finch.dat$`diff`,
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
boots <- boot(data = finch.dat$`diff`,
              statistic = boot.mean,
              R = R)
boot.ci(boots, type="bca")

# shift so H0 is true
mean(resamples$xbars)
(delta <- mean(resamples$xbars) - 2.9)


resamples <- resamples |>
  mutate(xbars.shifted = xbars - delta)

low <- 2.9 - delta
high <- 2.9 + delta

resamples |>
  summarize(mean = mean(xbars.shifted),
            p.low = mean(xbars.shifted <= low),
            p.high = mean(xbars.shifted >= high))|>
  mutate(p = p.low + p.high) |>
  view()


boot.pval(boots, theta_null = mu0) # this does something slightly different




####Further####
R <- 10000
resamples <- tibble(xbars = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finch.dat$`further`,
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
boots <- boot(data = finch.dat$`further`,
              statistic = boot.mean,
              R = R)
boot.ci(boots, type="bca")

# shift so H0 is true
mean(resamples$xbars)
(delta <- mean(resamples$xbars) - 2.9)


resamples <- resamples |>
  mutate(xbars.shifted = xbars - delta)

low <- 2.9 - delta
high <- 2.9 + delta

resamples |>
  summarize(mean = mean(xbars.shifted),
            p.low = mean(xbars.shifted <= low),
            p.high = mean(xbars.shifted >= high))|>
  mutate(p = p.low + p.high) |>
  view()


boot.pval(boots, theta_null = mu0) # this does something slightly different


####Diff####
R <- 10000
resamples <- tibble(xbars = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finch.dat$`diff`,
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
boots <- boot(data = finch.dat$`diff`,
              statistic = boot.mean,
              R = R)
boot.ci(boots, type="bca")

# shift so H0 is true
mean(resamples$xbars)
(delta <- mean(resamples$xbars) - 2.9)


resamples <- resamples |>
  mutate(xbars.shifted = xbars - delta)

low <- 2.9 - delta
high <- 2.9 + delta

resamples |>
  summarize(mean = mean(xbars.shifted),
            p.low = mean(xbars.shifted <= low),
            p.high = mean(xbars.shifted >= high))|>
  mutate(p = p.low + p.high) |>
  view()


boot.pval(boots, theta_null = mu0) # this does something slightly different

###########################################
# Question 3: Randomization Test
###########################################
mu0 = 0
#Closer
R <- 10000
rand <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- finch.dat$`closer` - mu0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  
  rand$xbars[i] <- mean(curr.rand)
}
# Thinking is hard
rand <- rand |>
  mutate(xbars = xbars + mu0) # shifting back

# p-value 
(delta <- abs(mean(finch.dat$`closer`) - mu0))
(low <- mu0 - delta) # mirror
(high<- mu0 + delta)   # xbar

mean(rand$xbars <= low) +
  mean(rand$xbars >= high)

#Confidence Interval
mu0.iterate <- 0.01
starting.point <- mean(finch.dat$`closer`)

mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- finch.dat$`closer` - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(finch.dat$`closer`) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- finch.dat$`closer` - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(finch.dat$`closer`) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

c(mu.lower, mu.upper)







#Further
R <- 10000
rand <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- finch.dat$`further` - mu0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  
  rand$xbars[i] <- mean(curr.rand)
}
# Thinking is hard
rand <- rand |>
  mutate(xbars = xbars + mu0) # shifting back

# p-value 
(delta <- abs(mean(finch.dat$`further`) - mu0))
(low <- mu0 - delta) # mirror
(high<- mu0 + delta)   # xbar

mean(rand$xbars <= low) +
  mean(rand$xbars >= high)

#Confidence Interval
mu0.iterate <- 0.01
starting.point <- mean(finch.dat$`further`)

mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- finch.dat$`further` - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(finch.dat$`further`) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- finch.dat$`further` - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(finch.dat$`further`) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

c(mu.lower, mu.upper)







#Difference
R <- 10000
rand <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- finch.dat$`diff` - mu0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  
  rand$xbars[i] <- mean(curr.rand)
}
# Thinking is hard
rand <- rand |>
  mutate(xbars = xbars + mu0) # shifting back

# p-value 
(diff.delta <- abs(mean(finch.dat$`diff`) - mu0))
(diff.low <- mu0 - delta) # mirror
(diff.high<- mu0 + delta)   # xbar

diff.mean = mean(rand$xbars <= low) +
  mean(rand$xbars >= high)
diff.mean  

# p-value 
(delta <- abs(mean(finch.dat$`diff`) - mu0))
(low <- mu0 - delta) # mirror
(high<- mu0 + delta)   # xbar

mean(rand$xbars <= low) +
  mean(rand$xbars >= high)

#Confidence Interval
mu0.iterate <- 0.01
starting.point <- mean(finch.dat$`diff`)

mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- finch.dat$`diff` - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(finch.dat$`diff`) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- finch.dat$`diff` - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(finch.dat$`diff`) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

c(mu.lower, mu.upper)

t.test(finch.dat$closer)
