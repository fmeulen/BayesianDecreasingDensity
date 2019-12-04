setwd("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodA")

library(plyr)
library(ggplot2)
library(ggthemes)
library(scales)
library(gridExtra)
library(tidyverse)


theme_set(theme_minimal())

d <- read_csv("postmean_atzero10000.csv")
IT <- nrow(d)
BI <- round(nrow(d)/2,0)

f0true <- 1
f0post <- colMeans(d[BI:IT,])

PL <- TRUE

if (PL==TRUE)
{
  d <- d %>% mutate(iter=1:IT)
  colnames(d)[1] <- "postmean" # then plot the result for the first simulation run
  
  postmean_minBI <- d$postmean[BI:nrow(d)]
  autocor <- acf(postmean_minBI,plot=FALSE)
  d_acf <- data.frame(lag=autocor$lag, acf=autocor$acf)
  
  p1 <- d %>% ggplot(aes(x=iter,y=postmean)) + 
    geom_line()+geom_hline(yintercept = mean(postmean_minBI),colour='yellow')+
    xlab("iteration")+ylab("posterior mean")
  
  p2 <- ggplot(data = d_acf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0))
  
  grid.arrange(p1,p2,ncol=2)
}



print(f0post)
print(log(sqrt(mean(f0post-f0true)^2)))

