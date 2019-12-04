setwd("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/data4figs for lixue")

library(plyr)
library(ggplot2)
library(ggthemes)
library(scales)
library(gridExtra)
library(tidyverse)

# Exponential
e50 <- read_csv("exponential/postmean_atzero50.csv")
e200 <- read_csv("exponential/postmean_atzero200.csv")
e10000 <- read_csv("exponential/postmean_atzero10000.csv")

em50 <- colMeans(e50)
em200 <- colMeans(e200)
em10000 <- colMeans(e10000)[1:50]   # I originally did 100 runs, we only need 50

boxplot(em50,em200,em10000)

# half-Cauchy
hc50 <- read_csv("half cauchy/postmean_atzero50.csv")
hc200 <- read_csv("half cauchy/postmean_atzero200.csv")
hc10000 <- read_csv("half cauchy/postmean_atzero10000.csv")

hcm50 <- colMeans(hc50)
hcm200 <- colMeans(hc200)
hcm10000 <- colMeans(hc10000)

boxplot(hcm50,hcm200,hcm10000)
