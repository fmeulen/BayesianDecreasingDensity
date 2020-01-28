setwd("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/")

library(plyr)
library(ggplot2)
library(ggthemes)
library(scales)
library(gridExtra)
library(tidyverse)
library(latex2exp)


theme_set(theme_minimal())


# read data
A1000 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodA/postmean_atzero1000.csv", 
                  col_names = FALSE, skip = 1)
A1500 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodA/postmean_atzero1000.csv", # adjust
                  col_names = FALSE, skip = 1)
A2500 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodA/postmean_atzero2500.csv", 
                  col_names = FALSE, skip = 1)
A5000 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodA/postmean_atzero5000.csv", 
                  col_names = FALSE, skip = 1)
A10000 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodA/postmean_atzero10000.csv", 
                   col_names = FALSE, skip = 1)
A12500 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodA/postmean_atzero12500.csv", 
                   col_names = FALSE, skip = 1)
A20000 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodA/postmean_atzero20000.csv", 
                   col_names = FALSE, skip = 1)


B1000 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodB/postmean_atzero1000.csv", 
                                     col_names = FALSE, skip = 1)
B1500 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodB/postmean_atzero1500.csv", 
                  col_names = FALSE, skip = 1)
B2500 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodB/postmean_atzero2500.csv", 
                  col_names = FALSE, skip = 1)
B5000 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodB/postmean_atzero5000.csv", 
                  col_names = FALSE, skip = 1)
B10000 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodB/postmean_atzero10000.csv", 
                  col_names = FALSE, skip = 1)
B12500 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodB/postmean_atzero12500.csv", 
                   col_names = FALSE, skip = 1)
B20000 <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/methodB/postmean_atzero20000.csv", 
                   col_names = FALSE, skip = 1)


RMSE <- function(x) sqrt(mean((x-1)^2))  # true value at zero equals 1

As1000 <- RMSE(colMeans(A1000))
As1500 <- RMSE(colMeans(A1500))
As2500 <- RMSE(colMeans(A2500))
As5000 <- RMSE(colMeans(A5000))
As10000 <- RMSE(colMeans(A10000))
As12500 <- RMSE(colMeans(A12500))
As20000 <- RMSE(colMeans(A20000))

Bs1000 <- RMSE(colMeans(B1000))
Bs1500 <- RMSE(colMeans(B1500))
Bs2500 <- RMSE(colMeans(B2500))
Bs5000 <- RMSE(colMeans(B5000))
Bs10000 <- RMSE(colMeans(B10000))
Bs12500 <- RMSE(colMeans(B12500))
Bs20000 <- RMSE(colMeans(B20000))
# construct dataframe
summ <- data.frame(n=rep(c(1000,1500,2500,5000,10000,12500,20000),each=2),
                    RMSE=c(As1000,Bs1000,As1500,Bs1500,As2500,Bs2500,
                           As5000,Bs5000,As10000,Bs10000,As12500,Bs12500,
                           As20000,Bs20000),
                    method=rep(c("A","B"),7))
summ <- summ %>% mutate(logn=log(n),logRMSE=log(RMSE))

# compute slopes
Asumm <-  summ %>% filter(method=="A")
Afit <- lm(logRMSE ~ logn, data=Asumm)
Afitlarge <- lm(logRMSE ~ logn, data=Asumm[4:7,])
Aslope <- Afit$coef[2]
Aslopelarge <- Afitlarge$coef[2]
summary(Afit)
summary(Afitlarge)

Bsumm <-  summ %>% filter(method=="B")
Bfit <- lm(logRMSE ~ logn, data=Bsumm)
Bfitlarge <- lm(logRMSE ~ logn, data=Bsumm[4:7,])
Bslope <- Bfit$coef[2]
Bslopelarge <- Bfitlarge$coef[2]
summary(Bfit)
summary(Bfitlarge)

# make figs

Btitel=TeX("$g(\\theta)  \\propto \\theta \\, exp(- \\theta)$")
# make fig
# pdf('empirical_rate_gamma.pdf',width=5,height=5)
# ggplot(data=Bsumm,aes(x=logn, y=logRMSE)) + geom_point(size=3)+   
#   geom_smooth(method='lm',se=FALSE) + xlab("log(n)") + ylab("log(RMSE)")+
#   geom_smooth(data=Bsumm[4:6,],aes(x=logn, y=logRMSE),method='lm',se=FALSE,colour='green') + xlab("log(n)") + ylab("log(RMSE)")+
#   annotate("text", x = 8, y = -1.9, label = paste('slope = ',round(Bslope,2),sep=''),size=5)+ggtitle(Btitel)+
#   labs(title=Btitel, subtitle=paste('slope = ',round(Bslope,3),sep=''))+ 
#       theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
# dev.off()
Atitel=TeX("$g(\\theta) \\propto exp(- \\theta - 1 / \\theta)$")

# joint fig
titel <- TeX("method A: $g(\\theta) \\propto exp(- \\theta - 1 / \\theta)$,  method B: $g(\\theta)  \\propto \\theta \\, exp(- \\theta)$")
p <- ggplot(data=summ,aes(x=logn, y=logRMSE, colour=method)) + geom_point(size=3)+   
  geom_smooth(method='lm',se=FALSE) + xlab("log(n)") + ylab("log(RMSE)")+
  geom_smooth(data=Asumm[4:7,],aes(x=logn, y=logRMSE,colour=method),method='lm',se=FALSE,colour="#F8766D",linetype='dashed') + xlab("log(n)") + ylab("log(RMSE)")+
  geom_smooth(data=Bsumm[4:7,],aes(x=logn, y=logRMSE,colour=method),method='lm',se=FALSE,colour="#00BFC4",linetype='dashed') + xlab("log(n)") + ylab("log(RMSE)")+
  annotate("text", x = 9.5, y = -1.7, label = paste('slope = ',round(Aslope,2),'/',round(Aslopelarge,2),sep=''),size=5,colour="#F8766D")+
  annotate("text", x = 7.25, y = -2.7, label = paste('slope = ',round(Bslope,2),'/',round(Bslopelarge,2),sep=''),size=5,colour="#00BFC4")+
  labs(title=titel)+ 
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.position = 'bottom')

p
cat("Aslopelarge ",Aslopelarge,"\n")
cat("Bslopelarge ",Bslopelarge,"\n")

pdf('rate-comparison.pdf',width=8,height=5)
show(p)
dev.off()
