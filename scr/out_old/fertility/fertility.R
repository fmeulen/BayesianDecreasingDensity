setwd("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/fertility")

library(plyr)
library(ggplot2)
library(ggthemes)
library(scales)
library(gridExtra)
library(tidyverse)
library(latex2exp)

# read the current duration data
d_init <- read_csv("~/Sync/DOCUMENTS/onderzoek/LiXue/bayesdec_julia/fertility/postmean.csv",col_names = FALSE)
# columns correspond to different point on the horizontal axis,
# rows correspond to iterations
# first row is the grid

d <- d_init[-1,]
grid <- as.numeric(d_init[1,])

# trace plots at two locations
d %>% mutate(iterate=1:nrow(d)) %>% ggplot() + geom_line(aes(x=iterate,y=X1))
d %>% mutate(iterate=1:nrow(d)) %>% ggplot() + geom_line(aes(x=iterate,y=X40))

# density plot
BI <- round(nrow(d_init)/2,0)
dminBI <- d_init[-(1:BI),]
dmean = apply(dminBI,2,mean)
dquantl <- apply(dminBI,2,quantile,probs=0.025)
dquantu <- apply(dminBI,2,quantile,probs=0.975)

res <- data.frame(x=grid,m=dmean,ql=dquantl,qu=dquantu)
res


pdf('density.pdf',width=3.5,height=3.5)
res %>% ggplot() +  geom_ribbon(aes(x=x,ymin = ql, ymax = qu), fill = "grey80") +
  geom_line(aes(x=x,y=m)) +  ylab(TeX("f_0(x)")) + xlab("x")#+coord_cartesian(xlim = c(0, 5))+  #+ #+ geom_rug(data=x,aes(x=V1))
dev.off()

# transform to distribution function estimate
dminBI_df <- dminBI
for (i in 1:ncol(dminBI))   dminBI_df[,i] = 1-dminBI[,i]/dminBI[,1]

# distribution plot
dmean_df = apply(dminBI_df,2,mean)
dquantl_df <- apply(dminBI_df,2,quantile,probs=0.025)
dquantu_df <- apply(dminBI_df,2,quantile,probs=0.975)

res_df <- data.frame(x=grid,m=dmean_df,ql=dquantl_df,qu=dquantu_df)
res_df

pdf('implied_df.pdf',width=3.5,height=3.5)
res_df %>% ggplot() +  geom_ribbon(aes(x=x,ymin = ql, ymax = qu), fill = "grey80") +
  geom_line(aes(x=x,y=m)) + ylab(TeX("H_0(x)")) + xlab("x")# +coord_cartesian(xlim = c(0, 40))
#+xlim(0,5)
dev.off()

