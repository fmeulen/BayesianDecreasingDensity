# get directory of source script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(plyr)
library(ggplot2)
library(ggthemes)
library(scales)
library(gridExtra)
library(tidyverse)
library(latex2exp)
library(fdrtool) # for Grenander estimator

d_init <- read_csv("out/postmean.csv",  col_names = FALSE)

#d_init <- read_csv("./methodA/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodB/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodC-tau0-5/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodC-tau0-05/postmean.csv",col_names = FALSE)
d_init <- read_csv("./methodC-tau0-005/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodD/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodDtest/postmean.csv",col_names = FALSE)

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

true <- data.frame(x=grid,y=dexp(grid))

#titel=TeX("$g(\\theta) \\propto exp(- \\theta - 1 / \\theta)$")
#titel <- "Gamma(2,1)"
#titel <- "Pareto(1,0.5)"
#titel <- "Pareto(1,0.05)"
titel <- "Pareto(1,0.005)"
#titel <- TeX("Mixture Pareto(1,$\\tau$)")

pdf('postmean.pdf',width=3.5,height=2.5)
res %>% ggplot() +  geom_ribbon(aes(x=x,ymin = ql, ymax = qu), fill = "grey80") +
  geom_line(aes(x=x,y=m)) +ggtitle(titel)+coord_cartesian(xlim = c(0, 5),ylim=c(0,1.1))+
 geom_line(data=true,aes(x=x,y=y),colour='brown2',size=1.2,linetype='dashed')+
  xlab("")+ ylab("")+
  theme(plot.title = element_text(hjust = 0.5))+scale_y_continuous(breaks=seq(0,1.1,by=0.2))
dev.off()

# trace plot for tau (only interesting with mixture of Pareto (=method D))
post_tau <- read_csv("out/post_tau.csv", col_names = FALSE)
post_tau %>% mutate(iterate=1:nrow(post_tau)) %>% ggplot() + 
  geom_line(aes(x=iterate,y=X1)) + ylab(expression(tau))


post_tau_minBI <- post_tau[-(1:BI),]

post_tau_minBI %>% ggplot() + 
  geom_histogram(aes(x=X1),bins=50,fill='white',colour='black')


pareto_small <- FALSE
if (pareto_small==TRUE)
{

# add grenander estimator
dat <- read_csv("exp100.csv")$x # read the data

mle <- grenander(ecdf(dat),type="decreasing")
mle_df <- data.frame(x=mle$x.knots,y=mle$f.knots)

mle_df %>% ggplot() +
#  geom_step(mapping=aes(x=x, y=y), linetype=3) +
  geom_step(data=mle_df,mapping=aes(x=x, y=y), direction="vh") +
  geom_point( mapping=aes(x=x, y=y), color="red") 

# only in case of Pareta(0.005,1)
pdf('exp100-density.pdf',width=3.5,height=2.5)
res %>% ggplot() +  geom_ribbon(aes(x=x,ymin = ql, ymax = qu), fill = "grey80") +
  geom_line(aes(x=x,y=m)) +ggtitle(titel)+coord_cartesian(xlim = c(0, 4),ylim=c(0,3.5))+
  geom_line(data=true,aes(x=x,y=y),colour='brown2',size=1.2,linetype='dashed')+
    geom_step(data=mle_df,mapping=aes(x=x, y=y),
              direction="vh",colour='royalblue2',linetype=1,size=1.1)+ 
  xlab("")+ ylab("")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
}

