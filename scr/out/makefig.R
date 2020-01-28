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

d_init <- read_csv("out/postsummary.csv")
grid = d_init$x
true <- data.frame(x=grid,y=dexp(grid))
iterates <- read_csv("out/iterates.csv") %>% gather(key="variable", value="value", posttau, postmean0)


iterates %>% ggplot(aes(x=iteratenr,y=value)) + geom_path() + facet_wrap(~variable, scales='free')
iterates %>% ggplot(aes(x=value)) + geom_histogram(fill='white',colour='black') + facet_wrap(~variable, scales='free')

#titel=TeX("$g(\\theta) \\propto exp(- \\theta - 1 / \\theta)$")
#titel <- "Gamma(2,1)"
#titel <- "Pareto(1,0.5)"
#titel <- "Pareto(1,0.05)"
#titel <- "Pareto(1,0.005)"
titel <- TeX("Mixture Pareto(1,$\\tau$)")

p <- d_init %>% ggplot() +  geom_ribbon(aes(x=x,ymin = lower, ymax = upper), fill = "grey80") +
  geom_line(aes(x=x,y=ave)) +ggtitle(titel)+coord_cartesian(xlim = c(0, 5),ylim=c(0,1.1))+
  geom_line(data=true,aes(x=x,y=y),colour='brown2',size=1.2,linetype='dashed')+
  xlab("")+ ylab("")+
  theme(plot.title = element_text(hjust = 0.5))+scale_y_continuous(breaks=seq(0,1.1,by=0.2))
p

pdf('postmean.pdf',width=3.5,height=2.5)
  show(p)
dev.off()

f0post100 <- read_csv("out/f0post100.csv")
f0post100 %>% ggplot() + geom_histogram(aes(x=f0post), fill="white",colour="black")


#d_init <- read_csv("./methodA/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodB/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodC-tau0-5/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodC-tau0-05/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodC-tau0-005/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodD/postmean.csv",col_names = FALSE)
#d_init <- read_csv("./methodDtest/postmean.csv",col_names = FALSE)

# columns correspond to different point on the horizontal axis,
# rows correspond to iterations
# first row is the grid










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

