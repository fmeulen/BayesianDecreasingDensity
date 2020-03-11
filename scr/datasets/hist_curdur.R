setwd("~/.julia/dev/BayesianDecreasingDensity/scr/datasets")
library(tidyverse)
theme_set(theme_light())

curdur <- read_csv("curdur.csv") %>% arrange(V1) %>% slice(1:618)


p <- curdur %>% ggplot(aes(x=V1)) + geom_histogram(fill='white',colour='black') + xlab("")

pdf("curdur_histogram.pdf",width =7, height=2.5)
show(p)
dev.off()