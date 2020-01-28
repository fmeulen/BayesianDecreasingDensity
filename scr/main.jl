using Distributions
using DelimitedFiles
using QuadGK
using LinearAlgebra
using CSV
using DataFrames
using RCall

R"""
    library(plyr)
    library(ggplot2)
    library(ggthemes)
    library(scales)
    library(gridExtra)
    library(tidyverse)
    library(latex2exp)
    library(fdrtool) # for Grenander estimator
    theme_set(theme_light())
"""

workdir = @__DIR__
include("bd_funcdefs.jl")

# make directory for figures
figdir = workdir*"/out/fig/"
mkpath(figdir)

#------------------ specify the prior ------------------
α = 1.0                       # con centration par

datasets = ["curdur", "exp100", "halfnormal100"]
methods = ["A", "B", "C", "D", "E", "F"]

# specify simulation settings
std_mhstep = 0.7 # sd for mh step updating θ
IT = 250#25_000

#------------------  read data ------------------
for dataset in datasets[[2,3]], method in methods
    println(dataset)
    println(method)
    println("----")


    if dataset=="curdur"
        y = readdlm("datasets/curdur.csv",',', header=true)[1][:,2]
        x = sort(y)[1:618]
        grid = range(0,36;length=300)
    elseif dataset=="exp100"
        x = readdlm("datasets/exp100.csv",',',header=true)[1][:,2]#[2:end,2]
        grid = range(0,5;length=100)
    elseif dataset=="halfnormal100"
        x = vec(readdlm("datasets/halfnormal100.csv",','))
        grid = range(0,5;length=100)
    end

    #------------------ mcmc by Neal ------------------
    start = time()
        out = mcmc(x, method, α, IT, std_mhstep, grid)
    elapsedtime = time() - start

    #------------------ postprocessing ------------------
    outdir =  workdir*"/out/"*dataset*"_"*method*"/"
    mkpath(outdir)

    # write results to csv file
    CSV.write(outdir*"postsummary.csv",out.dout)
    CSV.write(outdir*"iterates.csv",out.iterates)

    # save info to file
    facc = open(outdir*"info.txt","w")
    write(facc, "Average acceptance probability equals: ",string(round(out.mean_acc;digits=3)),"\n")
    write(facc, "Number of iterations: ",string(IT),"\n")
    write(facc, "std_mhstep = ",string(std_mhstep),"\n\n")
    write(facc, "elapsed time ",string(elapsedtime), "\n\n")
    write(facc, "---- Prior specification ----","\n")
    write(facc, "method: ", string(method),"\n")
    write(facc, "concentration parameter: ", string(α),"\n")
    write(facc, "extra_pars= ",string(out.ep),"\n")
    close(facc)

    d = out.dout
    if dataset=="curdur"
        newcol= 0.0 * d[:x]   # adjust later, true density is unknown
    elseif dataset=="exp100"
        newcol = pdf(Exponential(),d[:x])
    elseif dataset=="halfnormal100"
        newcol = 2.0*pdf(Normal(),d[:x])
    end
    insertcols!(d, 5, :truedens => newcol) # add true val of density

    #fn = outdir*"pm_"*dataset*"_"*method*".pdf"
    fn = figdir * "pm_"*dataset*"_"*method*".pdf"

    @rput fn # filename
    @rput d # dataframe with estimates
    @rput method
    @rput x # the data
    R"""
        if (method=="A")
                {titel=TeX("$g(\\theta) \\propto exp(- \\theta - 1 / \\theta)$")}
        if (method=="B")
            {titel <- "Gamma(2,1)"}
        if (method=="C")
            {titel <- "Pareto(1,0.5)"}
        if (method=="D")
            {titel <- "Pareto(1,0.05)"}
        if (method=="E")
            {titel <- "Pareto(1,0.005)"}
        if (method=="F")
            {titel <- TeX("Mixture Pareto(1,$\\tau$)")}

        p <- d %>% ggplot() +  geom_ribbon(aes(x=x,ymin = lower, ymax = upper), fill = "grey80") +
            geom_line(aes(x=x,y=ave)) +ggtitle(titel)+coord_cartesian(xlim = c(0, 5),ylim=c(0,1.1))+
            geom_line(aes(x=x,y=truedens),colour='brown2',size=1.2,linetype='dashed')+
            xlab("")+ ylab("")+
            theme(plot.title = element_text(hjust = 0.5))#+scale_y_continuous(breaks=seq(0,1.1,by=0.2))
        if (method=="E") # add Grenander estimator
        {
            mle <- grenander(ecdf(x),type="decreasing")
            mle_df <- data.frame(x=mle$x.knots,y=mle$f.knots)
            p <- p + geom_step(data=mle_df,mapping=aes(x=x, y=y),
                      direction="vh",colour='royalblue2',linetype=1,size=1.1)
        }
        pdf(fn,width=3.5,height=2.5)
            show(p)
        dev.off()

    """


end
