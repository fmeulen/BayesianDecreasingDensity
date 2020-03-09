# file for reproducing the experiments in section 5.1 of the paper. 

using Distributions
using DelimitedFiles
using QuadGK
using LinearAlgebra
using CSV
using DataFrames
using RCall

R"""
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
figdir = workdir*"/out/exp1/"
mkpath(figdir)

#------------------ specify the prior ------------------
α = 1.0                       # con centration par

datasets = ["curdur", "exp100", "halfnormal100", "exp1000", "halfnormal1000"]
methods = ["A", "B", "C", "D", "E", "F"]

# specify simulation settings
std_mhstep = 0.4 # sd for mh step updating θ
IT = 25_000


#------------------  read data ------------------
#for dataset in datasets[[2,3]], method in methods
for dataset in datasets[[4,5]], method in methods
    println(dataset)
    println(method)
    println("----")
    if dataset=="curdur"
        fn = joinpath(workdir, "datasets/curdur.csv")
        y = readdlm(fn,',', header=true)[1][:,2]
        x = sort(y)[1:618]
        grid = range(0,36;length=300)
    elseif dataset=="exp100"
        fn = joinpath(workdir, "datasets/exp100.csv")
        x = readdlm(fn,',',header=true)[1][:,2]#[2:end,2]
        grid = range(0,5;length=100)
    elseif dataset=="exp1000"
        fn = joinpath(workdir, "datasets/exp1000.csv")
        x = vec(readdlm(fn,','))
        grid = range(0,5;length=100)
    elseif dataset=="halfnormal100"
        fn = joinpath(workdir, "datasets/halfnormal100.csv")
        x = vec(readdlm(fn,','))
        grid = range(0,5;length=100)
    elseif dataset=="halfnormal1000"
        fn = joinpath(workdir, "datasets/halfnormal1000.csv")
        x = vec(readdlm(fn,','))
        grid = range(0,5;length=100)
    end

    #------------------ mcmc by Neal ------------------
    start = time()
        out = mcmc(x, method, α, IT, std_mhstep, grid)
    elapsedtime = time() - start

    #------------------ postprocessing ------------------
    outdir =  workdir*"/out/exp1/"*dataset*"_"*method*"/"
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
    # plotting
    d = out.dout
    if dataset=="curdur"
        newcol= 0.0 * d[:x]   # adjust later, true density is unknown
        ymax = 0.5  # check if reasonable
    elseif dataset in ["exp100", "exp1000"]
        newcol = pdf.(Exponential(),d[!,:x])
        ymax = 1.2
    elseif dataset in ["halfnormal100", "halfnormal1000"]
        newcol = 2.0*pdf.(Normal(),d[!,:x])
        ymax = 1.0
    end
    insertcols!(d, 5, :truedens => newcol) # add true val of density

    #fn = outdir*"pm_"*dataset*"_"*method*".pdf"
    fn = figdir * "pm_"*dataset*"_"*method*".pdf"

    @rput fn # filename
    @rput d # dataframe with estimates
    @rput method
    @rput x # the data
    @rput ymax # max value on y-axis
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
            geom_line(aes(x=x,y=ave)) +ggtitle(titel)+
            coord_cartesian(xlim = c(0, 5),ylim=c(0,ymax))+
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
