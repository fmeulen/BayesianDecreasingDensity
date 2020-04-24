using Distributions
using DelimitedFiles
using QuadGK
using LinearAlgebra
using CSV
using DataFrames
using RCall
using JLD2
using FileIO

R"""
    library(ggplot2)
    library(ggthemes)
    library(scales)
    library(gridExtra)
    library(tidyverse)
    library(latex2exp)
    library(hablar)
    theme_set(theme_light())
"""


workdir = @__DIR__
cd(workdir)
include("bd_funcdefs.jl")

#------------------ specify the prior ------------------
α = 1.0                       # concentration par

# specify simulation settings
std_mhstep = 0.4 # sd for mh step updating θ

rmse(x, truedist) = sqrt( mean((x .- pdf(truedist,0.0)).^2) )

IT = 10_000 # nr of iterations
NSIM = 500

methods = ["A", "B", "F"]
samplesizes = [50, 250]
truedistributions = [Exponential(), HalfNormal()]
function label(truedist) # only needed for plotting
    if truedist == Exponential()
        return("Exponential")
    end
    if truedist == HalfNormal()
      return("HalfNormal")
    end
end

grid0 =[0.0]
f0post = zeros(NSIM)
mean_acc = zeros(NSIM)



for truedistribution in truedistributions
  res = []
  for n in samplesizes
      for method in methods
        println(truedistribution, "   ", n, "    ", method)
        # do one preliminary run to obtain initial configuration
        x = rand(truedistribution, n)     #simdata(n, truedistribution)
        outprelim = mcmc(x, method, α, IT, std_mhstep, grid0)
        #println(outprelim.mean_acc)

        cpre = outprelim.config_end
        for k in 1:NSIM
          println(k)
          x = rand(truedistribution, n)
          sim = mcmc(x, method, α, IT, std_mhstep, grid0, config_init=cpre)
          f0post[k], mean_acc[k] = sim.dout.ave[1], sim.mean_acc
        end
        save("./out/exp2/"*label(truedistribution)*"_"*string(n)*"_"*string(method)*".jld2", "f0post", f0post)
    end
  end
end


############# postprocessing ####################
df = DataFrame(simnr = Float64[], f0=Float64[], method=String[], n=Int64[], true_distribution=String[], true0=Float64[])
for truedistribution in truedistributions
  for n in samplesizes
      for method in methods
          global df
          simrun = load("./out/exp2/"*label(truedistribution)*"_"*string(n)*"_"*string(method)*".jld2")["f0post"]
          df = vcat(df,DataFrame(simnr=1:NSIM, f0=simrun,
                                  method=fill(method,NSIM),
                                  n= fill(n,NSIM),
                                  true_distribution = fill(label(truedistribution),NSIM),
                                  true0 = fill(pdf(truedistribution,0.0),NSIM)
                                  ))
      end
  end
end


@rput df
R"""
df1 <- df %>% convert(chr(n,method,true_distribution))%>% mutate(method=fct_recode(method, "D" = "F")) %>% # recode F to D to align with article
    mutate(n=fct_inseq(n))

p1 <- df1 %>% ggplot() +
   geom_density(aes(x=f0,group=method,fill=method), alpha=0.2,colour='grey') +
      geom_vline(aes(xintercept= true0),colour="gray",linetype = "dashed") +
      facet_grid(n ~ true_distribution, scales="free")+xlab("")


pdf("res_exp2A.pdf",width=7.5,height=3.8)
        show(p1)
dev.off()

#p2 <- df1 %>%  ggplot() +
#   geom_boxplot(aes(f0,method, group=method), alpha=0.2) +
      #geom_vline(aes(xintercept= true0),colour="gray",linetype = "dashed")# +
      #facet_grid(n ~ true_distribution, scales="free")+xlab("")
 #pdf("res_exp2B.pdf",width=7.5,height=4.5)
#              show(p2)
#      dev.off()
"""
