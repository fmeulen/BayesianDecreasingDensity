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
std_mhstep = 0.2 # sd for mh step updating θ

rmse(x, truedist) = sqrt( mean((x .- pdf(truedist,0.0)).^2) )

IT = 5#_000 # nr of iterations
NSIM = 10#00

methods = ["A", "B", "F"]
samplesizes = [50, 500, 5000]
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
df = DataFrame(iter = Float64[], f0=Float64[], method=String[], n=Int64[], true_distribution=String[], true0=Float64[])
for truedistribution in truedistributions
  for n in samplesizes
      for method in methods
          global df
          simrun = load("./out/exp2/"*label(truedistribution)*"_"*string(n)*"_"*string(method)*".jld2")["f0post"]
          n_sr = length(simrun)
          df = vcat(df,DataFrame(iter=1:n_sr, f0=simrun,
                                  method=fill(method,n_sr),
                                  n= fill(n,n_sr),
                                  true_distribution = fill(label(truedistribution),n_sr),
                                  true0 = fill(pdf(truedistribution,0.0),n_sr)
                                  ))
      end
  end
end

@rput IT
@rput df
R"""
print(class(df))
df %>% convert(chr(n,method,true_distribution))%>% mutate(method=fct_recode(method, "D" = "F")) %>% # recode F to D to align with article
 dplyr::filter(iter >IT/2) %>%
   ggplot() +
   geom_density(aes(x=f0,group=method,fill=method), alpha=0.4) +
      geom_vline(aes(xintercept= true0),colour="gray",linetype = "dashed")+facet_grid(true_distribution ~n)
"""
