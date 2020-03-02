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

IT = 2500 # nr of iterations
NSIM = 50 # 100

methods = ["A", "B", "F"]
samplesizes = [1000, 1500, 2500, 5000, 10_000, 12_500, 20_000]
#samplesizes = [1000, 5000, 10_000]
truedistributions = [Exponential(), HalfNormal()]
function label(truedist) # only needed for plotting
    if truedist == Exponential()
        return("Exp")
    end
    if truedist == HalfNormal()
      return("HalfNormal")
    end
end

grid0 =[0.0]
f0post = zeros(NSIM)
mean_acc = zeros(NSIM)



for truedistribution in [truedistributions[2]]
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

        # write results to csv file
        #CSV.write("./out/ratecomparison/f0post$n.csv",DataFrame(f0post=f0post))
        err = rmse(f0post, truedistribution)
        push!(res, [label(truedistribution), n, method, err, mean(mean_acc)])

        # just to be sure, write intermediate results
        # save("./out/ratecomparison/res_ratecomparison.jld", "res", res)
    end
  end
  # just to be sure, write intermediate results
  save("./out/ratecomparison/res_ratecomparison"*string(truedistribution)*".jld2", "res", res)
end


############# postprocessing ####################

res1 = load("./out/ratecomparison/res_ratecomparison"*string(truedistributions[1])*".jld2")["res"]
res2 = load("./out/ratecomparison/res_ratecomparison"*string(truedistributions[2])*".jld2")["res"]
res = vcat(res1,res2)
# combine res1 and res2 res =



# plotting
ec(X,i) = map(x->x[i],X)
df = DataFrame(true_distribution = ec(res,1), n=ec(res,2),
  method=ec(res,3), rmse=ec(res,4), log10n = log10.(ec(res,2)), log10_rmse=log10.(ec(res,4)),
   logn = log.(ec(res,2)), log_rmse=log.(ec(res,4)),
  meanacc=ec(res,5))
@rput df
R"""
p <- df %>% ggplot(aes(x=log10n, y = log10_rmse, colour=method)) + geom_point() +
        stat_smooth(method = "lm", se=FALSE) +
          facet_wrap(~true_distribution, scales='free') +
        theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),legend.position = 'bottom') +
        xlab(TeX("$\\log_{10}(n)$")) + ylab("RMSE(n)")
pdf("rate_at0_comparison.pdf",width=7,height=4)
  show(p)
dev.off()
"""

println("Exponential distribution:")
@rput df
R"""
dfAexp <- df %>% filter(method=="A", true_distribution=="Exp")
fitA <- lm(log10_rmse ~ log10n, data=dfAexp)
print(coef(fitA))
"""

R"""
dfBexp <- df %>% filter(method=="B", true_distribution=="Exp")
fitB <- lm(log10_rmse ~ log10n, data=dfBexp)
print(coef(fitB))
"""

R"""
dfFexp <- df %>% filter(method=="F", true_distribution=="Exp")
fitF <- lm(log10_rmse ~ log10n, data=dfFexp)
print(coef(fitF))
"""
println("HalfNormal distribution:")
R"""
dfA_HN <- df %>% filter(method=="A", true_distribution=="HalfNormal")
fitA_HN <- lm(log10_rmse ~ log10n, data=dfA_HN)
print(coef(fitA_HN))
"""

R"""
dfB_HN <- df %>% filter(method=="B", true_distribution=="HalfNormal")
fitB_HN <- lm(log10_rmse ~ log10n, data=dfB_HN)
print(coef(fitB_HN))
"""

R"""
dfF_HN <- df %>% filter(method=="F", true_distribution=="HalfNormal")
fitF_HN <- lm(log10_rmse ~ log10n, data=dfF_HN)
print(coef(fitF_HN))
"""
