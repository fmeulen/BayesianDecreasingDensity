using Distributions
using DelimitedFiles
using QuadGK
using LinearAlgebra
using CSV
using DataFrames

workdir = @__DIR__
cd(workdir)
include("bd_funcdefs.jl")

#------------------ specify the prior ------------------
α = 1.0                       # concentration par
method = "D"

# specify simulation settings
std_mhstep = 0.2 # sd for mh step updating θ


start = time()

IT = 2500 # nr of iterations
NSIM = 50 # 100

# specify data generation
n = 100
truedistribution = ["Exponential(1)", "HalfNormal"][1]

f0post = zeros(NSIM)
mean_acc = zeros(NSIM)

grid0 =[0.0]

# do one preliminary run to obtain initial configuration
x = simdata(n, truedistribution)
outprelim = mcmc(x, method, α, IT, std_mhstep, grid0)
println(outprelim.mean_acc)

cpre = outprelim.config_end

println("Start Monte-Carlo study")
for k in 1:NSIM
  print(k)
  x = simdata(n, truedistribution)
  sim = mcmc(x, method, α, IT, std_mhstep, grid0, config_init=cpre)
  f0post[k], mean_acc[k] = sim.dout.ave[1], sim.mean_acc
end

# write results to csv file
CSV.write("./out/f0post$n.csv",DataFrame(f0post=f0post))

elapsedtime = time() - start
