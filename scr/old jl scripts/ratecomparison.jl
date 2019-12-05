using Distributions
using DelimitedFiles
using QuadGK
using LinearAlgebra
using CSV
using DataFrames

workdir = @__DIR__
println(workdir)
cd(workdir)

include("bd_funcdefs.jl")
# At each iteration I keep track of the configuration, for that I define the
# following data-structure


# specify simulation settings
mh_step = 0.06 # sd for mh step updating θ
IT = 2500 # nr of iterations
NSIM = 50 # 100

# specify data generation
true_distribution = Exponential(1)
n= 10000

# do the simulations
f0_post = zeros(IT,NSIM)
mean_acc = zeros(NSIM)
p0 = α * quadgk(Ψ(0), 0, Inf)[1]

# do one preliminary run to obtain initial configuration
tic()
#x = rand(true_distribution,n)
x = abs.(randn(n))  # half normal distr
config_simple = Config(ones(n),[1], [n],1,[maximum(x)+1])
f0_post_simple, mean_acc_simple, config_init = postmean0_simulation(x,20000, mh_step, α, p0, method, config_simple)
println(mean_acc_simple)

println("Start Monte-Carlo study")
for k in 1:NSIM
  print(k)
  #x = rand(true_distribution,n)
  x = abs.(randn(n))
  f0_post[:,k], mean_acc[k], config_last = postmean0_simulation(x,IT, mh_step, α, p0, method, config_init)
end
elapsed_time = toc()

# write results to csv file
writecsv("postmean_atzero$n.csv",f0_post)

print(mean_acc)

println("True value of density at zero equals: ",pdf(true_distribution,0))
