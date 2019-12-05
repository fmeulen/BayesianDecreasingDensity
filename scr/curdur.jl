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

#------------------ specify the prior ------------------
α = 1.0                       # concentration par
method = "D"

# specify simulation settings
std_mhstep = 0.2 # sd for mh step updating θ
IT = 25_000

#------------------  read data ------------------
dataset = ["curdur", "exp100"][2]
if dataset=="curdur"
    y = readdlm("curdur.csv",',', header=true)[1][:,2]
    x = sort(y)[1:618]
    grid = range(0,36;length=300)
elseif dataset=="exp100"
    x = readdlm("exp100.csv",',',header=true)[1][:,2]#[2:end,2]
    grid = range(0,5;length=100)
end

#------------------ mcmc by Neal ------------------
start = time()
    out = mcmc(x, method, α, IT, std_mhstep, grid)
elapsedtime = time() - start

#------------------ postprocessing ------------------
# write results to csv file
CSV.write("./out/postsummary.csv",out.dout)
CSV.write("./out/iterates.csv",out.iterates)

# save info to file
facc = open("./out/info.txt","w")
write(facc, "Average acceptance probability equals: ",string(round(out.mean_acc;digits=3)),"\n")
write(facc, "Number of iterations: ",string(IT),"\n")
write(facc, "std_mhstep = ",string(std_mhstep),"\n\n")
write(facc, "elapsed time ",string(elapsedtime), "\n\n")
write(facc, "---- Prior specification ----","\n")
write(facc, "method: ", string(method),"\n")
write(facc, "concentration parameter: ", string(α),"\n")
write(facc, "extra_pars= ",string(out.ep),"\n")
close(facc)

ratecomparison = false

start = time()
if ratecomparison
    IT = 2500 # nr of iterations
    NSIM = 50 # 100

    # specify data generation
    true_distribution = Exponential(1)
    n= 100

    # do the simulations
    f0post = zeros(NSIM)
    mean_acc = zeros(NSIM)
#    ρ0 =  quadgk(Ψ(0, base_density), 0, Inf)[1]
    grid0 =[0.0]
    # do one preliminary run to obtain initial configuration
    #x = rand(true_distribution,n)
    x = abs.(randn(n))  # half normal distr
    outprelim = mcmc(x, method, α, IT, std_mhstep, grid0)
    println(outprelim.mean_acc)

    cpre = outprelim.config_end

    println("Start Monte-Carlo study")
    for k in 1:NSIM
      print(k)
      #x = rand(true_distribution,n)
      x = abs.(randn(n))
      sim = mcmc(x, method, α, IT, std_mhstep, grid0, config_init=cpre)
      f0post[it] = sim.dout.ave[1]
    end


    # write results to csv file
    CSV.write("f0post$n.csv",f0_post)
end
elapsedtime = time() - start
