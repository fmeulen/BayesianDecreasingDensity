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
method = "A"

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
    dout, iterates, mean_acc, ep = mcmc(x, method, α, IT, std_mhstep, grid)
elapsedtime = time() - start

#------------------ postprocessing ------------------
# write results to csv file
CSV.write("./out/postsummary.csv",dout)
CSV.write("./out/iterates.csv",iterates)

# save info to file
facc = open("./out/info.txt","w")
write(facc, "Average acceptance probability equals: ",string(round(mean_acc;digits=3)),"\n")
write(facc, "Number of iterations: ",string(IT),"\n")
write(facc, "std_mhstep = ",string(std_mhstep),"\n\n")
write(facc, "elapsed time ",string(elapsedtime), "\n\n")
write(facc, "---- Prior specification ----","\n")
write(facc, "method: ", string(method),"\n")
write(facc, "concentration parameter: ", string(α),"\n")
write(facc, "extra_pars= ",string(ep),"\n")
close(facc)
