using Distributions
using CSV
using QuadGK

workdir = @__DIR__
println(workdir)
cd(workdir)


#using Plots
include("bd_funcdefs.jl")

# specify the prior
α = 1 # concentration par
method = "D"
base_density, base_measure, extra_pars = setprior(method)

# specify simulation settings
mh_step = 0.4 # sd for mh step updating θ
IT = 50000

# read data
y = readdlm("curdur.csv",',', header=true)[1][:,2]
x = sort(y)[1:618]

# sample of size 100 from Exp(1) distribution
#x =readcsv("exp100.csv")[2:end,2]

n = length(x)

#grid = linspace(0,maximum(x)+1,50)  # compute estimates on this grid
grid = range(0,36;length=300)#linspace(0,36,300)
#grid = [linspace(0,0.01,11);linspace(0.05,5,100)]  # for exp100

# initialisation of the configuration (simply take one cluster)
config1 = Config(ones(n),[1], [n],1,[maximum(x)+1])
configs = [deepcopy(config1) for j in 1:IT]  #Array{Config}(IT)


# compute prior constants
lg = length(grid)
p = zeros(lg)
p[lg] = α * quadgk(Ψ(grid[lg]), grid[lg], Inf)[1]
for k in lg:-1:2
    p[k-1] = p[k] + α * quadgk(Ψ(grid[k-1]), grid[k-1], grid[k])[1]
end
ψ0 = calc_ψ0(x,α)

sum_acc = 0
postmean = zeros(IT,lg)

postτ = ones(IT) # only relevant with mixture of Pareto base measure
postτ[1] =10

# MH algorithm by Neal
for it in 2:IT
  if method=="D"
      extra_pars[2] = postτ[it-1]
  end
  configs[it] = update_config(configs[it-1],n,x, ψ0,method,extra_pars)
  configs[it], sum_acc = update_θ(configs[it],x,mh_step,sum_acc,method,extra_pars)
  if method=="D"
      postτ[it] = rand_trunc_gamma(λ + configs[it].n_tables * αα,β,minimum(configs[it].θ))
  end

  if it%100==0   println(it) end

  for k in 1:lg
    gk = grid[k]
    postmean[it,k] = (p[k] + dot(configs[it].counts,ψ(gk).(configs[it].θ)))/(α+n)
  end
end

elapsed_time = toc()

# compute average acc prob for mh updates for theta
nθupdates = sum([configs[it].n_tables for it in 2:IT])

mean_acc = sum_acc/nθupdates
println(mean_acc)

# write results to csv file
postmean[1,:] = grid
writecsv("./out/postmean.csv",postmean)
writecsv("./out/post_tau.csv",postτ)

# save info to file
facc = open("./out/info.txt","w")
write(facc, "Average acceptance probability equals: ",string(round(mean_acc,3)),"\n")
write(facc, "Number of iterations: ",string(IT),"\n")
write(facc, "mh_step = ",string(mh_step),"\n\n")
write(facc, "elapsed time ",string(elapsed_time), "\n\n")
write(facc, "---- Prior specification ----","\n")
write(facc, "method: ", string(method),"\n")
write(facc, "concentration parameter: ", string(α),"\n")
write(facc, "extra_pars= ",string(extra_pars),"\n")
close(facc)
