"""
    Set the prior used in Bayesian estimation of a decreasing density

    function call is
        base_density, base_measure, extra_pars = setprior(method)
"""
function setprior(method)
    if method=="A"
        base_measure = 0 # apparently not needed
        base_density = (θ) -> exp(-1/θ-θ) * (θ>0) / 0.27973176363304486
        extra_pars = 0 # no extra pars
    elseif method=="B"
        base_measure = Gamma(2.0,1.0)
        base_density = (θ) -> pdf(base_measure,θ)
        extra_pars = 0 # no extra pars
    elseif method=="C"
        αα = 1.0
        τ = 0.005
        extra_pars = [αα, τ] # no extra pars
        base_measure = Pareto(αα,τ) # second argument is the support parameter
        base_density = (θ) ->  pdf(base_measure,θ)
    elseif method=="D"
    # add mixture of Pareto
        αα = 1.0
        λ = 2.0 # prior on τ (which is the Pareto threshold) is assumed Ga(λ,β)
        β = 1.0
        τ = 0.1
        extra_pars = [αα,τ, αα, λ, β] # no extra pars
        base_measure = Pareto(αα,τ) # second argument is the support parameter
        base_density = (θ) ->  pdf(base_measure,θ)
    end
    base_density, base_measure, extra_pars
end


# Keep track of the configuration at each iteration using Config
mutable struct Config
    labels::Vector{Int64}               # (z_1,...,z_n)
    tableIDs::Vector{Int64}             # table indices
    counts::Vector{Int64}               # table counts
    n_tables::Int64                     # nr of tables
    θ::Vector{Float64}                  # parameters
end


function postmean0_simulation(x,IT, std_mhstep, α, p0, method, config_init,extra_pars)
    n = length(x)
    ψ0 = calc_ψ0(x,α)

    # initialisation of the configuration (simply take one cluster)
    configs = Array{Config}(IT)
    configs[1] = config_init

    sum_acc = 0
    postmean_atzero = zeros(IT)

    # MH algorithm by Neal
    for it in 2:IT
      configs[it] = updateconfig(configs[it-1],n,x, ψ0,method)
      configs[it], sum_acc = updateθ(configs[it],x,std_mhstep,sum_acc,method,extra_pars)
      postmean_atzero[it] = (p0 + dot(configs[it].counts,ψ(0).(configs[it].θ)))/(α+n)
    end

    # compute average acc prob for mh updates for theta
    nθupdates = 0
    for it in 2:IT
      nθupdates += configs[it].n_tables
    end

    postmean_atzero, sum_acc/nθupdates, configs[IT]
end

# define uniform density
ψ(x,θ) = pdf(Uniform(0,θ),x)
ψ(x) = (θ) -> ψ(x,θ) # define for fixed x, as a function of theta

# compute constants by numerical integration
Ψ(x) = (θ) -> ψ(x,θ) * base_density(θ)  # define as a function of θ

function calc_ψ0(x,α)
  # ψ0 = zeros(n)
  # for i=1:n
  #   ψ0[i] = quadgk(Ψ(x[i]), x[i], Inf)[1]
  # end
  # α * ψ0
  α * [quadgk(Ψ(x[i]), x[i], Inf)[1] for i in eachindex(x)]
end

function updateconfig(config,n,x,ψ0,method,extra_pars)
    configN = deepcopy(config)
    for i in 1:n # loop over all labels
        # find table of the i-th customer
        #ind = find(isequal(configN.labels[i]),configN.tableIDs) # old
        ind = findall(x->x==configN.labels[i], configN.tableIDs)[1] # new dec 2019
        # remove i-th customer from the table counts
        configN.counts[ind] += -1
        a = configN.counts[ind]

      # check for empty tables as a result of replacing the i-th label
      # in that case: remove the table from tableIDs, counts and θ, and reduce n_tables by one
        if a[1]==0  #true if removing customer i results in an empty table
            deleteat!(configN.counts, ind) # remove counts of table
            deleteat!(configN.tableIDs,ind) # remove table
            deleteat!(configN.θ,ind) # remove element ind from vector of θs
            configN.n_tables += -1
        end

        # construct weights to decide new table for customer i

        # original implementation
        # J = configN.n_tables  # total number of tables
        # w = zeros(J+1)
        # for j in 1:J
        #   w[j] = ψ(x[i],configN.θ[j]) * configN.counts[j]
        # end
        # w[J+1] = ψ0[i]
        w =[ψ(x[i],configN.θ[j]) * configN.counts[j] for j in eachindex(configN.n_tables)]
        append!(w,ψ0[i])

        # choose new table for customer i
        id_new = maximum(configN.tableIDs) + 1 # the label assigned, in case a new table is to be chosen
        ind = rand(Categorical(w/sum(w))) # index for chosen table
        newtable = [configN.tableIDs;id_new][ind]

        # add customer i back in, on its new table
        configN.labels[i] = newtable

        # bookkeep new table seating
        if newtable==id_new # open new table
            push!(configN.tableIDs,id_new) #add at the end
            push!(configN.counts,1)
            configN.n_tables += 1
            push!(configN.θ,sample_θ1(x[i],method,extra_pars))  # sample new θ and add it
        else
            configN.counts[ind] +=1
        end
    end
    configN
end


function updateθ!(config,sum_acc,x,std_mhstep,method,extra_pars)
    for k in 1:config.n_tables
        #ind =find(isequal(config.tableIDs[k]),config.labels)  # OLD find observation indices in table keys
        ind = findall(x->x==config.tableIDs[k], config.labels)[1] # new dec 2019
        config.θ[k], acc = sampleθ(config.θ[k], x[ind],std_mhstep,method,extra_pars)
        sum_acc += acc
    end
    #config, sum_acc
end

logtarget(θ,x) = log(base_density(θ)) - length(x) * log(θ)



# target = function(θ,x)
#     pdf(base_measure,θ) * θ^(-length(x)) *(θ>maximum(x))
# end

"""
    Get one (independent) realisation for θ, assuming x is a float

    i suspect this is just a sample from the prior on θ
"""
function sample_θ1(x,method,extra_pars)
    if method=="A"
        accep = false
        while !accep
            y = rand(Uniform(0,1/x))
            upper_bound = exp(-1/y-y)/(0.18*y)
            if (rand()) < upper_bound
                accep = true
                out = 1/y
            end
        end
    elseif method=="B"
        out = x - log(rand())
    elseif method in ["C","D"]
        out = rand(Pareto(1 + extra_pars[1],max(extra_pars[2],x)))
    end
    out
end

"""
    Draw from the posterior, given observations in x
"""
function sampleθ(θold,x,std_mhstep,method,extra_pars)
    if method in ["C","D"]
        acc = 1
        out = rand(Pareto(length(x)+extra_pars[1],maximum([x;extra_pars[2]])))
    else
        θ = θold + std_mhstep * randn()
        xmax = maximum(x)
        acc = 0
        if θ < xmax # reject right-away θs out of the support of the target density
           out = θold
        else
           a = logtarget(θ,x) - logtarget(θold,x) - logccdf(Normal(θ,std_mhstep),xmax) + logccdf(Normal(θold,std_mhstep),xmax)
           if log(rand()) < a
               out = θ
               acc = 1
            else
               out = θold
           end
       end
    end
    out, acc
end

"""
    sample once from Ga(shape=α, rate=β)-distribution, truncated to be less than thresh
"""
function randtruncgamma(α,β,thresh)
    y = cdf(Gamma(α,1/β),thresh) * rand()
    quantile(Gamma(α,1/β),y)
end


function priorconstants(grid, x, α)
    lg = length(grid)
    ρ = zeros(lg)
    ρ[lg] = α * quadgk(Ψ(grid[lg]), grid[lg], Inf)[1]
    for k in lg:-1:2
        ρ[k-1] = ρ[k] + α * quadgk(Ψ(grid[k-1]), grid[k-1], grid[k])[1]
    end
    ψ0 = calc_ψ0(x,α)
    ρ, ψ0
end


"""
    Do mcmc algorithm of Neal

    postτ, postmean, mean_acc = mcmc(IT, grid, method)

    x: data
    α: concentration parameter for Dirichlet process
    IT: number of iterations is then IT-1

    first row of postmean contains the grid on which we compute the posterior mean
"""
function mcmc(x, method, α, IT, std_mhstep, grid)
    base_density, base_measure, extra_pars = setprior(method)
    # compute prior constants
    ρ, ψ0 = priorconstants(grid, x, α)
    n = length(x)

    # initialisation of the configuration (simply take one cluster)
    config_init = Config(ones(n),[1], [n],1,[maximum(x)+1])
    configs = [deepcopy(config_init) for j in 1:IT]  #Array{Config}(IT)

    postmean = zeros(IT,length(grid))
    postτ = ones(IT) # only relevant with mixture of Pareto base measure
    postτ[1] = 10.0

    md = method=="D"

    if md
        αα, λ, β = extra_pars[3], extra_pars[4], extra_pars[5]
        extrapars[2] = postτ[1]
    end

    # MH algorithm by Neal
    sum_acc = 0
    for it in 2:IT
        configs[it] = updateconfig(configs[it-1],n,x, ψ0,method,extra_pars)
        updateθ!(configs[it],sum_acc,x,std_mhstep,method,extra_pars)
        if md
            postτ[it] = randtruncgamma(λ + configs[it].n_tables * αα, β, minimum(configs[it].θ))
            extra_pars[2] = postτ[it]
        end

        if it%100==0   println(it) end

        postmean[it,:] = [(p[k] + dot(configs[it].counts,ψ(grid[k]).(configs[it].θ)))/(α+n) for k in eachindex(grid)]
    end

    # compute average acc prob for mh updates for theta
    nθupdates = sum([configs[it].n_tables for it in 2:IT])
    mean_acc = sum_acc/nθupdates

    println("Average acceptance rate on updating θ: ", round(mean_acc;digits=3))
    postmean[1,:] = grid
    postτ, postmean, mean_acc, extra_pars
end
