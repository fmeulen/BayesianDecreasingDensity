# Keep track of the configuration at each iteration using Config
mutable struct Config
    labels::Vector{Int64}               # (z_1,...,z_n)
    tableIDs::Vector{Int64}             # table indices
    counts::Vector{Int64}               # table counts
    n_tables::Int64                     # nr of tables
    θ::Vector{Float64}                  # parameters
end

mutable struct ExtraPar
    αα::Float64
    λ::Float64
    β::Float64
    τ::Float64
end


#-------------------------- set prior --------------------------
"""
    Set the prior used in Bayesian estimation of a decreasing density

    function call is
        base_density, base_measure, ep = setprior(method)
"""
function setprior(method)
    if method=="A"
        base_density = (θ) -> exp(-1/θ-θ) * (θ>0) / 0.27973176363304486
        ep = ExtraPar(0.0, 0.0, 0.0, 0.0) # no "ep = extra pars"
    elseif method=="B"
        base_measure = Gamma(2.0,1.0)
        base_density = (θ) -> pdf(base_measure,θ)
        ep = ExtraPar(0.0, 0.0, 0.0, 0.0) # no extra pars
    elseif method=="C"
        αα = 1.0
        τ = 0.005
        ep = ExtraPar(αα, 0.0, 0.0,τ) # no extra pars
        base_measure = Pareto(αα,τ) # second argument is the support parameter
        base_density = (θ) ->  pdf(base_measure,θ)
    elseif method=="D"
    # add mixture of Pareto # prior on τ (which is the Pareto threshold) is assumed Ga(λ,β)
        αα = 1.0
        λ = 2.0
        β = 1.0
        τ = 0.1
        ep = ExtraPar(αα, λ, β, τ) # no extra pars
        base_measure = Pareto(αα,τ) # second argument is the support parameter
        base_density = (θ) ->  pdf(base_measure,θ)
    end
    base_density,  ep
end



#-------------------------- define uniform density --------------------------
ψ(x,θ) = pdf(Uniform(0,θ),x)
ψ(x) = (θ) -> ψ(x,θ)
Ψ(x, base_density) = (θ) -> ψ(x,θ) * base_density(θ)

#-------------------------- compute some numerical constants for the dataset x --------------------------
"""
    compψ0(x,α,base_density)
        x: data
        α: concentration par of Dirichlet process
        base_density: density of base meaure of Dirichlet process

    Integrate  ψ(x,θ) * base_density(θ) over x=x[i]..infinity for all i=1..n
    Return the result, multiplied by α

    These numbers need to be calculated once, and are needed for updating the labels.
"""
compψ0(x,α,base_density) = α * [quadgk(Ψ(x[i],base_density), x[i], Inf)[1] for i in eachindex(x)]

"""
    ρ: compute    int ψ_x d G_0(x)     for x in grid
    also call compψ0 and returns both ρ and ψ0
"""
function priorconstants(grid, x, α, base_density)
    lg = length(grid)
    ρ = zeros(lg)
    ρ[lg] = quadgk(Ψ(grid[lg], base_density), grid[lg], Inf)[1]
    for k in lg:-1:2
        ρ[k-1] = ρ[k] + quadgk(Ψ(grid[k-1], base_density), grid[k-1], grid[k])[1]
    end
    ψ0 = compψ0(x,α,base_density)
    ρ, ψ0
end

#-------------------------- update labels --------------------------
function updateconfig(config, x, ψ0, method, ep)
    configN = deepcopy(config)
    for i in eachindex(x) # loop over all labels
        # find tableID of the i-th customer
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

        w =[ψ(x[i],configN.θ[j]) * configN.counts[j] for j in 1:configN.n_tables]
        append!(w,ψ0[i]) # append for opening a new table

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
            push!(configN.θ,sampleθnewtable(x[i],method,ep))  # sample new θ and add it
        else
            configN.counts[ind] +=1
        end
    end
    configN
end

"""
    Get one (independent) realisation for θ, assuming x is a float
    sample θ for a newly opened table with observation x::Float64
    i suspect this is just a sample from the prior on θ
"""
function sampleθnewtable(x,method,ep)
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
        out = rand(Pareto(1 + ep.αα,max(ep.τ,x)))
    end
    out
end

#-------------------------- update parameter θ for each table --------------------------
"""
    log target density for updating the parameter at one particular table; x is a vector of
    observations that sit at this particular table
"""
logtarget(θ,x,base_density) = log(base_density(θ)) - length(x) * log(θ)


"""
    Draw from the posterior at a particular table,
    Here x is a vector with those observations that belong to this table
"""
function updateθtable!(θ,x,std_mhstep,method,ep,base_density)
    acc = 0
    if method in ["C","D"]
        acc = 1
        θ = rand(Pareto(length(x)+ep.αα,maximum([x;ep.τ])))
    else
        θᵒ = θ + std_mhstep * randn()
        xmax = maximum(x)
        if θᵒ >= xmax # else reject right-away θᵒ out of the support of the target density
           a = logtarget(θᵒ,x,base_density) - logtarget(θ,x,base_density) - logccdf(Normal(θᵒ,std_mhstep),xmax) + logccdf(Normal(θ,std_mhstep),xmax)
           if log(rand()) < a
               θ = θᵒ
               acc = 1
           end
       end
   end
   acc
end

function updateθ!(config,sum_acc,x,std_mhstep,method,ep,base_density)
    for k in 1:config.n_tables
        # for the k-th table, find all indices of customers on that table
        ind = findall(x->x==config.tableIDs[k], config.labels) # new dec 2019
        acc = updateθtable!(config.θ[k], x[ind],std_mhstep,method,ep,base_density)
        sum_acc += acc
    end
    sum_acc
end

#-------------------------- sim from truncated Gamma distribution (only relevant for mixture prior (method D))--------------------------
"""
    sample once from Ga(shape=α, rate=β)-distribution, truncated to be less than thresh
"""
function randtruncgamma(α,β,thresh)
    y = cdf(Gamma(α,1/β),thresh) * rand()
    quantile(Gamma(α,1/β),y)
end


#-------------------------- mcmc by Neal --------------------------
"""
    Do mcmc algorithm of Neal

    dout, iterates, mean_acc, ep = mcmc(x, method, α, IT, std_mhstep, grid; p=0.05, BI=div(IT,2))

    Arguments:
        x: data
        method: type of prior
        α: concentration parameter for Dirichlet process
        std_mhstep: standard deviation of Normal distribution for Random-Walk Metropolis-Hastings step
        IT: number of iterations
        grid: grid on which to compute the posterior mean
        p: get pointwise (1-p)100% credible intervals
        BI: number of burnin samples

    Returns:
        dout: a dataframe with as columns grid, average of postmean, quantile(1-p/2) of postmean, quantile(p/2) of postmean
        iterates: dataframe containing 3 columns:
            (1) posttau: all iterates for τ (only relevant if method==D)
            (2) postmean0: all iterates at grid[1]
            (3) iteratenr
        mean_acc: average acceptance rate on MH-step for θ
        ep: extra parameters (only relevant if method==D)
        final configuration

"""
function mcmc(x, method, α, IT, std_mhstep, grid;
                p=0.05, BI=div(IT,2),
                config_init=Config(ones(n),[1], [n],1,[maximum(x)+1]) )

    base_density, ep = setprior(method)
    # compute prior constants
    ρ, ψ0 = priorconstants(grid, x, α, base_density)
    n = length(x)

    # initialisation of the configuration (simply take one cluster)
    #config_init = Config(ones(n),[1], [n],1,[maximum(x)+1])
    configs = [deepcopy(config_init) for j in 1:IT]  #Array{Config}(IT)

    postmean = zeros(IT,length(grid))
    postτ = fill(10.0, IT) # only relevant with mixture of Pareto base measure
    postmean[1,:] = [(α*ρ[k] + dot(configs[1].counts,ψ(grid[k]).(configs[1].θ)))/(α + n) for k in eachindex(grid)]

    md = method=="D"
    if md    ep.τ = postτ[1]    end

    # MH algorithm by Neal
    sum_acc = 0
    for it in 2:IT
        configs[it] = updateconfig(configs[it-1], x, ψ0, method, ep)
        sum_acc = updateθ!(configs[it],sum_acc, x, std_mhstep, method, ep, base_density)
        if md
            postτ[it] = randtruncgamma(ep.λ + configs[it].n_tables * ep.αα, ep.β, minimum(configs[it].θ))
            ep.τ = postτ[it]
        end

        if it%250==0   println(it) end

        conf = configs[it]
        postmean[it,:] = [(α*ρ[k] + dot(conf.counts, ψ(grid[k]).(conf.θ)))/(α + n) for k in eachindex(grid)]
    end

    # compute average acc prob for mh updates for theta
    nθupdates = sum([configs[it].n_tables for it in 2:IT])
    mean_acc = sum_acc/nθupdates

    postmean_BI = postmean[BI:IT,:]
    upper = vec(mapslices(v-> quantile(v, 1 - p/2), postmean_BI; dims= 1))
    ave = vec(mean(postmean_BI,dims=1))
    lower = vec(mapslices(v-> quantile(v,  p/2), postmean_BI; dims= 1))
    dout = DataFrame(lower=lower,upper=upper, ave=ave,x=grid)

    println("Average acceptance rate on updating θ: ", round(mean_acc;digits=3))

    iterates = DataFrame(posttau=postτ, postmean0=postmean[:,1],iteratenr=1:IT)
    (dout=dout, iterates=iterates, mean_acc=mean_acc, ep=ep, config_end=configs[end])
end


#------------------------------- script for rate comparison at zero --------------------------
# probably obsolete
function postmean0_simulation(x,IT, std_mhstep, α, p0, method, config_init,ep)
    n = length(x)
    ψ0 = compψ0(x,α)

    # initialisation of the configuration (simply take one cluster)
    configs = Array{Config}(IT)
    configs[1] = Config(ones(n),[1], [n],1,[maximum(x)+1])

    sum_acc = 0
    postmean_atzero = zeros(IT)

    # MH algorithm by Neal
    for it in 2:IT
      configs[it] = updateconfig(configs[it-1],x, ψ0,method)
      configs[it], sum_acc = updateθ(configs[it],x,std_mhstep,sum_acc,method,ep)
      postmean_atzero[it] = (α*p0 + dot(configs[it].counts,ψ(0).(configs[it].θ)))/(α+n)
    end

    # compute average acc prob for mh updates for theta
    nθupdates = 0
    for it in 2:IT
      nθupdates += configs[it].n_tables
    end

    postmean_atzero, sum_acc/nθupdates, configs[IT]
end



    ### debugging  example
    if false

        configN = Config([1,2,1,3,1,4],[1,3,4,2],[3,1,1,1],4,[0.1, 0.3, 0.5, 0.7])
        i=2
        ind = findall(x->x==configN.labels[i], configN.tableIDs)[1]


        # labels::Vector{Int64}               # (z_1,...,z_n)
        # tableIDs::Vector{Int64}             # table indices
        # counts::Vector{Int64}               # table counts
        # n_tables::Int64                     # nr of tables
        # θ::Vector{F
    end
