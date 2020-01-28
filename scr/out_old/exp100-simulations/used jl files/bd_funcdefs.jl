postmean0_simulation = function(x,IT, mh_step, α, p0, method, config_init,extra_pars)    
    n = length(x)
    ψ0 = calc_ψ0(x,α)
  
    # initialisation of the configuration (simply take one cluster)
    configs = Array{Config}(IT)
    configs[1] = config_init
  
    sum_acc = 0 
    postmean_atzero = zeros(IT)
  
    # MH algorithm by Neal
    for it in 2:IT
      configs[it] = update_config(configs[it-1],n,x, ψ0,method)
      configs[it], sum_acc = update_θ(configs[it],x,mh_step,sum_acc,method,extra_pars)
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

calc_ψ0 = function(x,α)
  ψ0 = zeros(n)
  for i=1:n
    ψ0[i] = quadgk(Ψ(x[i]), x[i], Inf)[1]
  end
  α * ψ0
end

update_config = function(config,n,x,ψ0,method,extra_pars)
    configN = deepcopy(config)
    for i in 1:n # loop over all labels
        # find table of the i-th customer
        ind = find(isequal(configN.labels[i]),configN.tableIDs)  
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
        J = configN.n_tables  # total number of tables
        w = zeros(J+1)
        for j in 1:J
          w[j] = ψ(x[i],configN.θ[j]) * configN.counts[j]
        end
        w[J+1] = ψ0[i]

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


update_θ = function (config,x,mh_step,sum_acc,method,extra_pars)
for k in 1:config.n_tables
    ind =find(isequal(config.tableIDs[k]),config.labels)  # find observation indices in table keys
    config.θ[k], acc =sample_θ(config.θ[k], x[ind],mh_step,method,extra_pars)
    sum_acc += acc
end
    config, sum_acc
end

logtarget = function(θ,x)
  log(base_density(θ)) - length(x) * log(θ) 
 end 

# target = function(θ,x)
#     pdf(base_measure,θ) * θ^(-length(x)) *(θ>maximum(x))
# end 


sample_θ1 = function(x,method,extra_pars)  # get one (independent) realisation for θ, assuming x is a float
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
end
if method=="B"
        out = x - log(rand())
end
if method in ["C","D"]
    out = rand(Pareto(1 + extra_pars[1],max(extra_pars[2],x)))
end

out    
end


sample_θ = function(θold,x,mh_step,method,extra_pars) # draw from the posterior, given observations in x
    if method in ["C","D"]
      acc = 1
      out = rand(Pareto(length(x)+extra_pars[1],maximum([x;extra_pars[2]])))
    else
      θ = θold + mh_step * randn()
      xmax = maximum(x)
      acc = 0 
      if θ < xmax # reject right-away θs out of the support of the target density
         out = θold
      else
        a = logtarget(θ,x) - logtarget(θold,x) - logccdf(Normal(θ,mh_step),xmax) + logccdf(Normal(θold,mh_step),xmax)
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

rand_trunc_gamma = function(α,β,thresh) 
    # sample once from Ga(shape=α, rate=β) -distribution, truncated to be less than thresh 
    y = cdf(Gamma(α,1/β),thresh) * rand()
    quantile(Gamma(α,1/β),y)
end