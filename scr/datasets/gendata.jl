x = abs.(randn(1000))
writedlm("halfnormal1000.csv",x,',')

x = rand(Exponential(1.0),1000)
writedlm("exp1000.csv",x,',')
