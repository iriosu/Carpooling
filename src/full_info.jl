Pkg.add("Distributions")
using JuMP
using Gurobi
import Distributions: Uniform

srand(1)

# Parameters of the problem
N = 3 # number of users who submitted a request

kappa = rand(0:2, N)

c = rand(Uniform(0,0.1), N, N)
w = rand(Uniform(0,0.25), N) # cost of alternative transportation method
println("=======Parameters=========")
println(kappa)
println(c)
println(w)
println("==========================")
w = [min(c[i,i], w[i]) for i=1:N] # driving solo is always an option
println(w)
println("==========================")
m = Model(solver=GurobiSolver(MIPGap = 1e-12))

# @variable(m, x[i=1:N,j=1:N], Bin)
# @variable(m, t[i=1:N,j=1:N]>=0)
# @variable(m, p[i=1:N,j=1:N]>=0)
#
# @constraint(m, cap[i=1:N], sum(x[i,j] for j=1:N) <= kappa[i]*x[i,i])
# @constraint(m, spu[i=1:N], sum(x[k,i] for i=1:N) <= 1)
# # @constraint(m, tra[i=1:N], sum(t[i,j] for j=1:N) <= x[i,i])
# @constraint(m, pri[i=1:N,j=1:N], p[i,j] <= x[i,j])
# @constraint(m, fea1[i=1:N], p[i,i] == 0)
# # @constraint(m, fea2[i=1:N], t[i,i] == 0)
# @constraint(m, IR[i=1:N], sum(t[i,j] for j=1:N) - sum(p[k,i] for k=1:N) - sum(c[i,j]*x[i,j] for j=1:N) >= w[i] ) #*sum(x[k,i] for k=1:N)
#
#
# @objective(m, Max, sum(p[i,j] for i=1:N for j=1:N) - sum(t[i,j] for i=1:N for j=1:N))


@variable(m, x[i=1:N,j=1:N], Bin)
@variable(m, t[i=1:N]>=0)
@variable(m, p[i=1:N]>=0)

@constraint(m, cap[i=1:N], sum(x[i,j] for j=1:N) <= kappa[i]*x[i,i])
@constraint(m, spu[i=1:N], sum(x[k,i] for k=1:N) <= 1)

# Incentives
# @constraint(m, pri[i=1:N], p[i] <= sum(x[k,i] for k=1:N if k!=i))
# @constraint(m, IR[i=1:N], t[i]-p[i] - sum(c[i,j]*x[i,j] for j=1:N) >= -w[i] ) #*sum(x[k,i] for k=1:N)

# @objective(m, Max, sum(p[i] for i=1:N) - sum(t[i] for i=1:N))
@objective(m, Max, sum(x[i,j]*(w[j] - c[i,j]) for i=1:N for j=1:N))

print(m)

status = solve(m)
obj = getobjectivevalue(m)
x_vals = getvalue(x)
t_vals = getvalue(t)
p_vals = getvalue(p)

println(x_vals)
println(t_vals)
println(p_vals)


# @constraint(m, tra[i=1:N], sum(t[i,j] for j=1:N) <= x[i,i])
# @constraint(m, fea2[i=1:N], t[i,i] == 0)

# [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 1.0 1.0]
# [0.0, 0.0, 0.0251662]
# [0.0, 0.0999905, 0.0]
