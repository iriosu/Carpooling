using JuMP
using PyPlot
using Gurobi, KNITRO
include("utilities.jl")

# ==============================
# Formulates the problem and solves it
# ==============================
function GenerateInputs(vt, ct, fm, gm)
    # fm = marginal distributions of valuations
    # gm = marginal distributions of costs
    if length(vt) != length(ct)
        println("***ERROR: current implementation does not allow unequal number of value and cost types")
        exit()
    end

    nvt = length(vt)
    nct = length(ct)
    nvct = nvt*nct
    ntypes = length(vt)
    N = length(fm)
    Theta = combwithrep(N, nvct) #2* because we have two parameters for each agent
    Theta_i = combtwovec(vt,ct)
    if length(Theta_i) != nvct
        println("***ERROR: Theta_i and number of types do not match")
        exit()
    end

    # each vector of types is a a vector (v_1, v_2, ..., v_n, c_1, c_2, ..., c_n)
    # joint distribution based on marginals
    Idx = Dict(Theta[i]=>i for i=1:length(Theta))





    f = Dict()
    for perm in Theta
        # perm[i] is in 1:nvtc; need to convert it to tuple (v,c)
        # to do that, we have Theta_i: Theta_i[perm[i]] = (v,c)
        f[perm] = prod([fm[i][Theta_i[perm[i]][1]] for i=1:N])*prod([gm[i][Theta_i[perm[i]][2]] for i=1:N])
    end



    # ------------------------------------------------------------------------------
    # CENTRALIZED WITH IC AND IR
    # ------------------------------------------------------------------------------
    nvars = length(Theta)*N # one variable for each supplier and type vector
    sts = length(Theta) # size of type space

    # this dict tell us the probability of all other subjects being of their type $f_{-i}(\theta_{-i})$
    f_woi = Dict(i=>Dict(j=>0.0 for j in Theta) for i in 1:N)

    for i in keys(f_woi)
        supp_woi = [j for j in 1:N if j!=i]
        for perm in Theta
            f_woi[i][perm] = prod([fm[j][Theta_i[perm[j]][1]] for j in supp_woi])*prod([gm[j][Theta_i[perm[j]][2]] for j in supp_woi])
        end
    end


    # perm[i] is in 1:nvtc; need to convert it to tuple (v,c)
    # to do that, we have Theta_i: Theta_i[perm[i]] = (v,c)
    # individual rationality constraints
    IR_y = zeros((nvct*N, nvars))
    IR_z = zeros((nvct*N, nvars))
    IR_t = zeros((nvct*N, nvars))
    row = 1
    for i in 1:N
        for s=1:nvct
            v_i, c_i = Theta_i[s][1], Theta_i[s][2]
            idxs_tt = [r for r in 1:length(Theta) if Theta_i[Theta[r][i]][1] == v_i && Theta_i[Theta[r][i]][2] == c_i]
            println(v_i, " ", c_i, " ", idxs_tt)
            for l in idxs_tt
                IR_y[row, N*(l-1)+i] = v_i*f_woi[i][Theta[l]]
                IR_z[row, N*(l-1)+i] = -c_i*f_woi[i][Theta[l]]
                IR_t[row, N*(l-1)+i] = f_woi[i][Theta[l]]
            end
            row=row+1
        end
    end

    IC_y = zeros((nvct*(nvct-1)*N,nvars))
    IC_z = zeros((nvct*(nvct-1)*N,nvars))
    IC_t = zeros((nvct*(nvct-1)*N,nvars))
    row = 1
    for i in 1:N
        for s=1:nvct
            v_i, c_i = Theta_i[s][1], Theta_i[s][2]
            ot = [l for l=1:nvct if l!=s]
            for l in ot
                v_ip, c_ip = Theta_i[l][1], Theta_i[l][2]
                idxs_tt = [r for r in 1:length(Theta) if Theta_i[Theta[r][i]][1] == v_i && Theta_i[Theta[r][i]][2] == c_i]
                idxs_tp = [r for r in 1:length(Theta) if Theta_i[Theta[r][i]][1] == v_ip && Theta_i[Theta[r][i]][2] == c_ip]
                for r in idxs_tt
                    IC_y[row, N*(r-1)+i] = v_i*f_woi[i][Theta[r]]
                    IC_z[row, N*(r-1)+i] = -c_i*f_woi[i][Theta[r]]
                    IC_t[row, N*(r-1)+i] = f_woi[i][Theta[r]]
                end
                for r in idxs_tp
                    IC_y[row, N*(r-1)+i] = -v_i*f_woi[i][Theta[r]]
                    IC_z[row, N*(r-1)+i] = c_i*f_woi[i][Theta[r]]
                    IC_t[row, N*(r-1)+i] = -f_woi[i][Theta[r]]
                end
                row+=1
            end
        end
    end



    # if version == "centralized"
    m = Model(solver=GurobiSolver(MIPGap = 1e-12))
    # elseif version == "decentralized"
    #     m = Model(solver=KnitroSolver(mip_method = KTR_MIP_METHOD_BB, honorbnds=0,
    #                                   ms_enable = 1, ms_maxsolves = 5000,
    #                                   algorithm = KTR_ALG_ACT_CG,
    #                                   outmode = KTR_OUTMODE_SCREEN,
    #                                   KTR_PARAM_OUTLEV = KTR_OUTLEV_ALL,
    #                                   KTR_PARAM_MIP_OUTINTERVAL = 1,
    #                                   KTR_PARAM_MIP_MAXNODES = 10000,
    #                                   KTR_PARAM_HESSIAN_NO_F = KTR_HESSIAN_NO_F_ALLOW))
    # else
    #     println("***ERROR: unknown version")
    #     exit()
    # end


    @variable(m, x[s=1:sts,i=1:N,j=1:N], Bin)
    @variable(m, y[s=1:nvars])
    @variable(m, z[s=1:nvars], Bin)
    @variable(m, t[s=1:nvars])

    @constraint(m, IR, IR_y*y + IR_t*t .>= 0)
    @constraint(m, IC, IC_y*y + IC_t*t .>= 0)

    @constraint(m, feas1[s=1:sts, i=1:N], sum(x[s,k,i] for k=1:N) - y[(s-1)*N+i] == 0 )
    @constraint(m, feas2[s=1:sts, i=1:N], sum(x[s,k,i] for k=1:N) <= 1 )
    @constraint(m, feas3[s=1:sts, i=1:N], sum(x[s,i,j] for j=1:N) <= 1+z[(s-1)*N+i] )




    @objective(m, Max, sum(f[Theta[s]]*( sum(x[s,i,j]*Theta_i[Theta[s][j]][1] for i=1:N for j=1:N) - sum(z[(s-1)*N+i]*Theta_i[Theta[s][i]][2] for i=1:N) ) for s=1:sts ))
    print(m)

    status = solve(m)
    obj = -getobjectivevalue(m)
    x_vals = getvalue(x)
    t_vals = getvalue(t)
    # transfers = wq_t'*t_vals
    # # println("Objective value: ", getobjectivevalue(m))
    # #
    println("===============================")
    println("Objective value: ", getobjectivevalue(m))
    println("Allocations: ", getvalue(x))
    println("Transfers: ", getvalue(t))
    println("===============================")
    # return obj, transfers, x_vals, t_vals
end




# IMPORTANT: types must be sorted in increasing order
### INPUTS ###
vt = Dict(1=>8, 2=>10)
ct = Dict(1=>10, 2=>12)
fm = Dict(1=>[0.5,0.5],2=>[0.5,0.5])
gm = Dict(1=>[0.5,0.5],2=>[0.5,0.5])

GenerateInputs(vt, ct, fm, gm)
