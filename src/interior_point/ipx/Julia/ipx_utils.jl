using MAT

# pres, dres, compl = CheckInteriorSolution(model, solution)
function CheckInteriorSolution(model::Model, solution::InteriorSolution)
    m,n = size(model.A)
    lb = model.lb
    ub = model.ub
    x = solution.x
    slack = solution.slack
    xl = solution.xl
    xu = solution.xu
    y = solution.y
    zl = solution.zl
    zu = solution.zu
    @assert all(isfinite.(x))
    @assert all(isfinite.(slack))
    @assert all(isfinite.(y))
    @assert all(isfinite.(zl))
    @assert all(isfinite.(zu))
    @assert all(xl .>= 0.0)
    @assert all(xu .>= 0.0)
    @assert all(zl .>= 0.0)
    @assert all(zu .>= 0.0)
    compl = 0.0                 # sum of pairwise complementarity products
    for i = 1:m
        if model.constr_type[i] == '<'
            @assert slack[i] >= 0.0
            @assert y[i] <= 0.0
            compl -= y[i]*slack[i]
        end
        if model.constr_type[i] == '>'
            @assert slack[i] <= 0.0
            @assert y[i] >= 0.0
            compl -= y[i]*slack[i]
        end
        if model.constr_type[i] == '='
            @assert slack[i] == 0
        end
    end
    rb = model.rhs - model.A*x  - slack
    rc = model.obj - model.A'*y - zl + zu
    pres = norm(rb, Inf)
    dres = norm(rc, Inf)
    for j = 1:n
        if isfinite(lb[j])
            pres = max(pres, abs(lb[j]-x[j]+xl[j]))
            compl += xl[j]*zl[j]
        else
            @assert isinf(xl[j])
        end
        if isfinite(ub[j])
            pres = max(pres, abs(ub[j]-x[j]-xu[j]))
            compl += xu[j]*zu[j]
        else
            @assert isinf(xu[j])
        end
    end
    return (pres, dres, compl)
end

# pres, dres, pinf, dinf = CheckBasicSolution(model, solution)
function CheckBasicSolution(model::Model, solution::BasicSolution)
    m,n = size(model.A)
    lb = model.lb
    ub = model.ub
    x = solution.x
    slack = solution.slack
    y = solution.y
    z = solution.z
    cbasis = solution.cbasis
    vbasis = solution.vbasis
    @assert all(isfinite.(x))
    @assert all(isfinite.(slack))
    @assert all(isfinite.(y))
    @assert all(isfinite.(z))
    pinf = 0.0
    dinf = 0.0
    for j = 1:n
        @assert vbasis[j] >= -3 && vbasis[j] <= 0
        if vbasis[j] == 0
            @assert z[j] == 0.0
            pinf = max(pinf, lb[j]-x[j])
            pinf = max(pinf, x[j]-ub[j])
        end
        if vbasis[j] == -1
            @assert x[j] == lb[j]
            dinf = max(dinf, -z[j])
        end
        if vbasis[j] == -2
            @assert x[j] == ub[j]
            dinf = max(dinf, z[j])
        end
        if vbasis[j] == -3
            @assert lb[j] == -Inf
            @assert ub[j] == +Inf
            @assert x[j] == 0.0
            dinf = max(dinf, abs(z[j]))
        end
    end
    for i = 1:m
        @assert cbasis[i] == -1 || cbasis[i] == 0
        if cbasis[i] == -1
            @assert slack[i] == 0.0
            if model.constr_type[i] == '<'
                dinf = max(dinf, y[i])
            end 
            if model.constr_type[i] == '>'
                dinf = max(dinf, -y[i])
            end
       end
        if cbasis[i] == 0
            @assert y[i] == 0.0
            if model.constr_type[i] == '<'
                pinf = max(pinf, -slack[i])
            end 
            if model.constr_type[i] == '>'
                pinf = max(pinf, slack[i])
            end
            if model.constr_type[i] == '='
                pinf = max(pinf, abs(slack[i]))
            end
        end
    end
    rb = model.rhs - model.A*x  - slack
    rc = model.obj - model.A'*y - z
    pres = norm(rb, Inf)
    dres = norm(rc, Inf)
    return (pres, dres, pinf, dinf)
end

function Check(solver::LPSolver, model::Model)
    ipxinfo = GetInfo(solver)
    if ipxinfo.status_ipm == 1 || ipxinfo.status_ipm == 2
        solution = GetInteriorSolution(solver)
        pres,dres,compl = CheckInteriorSolution(model, solution)
        @printf("Interior solution: pres = %.2e, dres = %.2e", pres, dres)
        @printf(", cmpl = %.2e\n", compl)
    end
    if ipxinfo.status_crossover == 1 || ipxinfo.status_crossover == 2
        solution = GetBasicSolution(solver)
        pres,dres,pinf,dinf = CheckBasicSolution(model, solution)
        @printf("Basic solution   : pres = %.2e, dres = %.2e", pres, dres)
        @printf(", pinf = %.2e, dinf = %.2e\n", pinf, dinf)
    end
end

function WriteBasis(solver::LPSolver, model::Model, filename::String)
    solution = GetBasicSolution(solver)
    pres,dres,pinf,dinf = CheckBasicSolution(model, solution)
    method = dinf <= pinf ? 1 : 0 # dual or primal simplex
    f = matopen(filename, "w")
    write(f, "cbasis", solution.cbasis)
    write(f, "vbasis", solution.vbasis)
    write(f, "method", method)
    close(f)
end
