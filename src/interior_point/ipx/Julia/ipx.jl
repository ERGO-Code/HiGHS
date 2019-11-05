module ipx

const ipxint = Int64
const spmatrix = SparseMatrixCSC{Cdouble,ipxint}
const IPX_STATUS_not_run                  = 0
const IPX_STATUS_solved                   = 1000
const IPX_STATUS_stopped                  = 1005
const IPX_STATUS_invalid_input            = 1002
const IPX_STATUS_out_of_memory            = 1003
const IPX_STATUS_internal_error           = 1004
const IPX_STATUS_optimal                  = 1
const IPX_STATUS_imprecise                = 2
const IPX_STATUS_primal_infeas            = 3
const IPX_STATUS_dual_infeas              = 4
const IPX_STATUS_time_limit               = 5
const IPX_STATUS_iter_limit               = 6
const IPX_STATUS_no_progress              = 7
const IPX_STATUS_failed                   = 8
const IPX_STATUS_debug                    = 9
const IPX_ERROR_argument_null             = 102
const IPX_ERROR_invalid_dimension         = 103
const IPX_ERROR_invalid_matrix            = 104
const IPX_ERROR_invalid_vector            = 105
const IPX_ERROR_invalid_bound             = 106
const IPX_ERROR_invalid_basis             = 107
const IPX_ERROR_cr_iter_limit             = 201
const IPX_ERROR_cr_matrix_not_posdef      = 202
const IPX_ERROR_cr_precond_not_posdef     = 203
const IPX_ERROR_cr_no_progress            = 204
const IPX_ERROR_cr_inf_or_nan             = 205
const IPX_ERROR_basis_singular            = 301
const IPX_ERROR_basis_almost_singular     = 302
const IPX_ERROR_basis_update_singular     = 303
const IPX_ERROR_basis_repair_overflow     = 304
const IPX_ERROR_basis_repair_search       = 305
const IPX_ERROR_basis_too_ill_conditioned = 306
const IPX_ERROR_interrupt_time            = 999
const IPX_basic                           = 0
const IPX_nonbasic                        = -1
const IPX_nonbasic_lb                     = -1
const IPX_nonbasic_ub                     = -2
const IPX_superbasic                      = -3

# =============================================================================
# LP model in the form
#
#  minimize     obj'*x
#  subject to   A*x {=,<,>} rhs,
#               lb <= x <= ub
#
# lb,ub  are -/+Inf when a variable has no lower/upper bound
# =============================================================================

type Model
    obj::Vector{Cdouble}
    A::spmatrix
    rhs::Vector{Cdouble}
    constr_type::Vector{Char}
    lb::Vector{Cdouble}
    ub::Vector{Cdouble}
end

type Parameters
    display::ipxint
    logfile_ptr::Ptr{Cchar}
    print_interval::Cdouble
    time_limit::Cdouble
    dualize::ipxint
    scale::ipxint
    ipm_maxiter::ipxint
    ipm_feasibility_tol::Cdouble
    ipm_optimality_tol::Cdouble
    ipm_drop_primal::Cdouble
    ipm_drop_dual::Cdouble
    kkt_tol::Cdouble
    crash_basis::ipxint
    dependency_tol::Cdouble
    volume_tol::Cdouble
    rows_per_slice::ipxint
    maxskip_updates::ipxint
    lu_kernel::ipxint
    lu_pivottol::Cdouble
    crossover::ipxint
    crossover_start::Cdouble
    pfeasibility_tol::Cdouble
    dfeasibility_tol::Cdouble
    debug::ipxint
    switchiter::ipxint
    stop_at_switch::ipxint
    update_heuristic::ipxint
    maxpasses::ipxint
    function Parameters()
        return ccall((:ipx_default_parameters, "libipx.so"), Parameters, ())
    end
end

type Info
    status::ipxint
    status_ipm::ipxint
    status_crossover::ipxint
    errflag::ipxint
    num_var::ipxint
    num_constr::ipxint
    num_entries::ipxint
    num_rows_solver::ipxint
    num_cols_solver::ipxint
    num_entries_solver::ipxint
    dualized::ipxint
    dense_cols::ipxint
    dependent_rows::ipxint
    dependent_cols::ipxint
    rows_inconsistent::ipxint
    cols_inconsistent::ipxint
    primal_dropped::ipxint
    dual_dropped::ipxint
    abs_presidual::Cdouble
    abs_dresidual::Cdouble
    rel_presidual::Cdouble
    rel_dresidual::Cdouble
    pobjval::Cdouble
    dobjval::Cdouble
    rel_objgap::Cdouble
    complementarity::Cdouble
    normx::Cdouble
    normy::Cdouble
    normz::Cdouble
    objval::Cdouble
    primal_infeas::Cdouble
    dual_infeas::Cdouble
    iter::ipxint
    kktiter1::ipxint
    kktiter2::ipxint
    basis_repairs::ipxint
    updates_start::ipxint
    updates_ipm::ipxint
    updates_crossover::ipxint
    time_total::Cdouble
    time_ipm1::Cdouble
    time_ipm2::Cdouble
    time_starting_basis::Cdouble
    time_crossover::Cdouble
    time_kkt_factorize::Cdouble
    time_kkt_solve::Cdouble
    time_maxvol::Cdouble
    time_cr1::Cdouble
    time_cr1_AAt::Cdouble
    time_cr1_pre::Cdouble
    time_cr2::Cdouble
    time_cr2_NNt::Cdouble
    time_cr2_B::Cdouble
    time_cr2_Bt::Cdouble
    ftran_sparse::Cdouble
    btran_sparse::Cdouble
    time_ftran::Cdouble
    time_btran::Cdouble
    time_lu_invert::Cdouble
    time_lu_update::Cdouble
    mean_fill::Cdouble
    max_fill::Cdouble
    time_symb_invert::Cdouble
    maxvol_updates::ipxint
    maxvol_skipped::ipxint
    maxvol_passes::ipxint
    tbl_nnz::ipxint
    tbl_max::Cdouble
    frobnorm_squared::Cdouble
    lambdamax::Cdouble
    volume_increase::Cdouble
    # constructs an Info object with all fields set to 0
    function Info()
        self = new()
        for f = fieldnames(Info)
            setfield!(self, f, convert(fieldtype(Info,f), 0))
        end
        return self
    end
end

# parses all lines of the form "info.fieldname value" where fieldname is
# alphanumeric and value is a number
function ParseInfo(lines::Array{String})
    ipxinfo = Info();
    for line = lines
        m = match(r"^\s*info\.([a-zA-Z0-9_]+)\s+([e\+\-\.0-9]+)$", line)
        if m != nothing && length(m.captures) == 2
            f = Symbol(m.captures[1]) # fieldname
            value = m.captures[2]     # value as string
            setfield!(ipxinfo, f, parse(fieldtype(Info, f), value))
        end
    end
    return ipxinfo
end

function ParseInfo(filename::String)
    f = open(filename, "r")
    lines = readlines(f)
    close(f)
    return ParseInfo(lines)
end

type InteriorSolution
    x::Vector{Cdouble}
    xl::Vector{Cdouble}
    xu::Vector{Cdouble}
    slack::Vector{Cdouble}
    y::Vector{Cdouble}
    zl::Vector{Cdouble}
    zu::Vector{Cdouble}
    function InteriorSolution(m::ipxint, n::ipxint)
        self = new()
        self.x = zeros(Cdouble,n)
        self.xl = zeros(Cdouble,n)
        self.xu = zeros(Cdouble,n)
        self.slack = zeros(Cdouble,m)
        self.y = zeros(Cdouble,m)
        self.zl = zeros(Cdouble,n)
        self.zu = zeros(Cdouble,n)
        return self
    end
end

type BasicSolution
    x::Vector{Cdouble}
    slack::Vector{Cdouble}
    y::Vector{Cdouble}
    z::Vector{Cdouble}
    cbasis::Vector{ipxint}
    vbasis::Vector{ipxint}
    function BasicSolution(m::ipxint, n::ipxint)
        self = new()
        self.x = zeros(Cdouble,n)
        self.slack = zeros(Cdouble,m)
        self.y = zeros(Cdouble,m)
        self.z = zeros(Cdouble,n)
        self.cbasis = zeros(ipxint,m)
        self.vbasis = zeros(ipxint,n)
        return self
    end
end

type LPSolver
    solver::Ptr{Void}
    function LPSolver()
        p_solver = Ref{Ptr{Void}}(C_NULL)
        ccall((:ipx_new, "libipx.so"), Void, (Ptr{Ptr{Void}},), p_solver)
        if p_solver[] == C_NULL
            error("ipx_new failed")
        end
        this = new(p_solver[])
        finalizer(this, ipx_free)
        return this
    end
end

function ipx_free(this::LPSolver)
    p_solver = Ref{Ptr{Void}}(this.solver)
    ccall((:ipx_free, "libipx.so"), Void, (Ptr{Ptr{Void}},), p_solver)
    this.solver = p_solver[]
end

function Solve(this::LPSolver, model::Model, logfile::String="")
    m,n = size(model.A)
    Ap = model.A.colptr - 1
    Ai = model.A.rowval - 1
    Ax = model.A.nzval
    rhs = model.rhs
    constr_type = convert(Array{Cchar,1}, model.constr_type)
    obj = model.obj
    lb = model.lb
    ub = model.ub
    logfile = convert(Vector{UInt8}, logfile)
    logfile = convert(Vector{Int8}, logfile)
    append!(logfile, 0)
    logfile_ptr = GetParameter(this, :logfile_ptr)
    SetParameter(this, :logfile_ptr, pointer(logfile))
    status = ccall((:ipx_solve, "libipx.so"), ipxint,
                   (Ptr{Void}, ipxint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                    ipxint, Ptr{ipxint}, Ptr{ipxint}, Ptr{Cdouble},
                    Ptr{Cdouble}, Ptr{Cchar}),
                   this.solver, n, obj, lb, ub,
                   m, Ap, Ai, Ax, rhs, constr_type)
    SetParameter(this, :logfile_ptr, logfile_ptr)
    return status
end

function Clear(this::LPSolver)
    ccall((:ipx_clear_model, "libipx.so"), Void,
          (Ptr{Void},), this.solver)
    nothing
end

function GetInteriorSolution(this::LPSolver)
    solverinfo = GetInfo(this)
    m = solverinfo.num_constr
    n = solverinfo.num_var
    solution = InteriorSolution(m,n)
    err = ccall((:ipx_get_interior_solution, "libipx.so"), ipxint,
                (Ptr{Void}, Ptr{Cdouble}, Ptr{Cdouble},
                 Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                 Ptr{Cdouble}, Ptr{Cdouble}),
                this.solver, solution.x, solution.xl, solution.xu,
                solution.slack, solution.y, solution.zl, solution.zu)
    if err != 0
        error("ipx_get_interior_solution failed")
    end
    return solution
end

function GetBasicSolution(this::LPSolver)
    solverinfo = GetInfo(this)
    m = solverinfo.num_constr
    n = solverinfo.num_var
    solution = BasicSolution(m,n)
    err = ccall((:ipx_get_basic_solution, "libipx.so"), ipxint,
                (Ptr{Void}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                 Ptr{Cdouble}, Ptr{ipxint}, Ptr{ipxint}),
                this.solver, solution.x, solution.slack, solution.y,
                solution.z, solution.cbasis, solution.vbasis)
    if err != 0
        error("ipx_get_basic_solution failed")
    end
    return solution
end

function GetBasis(this::LPSolver)
    solverinfo = GetInfo(this)
    m = solverinfo.num_constr
    n = solverinfo.num_var
    cbasis = Array{ipxint}(m)
    vbasis = Array{ipxint}(n)
    err = ccall((:ipx_get_basis, "libipx.so"), ipxint,
                (Ptr{Void}, Ptr{ipxint}, Ptr{ipxint}),
                this.solver, cbasis, vbasis)
    if err != 0
        error("ipx_get_basis failed")
    end
    return (cbasis,vbasis)
end

function GetIterate(this::LPSolver)
    solverinfo = GetInfo(this)
    m = solverinfo.num_rows_solver
    n = solverinfo.num_cols_solver
    x = Array{Cdouble}(n)
    y = Array{Cdouble}(m)
    zl = Array{Cdouble}(n)
    zu = Array{Cdouble}(n)
    xl = Array{Cdouble}(n)
    xu = Array{Cdouble}(n)
    err = ccall((:ipx_get_iterate, "libipx.so"), ipxint,
                (Ptr{Void}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                 Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                this.solver, x, y, zl, zu, xl, xu)
    if err != 0
        error("ipx_get_iterate failed")
    end
    return (x,y,zl,zu,xl,xu)
end

function GetKKTMatrix(this::LPSolver)
    solverinfo = GetInfo(this)
    m = solverinfo.num_rows_solver
    n = solverinfo.num_cols_solver
    nz = solverinfo.num_entries_solver
    AIp = Array{ipxint}(n+1)
    AIi = Array{ipxint}(nz)
    AIx = Array{Cdouble}(nz)
    g = Array{Cdouble}(n)
    err = ccall((:ipx_get_kktmatrix, "libipx.so"), ipxint,
                (Ptr{Void}, Ptr{ipxint}, Ptr{ipxint}, Ptr{Cdouble},
                 Ptr{Cdouble}),
                this.solver, AIp, AIi, AIx, g)
    if err != 0
        error("ipx_get_kktmatrix failed")
    end
    AIp[:] += 1
    AIi[:] += 1
    AI = spmatrix(m, n, AIp, AIi, AIx)
    return (AI, g)
end

function SymbolicInvert(this::LPSolver)
    solverinfo = GetInfo(this)
    m = solverinfo.num_rows_solver
    rowcounts = Array{ipxint}(m)
    colcounts = Array{ipxint}(m)
    err = ccall((:ipx_symbolic_invert, "libipx.so"), ipxint,
                (Ptr{Void}, Ptr{ipxint}, Ptr{ipxint}),
                this.solver, rowcounts, colcounts)
    if err != 0
        error("ipx_symbolic_invert failed")
    end
    @assert sum(rowcounts) == sum(colcounts)
    return (rowcounts,colcounts)
end

function GetParameters(this::LPSolver)
    return ccall((:ipx_get_parameters, "libipx.so"),
                 Parameters, (Ptr{Void},), this.solver)
end

function SetParameters(this::LPSolver, params::Parameters)
    ccall((:ipx_set_parameters, "libipx.so"), Void,
          (Ptr{Void}, Parameters), this.solver, params)
    nothing
end

# value = GetParameter(solver, :name)
function GetParameter(this::LPSolver, name::Symbol)
    params = GetParameters(this)
    return getfield(params, name)
end

# SetParameter(solver, :name, value)
function SetParameter(this::LPSolver, name::Symbol, value)
    params = GetParameters(this)
    setfield!(params, name, value)
    SetParameters(this, params)
end

function GetInfo(this::LPSolver)
    return ccall((:ipx_get_info, "libipx.so"),
                 Info, (Ptr{Void},), this.solver)
end

include("ipx_utils.jl")

end # module ipx
