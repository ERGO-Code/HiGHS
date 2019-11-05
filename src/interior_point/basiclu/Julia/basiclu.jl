#
# Julia interface to BASICLU
#
# BASICLU must have been compiled with the integer type matching cint (see
# below). The shared library must be in the load path.
#

module basiclu

const cint = Int64
const cdbl = Cdouble
const cint_ptr = Ptr{cint}
const cdbl_ptr = Ptr{cdbl}
const cvec = Array{cdbl,1}
const spvector = SparseVector{cdbl,cint}
const spmatrix = SparseMatrixCSC{cdbl,cint}

# comment out if the MAT module is installed
# include("test.jl")

# =============================================================================
# BASICLU defines
# =============================================================================

# size of istore, xstore
const BASICLU_SIZE_ISTORE_1 = 1024
const BASICLU_SIZE_ISTORE_M = 21
const BASICLU_SIZE_XSTORE_1 = 1024
const BASICLU_SIZE_XSTORE_M = 4

# status codes
const BASICLU_OK = 0
const BASICLU_REALLOCATE = 1
const BASICLU_WARNING_singular_matrix = 2
const BASICLU_ERROR_invalid_store = -1
const BASICLU_ERROR_invalid_call = -2
const BASICLU_ERROR_argument_missing = -3
const BASICLU_ERROR_invalid_argument = -4
const BASICLU_ERROR_maximum_updates = -5
const BASICLU_ERROR_singular_update = -6

# user parameters in xstore (offset by 1)
const BASICLU_MEMORYL = 2
const BASICLU_MEMORYU = 3
const BASICLU_MEMORYW = 4
const BASICLU_DROP_TOLERANCE = 5
const BASICLU_ABS_PIVOT_TOLERANCE = 6
const BASICLU_REL_PIVOT_TOLERANCE = 7
const BASICLU_BIAS_NONZEROS = 8
const BASICLU_MAXN_SEARCH_PIVOT = 9
const BASICLU_PAD = 10
const BASICLU_STRETCH = 11
const BASICLU_COMPRESSION_THRESHOLD = 12
const BASICLU_SPARSE_THRESHOLD = 13
const BASICLU_REMOVE_COLUMNS = 14
const BASICLU_SEARCH_ROWS = 15

# user readable from xstore (offset by 1)
const BASICLU_DIM = 65
const BASICLU_STATUS = 66
const BASICLU_ADD_MEMORYL = 67
const BASICLU_ADD_MEMORYU = 68
const BASICLU_ADD_MEMORYW = 69

const BASICLU_NUPDATE = 71
const BASICLU_NFORREST = 72
const BASICLU_NFACTORIZE = 73
const BASICLU_NUPDATE_TOTAL = 74
const BASICLU_NFORREST_TOTAL = 75
const BASICLU_NSYMPERM_TOTAL = 76
const BASICLU_LNZ = 77
const BASICLU_UNZ = 78
const BASICLU_RNZ = 79
const BASICLU_MIN_PIVOT = 80
const BASICLU_MAX_PIVOT = 81
const BASICLU_UPDATE_COST = 82
const BASICLU_TIME_FACTORIZE = 83
const BASICLU_TIME_SOLVE = 84
const BASICLU_TIME_UPDATE = 85
const BASICLU_TIME_FACTORIZE_TOTAL = 86
const BASICLU_TIME_SOLVE_TOTAL = 87
const BASICLU_TIME_UPDATE_TOTAL = 88
const BASICLU_LFLOPS = 89
const BASICLU_UFLOPS = 90
const BASICLU_RFLOPS = 91
const BASICLU_CONDEST_L = 92
const BASICLU_CONDEST_U = 93
const BASICLU_NORM_L = 95
const BASICLU_NORM_U = 96
const BASICLU_NORMEST_LINV = 97
const BASICLU_NORMEST_UINV = 98
const BASICLU_MATRIX_ONENORM = 99
const BASICLU_MATRIX_INFNORM = 100
const BASICLU_RESIDUAL_TEST = 112

const BASICLU_MATRIX_NZ = 101
const BASICLU_RANK = 102
const BASICLU_BUMP_SIZE = 103
const BASICLU_BUMP_NZ = 104
const BASICLU_NSEARCH_PIVOT = 105
const BASICLU_NEXPAND = 106
const BASICLU_NGARBAGE = 107
const BASICLU_FACTOR_FLOPS = 108
const BASICLU_TIME_SINGLETONS = 109
const BASICLU_TIME_SEARCH_PIVOT = 110
const BASICLU_TIME_ELIM_PIVOT = 111

const BASICLU_PIVOT_ERROR = 121

# parameters for Julia driver function
const realloc_factor = 1.5

# =============================================================================
# type BLU
# =============================================================================
"""
    type BLU

Hold BASICLU object. You can access @xstore to set parameters and get output
values from BASICLU routines.
"""
type BLU
    istore::Array{cint,1}
    xstore::Array{cdbl,1}
    Li::Array{cint,1}
    Lx::Array{cdbl,1}
    Ui::Array{cint,1}
    Ux::Array{cdbl,1}
    Wi::Array{cint,1}
    Wx::Array{cdbl,1}
    work0::Array{cdbl,1}        # workspace zeroed
    iwork::Array{cint,1}        # workspace uninitialized
    function BLU(m::cint)
        isize = BASICLU_SIZE_ISTORE_1 + BASICLU_SIZE_ISTORE_M * m
        xsize = BASICLU_SIZE_XSTORE_1 + BASICLU_SIZE_XSTORE_M * m
        fsize = m               # initial length of Li, Lx, Ui, Ux, Wi, Wx
        istore = Array{cint}(isize)
        xstore = Array{cdbl}(xsize)
        Li = Array{cint}(fsize)
        Lx = Array{cdbl}(fsize)
        Ui = Array{cint}(fsize)
        Ux = Array{cdbl}(fsize)
        Wi = Array{cint}(fsize)
        Wx = Array{cdbl}(fsize)
        work0 = zeros(cdbl, m)
        iwork = Array{cint}(m)
        err = ccall((:basiclu_initialize, "libbasiclu.so"), cint,
                    (cint, cint_ptr, cdbl_ptr), m, istore, xstore)
        if err != BASICLU_OK
            msg = @sprintf("basiclu_initialize() status code: %d", err)
            error(msg)
        end
        xstore[BASICLU_MEMORYL] = fsize
        xstore[BASICLU_MEMORYU] = fsize
        xstore[BASICLU_MEMORYW] = fsize
        return new(istore, xstore, Li, Lx, Ui, Ux, Wi, Wx, work0, iwork)
    end
end

# =============================================================================
# realloc
# =============================================================================

function realloc(this::BLU)
        addmemL = convert(cint, this.xstore[BASICLU_ADD_MEMORYL])
        addmemU = convert(cint, this.xstore[BASICLU_ADD_MEMORYU])
        addmemW = convert(cint, this.xstore[BASICLU_ADD_MEMORYW])
        if addmemL > 0
            memsize = length(this.Li) + addmemL
            memsize = ceil(realloc_factor*memsize)
            memsize = convert(cint, memsize)
            @assert memsize >= length(this.Li) + addmemL
            resize!(this.Li, memsize)
            resize!(this.Lx, memsize)
            this.xstore[BASICLU_MEMORYL] = memsize
        end
        if addmemU > 0
            memsize = length(this.Ui) + addmemU
            memsize = ceil(realloc_factor*memsize)
            memsize = convert(cint, memsize)
            @assert memsize >= length(this.Ui) + addmemU
            resize!(this.Ui, memsize)
            resize!(this.Ux, memsize)
            this.xstore[BASICLU_MEMORYU] = memsize
        end
        if addmemW > 0
            memsize = length(this.Wi) + addmemW
            memsize = ceil(realloc_factor*memsize)
            memsize = convert(cint, memsize)
            @assert memsize >= length(this.Wi) + addmemW
            resize!(this.Wi, memsize)
            resize!(this.Wx, memsize)
            this.xstore[BASICLU_MEMORYW] = memsize
        end
end

# =============================================================================
# gather_work
# =============================================================================

function gather_work(this::BLU, nz::cint)
    m = convert(cint, this.xstore[BASICLU_DIM])
    nzind = this.iwork[1:nz] + 1
    sort!(nzind)
    lhs = spvector(m, nzind, this.work0[nzind])
    this.work0[nzind] = 0
    return lhs
end

# =============================================================================
# initialize
# =============================================================================
"""
    initialize(m)

Return a new BASICLU object for matrices of dimension m.
"""
function initialize(m::cint)
    return BLU(m)
end

# =============================================================================
# factorize
# =============================================================================
"""
    factorize(this, B)

Load sparse matrix into BASICLU object and factorize it. @B must be square
and have the same dimension for which @this was initialized.

Return status code.
"""
function factorize(this::BLU, B::spmatrix)
    m = convert(cint, this.xstore[BASICLU_DIM])
    nrow, ncol = size(B)
    @assert nrow == ncol
    @assert nrow == m
    Bp = B.colptr-1
    Bi = B.rowval-1
    Bx = B.nzval                # don't need a copy
    c0ntinue = 0
    err = BASICLU_OK
    while true
        err = ccall((:basiclu_factorize, "libbasiclu.so"), cint,
                    (cint_ptr, cdbl_ptr,
                     cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr,
                     cint_ptr, cint_ptr, cint_ptr, cdbl_ptr, cint),
                    this.istore, this.xstore,
                    this.Li, this.Lx, this.Ui, this.Ux, this.Wi, this.Wx,
                    Bp, pointer(Bp, 2), Bi, Bx, c0ntinue)
        if err != BASICLU_REALLOCATE break; end
        realloc(this)
        c0ntinue = 1
    end
    return err
end

# =============================================================================
# get_factors
# =============================================================================
"""
    get_factors(this)

Export LU factors after fresh factorization.

# Example
```jldoctest
julia> m = 100;
julia> B = spdiagm((ones(m-1), 4*ones(m), ones(m-1)), [-1 0 1]);
julia> this = basiclu.initialize(m);
julia> basiclu.factorize(this, B)
0
julia> (L,U,p,q) = basiclu.get_factors(this);
julia> norm(L*U-B[p,q], Inf)
1.1102230246251565e-16
```
"""
function get_factors(this::BLU)
    m = convert(cint, this.xstore[BASICLU_DIM])
    Lnz = convert(cint, this.xstore[BASICLU_LNZ])
    Unz = convert(cint, this.xstore[BASICLU_UNZ])
    rowperm = Array{cint}(m)
    colperm = Array{cint}(m)
    Lp = Array{cint}(m+1)
    Li = Array{cint}(Lnz+m)
    Lx = Array{cdbl}(Lnz+m)
    Up = Array{cint}(m+1)
    Ui = Array{cint}(Unz+m)
    Ux = Array{cdbl}(Unz+m)
    err = ccall((:basiclu_get_factors, "libbasiclu.so"), cint,
                (cint_ptr, cdbl_ptr,
                 cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr,
                 cint_ptr, cint_ptr,
                 cint_ptr, cint_ptr, cdbl_ptr,
                 cint_ptr, cint_ptr, cdbl_ptr),
                this.istore, this.xstore,
                this.Li, this.Lx, this.Ui, this.Ux, this.Wi, this.Wx,
                rowperm, colperm, Lp, Li, Lx, Up, Ui, Ux)
    if err != BASICLU_OK
        msg = @sprintf("basiclu_get_factors() status code: %d", err)
        error(msg)
    end
    rowperm += 1
    colperm += 1
    L = spmatrix(m, m, Lp+1, Li+1, Lx)
    U = spmatrix(m, m, Up+1, Ui+1, Ux)
    @assert istril(L)
    @assert istriu(U)
    return (L, U, rowperm, colperm)
end

# =============================================================================
# solve
# =============================================================================
"""
    solve(this, rhs, trans)

Solve linear system with factorized matrix.

@rhs must be a dense or sparse vector. The solution is returned sparse when @rhs
is sparse and dense otherwise.

@trans must be 't' or 'T' for transposed solve, and any other character for
forward solve.
"""
function solve(this::BLU, rhs::cvec, trans::Char) # dense solve
    m = convert(cint, this.xstore[BASICLU_DIM])
    @assert length(rhs) == m
    lhs = cvec(m)
    err = ccall((:basiclu_solve_dense, "libbasiclu.so"), cint,
                (cint_ptr, cdbl_ptr,
                 cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr,
                 cdbl_ptr, cdbl_ptr, Char),
                this.istore, this.xstore,
                this.Li, this.Lx, this.Ui, this.Ux, this.Wi, this.Wx,
                rhs, lhs, trans)
    if err != BASICLU_OK
        msg = @sprintf("basiclu_solve_dense() status code: %d", err)
        error(msg)
    end
    return lhs
end

function solve(this::BLU, rhs::spvector, trans::Char) # sparse solve
    m = convert(cint, this.xstore[BASICLU_DIM])
    @assert length(rhs) == m
    nzlhs = Ref{cint}(0)
    err = ccall((:basiclu_solve_sparse, "libbasiclu.so"), cint,
                (cint_ptr, cdbl_ptr,
                 cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr,
                 cint, cint_ptr, cdbl_ptr,
                 Ref{cint}, cint_ptr, cdbl_ptr, Char),
                this.istore, this.xstore,
                this.Li, this.Lx, this.Ui, this.Ux, this.Wi, this.Wx,
                nnz(rhs), rhs.nzind-1, rhs.nzval,
                nzlhs, this.iwork, this.work0, trans)
    nzlhs = nzlhs[]
    lhs = gather_work(this, nzlhs)
    if err != BASICLU_OK
        msg = @sprintf("basiclu_solve_sparse() status code: %d", err)
        error(msg)
    end
    return lhs
end

# =============================================================================
# solve4update
# =============================================================================
"""
    solve4update(this, rhs, getsolution=false)

Solve linear system in preparation to update the factorization.

@rhs must be a sparse vector or a column index. When @rhs is a vector, then it
is the column to be inserted into the factorized matrix. When @rhs is an index,
than it is the column of the factorized matrix to be replaced.

@getsolution indicates if the solution is to be returned, or only the update is
to be prepared.
"""
function solve4update(this::BLU, rhs, getsolution::Bool=false)
    if typeof(rhs) != spvector && !(length(rhs) == 1 && isinteger(rhs))
        error("rhs must be a sparse vector or a column index")
    end
    m = convert(cint, this.xstore[BASICLU_DIM])
    trans = typeof(rhs) != spvector
    if trans
        @assert 1 <= rhs && rhs <= m
    else
        @assert length(rhs) == m
    end
    nzrhs = trans ? 0 : nnz(rhs)
    irhs = trans ? Ref{cint}(rhs-1) : rhs.nzind-1
    xrhs = trans ? C_NULL : rhs.nzval
    nzlhs = getsolution ? Ref{cint}(0) : C_NULL
    err = BASICLU_OK
    while true
        err = ccall((:basiclu_solve_for_update, "libbasiclu.so"), cint,
                    (cint_ptr, cdbl_ptr,
                     cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr,
                     cint, cint_ptr, cdbl_ptr,
                     cint_ptr, cint_ptr, cdbl_ptr, Char),
                    this.istore, this.xstore,
                    this.Li, this.Lx, this.Ui, this.Ux, this.Wi, this.Wx,
                    nzrhs, irhs, xrhs,
                    nzlhs, this.iwork, this.work0, trans ? 'T' : 'N')
        if err != BASICLU_REALLOCATE break; end
        realloc(this)
    end
    lhs = nothing
    if getsolution
        nzlhs = nzlhs[]
        lhs = gather_work(this, nzlhs)
    end
    if err != BASICLU_OK
        msg = @sprintf("basiclu_solve_for_update() status code: %d", err)
        error(msg)
    end
    return lhs
end

# =============================================================================
# update
# =============================================================================
"""
    update(this, xtbl)

Update the factorization after a column modification. The column position and
the new column must have been set in previous calls to @solve4update.

@xtbl is the pivot element in the simplex tableau. Used only to compute the
pivot error.

Return estimated error of the new pivot element.
"""
function update(this::BLU, xtbl::cdbl)
    err = BASICLU_OK
    while true
        err = ccall((:basiclu_update, "libbasiclu.so"), cint,
                    (cint_ptr, cdbl_ptr,
                     cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr, cint_ptr, cdbl_ptr,
                     cdbl),
                    this.istore, this.xstore,
                    this.Li, this.Lx, this.Ui, this.Ux, this.Wi, this.Wx,
                    xtbl)
        if err != BASICLU_REALLOCATE break; end
        realloc(this)
    end        
    if err != BASICLU_OK
        msg = @sprintf("basiclu_update() status code: %d", err)
        error(msg)
    end
    return this.xstore[BASICLU_PIVOT_ERROR]
end

end
