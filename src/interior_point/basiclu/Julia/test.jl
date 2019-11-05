using MAT

# =============================================================================
# readmat: read variables from MATLAB file
# =============================================================================

function readmat(file, names...)
    n = length(names)
    objects = Array{Any}(n)
    f1 = matopen(file)
    for k = 1:n
        objects[k] = read(f1, names[k])
    end
    close(f1)
    return tuple(objects...)
end

# =============================================================================
# test_factorize
# =============================================================================
"""
    test_factorize(testdir, trans=false)

For all *.mat files in @testdir read matrix B and factorize it. Monitor residual
of factorization and forward/backward solves.

@trans specifies if B or its transposed are factorized.
"""
function test_factorize(testdir::String, trans=false)
    files = readdir(testdir)
    for f in files
        if length(f) < 4 || f[end-3:end] != ".mat"
            continue
        end
        @printf(" %-24s", f)
        B, = readmat(string(testdir, f), "B")
        Bt = transpose(B)
        if trans
            (B,Bt) = (Bt,B)
        end
        m = size(B,1)
        blu = initialize(m)
        err = factorize(blu, B)
        if err != BASICLU_OK
            @printf(" failed (%d)\n", err)
            continue
        end
        (L,U,p,q) = get_factors(blu)
        res = norm(L*U-B[p,q], Inf)
        rhs = ones(m)
        lhs = solve(blu, rhs, 'N')
        res = max(res, norm(B*lhs-rhs, Inf))
        lhs = solve(blu, rhs, 'T')
        res = max(res, norm(Bt*lhs-rhs, Inf))
        rhs = spvector(m, [1;m], [1.0;1.0])
        lhs = solve(blu, rhs, 'N')
        res = max(res, norm(B*lhs-rhs, Inf))
        lhs = solve(blu, rhs, 'T')
        res = max(res, norm(Bt*lhs-rhs, Inf))
        lhs = solve4update(blu, rhs, true)
        res = max(res, norm(B*lhs-rhs, Inf))
        lhs = solve4update(blu, m-1, true)
        rhs = zeros(m); rhs[m-1] = 1
        res = max(res, norm(Bt*lhs-rhs, Inf))
        @printf("%.2e\n", res)
    end
    nothing
end

# =============================================================================
# test_update
# =============================================================================
"""
    test_update(testdir)

For all *.mat files in @testdir read matrix A and the pivot sequence. Starting
from the slack basis apply pivot operations. Monitor residuals to
forward/backward solves after each 100 updates.
"""
function test_update(testdir::String)
    files = readdir(testdir)
    for f in files
        if length(f) < 4 || f[end-3:end] != ".mat"
            continue
        end
        @printf(" %-24s", f)
        A,invar,outvar = readmat(string(testdir, f), "A", "invar", "outvar")
        invar = invar[:]
        outvar = outvar[:]
        m,n = size(A)
        A1 = [A speye(m)]
        basis = Array{cint}(m)
        basis[:] = collect(1:m) + n # slack basis
        res,nfactor,nforrest,nperm = test_update(A1, basis, invar, outvar)
        @printf("%.2e      %5d factor, %6d forrest, %6d perm\n",
                res, nfactor, nforrest, nperm)
    end
    nothing
end

function test_update(A::spmatrix, basis::Array{cint,1}, invar::Array{cint,1},
                     outvar::Array{cint,1})
    m,n = size(A)
    niter = length(invar)
    map2basis = zeros(cint, n)
    map2basis[basis] = 1:m
    blu = initialize(m)
    # blu.xstore[BASICLU_REL_PIVOT_TOLERANCE] = 0.5
    err = factorize(blu, A[:,basis])
    if err != BASICLU_OK
        return (NaN, 0, 0, 0)
    end
    res = 0.0
    for k = 1:niter
        if invar[k] == outvar[k] # bound flip
            continue
        end
        j = invar[k]
        p = map2basis[outvar[k]]
        @assert p > 0
        lhs = solve4update(blu, A[:,j], true)
        xtbl = lhs[p]
        solve4update(blu, p)
        piverr = update(blu, xtbl)
        basis[p] = j
        map2basis[j] = p
        map2basis[outvar[k]] = 0
        cost = blu.xstore[BASICLU_UPDATE_COST]
        if piverr > 1e-10 || cost > 1.0
            err = factorize(blu, A[:,basis])
            if err != BASICLU_OK
                return (NaN, 0, 0, 0)
            end
        end
        if k%100 == 0
            B = A[:,basis]
            rhs = spvector(m, [1;m], [1.0;1.0])
            lhs = solve(blu, rhs, 'N')
            res = max(res, norm(B*lhs-rhs, Inf))
            lhs = solve(blu, rhs, 'T')
            res = max(res, norm(B'*lhs-rhs, Inf))
        end
    end
    nfactor = convert(cint, blu.xstore[BASICLU_NFACTORIZE])
    nforrest = convert(cint, blu.xstore[BASICLU_NFORREST_TOTAL])
    nperm = convert(cint, blu.xstore[BASICLU_NUPDATE_TOTAL]) - nforrest
    return (res, nfactor, nforrest, nperm)
end
