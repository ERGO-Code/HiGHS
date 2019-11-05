/*
 * maxvolume.c
 *
 * Copyright (C) 2016-2017  Lukas Schork
 *
 * Command line program to compute a maximum volume basis [1,2]
 *
 * Purpose:
 *
 * Given an m-by-n matrix A, usually m < n, find m columns of A that form a
 * nonsingular matrix B with maximum volume in A. When A does not have m
 * linearly independent columns (in particular when m > n), then the remaining
 * columns of B are columns of the identity matrix. The number of such columns
 * is the row rank deficiency of A.
 *
 * The basis B has (local) maximum volume in A iff for every column a_j of A
 * all entries in B^{-1}*a_j are less than or equal to 1 in absolute value.
 *
 * Algorithm:
 *
 * Append 1e-8 * identity matrix to A and use these columns as initial basis.
 * Then repeatedly pass over the columns of A. If column a_j is not in B and
 * B^{-1}*a_j has an entry larger than @volumetol in absolute value, then a_j
 * replaces one column in B. The algorithm stops when one pass over the columns
 * of A did not change B, or when @maxpass passes are done.
 *
 * Usually the major cost of the algorithm is computing B^{-1}*a_j. When these
 * vectors are dense, then one pass over the columns of A has time complexity
 * Omega(m*n) and the algorithm is impractical for large scale matrices. When
 * these vectors are sparse (as they frequently are in linear programming
 * problems), then the BASICLU routines will exploit the sparsity and the
 * computational cost for one pass is proportional to the number of nonzeros in
 * B^{-1}*A. See the LP problems in data/.
 *
 * Note:
 *
 * On Maragal_4 the algorithm fails with default parameters because a
 * refactorization reports a singularity in the basis. Debugging has shown
 * that the final entry in the active submatrix is 2.03e-15, which is below
 * the default absolute pivot tolerance (1e-14). Changing the absolute pivot
 * tolerance to 1e-15 makes the algorithm complete, as does tightening the
 * relative pivot tolerance to 0.5. In the latter case I'm not sure if that
 * results accidently from a different refactorization point or if it fixes
 * numerical instability in the factorization.
 *
 * Parameters:
 *
 * maxvolume <matrix.mtx> [volumetol [maxpass]]
 *
 * <matrix.mtx>: name of a file containing the matrix A in matrix market format
 *               (real, sparse, general)
 * volumetol:    tolerance >= 1 on the absolute value of B^{-1}*a_j
 *               default: 1.1
 * maxpass:      maximum # passes over the columns of A
 *               default: 2
 *
 * [1] C. T. Pan, "On the existence and computation of rank-revealing LU
 *     factorizations", Linear Algebra Appl., 2000
 * [2] S. A. Goreinov, I. V. Oseledets, D. V. Savostyanov, E. E. Tyrtyshnikov,
 *     N. L. Zamarashkin, "How to find a good submatrix", technical report, 2008
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "basiclu.h"
#include "mmio.h"

/* internal error codes, supplementing BASICLU status codes */
#define IO_ERROR 102

static lu_int maxvolume(struct basiclu_object *factor, lu_int ncol,
                        const lu_int *Ap, const lu_int *Ai, const double *Ax,
                        lu_int *basis, lu_int *isbasic, double volumetol,
                        lu_int *p_nupdate);

int main(int argc, const char *argv[])
{
    lu_int *Ap      = NULL;
    lu_int *Ai      = NULL;
    double *Ax      = NULL;
    lu_int *basis   = NULL;
    lu_int *isbasic = NULL;
    lu_int *count   = NULL;
    int *mm_I       = NULL;
    int *mm_J       = NULL;
    double *mm_val  = NULL;
    int mm_M, mm_N, mm_nz;
    lu_int m, n, nz, i, j, l, put, err = 0;
    struct basiclu_object factor;

    double volumetol = 1.1;     /* command line parameters */
    long maxpass = 2;

    basiclu_obj_initialize(&factor, 0); /* nullify members */

    if (argc < 2 || argc > 4) {
        printf(" usage: maxvolume <matrix.mtx> [volumetol [maxpass]]\n");
        return 101;
    }
    if (argc >= 3)
        volumetol = atof(argv[2]);
    if (argc >= 4)
        maxpass = atol(argv[3]);

    err = mm_read_unsymmetric_sparse(argv[1], &mm_M, &mm_N, &mm_nz,
                                     &mm_val, &mm_I, &mm_J);
    if (err) {
        err = IO_ERROR;
        goto cleanup;
    }
    printf(" matrix: %d rows, %d columns, %d nonzeros\n", mm_M, mm_N, mm_nz);
    printf(" parameters: volumetol = %.2f, maxpass = %ld\n",
           volumetol, maxpass);

    m = mm_M;                   /* convert to lu_int */
    n = mm_N;
    nz = mm_nz;

    /*
     * Convert coordinate format to compressed column format and append 1e-8 *
     * identity matrix. Use work array @count to count the number of nonzeros
     * per column. While filling the matrix, @count[j] is the next unused
     * position in column j.
     */
    Ap = malloc((m+n+1)*sizeof(lu_int));
    Ai = malloc((nz+m)*sizeof(lu_int));
    Ax = malloc((nz+m)*sizeof(double));
    count = calloc(n, sizeof(lu_int));
    if (!Ap || !Ai || !Ax || !count) {
        err = BASICLU_ERROR_out_of_memory;
        goto cleanup;
    }

    for (l = 0; l < nz; l++)
        count[mm_J[l]]++;
    put = 0;
    for (j = 0; j < n; j++)
    {
        Ap[j] = put;
        put += count[j];
        count[j] = Ap[j];
    }
    Ap[n] = put;
    for (l = 0; l < nz; l++)
    {
        j = mm_J[l];
        put = count[j]++;
        Ai[put] = mm_I[l];
        Ax[put] = mm_val[l];
    }
    for (i = 0; i < m; i++)
    {
        Ai[Ap[n+i]] = i;
        Ax[Ap[n+i]] = 1e-8;
        Ap[n+i+1] = Ap[n+i] + 1;
    }
    free(mm_I);
    free(mm_J);
    free(mm_val);
    free(count);
    mm_I = NULL;
    mm_J = NULL;
    mm_val = NULL;
    count = NULL;

    /*
     * Initialize @factor. Initialize @basis, @isbasic to logical basis and run
     * maxvolume() until the basis does not change any more or @maxpass passes
     * through the matrix are done.
     */
    err = basiclu_obj_initialize(&factor, m);
    if (err != BASICLU_OK)
        goto cleanup;

    basis = malloc(m*sizeof(lu_int));
    isbasic = calloc(m+n, sizeof(lu_int));
    if (!basis || !isbasic) {
        err = BASICLU_ERROR_out_of_memory;
        goto cleanup;
    }
    for (i = 0; i < m; i++)
        basis[i] = n+i;
    for (j = n; j < m+n; j++)
        isbasic[j] = 1;

    long pass, changed = 1;
    for (pass = 0; pass < maxpass && changed; pass++)
    {
        lu_int nupdate;
        err = maxvolume(&factor, m+n, Ap, Ai, Ax, basis, isbasic, volumetol,
                        &nupdate);
        if (err != BASICLU_OK)
            goto cleanup;
        changed = nupdate > 0;
    }

    /*
     * The number of logical columns in the basis is the number of row rank
     * deficiencies in the input matrix. Print total statistics.
     */
    long rankdef = 0;
    for (j = n; j < m+n; j++)
        rankdef += isbasic[j] != 0;

    long nupdate = factor.xstore[BASICLU_NUPDATE_TOTAL];
    long nforrest = factor.xstore[BASICLU_NFORREST_TOTAL];
    long nperm = nupdate-nforrest;
    long nfactorize = factor.xstore[BASICLU_NFACTORIZE];
    double time_factorize = factor.xstore[BASICLU_TIME_FACTORIZE_TOTAL];
    double time_solve = factor.xstore[BASICLU_TIME_SOLVE_TOTAL];
    double time_update = factor.xstore[BASICLU_TIME_UPDATE_TOTAL];

    printf(" status:               %s\n",
           changed ? "max # passes done" : "optimal basis");
    printf(" # passes:             %ld\n", pass);
    printf(" row rank deficiency:  %ld\n", rankdef);
    printf(" updates [perm + FT]:  %ld [%ld + %ld]\n",
           nupdate, nperm, nforrest);
    printf(" # factorizations:     %ld\n", nfactorize);
    printf(" time factorize:       %.2f sec\n", time_factorize);
    printf(" time solve:           %.2f sec\n", time_solve);
    printf(" time update:          %.2f sec\n", time_update);

cleanup:
    if (Ap) free(Ap);
    if (Ai) free(Ai);
    if (Ax) free(Ax);
    if (basis) free(basis);
    if (isbasic) free(isbasic);
    if (count) free(count);
    if (mm_I) free(mm_I);
    if (mm_J) free(mm_J);
    if (mm_val) free(mm_val);
    basiclu_obj_free(&factor);

    if (err != BASICLU_OK)
        printf(" error (%ld)\n", (long) err);
    return err;
}

/*
 * factorize() - factorize A[:,basis]
 */
static lu_int factorize(struct basiclu_object *factor, const lu_int *Ap,
                        const lu_int *Ai, const double *Ax, const lu_int *basis)
{
    double *xstore = factor->xstore;
    const lu_int m = xstore[BASICLU_DIM];
    lu_int *begin = NULL;
    lu_int *end = NULL;
    lu_int i, status = BASICLU_OK;

    begin = malloc(m*sizeof(lu_int));
    end = malloc(m*sizeof(lu_int));
    if (!begin || !end) {
        status = BASICLU_ERROR_out_of_memory;
        goto cleanup;
    }
    for (i = 0; i < m; i++)
    {
        begin[i] = Ap[basis[i]];
        end[i] = Ap[basis[i]+1];
    }

    status = basiclu_obj_factorize(factor, begin, end, Ai, Ax);

cleanup:
    if (begin) free(begin);
    if (end) free(end);

    return status;
}

/*
 * refactorize_if_needed() - refactorize the basis if required or favourable
 *
 * The basis matrix is refactorized if
 * - the maximum number of updates is reached, or
 * - the previous update had a large pivot error, or
 * - it is favourable for performance
 *
 * factorize() is called for the actual factorization.
 *
 * Note: refactorize_if_needed() will not do an initial factorization.
 */
static lu_int refactorize_if_needed(struct basiclu_object *factor,
                                    const lu_int *Ap,
                                    const lu_int *Ai,
                                    const double *Ax,
                                    const lu_int *basis)
{
    lu_int status = BASICLU_OK;
    const double piverr_tol = 1e-8;
    double *xstore = factor->xstore;

    if (xstore[BASICLU_NFORREST] == xstore[BASICLU_DIM] ||
        xstore[BASICLU_PIVOT_ERROR] > piverr_tol ||
        xstore[BASICLU_UPDATE_COST] > 1.0)
        status = factorize(factor, Ap, Ai, Ax, basis);
    return status;
}

/*
 * print_log() - print log message to stdout all 3 sec
 */
static void print_log(struct basiclu_object *factor)
{
    double *xstore = factor->xstore;
    long nfactor = xstore[BASICLU_NFACTORIZE];
    long nupdate = xstore[BASICLU_NUPDATE_TOTAL];
    double total_time = xstore[BASICLU_TIME_FACTORIZE_TOTAL] +
                        xstore[BASICLU_TIME_UPDATE_TOTAL] +
                        xstore[BASICLU_TIME_SOLVE_TOTAL];
    static double last_log = 0;

    if (total_time - last_log > 3.0) {
        printf(" %6.1fs  %6ld update  %3ld factor\n", total_time,
               nupdate, nfactor);
        last_log = total_time;
    }
}

/*
 * maxvolume() - one pass over columns of A doing basis updates
 *
 * For each column a_j not in B, compute lhs = B^{-1}*a_j and find the maximum
 * entry lhs[imax]. If it is bigger than @volumetol in absolute value, then
 * replace position imax of the basis by index j. On return *p_nupdate is the
 * number of basis updates done.
 */
static lu_int maxvolume(struct basiclu_object *factor, lu_int ncol,
                        const lu_int *Ap, const lu_int *Ai, const double *Ax,
                        lu_int *basis, lu_int *isbasic, double volumetol,
                        lu_int *p_nupdate)
{
    lu_int i, j, k;
    lu_int nzrhs, imax, begin, nupdate = 0;
    double xtbl, xmax;
    lu_int status = BASICLU_OK;

    print_log(factor);

    /* Compute initial factorization. */
    status = factorize(factor, Ap, Ai, Ax, basis);
    if (status != BASICLU_OK)
        goto cleanup;

    for (j = 0; j < ncol; j++)
    {
        if (isbasic[j])
            continue;

        print_log(factor);

        /* compute B^{-1}*a_j */
        nzrhs = Ap[j+1] - Ap[j];
        begin = Ap[j];
        status = basiclu_obj_solve_for_update(factor, nzrhs, Ai+begin, Ax+begin,
                                              'N', 1);
        if (status != BASICLU_OK)
            goto cleanup;

        /* Find the maximum entry. */
        xmax = 0.0;
        xtbl = 0.0;
        imax = 0;
        for (k = 0; k < factor->nzlhs; k++) {
            i = factor->ilhs[k];
            if (fabs(factor->lhs[i]) > xmax) {
                xtbl = factor->lhs[i];
                xmax = fabs(xtbl);
                imax = i;
            }
        }

        if (xmax <= volumetol)
            continue;

        /* Update basis. */
        isbasic[basis[imax]] = 0;
        isbasic[j] = 1;
        basis[imax] = j;
        nupdate++;

        /* Prepare to update factorization. */
        status = basiclu_obj_solve_for_update(factor, 0, &imax, NULL, 'T', 0);
        if (status != BASICLU_OK)
            goto cleanup;

        status = basiclu_obj_update(factor, xtbl);
        if (status != BASICLU_OK)
            goto cleanup;

        status = refactorize_if_needed(factor, Ap, Ai, Ax, basis);
        if (status != BASICLU_OK)
            goto cleanup;
    }

cleanup:
    *p_nupdate = nupdate;
    return status;
}
