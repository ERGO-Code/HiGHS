#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "ipx_c.h"

static void check_double_vector(const mxArray *pm, size_t dim, const char *name)
{
    if (!mxIsDouble(pm) || mxIsSparse(pm) || mxIsComplex(pm)) {
        char msg[200];
        sprintf(msg, "%s must be a dense double vector of real numbers.", name);
        mexErrMsgTxt(msg);
    }
    size_t m = mxGetM(pm);
    size_t n = mxGetN(pm);
    if ((m != 1 && n != 1) || m*n != dim) {
        char msg[200];
        sprintf(msg, "%s has wrong dimension.", name);
        mexErrMsgTxt(msg);
    }
}

static void check_char_vector(const mxArray *pm, size_t dim, const char *name)
{
    if (!mxIsChar(pm)) {
        char msg[200];
        sprintf(msg, "%s must be a char vector.", name);
        mexErrMsgTxt(msg);
    }
    size_t m = mxGetM(pm);
    size_t n = mxGetN(pm);
    if ((m != 1 && n != 1) || m*n != dim) {
        char msg[200];
        sprintf(msg, "%s has wrong dimension.", name);
        mexErrMsgTxt(msg);
    }
}

static ipxint parse_int(const mxArray *pp, int fieldnumber)
{
    const char *fieldname = mxGetFieldNameByNumber(pp, fieldnumber);
    const mxArray *pm = mxGetFieldByNumber(pp, 0, fieldnumber);
    if (!mxIsScalar(pm) || !mxIsDouble(pm) || mxIsComplex(pm) || mxIsSparse(pm))
    {
        char msg[200];
        sprintf(msg, "parameter '%s' must be a real double scalar.", fieldname);
        mexErrMsgTxt(msg);
    }
    double value = mxGetScalar(pm);
    if ((ipxint) value != value) {
        char msg[200];
        sprintf(msg, "parameter '%s' must have an integer value.", fieldname);
        mexErrMsgTxt(msg);
    }
    return value;
}

static double parse_double(const mxArray *pp, int fieldnumber)
{
    const char *fieldname = mxGetFieldNameByNumber(pp, fieldnumber);
    const mxArray *pm = mxGetFieldByNumber(pp, 0, fieldnumber);
    if (!mxIsScalar(pm) || !mxIsDouble(pm) || mxIsComplex(pm) || mxIsSparse(pm))
    {
        char msg[200];
        sprintf(msg, "parameter '%s' must be a real double scalar.", fieldname);
        mexErrMsgTxt(msg);
    }
    return mxGetScalar(pm);
}

static const char* parse_string(const mxArray *pp, int fieldnumber)
{
    const char *fieldname = mxGetFieldNameByNumber(pp, fieldnumber);
    const mxArray *pm = mxGetFieldByNumber(pp, 0, fieldnumber);
    if (!mxIsChar(pm) || mxGetM(pm) > 1) {
        char msg[200];
        sprintf(msg, "parameter '%s' must be a char array.", fieldname);
        mexErrMsgTxt(msg);
    }
    size_t n = mxGetN(pm);
    const mxChar *data = mxGetData(pm);
    char *string = mxMalloc((n+1)*sizeof(char));
    for (size_t i = 0; i < n; i++)
        string[i] = data[i];
    string[n] = 0;
    return string;
}

static double* add_vector_field(mxArray *pm, const char *name, mwSize dim)
{
    int fieldnumber = mxAddField(pm, name);
    if (fieldnumber < 0)
        mexErrMsgTxt("Creating field in output struct failed.\n");
    mxArray* value = mxCreateDoubleMatrix(dim, 1, mxREAL);
    mxSetFieldByNumber(pm, 0, fieldnumber, value);
    return mxGetData(value);
}

static void add_scalar_field(mxArray *pm, const char *name, double value)
{
    int fieldnumber = mxAddField(pm, name);
    if (fieldnumber < 0)
        mexErrMsgTxt("Creating field in output struct failed.\n");
    mxSetFieldByNumber(pm, 0, fieldnumber, mxCreateDoubleScalar(value));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1 || nrhs > 2)
        mexErrMsgTxt("Wrong number of input arguments.");
    if (nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Check type and fields of first argument. */
    const mxArray *pm = prhs[0];
    if (!mxIsStruct(pm) || !mxIsScalar(pm))
        mexErrMsgTxt("First argument must be a struct.");
    if (mxGetFieldNumber(pm, "A") < 0)
        mexErrMsgTxt("Field 'A' missing in first argument.");
    if (mxGetFieldNumber(pm, "rhs") < 0)
        mexErrMsgTxt("Field 'rhs' missing in first argument.");
    if (mxGetFieldNumber(pm, "obj") < 0)
        mexErrMsgTxt("Field 'obj' missing in first argument.");
    if (mxGetFieldNumber(pm, "lb") < 0)
        mexErrMsgTxt("Field 'lb' missing in first argument.");
    if (mxGetFieldNumber(pm, "ub") < 0)
        mexErrMsgTxt("Field 'ub' missing in first argument.");
    if (mxGetFieldNumber(pm, "sense") < 0)
        mexErrMsgTxt("Field 'sense' missing in first argument.");

    /* Check constraint matrix. */
    const mxArray *A = mxGetField(pm, 0, "A");
    if (!mxIsSparse(A) || !mxIsDouble(A) || mxIsComplex(A))
        mexErrMsgTxt("A must be a sparse matrix of real numbers.");
    size_t m = mxGetM(A);
    size_t n = mxGetN(A);
    const double *Ax = mxGetPr(A);
    const mwIndex *AIr = mxGetIr(A);
    const mwIndex *AJc = mxGetJc(A);
    mwIndex nz = AJc[n];

    /* Check vectors. */
    const mxArray *rhs = mxGetField(pm, 0, "rhs");
    const mxArray *obj = mxGetField(pm, 0, "obj");
    const mxArray *lb = mxGetField(pm, 0, "lb");
    const mxArray *ub = mxGetField(pm, 0, "ub");
    const mxArray *sense = mxGetField(pm, 0, "sense");
    check_double_vector(rhs, m, "rhs");
    check_double_vector(obj, n, "obj");
    check_double_vector(lb, n, "lb");
    check_double_vector(ub, n, "ub");
    check_char_vector(sense, m, "sense");

    /* Copy integer and char arrays to allow incompatible types. */
    ipxint *Ap = mxMalloc((n+1)*sizeof(ipxint));
    ipxint *Ai = mxMalloc(nz*sizeof(ipxint));
    char *constr_type = mxMalloc(m*sizeof(char));
    for (size_t j = 0; j <= n; j++)
        Ap[j] = AJc[j];
    for (mwIndex p = 0; p < nz; p++)
        Ai[p] = AIr[p];
    mxChar *sense_data = mxGetData(sense);
    for (size_t i = 0; i < m; i++)
        constr_type[i] = sense_data[i];

    /* Parse parameters. */
    struct ipx_parameters params = ipx_default_parameters();
    if (nrhs >= 2) {
        const mxArray *pp = prhs[1];
        if (!mxIsStruct(pp) || !mxIsScalar(pp))
            mexErrMsgTxt("Second argument must be a struct.");
        int num_params = mxGetNumberOfFields(pp);
        for (int i = 0; i < num_params; i++) {
            const char *name = mxGetFieldNameByNumber(pp, i);

            if (strcmp(name, "logfile") == 0)
                params.logfile = parse_string(pp, i);

            else if (strcmp(name, "display") == 0)
                params.display = parse_int(pp, i);
            else if (strcmp(name, "dualize") == 0)
                params.dualize = parse_int(pp, i);
            else if (strcmp(name, "scale") == 0)
                params.scale = parse_int(pp, i);
            else if (strcmp(name, "ipm_maxiter") == 0)
                params.ipm_maxiter = parse_int(pp, i);
            else if (strcmp(name, "precond_dense_cols") == 0)
                params.precond_dense_cols = parse_int(pp, i);
            else if (strcmp(name, "crash_basis") == 0)
                params.crash_basis = parse_int(pp, i);
            else if (strcmp(name, "rows_per_slice") == 0)
                params.rows_per_slice = parse_int(pp, i);
            else if (strcmp(name, "maxskip_updates") == 0)
                params.maxskip_updates = parse_int(pp, i);
            else if (strcmp(name, "lu_kernel") == 0)
                params.lu_kernel = parse_int(pp, i);
            else if (strcmp(name, "crossover") == 0)
                params.crossover = parse_int(pp, i);
            else if (strcmp(name, "debug") == 0)
                params.debug = parse_int(pp, i);
            else if (strcmp(name, "switchiter") == 0)
                params.switchiter = parse_int(pp, i);
            else if (strcmp(name, "stop_at_switch") == 0)
                params.stop_at_switch = parse_int(pp, i);
            else if (strcmp(name, "update_heuristic") == 0)
                params.update_heuristic = parse_int(pp, i);
            else if (strcmp(name, "maxpasses") == 0)
                params.maxpasses = parse_int(pp, i);

            else if (strcmp(name, "print_interval") == 0)
                params.print_interval = parse_double(pp, i);
            else if (strcmp(name, "time_limit") == 0)
                params.time_limit = parse_double(pp, i);
            else if (strcmp(name, "ipm_feasibility_tol") == 0)
                params.ipm_feasibility_tol = parse_double(pp, i);
            else if (strcmp(name, "ipm_optimality_tol") == 0)
                params.ipm_optimality_tol = parse_double(pp, i);
            else if (strcmp(name, "ipm_drop_primal") == 0)
                params.ipm_drop_primal = parse_double(pp, i);
            else if (strcmp(name, "ipm_drop_dual") == 0)
                params.ipm_drop_dual = parse_double(pp, i);
            else if (strcmp(name, "kkt_tol") == 0)
                params.kkt_tol = parse_double(pp, i);
            else if (strcmp(name, "dependency_tol") == 0)
                params.dependency_tol = parse_double(pp, i);
            else if (strcmp(name, "volume_tol") == 0)
                params.volume_tol = parse_double(pp, i);
            else if (strcmp(name, "lu_pivottol") == 0)
                params.lu_pivottol = parse_double(pp, i);
            else if (strcmp(name, "crossover_start") == 0)
                params.crossover_start = parse_double(pp, i);
            else if (strcmp(name, "pfeasibility_tol") == 0)
                params.pfeasibility_tol = parse_double(pp, i);
            else if (strcmp(name, "dfeasibility_tol") == 0)
                params.dfeasibility_tol = parse_double(pp, i);

            else
                mexPrintf("warning: unrecognized parameter name ('%s')\n",
                          name);
        }
    }

    /* Initialize IPX. */
    void *solver = NULL;
    ipx_new(&solver);
    if (!solver)
        mexErrMsgTxt("Initializing IPX solver instance failed.");
    ipx_set_parameters(solver, params);

    /* Run IPX. */
    ipxint status = ipx_solve(solver, n, mxGetPr(obj), mxGetPr(lb), mxGetPr(ub),
                              m, Ap, Ai, Ax, mxGetPr(rhs), constr_type);
    struct ipx_info info = ipx_get_info(solver);

    /* Free temporary memory before allocating solution. */
    mxFree(Ap);
    mxFree(Ai);
    mxFree(constr_type);

    /* Return solution. */
    if (nlhs >= 1) {
        plhs[0] = mxCreateStructMatrix(1, 1, 0, NULL);
        add_scalar_field(plhs[0], "status", info.status);
        add_scalar_field(plhs[0], "status_ipm", info.status_ipm);
        add_scalar_field(plhs[0], "status_crossover", info.status_crossover);
        add_scalar_field(plhs[0], "errflag", info.errflag);
        if (info.status_ipm == IPX_STATUS_optimal ||
            info.status_ipm == IPX_STATUS_imprecise) {
            double *x = add_vector_field(plhs[0], "x", n);
            double *xl = add_vector_field(plhs[0], "xl", n);
            double *xu = add_vector_field(plhs[0], "xu", n);
            double *slack = add_vector_field(plhs[0], "slack", m);
            double *y = add_vector_field(plhs[0], "y", m);
            double *zl = add_vector_field(plhs[0], "zl", n);
            double *zu = add_vector_field(plhs[0], "zu", n);
            ipx_get_interior_solution(solver, x, xl, xu, slack, y, zl, zu);
        }
        if (info.status_crossover == IPX_STATUS_optimal ||
            info.status_crossover == IPX_STATUS_imprecise) {
            double *basis_x = add_vector_field(plhs[0], "basis_x", n);
            double *basis_slack = add_vector_field(plhs[0], "basis_slack", m);
            double *basis_y = add_vector_field(plhs[0], "basis_y", m);
            double *basis_z = add_vector_field(plhs[0], "basis_z", n);
            double *cbasis = add_vector_field(plhs[0], "cbasis", m);
            double *vbasis = add_vector_field(plhs[0], "vbasis", n);
            ipxint *cbasis_temp = mxMalloc(m*sizeof(ipxint));
            ipxint *vbasis_temp = mxMalloc(n*sizeof(ipxint));
            ipx_get_basic_solution(solver, basis_x, basis_slack, basis_y,
                                   basis_z, cbasis_temp, vbasis_temp);
            for (size_t i = 0; i < m; i++)
                cbasis[i] = cbasis_temp[i];
            for (size_t j = 0; j < n; j++)
                vbasis[j] = vbasis_temp[j];
            mxFree(cbasis_temp);
            mxFree(vbasis_temp);
        }
    }

    /* Return statistics. */
    if (nlhs >= 2) {
        plhs[1] = mxCreateStructMatrix(1, 1, 0, NULL);
        add_scalar_field(plhs[1], "num_var", info.num_var);
        add_scalar_field(plhs[1], "num_constr", info.num_constr);
        add_scalar_field(plhs[1], "num_entries", info.num_entries);
        add_scalar_field(plhs[1], "num_rows_solver", info.num_rows_solver);
        add_scalar_field(plhs[1], "num_cols_solver", info.num_cols_solver);
        add_scalar_field(plhs[1], "num_entries_solver",
                         info.num_entries_solver);
        add_scalar_field(plhs[1], "dualized", info.dualized);
        add_scalar_field(plhs[1], "dense_cols", info.dense_cols);
        add_scalar_field(plhs[1], "dependent_rows", info.dependent_rows);
        add_scalar_field(plhs[1], "dependent_cols", info.dependent_cols);
        add_scalar_field(plhs[1], "rows_inconsistent", info.rows_inconsistent);
        add_scalar_field(plhs[1], "cols_inconsistent", info.cols_inconsistent);
        add_scalar_field(plhs[1], "primal_dropped", info.primal_dropped);
        add_scalar_field(plhs[1], "dual_dropped", info.dual_dropped);
        add_scalar_field(plhs[1], "abs_presidual", info.abs_presidual);
        add_scalar_field(plhs[1], "abs_dresidual", info.abs_dresidual);
        add_scalar_field(plhs[1], "rel_presidual", info.rel_presidual);
        add_scalar_field(plhs[1], "rel_dresidual", info.rel_dresidual);
        add_scalar_field(plhs[1], "pobjval", info.pobjval);
        add_scalar_field(plhs[1], "dobjval", info.dobjval);
        add_scalar_field(plhs[1], "rel_objgap", info.rel_objgap);
        add_scalar_field(plhs[1], "complementarity", info.complementarity);
        add_scalar_field(plhs[1], "normx", info.normx);
        add_scalar_field(plhs[1], "normy", info.normy);
        add_scalar_field(plhs[1], "normz", info.normz);
        add_scalar_field(plhs[1], "objval", info.objval);
        add_scalar_field(plhs[1], "primal_infeas", info.primal_infeas);
        add_scalar_field(plhs[1], "dual_infeas", info.dual_infeas);
        add_scalar_field(plhs[1], "iter", info.iter);
        add_scalar_field(plhs[1], "kktiter1", info.kktiter1);
        add_scalar_field(plhs[1], "kktiter2", info.kktiter2);
        add_scalar_field(plhs[1], "basis_repairs", info.basis_repairs);
        add_scalar_field(plhs[1], "updates_start", info.updates_start);
        add_scalar_field(plhs[1], "updates_ipm", info.updates_ipm);
        add_scalar_field(plhs[1], "updates_crossover", info.updates_crossover);
        add_scalar_field(plhs[1], "time_total", info.time_total);
        add_scalar_field(plhs[1], "time_ipm1", info.time_ipm1);
        add_scalar_field(plhs[1], "time_ipm2", info.time_ipm2);
        add_scalar_field(plhs[1], "time_starting_basis",
                         info.time_starting_basis);
        add_scalar_field(plhs[1], "time_crossover", info.time_crossover);
    }

    /* Clean up. */
    ipx_free(&solver);
}
