using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;

// mcs -out:highscslib.dll -t:library highs_csharp_api.cs -unsafe

namespace Highs {
public enum HighsStatus
{
    kError = -1,
    kOk,
    kWarning
}

public enum HighsMatrixFormat
{
    kColwise = 1,
    kRowwise
}

public enum HessianFormat
{
    kTriangular = 1,
    kSquare
}

public enum HighsBasisStatus
{
    kLower = 0,
    kBasic,
    kUpper,
    kZero,
    kNonbasic
}

public enum HighsObjectiveSense
{
    kMinimize = 1,
    kMaximize = -1
}

public enum HighsModelStatus
{
    kNotset = 0,
    kLoadError,
    kModelError,
    kPresolveError,
    kSolveError,
    kPostsolveError,
    kModelEmpty,
    kOptimal,
    kInfeasible,
    kUnboundedOrInfeasible,
    kUnbounded,
    kObjectiveBound,
    kObjectiveTarget,
    kTimeLimit,
    kIterationLimit,
    kUnknown,
    kSolutionLimit,
    kInterrupt,
    kMemoryLimit
}

public enum HighsIntegrality
{
    kContinuous = 0,
    kInteger = 1,
    kSemiContinuous = 2,
    kSemiInteger = 3,
    kImplicitInteger = 4,
}

/// <summary>A category of log message</summary>
public enum HighsLogType
{
    Info = 1,
    Detailed,
    Verbose,
    Warning,
    Error
}

/// <summary>A category of callback</summary>
internal enum HighsCallbackType
{
    Logging = 0,
    SimplexInterrupt = 1,
    IpmInterrupt = 2,
    MipSolution = 3,
    MipImprovingSolution = 4,
    MipLogging = 5,
    MipInterrupt = 6,
    MipGetCutPool = 7,
    MipDefineLazyConstraints = 8,
}


public class HighsModel
{
    public HighsObjectiveSense sense;
    public double[] colcost;
    public double offset;
    public double[] collower;
    public double[] colupper;
    public double[] rowlower;
    public double[] rowupper;
    public HighsMatrixFormat a_format;
    public int[] astart;
    public int[] aindex;
    public double[] avalue;
    public int[] highs_integrality;

    public HighsModel()
    {

    }

    public HighsModel(double[] colcost, double[] collower, double[] colupper, double[] rowlower, double[] rowupper,
    int[] astart, int[] aindex, double[] avalue, int[] highs_integrality = null, double offset = 0, HighsMatrixFormat a_format = HighsMatrixFormat.kColwise, HighsObjectiveSense sense = HighsObjectiveSense.kMinimize)
    {
        this.colcost = colcost;
        this.collower = collower;
        this.colupper = colupper;
        this.rowlower = rowlower;
        this.rowupper = rowupper;
        this.astart = astart;
        this.aindex = aindex;
        this.avalue = avalue;
        this.offset = offset;
        this.a_format = a_format;
        this.sense = sense;
        this.highs_integrality = highs_integrality;
    }
}

public class HighsHessian
{
    public HessianFormat q_format;
    public int dim;
    public int[] qstart;
    public int[] qindex;
    public double[] qvalue;

    public HighsHessian()
    {

    }

    public HighsHessian(int dim, int[] qstart, int[] qindex, double[] qvalue, HessianFormat q_format = HessianFormat.kTriangular)
    {
        this.dim = dim;
        this.qstart = qstart;
        this.qindex = qindex;
        this.qvalue = qvalue;
        this.q_format = q_format;
    }
}

public class HighsSolution
{
    public double[] colvalue;
    public double[] coldual;
    public double[] rowvalue;
    public double[] rowdual;

    public HighsSolution(int numcol, int numrow)
    {
        this.colvalue = new double[numcol];
        this.coldual = new double[numcol];
        this.rowvalue = new double[numrow];
        this.rowdual = new double[numrow];
    }

    public HighsSolution(double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual)
    {
        this.colvalue = colvalue;
        this.coldual = coldual;
        this.rowvalue = rowvalue;
        this.rowdual = rowdual;
    }
}

public class HighsBasis
{
    public HighsBasisStatus[] colbasisstatus;
    public HighsBasisStatus[] rowbasisstatus;

    public HighsBasis(int numcol, int numrow)
    {
        this.colbasisstatus = new HighsBasisStatus[numcol];
        this.rowbasisstatus = new HighsBasisStatus[numrow];
    }

    public HighsBasis(HighsBasisStatus[] colbasisstatus, HighsBasisStatus[] rowbasisstatus)
    {
        this.colbasisstatus = colbasisstatus;
        this.rowbasisstatus = rowbasisstatus;
    }
}

public class HighsLpSolver : IDisposable
{
    private IntPtr highs;

    /// <summary>Read-only access to Highs instance pointer</summary>
    /// <remarks>Allows sub-classes to do meaningful things using the HiGHS C API.</remarks>
    protected IntPtr HighsObject => this.highs;

    private bool _disposed;

    private const string highslibname = "highs";

    /// <summary>Signature of functions that can be called by HiGHS when callback events occur</summary>
    private delegate void CallbackDelegate(
        HighsCallbackType cbType, IntPtr messagePtr, [In] ref HighsCallbackDataOut cbDataOut,
        ref HighsCallbackDataIn cbDataIn, IntPtr cbUserData);

    /// <summary>Pointer to function that is called when HiGHS callbacks occur</summary>
    private CallbackDelegate _cbDelegate;

    /// <summary>A C function pointer to the event-generating callback delegate</summary>
    /// <remarks>This property's primary purpose is to improve ability to test callback-triggered events</remarks>
    protected IntPtr CallbackFunctionPtr => Marshal.GetFunctionPointerForDelegate(_cbDelegate);

    [DllImport(highslibname)]
    private static extern int Highs_call(
        Int32 numcol,
        Int32 numrow,
        Int32 numnz,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        double[] colvalue,
        double[] coldual,
        double[] rowvalue,
        double[] rowdual,
        int[] colbasisstatus,
        int[] rowbasisstatus,
        ref int modelstatus);

    [DllImport(highslibname)]
    private static extern IntPtr Highs_create();

    [DllImport(highslibname)]
    private static extern void Highs_destroy(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_run(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_readModel(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_writeModel(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_writePresolvedModel(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_writeSolutionPretty(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_getInfinity(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_passLp(
        IntPtr highs,
        int numcol,
        int numrow,
        int numnz,
        int aformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue);

    [DllImport(highslibname)]
    private static extern int Highs_passMip(
        IntPtr highs,
        int numcol,
        int numrow,
        int numnz,
        int aformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        int[] highs_integrality);

    [DllImport(highslibname)]
    private static extern int Highs_passModel(
        IntPtr highs,
        int numcol,
        int numrow,
        int numnz,
        int qnumnz,
        int aformat,
        int qformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        int[] qstart,
        int[] qindex,
        double[] qvalue,
        int[] highs_integrality);

    [DllImport(highslibname)]
    private static extern int Highs_passHessian(
        IntPtr highs,
        int dim,
        int numnz,
        int q_format,
        int[] qstart,
        int[] qindex,
        double[] qvalue);

    [DllImport(highslibname)]
    private static extern int Highs_setOptionValue(IntPtr highs, string option, string value);

    [DllImport(highslibname)]
    private static extern int Highs_setBoolOptionValue(IntPtr highs, string option, int value);

    [DllImport(highslibname)]
    private static extern int Highs_setIntOptionValue(IntPtr highs, string option, int value);

    [DllImport(highslibname)]
    private static extern int Highs_setDoubleOptionValue(IntPtr highs, string option, double value);

    [DllImport(highslibname)]
    private static extern int Highs_setStringOptionValue(IntPtr highs, string option, string value);

    [DllImport(highslibname)]
    private static extern int Highs_getBoolOptionValue(IntPtr highs, string option, out int value);

    [DllImport(highslibname)]
    private static extern int Highs_getIntOptionValue(IntPtr highs, string option, out int value);

    [DllImport(highslibname)]
    private static extern int Highs_getDoubleOptionValue(IntPtr highs, string option, out double value);

    [DllImport(highslibname)]
    private static extern int Highs_getStringOptionValue(IntPtr highs, string option, [Out] StringBuilder value);

    [DllImport(highslibname)]
    private static extern int Highs_getSolution(IntPtr highs, double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual);

    [DllImport(highslibname)]
    private static extern int Highs_getNumCol(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getNumRow(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getNumNz(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getHessianNumNz(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getBasis(IntPtr highs, int[] colstatus, int[] rowstatus);

    [DllImport(highslibname)]
    private static extern double Highs_getObjectiveValue(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getIterationCount(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_getModelStatus(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_addRow(IntPtr highs, double lower, double upper, int num_new_nz, int[] indices, double[] values);

    [DllImport(highslibname)]
    private static extern int Highs_addRows(
        IntPtr highs,
        int num_new_row,
        double[] lower,
        double[] upper,
        int num_new_nz,
        int[] starts,
        int[] indices,
        double[] values);

    [DllImport(highslibname)]
    private static extern int Highs_addCol(
        IntPtr highs,
        double cost,
        double lower,
        double upper,
        int num_new_nz,
        int[] indices,
        double[] values);

    [DllImport(highslibname)]
    private static extern int Highs_addCols(
        IntPtr highs,
        int num_new_col,
        double[] costs,
        double[] lower,
        double[] upper,
        int num_new_nz,
        int[] starts,
        int[] indices,
        double[] values);

    [DllImport(highslibname)]
    private static extern int Highs_changeObjectiveSense(IntPtr highs, int sense);

    [DllImport(highslibname)]
    private static extern int Highs_changeColCost(IntPtr highs, int col, double cost);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsCostBySet(IntPtr highs, int num_set_entries, int[] set, double[] cost);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsCostByMask(IntPtr highs, int[] mask, double[] cost);

    [DllImport(highslibname)]
    private static extern int Highs_changeColBounds(IntPtr highs, int col, double lower, double upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsBoundsByRange(IntPtr highs, int from_col, int to_col, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsBoundsBySet(IntPtr highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsBoundsByMask(IntPtr highs, int[] mask, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeRowBounds(IntPtr highs, int row, double lower, double upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeRowsBoundsByRange(IntPtr highs, int from_row, int to_row, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeRowsBoundsBySet(IntPtr highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeRowsBoundsByMask(IntPtr highs, int[] mask, double[] lower, double[] upper);

    [DllImport(highslibname)]
    private static extern int Highs_changeColsIntegralityByRange(IntPtr highs, int from_col, int to_col, int[] integrality);

    [DllImport(highslibname)]
    private static extern int Highs_changeCoeff(IntPtr highs, int row, int col, double value);

    [DllImport(highslibname)]
    private static extern int Highs_deleteColsByRange(IntPtr highs, int from_col, int to_col);

    [DllImport(highslibname)]
    private static extern int Highs_deleteColsBySet(IntPtr highs, int num_set_entries, int[] set);

    [DllImport(highslibname)]
    private static extern int Highs_deleteColsByMask(IntPtr highs, int[] mask);

    [DllImport(highslibname)]
    private static extern int Highs_deleteRowsByRange(IntPtr highs, int from_row, int to_row);

    [DllImport(highslibname)]
    private static extern int Highs_deleteRowsBySet(IntPtr highs, int num_set_entries, int[] set);

    [DllImport(highslibname)]
    private static extern int Highs_deleteRowsByMask(IntPtr highs, int[] mask);

    [DllImport(highslibname)]
    private static extern int Highs_getDoubleInfoValue(IntPtr highs, string info, out double value);

    [DllImport(highslibname)]
    private static extern int Highs_getIntInfoValue(IntPtr highs, string info, out int value);

    [DllImport(highslibname)]
    private static extern int Highs_getInt64InfoValue(IntPtr highs, string info, out long value);

    [DllImport(highslibname)]
    private static extern int Highs_setSolution(IntPtr highs, double[] col_value, double[] row_value, double[] col_dual, double[] row_dual);

    [DllImport(highslibname)]
    private static extern int Highs_setSparseSolution(IntPtr highs, int num_entries, int[] index, double[] value);

    [DllImport(highslibname)]
    private static extern int Highs_getColsByRange(
        IntPtr highs,
        int from_col,
        int to_col,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getColsBySet(
        IntPtr highs,
        int num_set_entries,
        int[] set,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getColsByMask(
        IntPtr highs,
        int[] mask,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getRowsByRange(
        IntPtr highs,
        int from_row,
        int to_row,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getRowsBySet(
        IntPtr highs,
        int num_set_entries,
        int[] set,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getRowsByMask(
        IntPtr highs,
        int[] mask,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [DllImport(highslibname)]
    private static extern int Highs_getBasicVariables(IntPtr highs, int[] basic_variables);

    [DllImport(highslibname)]
    private static extern int Highs_getBasisInverseRow(IntPtr highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getBasisInverseCol(IntPtr highs, int col, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getBasisSolve(
        IntPtr highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_num_nz,
        int[] solution_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getBasisTransposeSolve(
        IntPtr highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_nz,
        int[] solution_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getReducedRow(IntPtr highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [DllImport(highslibname)]
    private static extern int Highs_getReducedColumn(IntPtr highs, int col, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [DllImport(highslibname)]
    private static extern int Highs_clearModel(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_clearSolver(IntPtr highs);

    [DllImport(highslibname)]
    private static extern int Highs_passColName(IntPtr highs, int col, string name);

    [DllImport(highslibname)]
    private static extern int Highs_passRowName(IntPtr highs, int row, string name);

    [DllImport(highslibname)]
    private static extern int Highs_writeOptions(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_writeOptionsDeviations(IntPtr highs, string filename);

    [DllImport(highslibname)]
    private static extern int Highs_setCallback(IntPtr highs, IntPtr cbFuncPtr, IntPtr cbUserData);

    [DllImport(highslibname)]
    private static extern int Highs_startCallback(IntPtr highs, HighsCallbackType cbType);

    [DllImport(highslibname)]
    private static extern int Highs_stopCallback(IntPtr highs, HighsCallbackType cbType);

    public static HighsStatus call(HighsModel model, ref HighsSolution sol, ref HighsBasis bas, ref HighsModelStatus modelstatus)
    {
        int nc = model.colcost.Length;
        int nr = model.rowlower.Length;
        int nnz = model.avalue.Length;

        int[] colbasstat = new int[nc];
        int[] rowbasstat = new int[nr];

        int modelstate = 0;

        HighsStatus status = (HighsStatus)HighsLpSolver.Highs_call(
            nc,
            nr,
            nnz,
            model.colcost,
            model.collower,
            model.colupper,
            model.rowlower,
            model.rowupper,
            model.astart,
            model.aindex,
            model.avalue,
            sol.colvalue,
            sol.coldual,
            sol.rowvalue,
            sol.rowdual,
            colbasstat,
            rowbasstat,
            ref modelstate);

        modelstatus = (HighsModelStatus)modelstate;

        bas.colbasisstatus = colbasstat.Select(x => (HighsBasisStatus)x).ToArray();
        bas.rowbasisstatus = rowbasstat.Select(x => (HighsBasisStatus)x).ToArray();

        return status;
    }

    public HighsLpSolver()
    {
        this.highs = HighsLpSolver.Highs_create();
        _cbDelegate = this.callbackFunction;
        Highs_setCallback(this.highs, Marshal.GetFunctionPointerForDelegate(_cbDelegate), IntPtr.Zero);
    }

    ~HighsLpSolver()
    {
        this.Dispose(false);
    }

    public void Dispose()
    {
        this.Dispose(true);
        GC.SuppressFinalize(this);
    }

    protected virtual void Dispose(bool disposing)
    {
        if (this._disposed)
        {
            return;
        }

        HighsLpSolver.Highs_destroy(this.highs);
        this._cbDelegate = null;
        this._disposed = true;
    }

    public HighsStatus run()
    {
        return (HighsStatus)HighsLpSolver.Highs_run(this.highs);
    }

    public HighsStatus readModel(string filename)
    {
        return (HighsStatus)HighsLpSolver.Highs_readModel(this.highs, filename);
    }

    public HighsStatus writeModel(string filename)
    {
        return (HighsStatus)HighsLpSolver.Highs_writeModel(this.highs, filename);
    }

    public HighsStatus writePresolvedModel(string filename)
    {
        return (HighsStatus)HighsLpSolver.Highs_writePresolvedModel(this.highs, filename);
    }

    public HighsStatus writeSolutionPretty(string filename)
    {
        return (HighsStatus)HighsLpSolver.Highs_writeSolutionPretty(this.highs, filename);
    }

    public Double getInfinity()
    {
        return (Double)HighsLpSolver.Highs_getInfinity(this.highs);
    }

    public HighsStatus passLp(HighsModel model)
    {
        return (HighsStatus)HighsLpSolver.Highs_passLp(
            this.highs,
            model.colcost.Length,
            model.rowlower.Length,
            model.avalue.Length,
            (int)model.a_format,
            (int)model.sense,
            model.offset,
            model.colcost,
            model.collower,
            model.colupper,
            model.rowlower,
            model.rowupper,
            model.astart,
            model.aindex,
            model.avalue);
    }

    public HighsStatus passMip(HighsModel model)
    {
        return (HighsStatus)HighsLpSolver.Highs_passMip(
            this.highs,
            model.colcost.Length,
            model.rowlower.Length,
            model.avalue.Length,
            (int)model.a_format,
            (int)model.sense,
            model.offset,
            model.colcost,
            model.collower,
            model.colupper,
            model.rowlower,
            model.rowupper,
            model.astart,
            model.aindex,
            model.avalue,
            model.highs_integrality);
    }

    public HighsStatus passHessian(HighsHessian hessian)
    {
        return (HighsStatus)HighsLpSolver.Highs_passHessian(
            this.highs,
            hessian.dim,
            hessian.qvalue.Length,
            (int)hessian.q_format,
            hessian.qstart,
            hessian.qindex,
            hessian.qvalue);
    }

    public HighsStatus setOptionValue(string option, string value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setOptionValue(this.highs, option, value);
    }

    public HighsStatus setStringOptionValue(string option, string value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setStringOptionValue(this.highs, option, value);
    }

    public HighsStatus setBoolOptionValue(string option, int value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setBoolOptionValue(this.highs, option, value);
    }

    public HighsStatus setDoubleOptionValue(string option, double value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setDoubleOptionValue(this.highs, option, value);
    }

    public HighsStatus setIntOptionValue(string option, int value)
    {
        return (HighsStatus)HighsLpSolver.Highs_setIntOptionValue(this.highs, option, value);
    }

    public HighsStatus getStringOptionValue(string option, out string value)
    {
        var stringBuilder = new StringBuilder();
        var result = (HighsStatus)HighsLpSolver.Highs_getStringOptionValue(this.highs, option, stringBuilder);
        value = stringBuilder.ToString();
        return result;
    }

    public HighsStatus getBoolOptionValue(string option, out int value)
    {
        return (HighsStatus)HighsLpSolver.Highs_getBoolOptionValue(this.highs, option, out value);
    }

    public HighsStatus getDoubleOptionValue(string option, out double value)
    {
        return (HighsStatus)HighsLpSolver.Highs_getDoubleOptionValue(this.highs, option, out value);
    }

    public HighsStatus getIntOptionValue(string option, out int value)
    {
        return (HighsStatus)HighsLpSolver.Highs_getIntOptionValue(this.highs, option, out value);
    }

    public int getNumCol()
    {
        return HighsLpSolver.Highs_getNumCol(this.highs);
    }

    public int getNumRow()
    {
        return HighsLpSolver.Highs_getNumRow(this.highs);
    }

    public int getNumNz()
    {
        return HighsLpSolver.Highs_getNumNz(this.highs);
    }

    public HighsSolution getSolution()
    {
        int nc = this.getNumCol();
        int nr = this.getNumRow();

        HighsSolution sol = new HighsSolution(nc, nr);
        HighsLpSolver.Highs_getSolution(this.highs, sol.colvalue, sol.coldual, sol.rowvalue, sol.rowdual);

        return sol;
    }

    public HighsBasis getBasis()
    {
        int nc = this.getNumCol();
        int nr = this.getNumRow();

        int[] colbasstat = new int[nc];
        int[] rowbasstat = new int[nr];

        HighsLpSolver.Highs_getBasis(this.highs, colbasstat, rowbasstat);
        HighsBasis bas = new HighsBasis(
            colbasstat.Select(x => (HighsBasisStatus)x).ToArray(),
            rowbasstat.Select(x => (HighsBasisStatus)x).ToArray());

        return bas;
    }

    public double getObjectiveValue()
    {
        return HighsLpSolver.Highs_getObjectiveValue(this.highs);
    }

    public HighsModelStatus GetModelStatus()
    {
        return (HighsModelStatus)HighsLpSolver.Highs_getModelStatus(this.highs);
    }

    public int getIterationCount()
    {
        return HighsLpSolver.Highs_getIterationCount(this.highs);
    }

    public HighsStatus addRow(double lower, double upper, int[] indices, double[] values)
    {
        return (HighsStatus)HighsLpSolver.Highs_addRow(this.highs, lower, upper, indices.Length, indices, values);
    }

    public HighsStatus addRows(double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
    {
        return (HighsStatus)HighsLpSolver.Highs_addRows(this.highs, lower.Length, lower, upper, indices.Length, starts, indices, values);
    }

    public HighsStatus addCol(double cost, double lower, double upper, int[] indices, double[] values)
    {
        return (HighsStatus)HighsLpSolver.Highs_addCol(this.highs, cost, lower, upper, indices.Length, indices, values);
    }

    public HighsStatus addCols(double[] costs, double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
    {
        return (HighsStatus)HighsLpSolver.Highs_addCols(
            this.highs,
            costs.Length,
            costs,
            lower,
            upper,
            indices.Length,
            starts,
            indices,
            values);
    }

    public HighsStatus changeObjectiveSense(HighsObjectiveSense sense)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeObjectiveSense(this.highs, (int)sense);
    }

    public HighsStatus changeColCost(int col, double cost)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeColCost(this.highs, col, cost);
    }

    public HighsStatus changeColsCostBySet(int[] cols, double[] costs)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeColsCostBySet(this.highs, cols.Length, cols, costs);
    }

    public HighsStatus changeColsCostByMask(bool[] mask, double[] cost)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeColsCostByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), cost);
    }

    public HighsStatus changeColBounds(int col, double lower, double upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeColBounds(this.highs, col, lower, upper);
    }

    public HighsStatus changeColsBoundsByRange(int from, int to, double[] lower, double[] upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeColsBoundsByRange(this.highs, from, to, lower, upper);
    }

    public HighsStatus changeColsBoundsBySet(int[] cols, double[] lower, double[] upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeColsBoundsBySet(this.highs, cols.Length, cols, lower, upper);
    }

    public HighsStatus changeColsBoundsByMask(bool[] mask, double[] lower, double[] upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeColsBoundsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), lower, upper);
    }

    public HighsStatus changeRowBounds(int row, double lower, double upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeRowBounds(this.highs, row, lower, upper);
    }

    public HighsStatus changeRowsBoundsByRange(int from, int to, double[] lower, double[] upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeRowsBoundsByRange(this.highs, from, to, lower, upper);
    }

    public HighsStatus changeRowsBoundsBySet(int[] rows, double[] lower, double[] upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeRowsBoundsBySet(this.highs, rows.Length, rows, lower, upper);
    }

    public HighsStatus changeRowsBoundsByMask(bool[] mask, double[] lower, double[] upper)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeRowsBoundsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), lower, upper);
    }

    public HighsStatus changeColsIntegralityByRange(int from_col, int to_col, HighsIntegrality[] integrality)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeColsIntegralityByRange(this.highs, from_col, to_col, Array.ConvertAll(integrality, item => (int)item));
    }

    public HighsStatus changeCoeff(int row, int col, double value)
    {
        return (HighsStatus)HighsLpSolver.Highs_changeCoeff(this.highs, row, col, value);
    }

    public HighsStatus deleteColsByRange(int from, int to)
    {
        return (HighsStatus)HighsLpSolver.Highs_deleteColsByRange(this.highs, from, to);
    }

    public HighsStatus deleteColsBySet(int[] cols)
    {
        return (HighsStatus)HighsLpSolver.Highs_deleteColsBySet(this.highs, cols.Length, cols);
    }

    public HighsStatus deleteColsByMask(bool[] mask)
    {
        return (HighsStatus)HighsLpSolver.Highs_deleteColsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray());
    }

    public HighsStatus deleteRowsByRange(int from, int to)
    {
        return (HighsStatus)HighsLpSolver.Highs_deleteRowsByRange(this.highs, from, to);
    }

    public HighsStatus deleteRowsBySet(int[] rows)
    {
        return (HighsStatus)HighsLpSolver.Highs_deleteRowsBySet(this.highs, rows.Length, rows);
    }

    public HighsStatus deleteRowsByMask(bool[] mask)
    {
        return (HighsStatus)HighsLpSolver.Highs_deleteRowsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray());
    }

    delegate int HighsGetInfoDelegate<TValue>(IntPtr highs, string infoName, out TValue output);

    private TValue GetValueOrFallback<TValue>(HighsGetInfoDelegate<TValue> highsGetInfoDelegate, string infoName, TValue fallback)
    {
        try
        {
            var status = (HighsStatus)highsGetInfoDelegate(this.highs, infoName, out var value);
            if (status != HighsStatus.kOk)
            {
                return fallback;
            }

            return value;
        }
        catch
        {
            return fallback;
        }
    }

    /// <summary>
    /// Gets the current solution info.
    /// </summary>
    /// <returns>The <see cref="SolutionInfo"/>.</returns>
    public SolutionInfo getInfo()
    {
        // TODO: This object does not contian the "complete" info from the C api. Add further props, if you need them.
        var info = new SolutionInfo()
        {
            MipGap = this.GetValueOrFallback(HighsLpSolver.Highs_getDoubleInfoValue, "mip_gap", double.NaN),
            DualBound = this.GetValueOrFallback(HighsLpSolver.Highs_getDoubleInfoValue, "mip_dual_bound", double.NaN),
            ObjectiveValue = this.GetValueOrFallback(HighsLpSolver.Highs_getDoubleInfoValue, "objective_function_value", double.NaN),
            NodeCount = this.GetValueOrFallback(HighsLpSolver.Highs_getInt64InfoValue, "mip_node_count", 0L),
            IpmIterationCount = this.GetValueOrFallback(HighsLpSolver.Highs_getIntInfoValue, "ipm_iteration_count", 0),
            SimplexIterationCount = this.GetValueOrFallback(HighsLpSolver.Highs_getIntInfoValue, "simplex_iteration_count", 0),
            PdlpIterationCount = this.GetValueOrFallback(HighsLpSolver.Highs_getIntInfoValue, "pdlp_iteration_count", 0),
        };
        return info;
    }

    public HighsStatus setSolution(HighsSolution solution)
    {
        return (HighsStatus)HighsLpSolver.Highs_setSolution(this.highs, solution.colvalue, solution.rowvalue, solution.coldual, solution.rowdual);
    }

    /// <summary>Set a partial primal solution by passing values for a set of variables</summary>
    /// <param name="valuesByIndex">A dictionary that maps variable indices to variable values</param>
    /// <remarks>The sparse solution set by this function has values for a subset of the model's variables.
    /// For each entry in <paramref name="valuesByIndex"/>, the key identifies a variable by index, and
    /// the value indicates the variable's value in the sparse solution.</remarks>
    /// <returns>A <see cref="HighsStatus"/> constant indicating whether the call succeeded</returns>
    public HighsStatus setSparseSolution(IReadOnlyDictionary<int, double> valuesByIndex)
    {
        return (HighsStatus)Highs_setSparseSolution(this.highs, valuesByIndex.Count, valuesByIndex.Keys.ToArray(), valuesByIndex.Values.ToArray());
    }

    public HighsStatus getBasicVariables(ref int[] basic_variables)
    {
        return (HighsStatus)Highs_getBasicVariables(this.highs, basic_variables);
    }

    public HighsStatus getBasisInverseRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
    {
        return (HighsStatus)Highs_getBasisInverseRow(this.highs, row, row_vector, ref row_num_nz, row_indices);
    }

    public HighsStatus getBasisInverseCol(int col, double[] col_vector, ref int col_num_nz, int[] col_indices)
    {
        return (HighsStatus)Highs_getBasisInverseCol(this.highs, col, col_vector, ref col_num_nz, col_indices);
    }

    public HighsStatus getBasisSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
    {
        return (HighsStatus)Highs_getBasisSolve(this.highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
    }

    public HighsStatus getBasisTransposeSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
    {
        return (HighsStatus)Highs_getBasisTransposeSolve(this.highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
    }

    public HighsStatus getReducedRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
    {
        return (HighsStatus)Highs_getReducedRow(this.highs, row, row_vector, ref row_num_nz, row_indices);
    }

    public HighsStatus getReducedColumn(int col, double[] col_vector, ref int col_num_nz, int[] col_indices)
    {
        return (HighsStatus)Highs_getReducedColumn(this.highs, col, col_vector, ref col_num_nz, col_indices);
    }

    public HighsStatus clearModel()
    {
        return (HighsStatus)Highs_clearModel(this.highs);
    }

    public HighsStatus clearSolver()
    {
        return (HighsStatus)Highs_clearSolver(this.highs);
    }

    public HighsStatus passColName(int col, string name)
    {
        return (HighsStatus)Highs_passColName(this.highs, col, name);
    }

    public HighsStatus passRowName(int row, string name)
    {
        return (HighsStatus)Highs_passRowName(this.highs, row, name);
    }

    public HighsStatus writeOptions(string filename)
    {
        return (HighsStatus)Highs_writeOptions(this.highs, filename);
    }

    public HighsStatus writeOptionsDeviations(string filename)
    {
        return (HighsStatus)Highs_writeOptionsDeviations(this.highs, filename);
    }

#region "Callbacks as events"
    private HighsStatus startCallback(HighsCallbackType cbType)
    {
        return (HighsStatus)Highs_startCallback(this.highs, cbType);
    }

    private HighsStatus stopCallback(HighsCallbackType cbType)
    {
        return (HighsStatus)Highs_stopCallback(this.highs, cbType);
    }

    private void callbackFunction(HighsCallbackType cbType, IntPtr messagePtr,
                                  [In] ref HighsCallbackDataOut cbDataOut,
                                  [In, Out] ref HighsCallbackDataIn cbDataIn, IntPtr cbUserData)
    {
        switch (cbType)
        {
            case HighsCallbackType.Logging:
                // We receive the message as an IntPtr instead of a string so that the marshaller
                // doesn't attempt to free the C string.
                string message = Marshal.PtrToStringAnsi(messagePtr);

                var loggingEventData = new LoggingEventArgs(cbDataOut.log_type, message);
                _innerLogReceived?.Invoke(this, loggingEventData);
                break;

            case HighsCallbackType.MipImprovingSolution:
                var mipImpArgs = new MipEventArgs(cbDataOut);
                _innerMipImproving?.Invoke(this, mipImpArgs);
                break;

            case HighsCallbackType.MipLogging:
                var mipLogArgs = new MipEventArgs(cbDataOut);
                _innerMipLogging?.Invoke(this, mipLogArgs);
                break;

            case HighsCallbackType.MipInterrupt:
            case HighsCallbackType.IpmInterrupt:
            case HighsCallbackType.SimplexInterrupt:
                var interruptArgs = new InterruptCheckEventArgs(cbDataOut);
                var evnt = cbType == HighsCallbackType.MipInterrupt ? _innerMipInterrupt :
                           cbType == HighsCallbackType.SimplexInterrupt ? _innerSimplexInterrupt :
                           cbType == HighsCallbackType.IpmInterrupt ? _innerIpmInterrupt :
                           null;
                evnt?.Invoke(this, interruptArgs);
                if (interruptArgs.InterruptSolver)
                    cbDataIn.user_interrupt = 1;
                break;

            default:
                break;
        }
    }

    // Expose callbacks as .NET events.
    // Event declarations use custom add/remove accessors to automatically start and stop
    // the relevant callbacks when the number of listeners moves between 0 and 1.

    // kCallbackLogging as an event
    private readonly object _logReceivedLockObject = new object();
    private EventHandler<LoggingEventArgs> _innerLogReceived;
    /// <summary>Occurs when a log message is generated by HiGHS</summary>
    public event EventHandler<LoggingEventArgs> LogMessageReceived
    {
        add
        {
            lock (_logReceivedLockObject)
            {
                // If this is the first subscription to the event, start the callback
                if (_innerLogReceived == null)
                    this.startCallback(HighsCallbackType.Logging);
                _innerLogReceived += value;
            }
        }
        remove
        {
            lock (_logReceivedLockObject)
            {
                _innerLogReceived -= value;
                // If this was the last subscription to the event, stop the callback
                if (_innerLogReceived == null)
                    this.stopCallback(HighsCallbackType.Logging);
            }
        }
    }

    // kCallbackMipImprovingSolution as an event
    private readonly object _mipImprovingLockObject = new object();
    private EventHandler<MipEventArgs> _innerMipImproving;
    /// <summary>Occurs when the MIP solver identifies an improving integer feasible solution</summary>
    public event EventHandler<MipEventArgs> MipImprovingSolutionFound
    {
        add
        {
            lock (_mipImprovingLockObject)
            {
                if (_innerMipImproving == null)
                    this.startCallback(HighsCallbackType.MipImprovingSolution);
                _innerMipImproving += value;
            }
        }
        remove
        {
            lock (_mipImprovingLockObject)
            {
                _innerMipImproving -= value;
                if (_innerMipImproving == null)
                    this.stopCallback(HighsCallbackType.MipImprovingSolution);
            }
        }
    }

    // kCallbackMipLogging as an event
    private readonly object _mipLoggingLockObject = new object();
    private EventHandler<MipEventArgs> _innerMipLogging;
    /// <summary>Occurs when the MIP solver receives a MIP status report</summary>
    public event EventHandler<MipEventArgs> MipStatusReported
    {
        add
        {
            lock (_mipLoggingLockObject)
            {
                if (_innerMipLogging == null)
                    this.startCallback(HighsCallbackType.MipLogging);
                _innerMipLogging += value;
            }
        }
        remove
        {
            lock (_mipLoggingLockObject)
            {
                _innerMipLogging -= value;
                if (_innerMipLogging == null)
                    this.stopCallback(HighsCallbackType.MipLogging);
            }
        }
    }

    // kCallbackMipInterrupt as an event
    private readonly object _mipInterruptLockObject = new object();
    private EventHandler<InterruptCheckEventArgs> _innerMipInterrupt;
    /// <summary>Occurs when the solver checks whether MIP stopping criteria have been satisfied</summary>
    /// <remarks>If the client wishes to terminate the solve, set the event's user_interrupt to true</remarks>
    public event EventHandler<InterruptCheckEventArgs> MipInterruptCheck
    {
        add
        {
            lock (_mipInterruptLockObject)
            {
                if (_innerMipInterrupt == null)
                {
                    this.startCallback(HighsCallbackType.MipInterrupt);
                }
                _innerMipInterrupt += value;
            }
        }
        remove
        {
            lock (_mipInterruptLockObject)
            {
                _innerMipInterrupt -= value;
                if (_innerMipInterrupt == null)
                {
                    this.stopCallback(HighsCallbackType.MipInterrupt);
                }
            }
        }
    }

    // kCallbackIpmInterrupt as an event
    private readonly object _ipmInterruptLockObject = new object();
    private EventHandler<InterruptCheckEventArgs> _innerIpmInterrupt;
    /// <summary>Occurs when the solver checks whether MIP stopping criteria have been satisfied</summary>
    /// <remarks>If the client wishes to terminate the solve, set the event's user_interrupt to true</remarks>
    public event EventHandler<InterruptCheckEventArgs> IpmInterruptCheck
    {
        add
        {
            lock (_ipmInterruptLockObject)
            {
                if (_innerIpmInterrupt == null)
                {
                    this.startCallback(HighsCallbackType.IpmInterrupt);
                }
                _innerIpmInterrupt += value;
            }
        }
        remove
        {
            lock (_ipmInterruptLockObject)
            {
                _innerIpmInterrupt -= value;
                if (_innerIpmInterrupt == null)
                {
                    this.stopCallback(HighsCallbackType.IpmInterrupt);
                }
            }
        }
    }

    // kCallbackSimplexInterrupt as an event
    private readonly object _simplexInterruptLockObject = new object();
    private EventHandler<InterruptCheckEventArgs> _innerSimplexInterrupt;
    /// <summary>Occurs when the solver checks whether MIP stopping criteria have been satisfied</summary>
    /// <remarks>If the client wishes to terminate the solve, set the event's user_interrupt to true</remarks>
    public event EventHandler<InterruptCheckEventArgs> SimplexInterruptCheck
    {
        add
        {
            lock (_simplexInterruptLockObject)
            {
                if (_innerSimplexInterrupt == null)
                {
                    this.startCallback(HighsCallbackType.SimplexInterrupt);
                }
                _innerSimplexInterrupt += value;
            }
        }
        remove
        {
            lock (_simplexInterruptLockObject)
            {
                _innerSimplexInterrupt -= value;
                if (_innerSimplexInterrupt == null)
                {
                    this.stopCallback(HighsCallbackType.SimplexInterrupt);
                }
            }
        }
    }
#endregion
}

    /// <summary>Data passed to the callback function from HiGHS</summary>
    [StructLayout(LayoutKind.Sequential)]
    internal struct HighsCallbackDataOut
    {
        private IntPtr _ignore;
        public HighsLogType log_type;
        public double running_time;
        public int simplex_iteration_count;
        public int ipm_iteration_count;
        public int pdlp_iteration_count;
        public double objective_function_value;
        public long mip_node_count;
        public long mip_total_lp_iterations;
        public double mip_primal_bound;
        public double mip_dual_bound;
        public double mip_gap;
        // Additional fields omitted, .NET marshaller will just ignore any fields beyond this point
    }

    /// <summary>Data passed from the callback function to HiGHS</summary>
    [StructLayout(LayoutKind.Sequential)]
    internal struct HighsCallbackDataIn
    {
        public int user_interrupt;
    }

    /// <summary>Data for message logging events</summary>
    public class LoggingEventArgs : EventArgs
    {
        /// <summary>The type/level of log message</summary>
        public HighsLogType LogType { get; }
        /// <summary>The log message</summary>
        public string Message { get; }

        public LoggingEventArgs(HighsLogType log_type, string message)
        {
            this.LogType = log_type;
            this.Message = message;
        }
    }

    /// <summary>Data for MIP-related events</summary>
    public class MipEventArgs : EventArgs
    {
        /// <summary>The execution time in seconds</summary>
        public double RunningTime { get; }
        /// <summary>The objective function value of the best integer feasible solution found so far</summary>
        public double ObjectiveFunctionValue { get; }
        /// <summary>The number of MIP nodes explored so far</summary>
        public long MipNodeCount { get; }
        /// <summary>The primal bound</summary>
        public double MipPrimalBound { get; }
        /// <summary>The dual bound</summary>
        public double MipDualBound { get; }
        /// <summary>The relative difference between the primal and dual bounds</summary>
        public double MipGap { get; }

        internal MipEventArgs(HighsCallbackDataOut data)
        {
            this.RunningTime = data.running_time;
            this.ObjectiveFunctionValue = data.objective_function_value;
            this.MipNodeCount = data.mip_node_count;
            this.MipPrimalBound = data.mip_primal_bound;
            this.MipDualBound = data.mip_dual_bound;
            this.MipGap = data.mip_gap;
        }
    }

    public class InterruptCheckEventArgs : EventArgs
    {
        internal InterruptCheckEventArgs(HighsCallbackDataOut data)
        {}

        /// <summary>Whether to interrupt the solver operation currently in progress</summary>
        public bool InterruptSolver { get; set; } = false;
    }

/// <summary>
/// The solution info.
/// </summary>
public class SolutionInfo
{
    /// <summary>
    /// Gets or sets the simplex iteration count.
    /// </summary>
    public int SimplexIterationCount { get; set; }

    /// <summary>
    /// Gets or sets the Interior Point Method (IPM) iteration count.
    /// </summary>
    public int IpmIterationCount { get; set; }

    /// <summary>
    /// Gets or sets the PDLP iteration count.
    /// </summary>
    public int PdlpIterationCount { get; set; }

    /// <summary>
    /// Gets or sets the MIP gap.
    /// </summary>
    public double MipGap { get; set; }

    /// <summary>
    /// Gets or sets the best dual bound.
    /// </summary>
    public double DualBound { get; set; }

    /// <summary>
    /// Gets or sets the MIP node count.
    /// </summary>
    public long NodeCount { get; set; }

    /// <summary>
    /// Gets or sets the objective value.
    /// </summary>
    public double ObjectiveValue { get; set; }
}
}
