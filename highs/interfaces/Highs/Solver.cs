using System.Runtime.InteropServices;
using System.Text;
using Highs.Enums;
using Highs.Records;

namespace Highs;

/// <summary>
/// The Highs Solver interface.
/// </summary>
public class Solver : IDisposable
{
    /// <summary>
    /// The pointer to the _highs instance.
    /// </summary>
    private readonly IntPtr _highs;

    /// <summary>
    /// Indicates whether the instance has been disposed.
    /// </summary>
    private bool _disposed;

    /// <summary>
    /// The name of the Highs library.
    /// </summary>
    private const string HighsLibName = "highs";

    #region Library Imports
    [DllImport(HighsLibName)]
    private static extern int Highs_lpCall(
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
        double[] colvalue,
        double[] coldual,
        double[] rowvalue,
        double[] rowdual,
        int[] colbasisstatus,
        int[] rowbasisstatus,
        ref int modelstatus);

    [DllImport(HighsLibName)]
    private static extern int Highs_mipCall(
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
        int[] integrality,
        double[] colvalue,
        double[] rowvalue,
        ref int modelstatus);

    [DllImport(HighsLibName)]
    private static extern int Highs_qpCall(
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

    double[] colvalue,
    double[] coldual,
    double[] rowvalue,
    double[] rowdual,
    int[] colbasisstatus,
    int[] rowbasisstatus,
    ref int modelstatus);

    [DllImport(HighsLibName)]
    private static extern IntPtr Highs_create();

    [DllImport(HighsLibName)]
    private static extern void Highs_destroy(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_run(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_readModel(IntPtr highs, string filename);

    [DllImport(HighsLibName)]
    private static extern int Highs_writeModel(IntPtr highs, string filename);

    [DllImport(HighsLibName)]
    private static extern int Highs_writePresolvedModel(IntPtr highs, string filename);

    [DllImport(HighsLibName)]
    private static extern int Highs_writeSolutionPretty(IntPtr highs, string filename);

    [DllImport(HighsLibName)]
    private static extern int Highs_getInfinity(IntPtr highs);

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
    private static extern int Highs_passHessian(
        IntPtr highs,
        int dim,
        int numnz,
        int q_format,
        int[] qstart,
        int[] qindex,
        double[] qvalue);

    [DllImport(HighsLibName)]
    private static extern int Highs_setOptionValue(IntPtr highs, string option, string value);

    [DllImport(HighsLibName)]
    private static extern int Highs_setBoolOptionValue(IntPtr highs, string option, int value);

    [DllImport(HighsLibName)]
    private static extern int Highs_setIntOptionValue(IntPtr highs, string option, int value);

    [DllImport(HighsLibName)]
    private static extern int Highs_setDoubleOptionValue(IntPtr highs, string option, double value);

    [DllImport(HighsLibName)]
    private static extern int Highs_setStringOptionValue(IntPtr highs, string option, string value);

    [DllImport(HighsLibName)]
    private static extern int Highs_getBoolOptionValue(IntPtr highs, string option, out int value);

    [DllImport(HighsLibName)]
    private static extern int Highs_getIntOptionValue(IntPtr highs, string option, out int value);

    [DllImport(HighsLibName)]
    private static extern int Highs_getDoubleOptionValue(IntPtr highs, string option, out double value);

    [DllImport(HighsLibName)]
    private static extern int Highs_getStringOptionValue(IntPtr highs, string option, [Out] StringBuilder value);

    [DllImport(HighsLibName)]
    private static extern int Highs_getSolution(IntPtr highs, double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual);

    [DllImport(HighsLibName)]
    private static extern int Highs_getNumCol(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_getNumRow(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_getNumNz(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_getHessianNumNz(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_getBasis(IntPtr highs, int[] colstatus, int[] rowstatus);

    [DllImport(HighsLibName)]
    private static extern double Highs_getObjectiveValue(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_getIterationCount(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_getModelStatus(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_addRow(IntPtr highs, double lower, double upper, int num_new_nz, int[] indices, double[] values);

    [DllImport(HighsLibName)]
    private static extern int Highs_addRows(
        IntPtr highs,
        int num_new_row,
        double[] lower,
        double[] upper,
        int num_new_nz,
        int[] starts,
        int[] indices,
        double[] values);

    [DllImport(HighsLibName)]
    private static extern int Highs_addCol(
        IntPtr highs,
        double cost,
        double lower,
        double upper,
        int num_new_nz,
        int[] indices,
        double[] values);

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
    private static extern int Highs_changeObjectiveSense(IntPtr highs, int sense);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeColCost(IntPtr highs, int column, double cost);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeColsCostBySet(IntPtr highs, int num_set_entries, int[] set, double[] cost);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeColsCostByMask(IntPtr highs, int[] mask, double[] cost);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeColBounds(IntPtr highs, int column, double lower, double upper);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeColsBoundsByRange(IntPtr highs, int from_col, int to_col, double[] lower, double[] upper);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeColsBoundsBySet(IntPtr highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeColsBoundsByMask(IntPtr highs, int[] mask, double[] lower, double[] upper);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeRowBounds(IntPtr highs, int row, double lower, double upper);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeRowsBoundsByRange(IntPtr highs, int from_row, int to_row, double[] lower, double[] upper);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeRowsBoundsBySet(IntPtr highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeRowsBoundsByMask(IntPtr highs, int[] mask, double[] lower, double[] upper);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeColsIntegralityByRange(IntPtr highs, int from_col, int to_col, int[] integrality);

    [DllImport(HighsLibName)]
    private static extern int Highs_changeCoeff(IntPtr highs, int row, int column, double value);

    [DllImport(HighsLibName)]
    private static extern int Highs_deleteColsByRange(IntPtr highs, int from_col, int to_col);

    [DllImport(HighsLibName)]
    private static extern int Highs_deleteColsBySet(IntPtr highs, int num_set_entries, int[] set);

    [DllImport(HighsLibName)]
    private static extern int Highs_deleteColsByMask(IntPtr highs, int[] mask);

    [DllImport(HighsLibName)]
    private static extern int Highs_deleteRowsByRange(IntPtr highs, int from_row, int to_row);

    [DllImport(HighsLibName)]
    private static extern int Highs_deleteRowsBySet(IntPtr highs, int num_set_entries, int[] set);

    [DllImport(HighsLibName)]
    private static extern int Highs_deleteRowsByMask(IntPtr highs, int[] mask);

    [DllImport(HighsLibName)]
    private static extern int Highs_getDoubleInfoValue(IntPtr highs, string info, out double value);

    [DllImport(HighsLibName)]
    private static extern int Highs_getIntInfoValue(IntPtr highs, string info, out int value);

    [DllImport(HighsLibName)]
    private static extern int Highs_getInt64InfoValue(IntPtr highs, string info, out long value);

    [DllImport(HighsLibName)]
    private static extern int Highs_setSolution(IntPtr highs, double[] col_value, double[] row_value, double[] col_dual, double[] row_dual);

    [DllImport(HighsLibName)]
    private static extern int Highs_setSparseSolution(IntPtr highs, int num_entries, int[] index, double[] value);

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
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

    [DllImport(HighsLibName)]
    private static extern int Highs_getBasicVariables(IntPtr highs, int[] basic_variables);

    [DllImport(HighsLibName)]
    private static extern int Highs_getBasisInverseRow(IntPtr highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [DllImport(HighsLibName)]
    private static extern int Highs_getBasisInverseCol(IntPtr highs, int column, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [DllImport(HighsLibName)]
    private static extern int Highs_getBasisSolve(
        IntPtr highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_num_nz,
        int[] solution_indices);

    [DllImport(HighsLibName)]
    private static extern int Highs_getBasisTransposeSolve(
        IntPtr highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_nz,
        int[] solution_indices);

    [DllImport(HighsLibName)]
    private static extern int Highs_getReducedRow(IntPtr highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [DllImport(HighsLibName)]
    private static extern int Highs_getReducedColumn(IntPtr highs, int column, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [DllImport(HighsLibName)]
    private static extern int Highs_clearModel(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_clearSolver(IntPtr highs);

    [DllImport(HighsLibName)]
    private static extern int Highs_passColName(IntPtr highs, int column, string name);

    [DllImport(HighsLibName)]
    private static extern int Highs_passRowName(IntPtr highs, int row, string name);

    [DllImport(HighsLibName)]
    private static extern int Highs_writeOptions(IntPtr highs, string filename);

    [DllImport(HighsLibName)]
    private static extern int Highs_writeOptionsDeviations(IntPtr highs, string filename);
    #endregion

    /// <summary>
    /// Calls the Highs solver in a single call of a LP.
    /// </summary>
    /// <param name="model"></param>
    /// <param name="solution"></param>
    /// <param name="basisInfo"></param>
    /// <param name="modelStatus"></param>
    /// <returns></returns>
    public static Status LpCall(Model model, ref Solution solution, out BasisInfo basisInfo, out ModelStatus modelStatus)
    {
        var numberOfColumns = model.ColumnCost.Length;
        var numberOfRows = model.RowLower.Length;
        var numberOfMatrixValues = model.MatrixValues.Length;

        var columnBasisStatus = new int[numberOfColumns];
        var rowBasisStatus = new int[numberOfRows];

        var modelstate = 0;

        var status = (Status)Highs_lpCall(
            numberOfColumns,
            numberOfRows,
            numberOfMatrixValues,
            (int)model.MatrixFormat,
            (int)model.ObjectiveSense,
            model.Offset,
            model.ColumnCost,
            model.ColumnLower,
            model.ColumnUpper,
            model.RowLower,
            model.RowUpper,
            model.MatrixStart,
            model.MatrixIndices,
            model.MatrixValues,
            solution.ColumnValue,
            solution.ColumnDual,
            solution.RowValue,
            solution.RowDual,
            columnBasisStatus,
            rowBasisStatus,
            ref modelstate);

        modelStatus = (ModelStatus)modelstate;
        basisInfo = new(Array.ConvertAll(columnBasisStatus, c => (BasisStatus)c), Array.ConvertAll(rowBasisStatus, r => (BasisStatus)r));

        return status;
    }

    /// <summary>
    /// Calls the Highs solver in a single call of a Mip.
    /// </summary>
    /// <param name="model"></param>
    /// <param name="solution"></param>
    /// <param name="modelStatus"></param>
    /// <returns></returns>
    public static Status MipCall(Model model, ref Solution solution, out ModelStatus modelStatus)
    {
        var numberOfColumns = model.ColumnCost.Length;
        var numberOfRows = model.RowLower.Length;
        var numberOfMatrixValues = model.MatrixValues.Length;

        var modelstate = 0;

        var status = (Status)Highs_mipCall(
            numberOfColumns,
            numberOfRows,
            numberOfMatrixValues,
            (int)model.MatrixFormat,
            (int)model.ObjectiveSense,
            model.Offset,
            model.ColumnCost,
            model.ColumnLower,
            model.ColumnUpper,
            model.RowLower,
            model.RowUpper,
            model.MatrixStart,
            model.MatrixIndices,
            model.MatrixValues,
            Array.ConvertAll(model.VariableTypes, v => (int)v),
            solution.ColumnValue,
            solution.RowValue,
            ref modelstate);

        modelStatus = (ModelStatus)modelstate;

        return status;
    }

    /// <summary>
    /// Calls the Highs solver in a single call of a QP.
    /// </summary>
    /// <param name="model"></param>
    /// <param name="solution"></param>
    /// <param name="modelStatus"></param>
    /// <returns></returns>
    public static Status QPCall(Model model, ref Solution solution, out BasisInfo basisInfo, out ModelStatus modelStatus)
    {
        var numberOfColumns = model.ColumnCost.Length;
        var numberOfRows = model.RowLower.Length;
        var numberOfMatrixValues = model.MatrixValues.Length;
        var numberOfHessianValues = model.Hessian.Values.Length;

        var columnBasisStatus = new int[numberOfColumns];
        var rowBasisStatus = new int[numberOfRows];

        var modelstate = 0;

        var status = (Status)Highs_qpCall(
            numberOfColumns,
            numberOfRows,
            numberOfMatrixValues,
            numberOfHessianValues,
            (int)model.MatrixFormat,
            (int)model.Hessian.HessianFormat,
            (int)model.ObjectiveSense,
            model.Offset,
            model.ColumnCost,
            model.ColumnLower,
            model.ColumnUpper,
            model.RowLower,
            model.RowUpper,
            model.MatrixStart,
            model.MatrixIndices,
            model.MatrixValues,
            model.Hessian.Start,
            model.Hessian.Index,
            model.Hessian.Values,
            solution.ColumnValue,
            solution.ColumnDual,
            solution.RowValue,
            solution.RowDual,
            columnBasisStatus,
            rowBasisStatus,
            ref modelstate);

        modelStatus = (ModelStatus)modelstate;
        basisInfo = new(Array.ConvertAll(columnBasisStatus, c => (BasisStatus)c), Array.ConvertAll(rowBasisStatus, r => (BasisStatus)r));

        return status;
    }

    /// <summary>
    /// The default constructor.
    /// </summary>
    public Solver() => _highs = Highs_create();

    /// <summary>
    /// The destructor.
    /// </summary>
    ~Solver() => Dispose(false);

    /// <summary>
    /// Disposes the instance.
    /// </summary>
    public void Dispose()
    {
        Dispose(true);
        GC.SuppressFinalize(this);
    }

    /// <summary>
    /// Disposes the instance.
    /// </summary>
    /// <param name="disposing"></param>
    protected virtual void Dispose(bool disposing)
    {
        if (_disposed)
        {
            return;
        }
        Highs_destroy(_highs);
        _disposed = true;
    }

    /// <summary>
    /// Runs the solver.
    /// </summary>
    /// <returns></returns>
    public Status Run() => (Status)Highs_run(_highs);

    /// <summary>
    /// Reads a model from file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status ReadModel(string filename) => (Status)Highs_readModel(_highs, filename);

    /// <summary>
    /// Writes the model to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WriteModel(string filename) => (Status)Highs_writeModel(_highs, filename);

    /// <summary>
    /// Writes the presolved model to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WritePresolvedModel(string filename) => (Status)Highs_writePresolvedModel(_highs, filename);

    /// <summary>
    /// Writes the solution to file in a pretty format.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WriteSolutionPretty(string filename) => (Status)Highs_writeSolutionPretty(_highs, filename);

    /// <summary>
    /// Gets the infinity value.
    /// </summary>
    /// <returns></returns>
    public double GetInfinity() => Highs_getInfinity(_highs);

    /// <summary>
    /// Passes the LP model to Highs.
    /// </summary>
    /// <param name="model"></param>
    /// <returns></returns>
    public Status PassLp(Model model)
    {
        return (Status)Highs_passLp(
            _highs,
            model.ColumnCost.Length,
            model.RowLower.Length,
            model.MatrixValues.Length,
            (int)model.MatrixFormat,
            (int)model.ObjectiveSense,
            model.Offset,
            model.ColumnCost,
            model.ColumnLower,
            model.ColumnUpper,
            model.RowLower,
            model.RowUpper,
            model.MatrixStart,
            model.MatrixIndices,
            model.MatrixValues);
    }

    /// <summary>
    /// Passes the MIP model to Highs.
    /// </summary>
    /// <param name="model"></param>
    /// <returns></returns>
    public Status PassMip(Model model)
    {
        return (Status)Highs_passMip(
            _highs,
            model.ColumnCost.Length,
            model.RowLower.Length,
            model.MatrixValues.Length,
            (int)model.MatrixFormat,
            (int)model.ObjectiveSense,
            model.Offset,
            model.ColumnCost,
            model.ColumnLower,
            model.ColumnUpper,
            model.RowLower,
            model.RowUpper,
            model.MatrixStart,
            model.MatrixIndices,
            model.MatrixValues,
            Array.ConvertAll(model.VariableTypes, v => (int)v));
    }

    /// <summary>
    /// Sets the option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status SetOptionValue(string option, string value) => (Status)Highs_setOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the string option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status SetStringOptionValue(string option, string value) => (Status)Highs_setStringOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the boolean option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status SetBoolOptionValue(string option, int value) => (Status)Highs_setBoolOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the double option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status SetDoubleOptionValue(string option, double value) => (Status)Highs_setDoubleOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the integer option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status SetIntOptionValue(string option, int value) => (Status)Highs_setIntOptionValue(_highs, option, value);

    /// <summary>
    /// Gets the string option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status GetStringOptionValue(string option, out string value)
    {
        var stringBuilder = new StringBuilder();
        var result = (Status)Highs_getStringOptionValue(_highs, option, stringBuilder);
        value = stringBuilder.ToString();
        return result;
    }

    /// <summary>
    /// Gets the boolean option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status GetBoolOptionValue(string option, out int value) => (Status)Highs_getBoolOptionValue(_highs, option, out value);

    /// <summary>
    /// Gets the double option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status GetDoubleOptionValue(string option, out double value) => (Status)Highs_getDoubleOptionValue(_highs, option, out value);

    /// <summary>
    /// Gets the integer option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status GetIntOptionValue(string option, out int value) => (Status)Highs_getIntOptionValue(_highs, option, out value);

    /// <summary>
    /// Gets the number of columns.
    /// </summary>
    /// <returns></returns>
    public int GetNumberOfColumns() => Highs_getNumCol(_highs);

    /// <summary>
    /// Gets the number of rows.
    /// </summary>
    /// <returns></returns>
    public int GetNumberOfRows() => Highs_getNumRow(_highs);

    /// <summary>
    /// Gets the number of non-zero entries.
    /// </summary>
    /// <returns></returns>
    public int GetNumberOfNonZeroEntries() => Highs_getNumNz(_highs);

    /// <summary>
    /// Gets the solution.
    /// </summary>
    /// <param name="solution"></param>
    /// <returns></returns>
    public Status GetSolution(out Solution solution)
    {
        var numberOfColumns = GetNumberOfColumns();
        var numberOfRows = GetNumberOfRows();

        solution = new Solution(numberOfColumns, numberOfRows);
        return (Status)Highs_getSolution(_highs, solution.ColumnValue, solution.ColumnDual, solution.RowValue, solution.RowDual);
    }

    /// <summary>
    /// Passes the Hessian to Highs.
    /// </summary>
    /// <param name="hessian"></param>
    /// <returns></returns>
    public Status PassHessian(Hessian hessian)
    {
        return (Status)Highs_passHessian(_highs,
                                         hessian.Dimension,
                                         hessian.Values.Length,
                                         (int)hessian.HessianFormat,
                                         hessian.Start,
                                         hessian.Index,
                                         hessian.Values);
    }

    /// <summary>
    /// Gets the basis information.
    /// </summary>
    /// <returns></returns>
    public Status GetBasis(out BasisInfo basisInfo)
    {
        var numberOfColumns = GetNumberOfColumns();
        var numberOfRows = GetNumberOfRows();

        var columnBasisStatus = new int[numberOfColumns];
        var rowBasisStatus = new int[numberOfRows];

        var status = (Status)Highs_getBasis(_highs, columnBasisStatus, rowBasisStatus);
        if (status == Status.Error)
        {
            basisInfo = null;
            return Status.Error;
        }

        basisInfo = new BasisInfo(Array.ConvertAll(columnBasisStatus, x => (BasisStatus)x), Array.ConvertAll(rowBasisStatus, x => (BasisStatus)x));

        return status;
    }

    /// <summary>
    /// Gets the objective value.
    /// </summary>
    /// <returns></returns>
    public double GetObjectiveValue() => Highs_getObjectiveValue(_highs);

    /// <summary>
    /// Gets the model status.
    /// </summary>
    /// <returns></returns>
    public ModelStatus GetModelStatus() => (ModelStatus)Highs_getModelStatus(_highs);

    /// <summary>
    /// Gets the iteration count.
    /// </summary>
    /// <returns></returns>
    public int GetIterationCount() => Highs_getIterationCount(_highs);

    /// <summary>
    /// Adds a row to the model.
    /// </summary>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <param name="indices"></param>
    /// <param name="values"></param>
    /// <returns></returns>
    public Status AddRow(double lower, double upper, int[] indices, double[] values)
    {
        return (Status)Highs_addRow(_highs, lower, upper, indices.Length, indices, values);
    }

    /// <summary>
    /// Adds multiple rows to the model.
    /// </summary>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <param name="starts"></param>
    /// <param name="indices"></param>
    /// <param name="values"></param>
    /// <returns></returns>
    public Status AddRows(double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
    {
        return (Status)Highs_addRows(_highs, lower.Length, lower, upper, indices.Length, starts, indices, values);
    }

    /// <summary>
    /// Adds a column to the model.
    /// </summary>
    /// <param name="cost"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <param name="indices"></param>
    /// <param name="values"></param>
    /// <returns></returns>
    public Status AddColumn(double cost, double lower, double upper, int[] indices, double[] values)
    {
        return (Status)Highs_addCol(_highs, cost, lower, upper, indices.Length, indices, values);
    }

    /// <summary>
    /// Adds multiple columns to the model.
    /// </summary>
    /// <param name="costs"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <param name="starts"></param>
    /// <param name="indices"></param>
    /// <param name="values"></param>
    /// <returns></returns>
    public Status AddColumns(double[] costs, double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
    {
        return (Status)Highs_addCols(
            _highs,
            costs.Length,
            costs,
            lower,
            upper,
            indices.Length,
            starts,
            indices,
            values);
    }

    /// <summary>
    /// Changes the objective sense.
    /// </summary>
    /// <param name="sense"></param>
    /// <returns></returns>
    public Status ChangeObjectiveSense(ObjectiveSense sense) => (Status)Highs_changeObjectiveSense(_highs, (int)sense);

    /// <summary>
    /// Changes the cost of a column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="cost"></param>
    /// <returns></returns>
    public Status ChangeColumnCost(int column, double cost) => (Status)Highs_changeColCost(_highs, column, cost);

    /// <summary>
    /// Changes the costs of multiple columns by set.
    /// </summary>
    /// <param name="columns"></param>
    /// <param name="costs"></param>
    /// <returns></returns>
    public Status ChangeColumnsCostBySet(int[] columns, double[] costs) => (Status)Highs_changeColsCostBySet(_highs, columns.Length, columns, costs);

    /// <summary>
    /// Changes the costs of multiple columns by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <param name="cost"></param>
    /// <returns></returns>
    public Status ChangeColumnsCostByMask(bool[] mask, double[] cost) => (Status)Highs_changeColsCostByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0), cost);

    /// <summary>
    /// Changes the bounds of a column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeColumnBounds(int column, double lower, double upper) => (Status)Highs_changeColBounds(_highs, column, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple columns by range.
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeColumnsBoundsByRange(int from, int to, double[] lower, double[] upper) => (Status)Highs_changeColsBoundsByRange(_highs, from, to, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple columns by set.
    /// </summary>
    /// <param name="columns"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeColumnsBoundsBySet(int[] columns, double[] lower, double[] upper)
    {
        return (Status)Highs_changeColsBoundsBySet(_highs, columns.Length, columns, lower, upper);
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="mask"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeColumnsBoundsByMask(bool[] mask, double[] lower, double[] upper)
    {
        return (Status)Highs_changeColsBoundsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0), lower, upper);
    }

    /// <summary>
    /// Changes the bounds of a row.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeRowBounds(int row, double lower, double upper) => (Status)Highs_changeRowBounds(_highs, row, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple rows by range.
    /// </summary>
    /// <param name="rows"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeRowsBoundsBySet(int[] rows, double[] lower, double[] upper) => (Status)Highs_changeRowsBoundsBySet(_highs, rows.Length, rows, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple rows by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeRowsBoundsByMask(bool[] mask, double[] lower, double[] upper) => (Status)Highs_changeRowsBoundsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0), lower, upper);

    /// <summary>
    /// Changes the integrality of multiple columns by range.
    /// </summary>
    /// <param name="from_col"></param>
    /// <param name="to_col"></param>
    /// <param name="integrality"></param>
    /// <returns></returns>
    public Status ChangeColumnsIntegralityByRange(int from_col, int to_col, VariableType[] integrality)
    {
        return (Status)Highs_changeColsIntegralityByRange(_highs, from_col, to_col, Array.ConvertAll(integrality, item => (int)item));
    }

    /// <summary>
    /// Changes a coefficient in the constraint matrix.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="column"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status ChangeCoefficient(int row, int column, double value) => (Status)Highs_changeCoeff(_highs, row, column, value);

    /// <summary>
    /// Deletes multiple columns by range.
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <returns></returns>
    public Status DeleteColumnsByRange(int from, int to) => (Status)Highs_deleteColsByRange(_highs, from, to);

    /// <summary>
    /// Deletes multiple columns by set.
    /// </summary>
    /// <param name="columns"></param>
    /// <returns></returns>
    public Status DeleteColumnsBySet(int[] columns) => (Status)Highs_deleteColsBySet(_highs, columns.Length, columns);

    /// <summary>
    /// Deletes multiple columns by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <returns></returns>
    public Status DeleteColumnsByMask(bool[] mask) => (Status)Highs_deleteColsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0));

    /// <summary>
    /// Deletes multiple rows by range.
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <returns></returns>
    public Status DeleteRowsByRange(int from, int to) => (Status)Highs_deleteRowsByRange(_highs, from, to);

    /// <summary>
    /// Deletes multiple rows by set.
    /// </summary>
    /// <param name="rows"></param>
    /// <returns></returns>
    public Status DeleteRowsBySet(int[] rows) => (Status)Highs_deleteRowsBySet(_highs, rows.Length, rows);

    /// <summary>
    /// Deletes multiple rows by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <returns></returns>
    public Status DeleteRowsByMask(bool[] mask) => (Status)Highs_deleteRowsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0));

    /// <summary>
    /// Delegate for getting info values.
    /// </summary>
    /// <typeparam name="TValue"></typeparam>
    /// <param name="highs"></param>
    /// <param name="infoName"></param>
    /// <param name="output"></param>
    /// <returns></returns>
    delegate int HighsGetInfoDelegate<TValue>(IntPtr highs, string infoName, out TValue output);

    /// <summary>
    /// Gets a value or a fallback if an error occurs.
    /// </summary>
    /// <typeparam name="TValue"></typeparam>
    /// <param name="highsGetInfoDelegate"></param>
    /// <param name="infoName"></param>
    /// <param name="fallback"></param>
    /// <returns></returns>
    private TValue GetValueOrFallback<TValue>(HighsGetInfoDelegate<TValue> highsGetInfoDelegate, string infoName, TValue fallback)
    {
        try
        {
            var status = (Status)highsGetInfoDelegate(_highs, infoName, out var value);
            return status != Status.Ok ? fallback : value;
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
    public SolutionInfo GetInfo()
    {
        // TODO: This object does not contian the "complete" info from the C api. Add further props, if you need them.
        return new SolutionInfo(GetValueOrFallback(Highs_getIntInfoValue, "simplex_iteration_count", 0),
                                GetValueOrFallback(Highs_getIntInfoValue, "ipm_iteration_count", 0),
                                GetValueOrFallback(Highs_getIntInfoValue, "pdlp_iteration_count", 0),
                                GetValueOrFallback(Highs_getDoubleInfoValue, "mip_gap", double.NaN),
                                GetValueOrFallback(Highs_getDoubleInfoValue, "mip_dual_bound", double.NaN),
                                GetValueOrFallback(Highs_getInt64InfoValue, "mip_node_count", 0L),
                                GetValueOrFallback(Highs_getDoubleInfoValue, "objective_function_value", double.NaN));
    }

    /// <summary>
    /// Sets the solution.
    /// </summary>
    /// <param name="solution"></param>
    /// <returns></returns>
    public Status SetSolution(Solution solution) => (Status)Highs_setSolution(_highs, solution.ColumnValue, solution.ColumnDual, solution.RowValue, solution.RowDual);

    /// <summary>
    /// Set a partial primal solution by passing values for a set of variables
    /// </summary>
    /// <remarks>
    /// The sparse solution set by this function has values for a subset of the model's variables.
    /// For each entry in <paramref name="valuesByIndex"/>, the key identifies a variable by index, and
    /// the value indicates the variable's value in the sparse solution.
    /// </remarks>
    /// <param name="valuesByIndex">A dictionary that maps variable indices to variable values</param>
    /// <returns></returns>
    public Status SetSparseSolution(IReadOnlyDictionary<int, double> valuesByIndex)
    {
        return (Status)Highs_setSparseSolution(_highs, valuesByIndex.Count, [.. valuesByIndex.Keys], [.. valuesByIndex.Values]);
    }

    /// <summary>
    /// Gets the basic variables.
    /// </summary>
    /// <param name="basic_variables"></param>
    /// <returns></returns>
    public Status GetBasicVariables(ref int[] basic_variables) => (Status)Highs_getBasicVariables(_highs, basic_variables);

    /// <summary>
    /// Gets a row of the basis inverse.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="row_vector"></param>
    /// <param name="row_num_nz"></param>
    /// <param name="row_indices"></param>
    /// <returns></returns>
    public Status GetBasisInverseRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
    {
        return (Status)Highs_getBasisInverseRow(_highs, row, row_vector, ref row_num_nz, row_indices);
    }

    /// <summary>
    /// Gets a column of the basis inverse.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="column_vector"></param>
    /// <param name="column_num_nz"></param>
    /// <param name="column_indices"></param>
    /// <returns></returns>
    public Status GetBasisInverseColumn(int column, double[] column_vector, ref int column_num_nz, int[] column_indices)
    {
        return (Status)Highs_getBasisInverseCol(_highs, column, column_vector, ref column_num_nz, column_indices);
    }

    /// <summary>
    /// Gets the basis solve.
    /// </summary>
    /// <param name="rhs"></param>
    /// <param name="solution_vector"></param>
    /// <param name="solution_num_nz"></param>
    /// <param name="solution_indices"></param>
    /// <returns></returns>
    public Status GetBasisSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
    {
        return (Status)Highs_getBasisSolve(_highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
    }

    /// <summary>
    /// Gets the basis transpose solve.
    /// </summary>
    /// <param name="rhs"></param>
    /// <param name="solution_vector"></param>
    /// <param name="solution_num_nz"></param>
    /// <param name="solution_indices"></param>
    /// <returns></returns>
    public Status GetBasisTransposeSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
    {
        return (Status)Highs_getBasisTransposeSolve(_highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
    }

    /// <summary>
    /// Gets a reduced row.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="row_vector"></param>
    /// <param name="row_num_nz"></param>
    /// <param name="row_indices"></param>
    /// <returns></returns>
    public Status GetReducedRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
    {
        return (Status)Highs_getReducedRow(_highs, row, row_vector, ref row_num_nz, row_indices);
    }

    /// <summary>
    /// Gets a reduced column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="column_vector"></param>
    /// <param name="column_num_nz"></param>
    /// <param name="column_indices"></param>
    /// <returns></returns>
    public Status GetReducedColumn(int column, double[] column_vector, ref int column_num_nz, int[] column_indices)
    {
        return (Status)Highs_getReducedColumn(_highs, column, column_vector, ref column_num_nz, column_indices);
    }

    /// <summary>
    /// Clears the model.
    /// </summary>
    /// <returns></returns>
    public Status ClearModel() => (Status)Highs_clearModel(_highs);

    /// <summary>
    /// Clears the solver.
    /// </summary>
    /// <returns></returns>
    public Status ClearSolver() => (Status)Highs_clearSolver(_highs);

    /// <summary>
    /// Passes the name of a column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="name"></param>
    /// <returns></returns>
    public Status PassColumnName(int column, string name) => (Status)Highs_passColName(_highs, column, name);

    /// <summary>
    /// Passes the name of a row.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="name"></param>
    /// <returns></returns>
    public Status PassRowName(int row, string name) => (Status)Highs_passRowName(_highs, row, name);

    /// <summary>
    /// Writes the options to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WriteOptions(string filename) => (Status)Highs_writeOptions(_highs, filename);

    /// <summary>
    /// Writes the options deviations to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WriteOptionsDeviations(string filename) => (Status)Highs_writeOptionsDeviations(_highs, filename);
}
