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
    private const string HighsLibName = "_highs";

    #region Library Imports
    [LibraryImport(HighsLibName)]
    private static extern int Highs_call(
        int numcol,
        int numrow,
        int numnz,
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

    [LibraryImport(HighsLibName)]
    private static extern IntPtr Highs_create();

    [LibraryImport(HighsLibName)]
    private static extern void Highs_destroy(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_run(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_readModel(IntPtr _highs, string filename);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_writeModel(IntPtr _highs, string filename);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_writePresolvedModel(IntPtr _highs, string filename);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_writeSolutionPretty(IntPtr _highs, string filename);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getInfinity(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_passLp(
        IntPtr _highs,
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

    [LibraryImport(HighsLibName)]
    private static extern int Highs_passMip(
        IntPtr _highs,
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

    [LibraryImport(HighsLibName)]
    private static extern int Highs_setOptionValue(IntPtr _highs, string option, string value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_setBoolOptionValue(IntPtr _highs, string option, int value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_setIntOptionValue(IntPtr _highs, string option, int value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_setDoubleOptionValue(IntPtr _highs, string option, double value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_setStringOptionValue(IntPtr _highs, string option, string value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getBoolOptionValue(IntPtr _highs, string option, out int value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getIntOptionValue(IntPtr _highs, string option, out int value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getDoubleOptionValue(IntPtr _highs, string option, out double value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getStringOptionValue(IntPtr _highs, string option, [Out] StringBuilder value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getSolution(IntPtr _highs, double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getNumCol(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getNumRow(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getNumNz(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getBasis(IntPtr _highs, int[] colstatus, int[] rowstatus);

    [LibraryImport(HighsLibName)]
    private static extern double Highs_getObjectiveValue(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getIterationCount(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getModelStatus(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_addRow(IntPtr _highs, double lower, double upper, int num_new_nz, int[] indices, double[] values);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_addRows(
        IntPtr _highs,
        int num_new_row,
        double[] lower,
        double[] upper,
        int num_new_nz,
        int[] starts,
        int[] indices,
        double[] values);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_addCol(
        IntPtr _highs,
        double cost,
        double lower,
        double upper,
        int num_new_nz,
        int[] indices,
        double[] values);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_addCols(
        IntPtr _highs,
        int num_new_col,
        double[] costs,
        double[] lower,
        double[] upper,
        int num_new_nz,
        int[] starts,
        int[] indices,
        double[] values);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeObjectiveSense(IntPtr _highs, int sense);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeColCost(IntPtr _highs, int col, double cost);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeColsCostBySet(IntPtr _highs, int num_set_entries, int[] set, double[] cost);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeColsCostByMask(IntPtr _highs, int[] mask, double[] cost);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeColBounds(IntPtr _highs, int col, double lower, double upper);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeColsBoundsByRange(IntPtr _highs, int from_col, int to_col, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeColsBoundsBySet(IntPtr _highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeColsBoundsByMask(IntPtr _highs, int[] mask, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeRowBounds(IntPtr _highs, int row, double lower, double upper);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeRowsBoundsBySet(IntPtr _highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeRowsBoundsByMask(IntPtr _highs, int[] mask, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeColsIntegralityByRange(IntPtr _highs, int from_col, int to_col, int[] integrality);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_changeCoeff(IntPtr _highs, int row, int col, double value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_deleteColsByRange(IntPtr _highs, int from_col, int to_col);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_deleteColsBySet(IntPtr _highs, int num_set_entries, int[] set);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_deleteColsByMask(IntPtr _highs, int[] mask);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_deleteRowsByRange(IntPtr _highs, int from_row, int to_row);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_deleteRowsBySet(IntPtr _highs, int num_set_entries, int[] set);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_deleteRowsByMask(IntPtr _highs, int[] mask);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getDoubleInfoValue(IntPtr _highs, string info, out double value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getIntInfoValue(IntPtr _highs, string info, out int value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getInt64InfoValue(IntPtr _highs, string info, out long value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_setSolution(IntPtr _highs, double[] col_value, double[] row_value, double[] col_dual, double[] row_dual);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getColsByRange(
        IntPtr _highs,
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

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getColsBySet(
        IntPtr _highs,
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

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getColsByMask(
        IntPtr _highs,
        int[] mask,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getRowsByRange(
        IntPtr _highs,
        int from_row,
        int to_row,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getRowsBySet(
        IntPtr _highs,
        int num_set_entries,
        int[] set,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getRowsByMask(
        IntPtr _highs,
        int[] mask,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getBasicVariables(IntPtr _highs, int[] basic_variables);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getBasisInverseRow(IntPtr _highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getBasisInverseCol(IntPtr _highs, int col, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getBasisSolve(
        IntPtr _highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_num_nz,
        int[] solution_indices);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getBasisTransposeSolve(
        IntPtr _highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_nz,
        int[] solution_indices);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getReducedRow(IntPtr _highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_getReducedColumn(IntPtr _highs, int col, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_clearModel(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_clearSolver(IntPtr _highs);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_passColName(IntPtr _highs, int col, string name);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_passRowName(IntPtr _highs, int row, string name);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_writeOptions(IntPtr _highs, string filename);

    [LibraryImport(HighsLibName)]
    private static extern int Highs_writeOptionsDeviations(IntPtr _highs, string filename);
    #endregion

    /// <summary>
    /// Calls the Highs solver in a single call.
    /// </summary>
    /// <param name="model"></param>
    /// <param name="solution"></param>
    /// <param name="basisInfo"></param>
    /// <param name="modelStatus"></param>
    /// <returns></returns>
    public static Status Call(Model model, ref Solution solution, out BasisInfo basisInfo, out ModelStatus modelStatus)
    {
        var numberOfColumns = model.ColumnCost.Length;
        var numberOfRows = model.RowLower.Length;
        var numberOfMatrixValues = model.MatrixValues.Length;

        var columnBasisStatus = new int[numberOfColumns];
        var rowBasisStatus = new int[numberOfRows];

        var modelstate = 0;

        var status = (Status)Highs_call(
            numberOfColumns,
            numberOfRows,
            numberOfMatrixValues,
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
        basisInfo = new([.. columnBasisStatus.Select(x => (BasisStatus)x)], [.. rowBasisStatus.Select(x => (BasisStatus)x)]);

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

    public Status passLp(HighsModel model)
    {
        return (Status)Highs_passLp(
            _highs,
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

    public Status passMip(HighsModel model)
    {
        return (Status)Highs_passMip(
            _highs,
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

    public Status getStringOptionValue(string option, out string value)
    {
        var stringBuilder = new StringBuilder();
        var result = (Status)Highs_getStringOptionValue(_highs, option, stringBuilder);
        value = stringBuilder.ToString();
        return result;
    }

    public Status getBoolOptionValue(string option, out int value)
    {
        return (Status)Highs_getBoolOptionValue(_highs, option, out value);
    }

    public Status getDoubleOptionValue(string option, out double value)
    {
        return (Status)Highs_getDoubleOptionValue(_highs, option, out value);
    }

    public Status getIntOptionValue(string option, out int value)
    {
        return (Status)Highs_getIntOptionValue(_highs, option, out value);
    }

    public int getNumCol()
    {
        return Highs_getNumCol(_highs);
    }

    public int getNumRow()
    {
        return Highs_getNumRow(_highs);
    }

    public int getNumNz()
    {
        return Highs_getNumNz(_highs);
    }

    public HighsSolution getSolution()
    {
        int numberOfColumns = getNumCol();
        int numberOfRows = getNumRow();

        HighsSolution sol = new HighsSolution(numberOfColumns, numberOfRows);
        Highs_getSolution(_highs, sol.colvalue, sol.coldual, sol.rowvalue, sol.rowdual);

        return sol;
    }

    public HighsBasis getBasis()
    {
        int numberOfColumns = getNumCol();
        int numberOfRows = getNumRow();

        int[] colbasstat = new int[numberOfColumns];
        int[] rowbasstat = new int[numberOfRows];

        Highs_getBasis(_highs, colbasstat, rowbasstat);
        HighsBasis bas = new HighsBasis(
            colbasstat.Select(x => (HighsBasisStatus)x).ToArray(),
            rowbasstat.Select(x => (HighsBasisStatus)x).ToArray());

        return bas;
    }

    public double getObjectiveValue()
    {
        return Highs_getObjectiveValue(_highs);
    }

    public HighsModelStatus GetModelStatus()
    {
        return (HighsModelStatus)Highs_getModelStatus(_highs);
    }

    public int getIterationCount()
    {
        return Highs_getIterationCount(_highs);
    }

    public Status addRow(double lower, double upper, int[] indices, double[] values)
    {
        return (Status)Highs_addRow(_highs, lower, upper, indices.Length, indices, values);
    }

    public Status addRows(double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
    {
        return (Status)Highs_addRows(_highs, lower.Length, lower, upper, indices.Length, starts, indices, values);
    }

    public Status addCol(double cost, double lower, double upper, int[] indices, double[] values)
    {
        return (Status)Highs_addCol(_highs, cost, lower, upper, indices.Length, indices, values);
    }

    public Status addCols(double[] costs, double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
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

    public Status changeObjectiveSense(HighsObjectiveSense sense)
    {
        return (Status)Highs_changeObjectiveSense(_highs, (int)sense);
    }

    public Status changeColCost(int col, double cost)
    {
        return (Status)Highs_changeColCost(_highs, col, cost);
    }

    public Status changeColsCostBySet(int[] cols, double[] costs)
    {
        return (Status)Highs_changeColsCostBySet(_highs, cols.Length, cols, costs);
    }

    public Status changeColsCostByMask(bool[] mask, double[] cost)
    {
        return (Status)Highs_changeColsCostByMask(_highs, mask.Select(x => x ? 1 : 0).ToArray(), cost);
    }

    public Status changeColBounds(int col, double lower, double upper)
    {
        return (Status)Highs_changeColBounds(_highs, col, lower, upper);
    }

    public Status changeColsBoundsByRange(int from, int to, double[] lower, double[] upper)
    {
        return (Status)Highs_changeColsBoundsByRange(_highs, from, to, lower, upper);
    }

    public Status changeColsBoundsBySet(int[] cols, double[] lower, double[] upper)
    {
        return (Status)Highs_changeColsBoundsBySet(_highs, cols.Length, cols, lower, upper);
    }

    public Status changeColsBoundsByMask(bool[] mask, double[] lower, double[] upper)
    {
        return (Status)Highs_changeColsBoundsByMask(_highs, mask.Select(x => x ? 1 : 0).ToArray(), lower, upper);
    }

    public Status changeRowBounds(int row, double lower, double upper)
    {
        return (Status)Highs_changeRowBounds(_highs, row, lower, upper);
    }

    public Status changeRowsBoundsBySet(int[] rows, double[] lower, double[] upper)
    {
        return (Status)Highs_changeRowsBoundsBySet(_highs, rows.Length, rows, lower, upper);
    }

    public Status changeRowsBoundsByMask(bool[] mask, double[] lower, double[] upper)
    {
        return (Status)Highs_changeRowsBoundsByMask(_highs, mask.Select(x => x ? 1 : 0).ToArray(), lower, upper);
    }

    public Status changeColsIntegralityByRange(int from_col, int to_col, HighsIntegrality[] integrality)
    {
        return (Status)Highs_changeColsIntegralityByRange(_highs, from_col, to_col, Array.ConvertAll(integrality, item => (int)item));
    }

    public Status changeCoeff(int row, int col, double value)
    {
        return (Status)Highs_changeCoeff(_highs, row, col, value);
    }

    public Status deleteColsByRange(int from, int to)
    {
        return (Status)Highs_deleteColsByRange(_highs, from, to);
    }

    public Status deleteColsBySet(int[] cols)
    {
        return (Status)Highs_deleteColsBySet(_highs, cols.Length, cols);
    }

    public Status deleteColsByMask(bool[] mask)
    {
        return (Status)Highs_deleteColsByMask(_highs, mask.Select(x => x ? 1 : 0).ToArray());
    }

    public Status deleteRowsByRange(int from, int to)
    {
        return (Status)Highs_deleteRowsByRange(_highs, from, to);
    }

    public Status deleteRowsBySet(int[] rows)
    {
        return (Status)Highs_deleteRowsBySet(_highs, rows.Length, rows);
    }

    public Status deleteRowsByMask(bool[] mask)
    {
        return (Status)Highs_deleteRowsByMask(_highs, mask.Select(x => x ? 1 : 0).ToArray());
    }

    delegate int HighsGetInfoDelegate<TValue>(IntPtr _highs, string infoName, out TValue output);

    private TValue GetValueOrFallback<TValue>(HighsGetInfoDelegate<TValue> highsGetInfoDelegate, string infoName, TValue fallback)
    {
        try
        {
            var status = (Status)highsGetInfoDelegate(_highs, infoName, out var value);
            if (status != Status.kOk)
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
            MipGap = GetValueOrFallback(Highs_getDoubleInfoValue, "mip_gap", double.NaN),
            DualBound = GetValueOrFallback(Highs_getDoubleInfoValue, "mip_dual_bound", double.NaN),
            ObjectiveValue = GetValueOrFallback(Highs_getDoubleInfoValue, "objective_function_value", double.NaN),
            NodeCount = GetValueOrFallback(Highs_getInt64InfoValue, "mip_node_count", 0L),
            IpmIterationCount = GetValueOrFallback(Highs_getIntInfoValue, "ipm_iteration_count", 0),
            SimplexIterationCount = GetValueOrFallback(Highs_getIntInfoValue, "simplex_iteration_count", 0),
            PdlpIterationCount = GetValueOrFallback(Highs_getIntInfoValue, "pdlp_iteration_count", 0),
        };
        return info;
    }

    public Status setSolution(HighsSolution solution)
    {
        return (Status)Highs_setSolution(_highs, solution.colvalue, solution.coldual, solution.rowvalue, solution.rowdual);
    }

    public Status getBasicVariables(ref int[] basic_variables)
    {
        return (Status)Highs_getBasicVariables(_highs, basic_variables);
    }

    public Status getBasisInverseRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
    {
        return (Status)Highs_getBasisInverseRow(_highs, row, row_vector, ref row_num_nz, row_indices);
    }

    public Status getBasisInverseCol(int col, double[] col_vector, ref int col_num_nz, int[] col_indices)
    {
        return (Status)Highs_getBasisInverseCol(_highs, col, col_vector, ref col_num_nz, col_indices);
    }

    public Status getBasisSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
    {
        return (Status)Highs_getBasisSolve(_highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
    }

    public Status getBasisTransposeSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
    {
        return (Status)Highs_getBasisTransposeSolve(_highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
    }

    public Status getReducedRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
    {
        return (Status)Highs_getReducedRow(_highs, row, row_vector, ref row_num_nz, row_indices);
    }

    public Status getReducedColumn(int col, double[] col_vector, ref int col_num_nz, int[] col_indices)
    {
        return (Status)Highs_getReducedColumn(_highs, col, col_vector, ref col_num_nz, col_indices);
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
    /// <param name="col"></param>
    /// <param name="name"></param>
    /// <returns></returns>
    public Status PassColumnName(int col, string name) => (Status)Highs_passColName(_highs, col, name);

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
