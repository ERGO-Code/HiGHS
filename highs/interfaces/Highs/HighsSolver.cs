using System.Text;
using Highs.Enums;
using Highs.Records;

namespace Highs;

/// <summary>
/// The Highs Solver interface.
/// </summary>
public class HighsSolver : IDisposable
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
    /// Calls the Highs solver in a single call of a LP.
    /// </summary>
    /// <param name="model"></param>
    /// <param name="solution"></param>
    /// <param name="basisInfo"></param>
    /// <param name="modelStatus"></param>
    /// <returns></returns>
    public static HighsStatus LpCall(HighsModel model, ref HighsSolution solution, out BasisInfo basisInfo, out ModelStatus modelStatus)
    {
        var numberOfColumns = model.ColumnCost.Length;
        var numberOfRows = model.RowLower.Length;
        var numberOfMatrixValues = model.MatrixValues.Length;

        var columnBasisStatus = new int[numberOfColumns];
        var rowBasisStatus = new int[numberOfRows];

        var modelstate = 0;

        var status = (HighsStatus)Imports.Highs_lpCall(
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
    public static HighsStatus MipCall(HighsModel model, ref HighsSolution solution, out ModelStatus modelStatus)
    {
        var numberOfColumns = model.ColumnCost.Length;
        var numberOfRows = model.RowLower.Length;
        var numberOfMatrixValues = model.MatrixValues.Length;

        var modelstate = 0;

        var status = (HighsStatus)Imports.Highs_mipCall(
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
    public static HighsStatus QPCall(HighsModel model, ref HighsSolution solution, out BasisInfo basisInfo, out ModelStatus modelStatus)
    {
        var numberOfColumns = model.ColumnCost.Length;
        var numberOfRows = model.RowLower.Length;
        var numberOfMatrixValues = model.MatrixValues.Length;
        var numberOfHessianValues = model.Hessian.Values.Length;

        var columnBasisStatus = new int[numberOfColumns];
        var rowBasisStatus = new int[numberOfRows];

        var modelstate = 0;

        var status = (HighsStatus)Imports.Highs_qpCall(
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
    public HighsSolver() => _highs = Imports.Highs_create();

    /// <summary>
    /// The destructor.
    /// </summary>
    ~HighsSolver() => Dispose(false);

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
        Imports.Highs_destroy(_highs);
        _disposed = true;
    }

    /// <summary>
    /// Runs the solver.
    /// </summary>
    /// <returns></returns>
    public HighsStatus Run() => (HighsStatus)Imports.Highs_run(_highs);

    /// <summary>
    /// Reads a model from file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public HighsStatus ReadModel(string filename) => (HighsStatus)Imports.Highs_readModel(_highs, filename);

    /// <summary>
    /// Writes the model to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public HighsStatus WriteModel(string filename) => (HighsStatus)Imports.Highs_writeModel(_highs, filename);

    /// <summary>
    /// Writes the presolved model to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public HighsStatus WritePresolvedModel(string filename) => (HighsStatus)Imports.Highs_writePresolvedModel(_highs, filename);

    /// <summary>
    /// Writes the solution to file in a pretty format.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public HighsStatus WriteSolutionPretty(string filename) => (HighsStatus)Imports.Highs_writeSolutionPretty(_highs, filename);

    /// <summary>
    /// Gets the infinity value.
    /// </summary>
    /// <returns></returns>
    public double GetInfinity() => Imports.Highs_getInfinity(_highs);

    /// <summary>
    /// Passes the LP model to Highs.
    /// </summary>
    /// <param name="model"></param>
    /// <returns></returns>
    public HighsStatus PassLp(HighsModel model)
    {
        return (HighsStatus)Imports.Highs_passLp(
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
    public HighsStatus PassMip(HighsModel model)
    {
        return (HighsStatus)Imports.Highs_passMip(
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
    public HighsStatus SetOptionValue(string option, string value) => (HighsStatus)Imports.Highs_setOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the string option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public HighsStatus SetStringOptionValue(string option, string value) => (HighsStatus)Imports.Highs_setStringOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the boolean option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public HighsStatus SetBoolOptionValue(string option, int value) => (HighsStatus)Imports.Highs_setBoolOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the double option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public HighsStatus SetDoubleOptionValue(string option, double value) => (HighsStatus)Imports.Highs_setDoubleOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the integer option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public HighsStatus SetIntOptionValue(string option, int value) => (HighsStatus)Imports.Highs_setIntOptionValue(_highs, option, value);

    /// <summary>
    /// Gets the string option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public HighsStatus GetStringOptionValue(string option, out string value)
    {
        var stringBuilder = new StringBuilder();
        var result = (HighsStatus)Imports.Highs_getStringOptionValue(_highs, option, stringBuilder);
        value = stringBuilder.ToString();
        return result;
    }

    /// <summary>
    /// Gets the boolean option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public HighsStatus GetBoolOptionValue(string option, out int value) => (HighsStatus)Imports.Highs_getBoolOptionValue(_highs, option, out value);

    /// <summary>
    /// Gets the double option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public HighsStatus GetDoubleOptionValue(string option, out double value) => (HighsStatus)Imports.Highs_getDoubleOptionValue(_highs, option, out value);

    /// <summary>
    /// Gets the integer option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public HighsStatus GetIntOptionValue(string option, out int value) => (HighsStatus)Imports.Highs_getIntOptionValue(_highs, option, out value);

    /// <summary>
    /// Gets the number of columns.
    /// </summary>
    /// <returns></returns>
    public int GetNumberOfColumns() => Imports.Highs_getNumCol(_highs);

    /// <summary>
    /// Gets the number of rows.
    /// </summary>
    /// <returns></returns>
    public int GetNumberOfRows() => Imports.Highs_getNumRow(_highs);

    /// <summary>
    /// Gets the number of non-zero entries.
    /// </summary>
    /// <returns></returns>
    public int GetNumberOfNonZeroEntries() => Imports.Highs_getNumNz(_highs);

    /// <summary>
    /// Gets the solution.
    /// </summary>
    /// <param name="solution"></param>
    /// <returns></returns>
    public HighsStatus GetSolution(out HighsSolution solution)
    {
        var numberOfColumns = GetNumberOfColumns();
        var numberOfRows = GetNumberOfRows();

        solution = new HighsSolution(numberOfColumns, numberOfRows);
        return (HighsStatus)Imports.Highs_getSolution(_highs, solution.ColumnValue, solution.ColumnDual, solution.RowValue, solution.RowDual);
    }

    /// <summary>
    /// Passes the Hessian to Highs.
    /// </summary>
    /// <param name="hessian"></param>
    /// <returns></returns>
    public HighsStatus PassHessian(Hessian hessian)
    {
        return (HighsStatus)Imports.Highs_passHessian(_highs,
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
    public HighsStatus GetBasis(out BasisInfo basisInfo)
    {
        var numberOfColumns = GetNumberOfColumns();
        var numberOfRows = GetNumberOfRows();

        var columnBasisStatus = new int[numberOfColumns];
        var rowBasisStatus = new int[numberOfRows];

        var status = (HighsStatus)Imports.Highs_getBasis(_highs, columnBasisStatus, rowBasisStatus);
        if (status == HighsStatus.Error)
        {
            basisInfo = null;
            return HighsStatus.Error;
        }

        basisInfo = new BasisInfo(Array.ConvertAll(columnBasisStatus, x => (BasisStatus)x), Array.ConvertAll(rowBasisStatus, x => (BasisStatus)x));

        return status;
    }

    /// <summary>
    /// Gets the objective value.
    /// </summary>
    /// <returns></returns>
    public double GetObjectiveValue() => Imports.Highs_getObjectiveValue(_highs);

    /// <summary>
    /// Gets the model status.
    /// </summary>
    /// <returns></returns>
    public ModelStatus GetModelStatus() => (ModelStatus)Imports.Highs_getModelStatus(_highs);

    /// <summary>
    /// Gets the iteration count.
    /// </summary>
    /// <returns></returns>
    public int GetIterationCount() => Imports.Highs_getIterationCount(_highs);

    /// <summary>
    /// Adds a row to the model.
    /// </summary>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <param name="indices"></param>
    /// <param name="values"></param>
    /// <returns></returns>
    public HighsStatus AddRow(double lower, double upper, int[] indices, double[] values)
    {
        return (HighsStatus)Imports.Highs_addRow(_highs, lower, upper, indices.Length, indices, values);
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
    public HighsStatus AddRows(double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
    {
        return (HighsStatus)Imports.Highs_addRows(_highs, lower.Length, lower, upper, indices.Length, starts, indices, values);
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
    public HighsStatus AddColumn(double cost, double lower, double upper, int[] indices, double[] values)
    {
        return (HighsStatus)Imports.Highs_addCol(_highs, cost, lower, upper, indices.Length, indices, values);
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
    public HighsStatus AddColumns(double[] costs, double[] lower, double[] upper, int[] starts, int[] indices, double[] values)
    {
        return (HighsStatus)Imports.Highs_addCols(
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
    public HighsStatus ChangeObjectiveSense(ObjectiveSense sense) => (HighsStatus)Imports.Highs_changeObjectiveSense(_highs, (int)sense);

    /// <summary>
    /// Changes the cost of a column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="cost"></param>
    /// <returns></returns>
    public HighsStatus ChangeColumnCost(int column, double cost) => (HighsStatus)Imports.Highs_changeColCost(_highs, column, cost);

    /// <summary>
    /// Changes the costs of multiple columns by set.
    /// </summary>
    /// <param name="columns"></param>
    /// <param name="costs"></param>
    /// <returns></returns>
    public HighsStatus ChangeColumnsCostBySet(int[] columns, double[] costs) => (HighsStatus)Imports.Highs_changeColsCostBySet(_highs, columns.Length, columns, costs);

    /// <summary>
    /// Changes the costs of multiple columns by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <param name="cost"></param>
    /// <returns></returns>
    public HighsStatus ChangeColumnsCostByMask(bool[] mask, double[] cost) => (HighsStatus)Imports.Highs_changeColsCostByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0), cost);

    /// <summary>
    /// Changes the bounds of a column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public HighsStatus ChangeColumnBounds(int column, double lower, double upper) => (HighsStatus)Imports.Highs_changeColBounds(_highs, column, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple columns by range.
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public HighsStatus ChangeColumnsBoundsByRange(int from, int to, double[] lower, double[] upper) => (HighsStatus)Imports.Highs_changeColsBoundsByRange(_highs, from, to, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple columns by set.
    /// </summary>
    /// <param name="columns"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public HighsStatus ChangeColumnsBoundsBySet(int[] columns, double[] lower, double[] upper)
    {
        return (HighsStatus)Imports.Highs_changeColsBoundsBySet(_highs, columns.Length, columns, lower, upper);
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="mask"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public HighsStatus ChangeColumnsBoundsByMask(bool[] mask, double[] lower, double[] upper)
    {
        return (HighsStatus)Imports.Highs_changeColsBoundsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0), lower, upper);
    }

    /// <summary>
    /// Changes the bounds of a row.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public HighsStatus ChangeRowBounds(int row, double lower, double upper) => (HighsStatus)Imports.Highs_changeRowBounds(_highs, row, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple rows by range.
    /// </summary>
    /// <param name="rows"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public HighsStatus ChangeRowsBoundsBySet(int[] rows, double[] lower, double[] upper) => (HighsStatus)Imports.Highs_changeRowsBoundsBySet(_highs, rows.Length, rows, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple rows by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public HighsStatus ChangeRowsBoundsByMask(bool[] mask, double[] lower, double[] upper) => (HighsStatus)Imports.Highs_changeRowsBoundsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0), lower, upper);

    /// <summary>
    /// Changes the integrality of multiple columns by range.
    /// </summary>
    /// <param name="from_col"></param>
    /// <param name="to_col"></param>
    /// <param name="integrality"></param>
    /// <returns></returns>
    public HighsStatus ChangeColumnsIntegralityByRange(int from_col, int to_col, VariableType[] integrality)
    {
        return (HighsStatus)Imports.Highs_changeColsIntegralityByRange(_highs, from_col, to_col, Array.ConvertAll(integrality, item => (int)item));
    }

    /// <summary>
    /// Changes a coefficient in the constraint matrix.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="column"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public HighsStatus ChangeCoefficient(int row, int column, double value) => (HighsStatus)Imports.Highs_changeCoeff(_highs, row, column, value);

    /// <summary>
    /// Deletes multiple columns by range.
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <returns></returns>
    public HighsStatus DeleteColumnsByRange(int from, int to) => (HighsStatus)Imports.Highs_deleteColsByRange(_highs, from, to);

    /// <summary>
    /// Deletes multiple columns by set.
    /// </summary>
    /// <param name="columns"></param>
    /// <returns></returns>
    public HighsStatus DeleteColumnsBySet(int[] columns) => (HighsStatus)Imports.Highs_deleteColsBySet(_highs, columns.Length, columns);

    /// <summary>
    /// Deletes multiple columns by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <returns></returns>
    public HighsStatus DeleteColumnsByMask(bool[] mask) => (HighsStatus)Imports.Highs_deleteColsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0));

    /// <summary>
    /// Deletes multiple rows by range.
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <returns></returns>
    public HighsStatus DeleteRowsByRange(int from, int to) => (HighsStatus)Imports.Highs_deleteRowsByRange(_highs, from, to);

    /// <summary>
    /// Deletes multiple rows by set.
    /// </summary>
    /// <param name="rows"></param>
    /// <returns></returns>
    public HighsStatus DeleteRowsBySet(int[] rows) => (HighsStatus)Imports.Highs_deleteRowsBySet(_highs, rows.Length, rows);

    /// <summary>
    /// Deletes multiple rows by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <returns></returns>
    public HighsStatus DeleteRowsByMask(bool[] mask) => (HighsStatus)Imports.Highs_deleteRowsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0));

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
            var status = (HighsStatus)highsGetInfoDelegate(_highs, infoName, out var value);
            return status != HighsStatus.Ok ? fallback : value;
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
        return new SolutionInfo(GetValueOrFallback(Imports.Highs_getIntInfoValue, "simplex_iteration_count", 0),
                                GetValueOrFallback(Imports.Highs_getIntInfoValue, "ipm_iteration_count", 0),
                                GetValueOrFallback(Imports.Highs_getIntInfoValue, "pdlp_iteration_count", 0),
                                GetValueOrFallback(Imports.Highs_getDoubleInfoValue, "mip_gap", double.NaN),
                                GetValueOrFallback(Imports.Highs_getDoubleInfoValue, "mip_dual_bound", double.NaN),
                                GetValueOrFallback(Imports.Highs_getInt64InfoValue, "mip_node_count", 0L),
                                GetValueOrFallback(Imports.Highs_getDoubleInfoValue, "objective_function_value", double.NaN));
    }

    /// <summary>
    /// Sets the solution.
    /// </summary>
    /// <param name="solution"></param>
    /// <returns></returns>
    public HighsStatus SetSolution(HighsSolution solution) => (HighsStatus)Imports.Highs_setSolution(_highs, solution.ColumnValue, solution.ColumnDual, solution.RowValue, solution.RowDual);

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
    public HighsStatus SetSparseSolution(IReadOnlyDictionary<int, double> valuesByIndex)
    {
        return (HighsStatus)Imports.Highs_setSparseSolution(_highs, valuesByIndex.Count, [.. valuesByIndex.Keys], [.. valuesByIndex.Values]);
    }

    /// <summary>
    /// Gets the basic variables.
    /// </summary>
    /// <param name="basic_variables"></param>
    /// <returns></returns>
    public HighsStatus GetBasicVariables(ref int[] basic_variables) => (HighsStatus)Imports.Highs_getBasicVariables(_highs, basic_variables);

    /// <summary>
    /// Gets a row of the basis inverse.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="row_vector"></param>
    /// <param name="row_num_nz"></param>
    /// <param name="row_indices"></param>
    /// <returns></returns>
    public HighsStatus GetBasisInverseRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
    {
        return (HighsStatus)Imports.Highs_getBasisInverseRow(_highs, row, row_vector, ref row_num_nz, row_indices);
    }

    /// <summary>
    /// Gets a column of the basis inverse.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="column_vector"></param>
    /// <param name="column_num_nz"></param>
    /// <param name="column_indices"></param>
    /// <returns></returns>
    public HighsStatus GetBasisInverseColumn(int column, double[] column_vector, ref int column_num_nz, int[] column_indices)
    {
        return (HighsStatus)Imports.Highs_getBasisInverseCol(_highs, column, column_vector, ref column_num_nz, column_indices);
    }

    /// <summary>
    /// Gets the basis solve.
    /// </summary>
    /// <param name="rhs"></param>
    /// <param name="solution_vector"></param>
    /// <param name="solution_num_nz"></param>
    /// <param name="solution_indices"></param>
    /// <returns></returns>
    public HighsStatus GetBasisSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
    {
        return (HighsStatus)Imports.Highs_getBasisSolve(_highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
    }

    /// <summary>
    /// Gets the basis transpose solve.
    /// </summary>
    /// <param name="rhs"></param>
    /// <param name="solution_vector"></param>
    /// <param name="solution_num_nz"></param>
    /// <param name="solution_indices"></param>
    /// <returns></returns>
    public HighsStatus GetBasisTransposeSolve(double[] rhs, double[] solution_vector, ref int solution_num_nz, int[] solution_indices)
    {
        return (HighsStatus)Imports.Highs_getBasisTransposeSolve(_highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
    }

    /// <summary>
    /// Gets a reduced row.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="row_vector"></param>
    /// <param name="row_num_nz"></param>
    /// <param name="row_indices"></param>
    /// <returns></returns>
    public HighsStatus GetReducedRow(int row, double[] row_vector, ref int row_num_nz, int[] row_indices)
    {
        return (HighsStatus)Imports.Highs_getReducedRow(_highs, row, row_vector, ref row_num_nz, row_indices);
    }

    /// <summary>
    /// Gets a reduced column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="column_vector"></param>
    /// <param name="column_num_nz"></param>
    /// <param name="column_indices"></param>
    /// <returns></returns>
    public HighsStatus GetReducedColumn(int column, double[] column_vector, ref int column_num_nz, int[] column_indices)
    {
        return (HighsStatus)Imports.Highs_getReducedColumn(_highs, column, column_vector, ref column_num_nz, column_indices);
    }

    /// <summary>
    /// Clears the model.
    /// </summary>
    /// <returns></returns>
    public HighsStatus ClearModel() => (HighsStatus)Imports.Highs_clearModel(_highs);

    /// <summary>
    /// Clears the solver.
    /// </summary>
    /// <returns></returns>
    public HighsStatus ClearSolver() => (HighsStatus)Imports.Highs_clearSolver(_highs);

    /// <summary>
    /// Passes the name of a column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="name"></param>
    /// <returns></returns>
    public HighsStatus PassColumnName(int column, string name) => (HighsStatus)Imports.Highs_passColName(_highs, column, name);

    /// <summary>
    /// Passes the name of a row.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="name"></param>
    /// <returns></returns>
    public HighsStatus PassRowName(int row, string name) => (HighsStatus)Imports.Highs_passRowName(_highs, row, name);

    /// <summary>
    /// Writes the options to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public HighsStatus WriteOptions(string filename) => (HighsStatus)Imports.Highs_writeOptions(_highs, filename);

    /// <summary>
    /// Writes the options deviations to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public HighsStatus WriteOptionsDeviations(string filename) => (HighsStatus)Imports.Highs_writeOptionsDeviations(_highs, filename);
}
