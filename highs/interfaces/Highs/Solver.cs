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

        var status = (Status)Imports.Highs_lpCall(
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

        var status = (Status)Imports.Highs_mipCall(
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

        var status = (Status)Imports.Highs_qpCall(
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
    public Solver() => _highs = Imports.Highs_create();

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
        Imports.Highs_destroy(_highs);
        _disposed = true;
    }

    /// <summary>
    /// Runs the solver.
    /// </summary>
    /// <returns></returns>
    public Status Run() => (Status)Imports.Highs_run(_highs);

    /// <summary>
    /// Reads a model from file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status ReadModel(string filename) => (Status)Imports.Highs_readModel(_highs, filename);

    /// <summary>
    /// Writes the model to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WriteModel(string filename) => (Status)Imports.Highs_writeModel(_highs, filename);

    /// <summary>
    /// Writes the presolved model to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WritePresolvedModel(string filename) => (Status)Imports.Highs_writePresolvedModel(_highs, filename);

    /// <summary>
    /// Writes the solution to file in a pretty format.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WriteSolutionPretty(string filename) => (Status)Imports.Highs_writeSolutionPretty(_highs, filename);

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
    public Status PassLp(Model model)
    {
        return (Status)Imports.Highs_passLp(
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
        return (Status)Imports.Highs_passMip(
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
    public Status SetOptionValue(string option, string value) => (Status)Imports.Highs_setOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the string option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status SetStringOptionValue(string option, string value) => (Status)Imports.Highs_setStringOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the boolean option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status SetBoolOptionValue(string option, int value) => (Status)Imports.Highs_setBoolOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the double option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status SetDoubleOptionValue(string option, double value) => (Status)Imports.Highs_setDoubleOptionValue(_highs, option, value);

    /// <summary>
    /// Sets the integer option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status SetIntOptionValue(string option, int value) => (Status)Imports.Highs_setIntOptionValue(_highs, option, value);

    /// <summary>
    /// Gets the string option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status GetStringOptionValue(string option, out string value)
    {
        var stringBuilder = new StringBuilder();
        var result = (Status)Imports.Highs_getStringOptionValue(_highs, option, stringBuilder);
        value = stringBuilder.ToString();
        return result;
    }

    /// <summary>
    /// Gets the boolean option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status GetBoolOptionValue(string option, out int value) => (Status)Imports.Highs_getBoolOptionValue(_highs, option, out value);

    /// <summary>
    /// Gets the double option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status GetDoubleOptionValue(string option, out double value) => (Status)Imports.Highs_getDoubleOptionValue(_highs, option, out value);

    /// <summary>
    /// Gets the integer option value.
    /// </summary>
    /// <param name="option"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status GetIntOptionValue(string option, out int value) => (Status)Imports.Highs_getIntOptionValue(_highs, option, out value);

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
    public Status GetSolution(out Solution solution)
    {
        var numberOfColumns = GetNumberOfColumns();
        var numberOfRows = GetNumberOfRows();

        solution = new Solution(numberOfColumns, numberOfRows);
        return (Status)Imports.Highs_getSolution(_highs, solution.ColumnValue, solution.ColumnDual, solution.RowValue, solution.RowDual);
    }

    /// <summary>
    /// Passes the Hessian to Highs.
    /// </summary>
    /// <param name="hessian"></param>
    /// <returns></returns>
    public Status PassHessian(Hessian hessian)
    {
        return (Status)Imports.Highs_passHessian(_highs,
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

        var status = (Status)Imports.Highs_getBasis(_highs, columnBasisStatus, rowBasisStatus);
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
    public Status AddRow(double lower, double upper, int[] indices, double[] values)
    {
        return (Status)Imports.Highs_addRow(_highs, lower, upper, indices.Length, indices, values);
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
        return (Status)Imports.Highs_addRows(_highs, lower.Length, lower, upper, indices.Length, starts, indices, values);
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
        return (Status)Imports.Highs_addCol(_highs, cost, lower, upper, indices.Length, indices, values);
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
        return (Status)Imports.Highs_addCols(
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
    public Status ChangeObjectiveSense(ObjectiveSense sense) => (Status)Imports.Highs_changeObjectiveSense(_highs, (int)sense);

    /// <summary>
    /// Changes the cost of a column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="cost"></param>
    /// <returns></returns>
    public Status ChangeColumnCost(int column, double cost) => (Status)Imports.Highs_changeColCost(_highs, column, cost);

    /// <summary>
    /// Changes the costs of multiple columns by set.
    /// </summary>
    /// <param name="columns"></param>
    /// <param name="costs"></param>
    /// <returns></returns>
    public Status ChangeColumnsCostBySet(int[] columns, double[] costs) => (Status)Imports.Highs_changeColsCostBySet(_highs, columns.Length, columns, costs);

    /// <summary>
    /// Changes the costs of multiple columns by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <param name="cost"></param>
    /// <returns></returns>
    public Status ChangeColumnsCostByMask(bool[] mask, double[] cost) => (Status)Imports.Highs_changeColsCostByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0), cost);

    /// <summary>
    /// Changes the bounds of a column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeColumnBounds(int column, double lower, double upper) => (Status)Imports.Highs_changeColBounds(_highs, column, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple columns by range.
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeColumnsBoundsByRange(int from, int to, double[] lower, double[] upper) => (Status)Imports.Highs_changeColsBoundsByRange(_highs, from, to, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple columns by set.
    /// </summary>
    /// <param name="columns"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeColumnsBoundsBySet(int[] columns, double[] lower, double[] upper)
    {
        return (Status)Imports.Highs_changeColsBoundsBySet(_highs, columns.Length, columns, lower, upper);
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
        return (Status)Imports.Highs_changeColsBoundsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0), lower, upper);
    }

    /// <summary>
    /// Changes the bounds of a row.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeRowBounds(int row, double lower, double upper) => (Status)Imports.Highs_changeRowBounds(_highs, row, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple rows by range.
    /// </summary>
    /// <param name="rows"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeRowsBoundsBySet(int[] rows, double[] lower, double[] upper) => (Status)Imports.Highs_changeRowsBoundsBySet(_highs, rows.Length, rows, lower, upper);

    /// <summary>
    /// Changes the bounds of multiple rows by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <param name="lower"></param>
    /// <param name="upper"></param>
    /// <returns></returns>
    public Status ChangeRowsBoundsByMask(bool[] mask, double[] lower, double[] upper) => (Status)Imports.Highs_changeRowsBoundsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0), lower, upper);

    /// <summary>
    /// Changes the integrality of multiple columns by range.
    /// </summary>
    /// <param name="from_col"></param>
    /// <param name="to_col"></param>
    /// <param name="integrality"></param>
    /// <returns></returns>
    public Status ChangeColumnsIntegralityByRange(int from_col, int to_col, VariableType[] integrality)
    {
        return (Status)Imports.Highs_changeColsIntegralityByRange(_highs, from_col, to_col, Array.ConvertAll(integrality, item => (int)item));
    }

    /// <summary>
    /// Changes a coefficient in the constraint matrix.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="column"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public Status ChangeCoefficient(int row, int column, double value) => (Status)Imports.Highs_changeCoeff(_highs, row, column, value);

    /// <summary>
    /// Deletes multiple columns by range.
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <returns></returns>
    public Status DeleteColumnsByRange(int from, int to) => (Status)Imports.Highs_deleteColsByRange(_highs, from, to);

    /// <summary>
    /// Deletes multiple columns by set.
    /// </summary>
    /// <param name="columns"></param>
    /// <returns></returns>
    public Status DeleteColumnsBySet(int[] columns) => (Status)Imports.Highs_deleteColsBySet(_highs, columns.Length, columns);

    /// <summary>
    /// Deletes multiple columns by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <returns></returns>
    public Status DeleteColumnsByMask(bool[] mask) => (Status)Imports.Highs_deleteColsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0));

    /// <summary>
    /// Deletes multiple rows by range.
    /// </summary>
    /// <param name="from"></param>
    /// <param name="to"></param>
    /// <returns></returns>
    public Status DeleteRowsByRange(int from, int to) => (Status)Imports.Highs_deleteRowsByRange(_highs, from, to);

    /// <summary>
    /// Deletes multiple rows by set.
    /// </summary>
    /// <param name="rows"></param>
    /// <returns></returns>
    public Status DeleteRowsBySet(int[] rows) => (Status)Imports.Highs_deleteRowsBySet(_highs, rows.Length, rows);

    /// <summary>
    /// Deletes multiple rows by mask.
    /// </summary>
    /// <param name="mask"></param>
    /// <returns></returns>
    public Status DeleteRowsByMask(bool[] mask) => (Status)Imports.Highs_deleteRowsByMask(_highs, Array.ConvertAll(mask, x => x ? 1 : 0));

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
    public Status SetSolution(Solution solution) => (Status)Imports.Highs_setSolution(_highs, solution.ColumnValue, solution.ColumnDual, solution.RowValue, solution.RowDual);

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
        return (Status)Imports.Highs_setSparseSolution(_highs, valuesByIndex.Count, [.. valuesByIndex.Keys], [.. valuesByIndex.Values]);
    }

    /// <summary>
    /// Gets the basic variables.
    /// </summary>
    /// <param name="basic_variables"></param>
    /// <returns></returns>
    public Status GetBasicVariables(ref int[] basic_variables) => (Status)Imports.Highs_getBasicVariables(_highs, basic_variables);

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
        return (Status)Imports.Highs_getBasisInverseRow(_highs, row, row_vector, ref row_num_nz, row_indices);
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
        return (Status)Imports.Highs_getBasisInverseCol(_highs, column, column_vector, ref column_num_nz, column_indices);
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
        return (Status)Imports.Highs_getBasisSolve(_highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
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
        return (Status)Imports.Highs_getBasisTransposeSolve(_highs, rhs, solution_vector, ref solution_num_nz, solution_indices);
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
        return (Status)Imports.Highs_getReducedRow(_highs, row, row_vector, ref row_num_nz, row_indices);
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
        return (Status)Imports.Highs_getReducedColumn(_highs, column, column_vector, ref column_num_nz, column_indices);
    }

    /// <summary>
    /// Clears the model.
    /// </summary>
    /// <returns></returns>
    public Status ClearModel() => (Status)Imports.Highs_clearModel(_highs);

    /// <summary>
    /// Clears the solver.
    /// </summary>
    /// <returns></returns>
    public Status ClearSolver() => (Status)Imports.Highs_clearSolver(_highs);

    /// <summary>
    /// Passes the name of a column.
    /// </summary>
    /// <param name="column"></param>
    /// <param name="name"></param>
    /// <returns></returns>
    public Status PassColumnName(int column, string name) => (Status)Imports.Highs_passColName(_highs, column, name);

    /// <summary>
    /// Passes the name of a row.
    /// </summary>
    /// <param name="row"></param>
    /// <param name="name"></param>
    /// <returns></returns>
    public Status PassRowName(int row, string name) => (Status)Imports.Highs_passRowName(_highs, row, name);

    /// <summary>
    /// Writes the options to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WriteOptions(string filename) => (Status)Imports.Highs_writeOptions(_highs, filename);

    /// <summary>
    /// Writes the options deviations to file.
    /// </summary>
    /// <param name="filename"></param>
    /// <returns></returns>
    public Status WriteOptionsDeviations(string filename) => (Status)Imports.Highs_writeOptionsDeviations(_highs, filename);
}
