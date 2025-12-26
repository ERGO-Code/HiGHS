using Highs.Enums;

namespace Highs.Records;

/// <summary>
/// This defines a model for Highs
/// </summary>
public record Model
{
    /// <summary>
    /// The objective sense
    /// </summary>
    public ObjectiveSense ObjectiveSense;
    /// <summary>
    /// The objective constant
    /// </summary>
    public double Offset;
    /// <summary>
    /// The column costs
    /// </summary>
    public double[] ColumnCost = [];
    /// <summary>
    /// The column lower bounds
    /// </summary>
    public double[] ColumnLower = [];
    /// <summary>
    /// The column upper bounds
    /// </summary>
    public double[] ColumnUpper = [];
    /// <summary>
    /// The row lower bounds
    /// </summary>
    public double[] RowLower = [];
    /// <summary>
    /// The row upper bounds
    /// </summary>
    public double[] RowUpper = [];
    /// <summary>
    /// The format of the constraint matrix.
    /// </summary>
    public MatrixFormat MatrixFormat;
    /// <summary>
    /// The starting index of each column (or row) in MatrixIndices.
    /// </summary>
    public int[] MatrixStart = [];
    /// <summary>
    /// The indices of matrix entries
    /// </summary>
    public int[] MatrixIndices = [];
    /// <summary>
    /// The values of matrix entries
    /// </summary>
    public double[] MatrixValues = [];
    /// <summary>
    /// The integrality of the variables
    /// </summary>
    public VariableType[] VariableTypes = [];
    /// <summary>
    /// The Hessian for quadratic objectives
    /// </summary>
    public Hessian Hessian = null;
}
