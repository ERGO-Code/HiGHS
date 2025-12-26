namespace Highs.Records;

/// <summary>
/// The solution.
/// </summary>
/// <param name="ColumnValue">The column value.</param>
/// <param name="ColumnDual">The column dual.</param>
/// <param name="RowValue">The row value.</param>
/// <param name="RowDual">The row dual.</param>
public record HighsSolution(double[] ColumnValue,
                       double[] ColumnDual,
                       double[] RowValue,
                       double[] RowDual)
{
    /// <summary>
    /// The default constructor creates empty arrays
    /// </summary>
    /// <param name="numberOfColumns">The number of columns</param>
    /// <param name="numberOfRows">The number of rows</param>
    public HighsSolution(int numberOfColumns, int numberOfRows) : this(new double[numberOfColumns],
                                                                  new double[numberOfColumns],
                                                                  new double[numberOfRows],
                                                                  new double[numberOfRows])
    {
    }
}
