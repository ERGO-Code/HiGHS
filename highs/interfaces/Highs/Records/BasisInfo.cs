using Highs.Enums;

namespace Highs.Records;

/// <summary>
/// This defines the basis status of the columns and rows in a basis
/// </summary>
/// <param name="ColumnBasisStatus">The column basis status</param>
/// <param name="RowBasisStatus">The row basis status</param>
public record BasisInfo(BasisStatus[] ColumnBasisStatus, BasisStatus[] RowBasisStatus)
{
    /// <summary>
    /// The default constructor creates empty arrays
    /// </summary>
    /// <param name="numberOfColumns">The number of columns</param>
    /// <param name="numberOfRows">The number of rows</param>
    public BasisInfo(int numberOfColumns, int numberOfRows) : this(new BasisStatus[numberOfColumns],
                                                                   new BasisStatus[numberOfRows])
    {
    }
}
