using Highs.Enums;

namespace Highs.Records;

/// <summary>
/// This defines the Hessian of the quadratic objective
/// </summary>
public record Hessian
{
    /// <summary>
    /// Format of the Hessian
    /// </summary>
    public HessianFormat HessianFormat;
    /// <summary>
    /// Dimension of the Hessian
    /// </summary>
    public int Dimension;
    /// <summary>
    /// Start of each compressed column in the Hessian
    /// </summary>
    public int[] Start = [];
    /// <summary>
    /// Indices of the nonzeros in the Hessian
    /// </summary>
    public int[] Index = [];
    /// <summary>
    /// Values of the nonzeros in the Hessian
    /// </summary>
    public double[] Values = [];
}
