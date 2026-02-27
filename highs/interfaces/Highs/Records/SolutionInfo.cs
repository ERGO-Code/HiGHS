namespace Highs.Records;

/// <summary>
/// The solution info.
/// </summary>
/// <param name="SimplexIterationCount">The simplex iteration count.</param>
/// <param name="IpmIterationCount">The Interior Point Method (IPM) iteration count.</param>
/// <param name="PdlpIterationCount">The PDLP iteration count.</param>
/// <param name="MipGap">The MIP gap.</param>
/// <param name="DualBound">The best dual bound.</param>
/// <param name="NodeCount">The MIP node count.</param>
/// <param name="ObjectiveValue">The objective value.</param>
public record SolutionInfo(int SimplexIterationCount,
                           int IpmIterationCount,
                           int PdlpIterationCount,
                           double MipGap,
                           double DualBound,
                           long NodeCount,
                           double ObjectiveValue)
{
}
