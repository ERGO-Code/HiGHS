namespace Highs.Enums;

/// <summary>
/// This defines optimization sense of a HighsLp
/// </summary>
public enum ObjectiveSense
{
    /// <summary>
    /// The objective is to be minimized
    /// </summary>
    Minimize = 1,
    /// <summary>
    /// The objective is to be maximized
    /// </summary>
    Maximize = -1
}
