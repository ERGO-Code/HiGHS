namespace Highs.Enums;

/// <summary>
/// This is (part of) the return value of most HiGHS methods
/// </summary>
public enum Status
{
    /// <summary>
    /// The method has exposed an error
    /// </summary>
    Error = -1,
    /// <summary>
    /// The method has completed successfully
    /// </summary>
    Ok,
    /// <summary>
    /// The method has recovered from an unusual event, or has terminated due to reaching a time or iteration limit
    /// </summary>
    Warning
}
