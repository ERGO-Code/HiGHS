namespace Highs.Enums;

/// <summary>
/// This defines the feasible values of a variable within a model
/// </summary>
public enum VariableType
{
    /// <summary>
    /// The variable can take continuous values between its bounds
    /// </summary>
    Continuous = 0,
    /// <summary>
    /// The variable must take integer values between its bounds
    /// </summary>
    Integer,
    /// <summary>
    /// The variable must be zero or take continuous values between its bounds
    /// </summary>
    SemiContinuous,
    /// <summary>
    /// The variable must be zero or take integer values between its bounds
    /// </summary>
    SemiInteger,
    /// <summary>
    /// The variable must take implicit integer values between its bounds
    /// </summary>
    ImplicitInteger,
}