namespace Highs.Enums;

/// <summary>
/// This defines the status of a variable (or slack variable for a constraint) in a basis
/// </summary>
public enum BasisStatus
{
    /// <summary>
    /// The variable is nonbasic at its lower bound (or fixed value)
    /// </summary>
    Lower = 0,
    /// <summary>
    /// The variable is basic
    /// </summary>
    Basic,
    /// <summary>
    /// he variable is at its upper bound
    /// </summary>
    Upper,
    /// <summary>
    /// A free variable is nonbasic and set to zero
    /// </summary>
    Zero,
    /// <summary>
    /// The variable is nonbasic
    /// </summary>
    Nonbasic
}
