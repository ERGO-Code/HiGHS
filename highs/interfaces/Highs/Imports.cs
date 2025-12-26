using System.Runtime.InteropServices;
using System.Text;

namespace Highs;

/// <summary>
/// This contains the imports for the Highs library
/// </summary>
public static class Imports
{
    /// <summary>
    /// The name of the Highs library.
    /// </summary>
    private const string HighsLibName = "highs";

    #region Library Imports
    [LibraryImport(HighsLibName)]
    internal static extern int Highs_lpCall(
        int numcol,
        int numrow,
        int numnz,
        int aformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        double[] colvalue,
        double[] coldual,
        double[] rowvalue,
        double[] rowdual,
        int[] colbasisstatus,
        int[] rowbasisstatus,
        ref int modelstatus);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_mipCall(
        int numcol,
        int numrow,
        int numnz,
        int aformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        int[] integrality,
        double[] colvalue,
        double[] rowvalue,
        ref int modelstatus);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_qpCall(
    int numcol,
    int numrow,
    int numnz,
    int qnumnz,
    int aformat,
    int qformat,
    int sense,
    double offset,
    double[] colcost,
    double[] collower,
    double[] colupper,
    double[] rowlower,
    double[] rowupper,
    int[] astart,
    int[] aindex,
    double[] avalue,
    int[] qstart,
    int[] qindex,
    double[] qvalue,

    double[] colvalue,
    double[] coldual,
    double[] rowvalue,
    double[] rowdual,
    int[] colbasisstatus,
    int[] rowbasisstatus,
    ref int modelstatus);

    [LibraryImport(HighsLibName)]
    internal static extern IntPtr Highs_create();

    [LibraryImport(HighsLibName)]
    internal static extern void Highs_destroy(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_run(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_readModel(IntPtr highs, string filename);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_writeModel(IntPtr highs, string filename);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_writePresolvedModel(IntPtr highs, string filename);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_writeSolutionPretty(IntPtr highs, string filename);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getInfinity(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_passLp(
        IntPtr highs,
        int numcol,
        int numrow,
        int numnz,
        int aformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_passMip(
        IntPtr highs,
        int numcol,
        int numrow,
        int numnz,
        int aformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        int[] highs_integrality);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_passModel(
        IntPtr highs,
        int numcol,
        int numrow,
        int numnz,
        int qnumnz,
        int aformat,
        int qformat,
        int sense,
        double offset,
        double[] colcost,
        double[] collower,
        double[] colupper,
        double[] rowlower,
        double[] rowupper,
        int[] astart,
        int[] aindex,
        double[] avalue,
        int[] qstart,
        int[] qindex,
        double[] qvalue,
        int[] highs_integrality);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_passHessian(
        IntPtr highs,
        int dim,
        int numnz,
        int q_format,
        int[] qstart,
        int[] qindex,
        double[] qvalue);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_setOptionValue(IntPtr highs, string option, string value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_setBoolOptionValue(IntPtr highs, string option, int value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_setIntOptionValue(IntPtr highs, string option, int value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_setDoubleOptionValue(IntPtr highs, string option, double value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_setStringOptionValue(IntPtr highs, string option, string value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getBoolOptionValue(IntPtr highs, string option, out int value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getIntOptionValue(IntPtr highs, string option, out int value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getDoubleOptionValue(IntPtr highs, string option, out double value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getStringOptionValue(IntPtr highs, string option, [Out] StringBuilder value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getSolution(IntPtr highs, double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getNumCol(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getNumRow(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getNumNz(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getHessianNumNz(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getBasis(IntPtr highs, int[] colstatus, int[] rowstatus);

    [LibraryImport(HighsLibName)]
    internal static extern double Highs_getObjectiveValue(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getIterationCount(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getModelStatus(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_addRow(IntPtr highs, double lower, double upper, int num_new_nz, int[] indices, double[] values);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_addRows(
        IntPtr highs,
        int num_new_row,
        double[] lower,
        double[] upper,
        int num_new_nz,
        int[] starts,
        int[] indices,
        double[] values);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_addCol(
        IntPtr highs,
        double cost,
        double lower,
        double upper,
        int num_new_nz,
        int[] indices,
        double[] values);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_addCols(
        IntPtr highs,
        int num_new_col,
        double[] costs,
        double[] lower,
        double[] upper,
        int num_new_nz,
        int[] starts,
        int[] indices,
        double[] values);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeObjectiveSense(IntPtr highs, int sense);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeColCost(IntPtr highs, int column, double cost);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeColsCostBySet(IntPtr highs, int num_set_entries, int[] set, double[] cost);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeColsCostByMask(IntPtr highs, int[] mask, double[] cost);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeColBounds(IntPtr highs, int column, double lower, double upper);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeColsBoundsByRange(IntPtr highs, int from_col, int to_col, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeColsBoundsBySet(IntPtr highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeColsBoundsByMask(IntPtr highs, int[] mask, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeRowBounds(IntPtr highs, int row, double lower, double upper);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeRowsBoundsByRange(IntPtr highs, int from_row, int to_row, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeRowsBoundsBySet(IntPtr highs, int num_set_entries, int[] set, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeRowsBoundsByMask(IntPtr highs, int[] mask, double[] lower, double[] upper);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeColsIntegralityByRange(IntPtr highs, int from_col, int to_col, int[] integrality);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_changeCoeff(IntPtr highs, int row, int column, double value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_deleteColsByRange(IntPtr highs, int from_col, int to_col);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_deleteColsBySet(IntPtr highs, int num_set_entries, int[] set);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_deleteColsByMask(IntPtr highs, int[] mask);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_deleteRowsByRange(IntPtr highs, int from_row, int to_row);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_deleteRowsBySet(IntPtr highs, int num_set_entries, int[] set);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_deleteRowsByMask(IntPtr highs, int[] mask);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getDoubleInfoValue(IntPtr highs, string info, out double value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getIntInfoValue(IntPtr highs, string info, out int value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getInt64InfoValue(IntPtr highs, string info, out long value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_setSolution(IntPtr highs, double[] col_value, double[] row_value, double[] col_dual, double[] row_dual);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_setSparseSolution(IntPtr highs, int num_entries, int[] index, double[] value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getColsByRange(
        IntPtr highs,
        int from_col,
        int to_col,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getColsBySet(
        IntPtr highs,
        int num_set_entries,
        int[] set,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getColsByMask(
        IntPtr highs,
        int[] mask,
        ref int num_col,
        double[] costs,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getRowsByRange(
        IntPtr highs,
        int from_row,
        int to_row,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getRowsBySet(
        IntPtr highs,
        int num_set_entries,
        int[] set,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getRowsByMask(
        IntPtr highs,
        int[] mask,
        ref int num_row,
        double[] lower,
        double[] upper,
        ref int num_nz,
        int[] matrix_start,
        int[] matrix_index,
        double[] matrix_value);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getBasicVariables(IntPtr highs, int[] basic_variables);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getBasisInverseRow(IntPtr highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getBasisInverseCol(IntPtr highs, int column, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getBasisSolve(
        IntPtr highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_num_nz,
        int[] solution_indices);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getBasisTransposeSolve(
        IntPtr highs,
        double[] rhs,
        double[] solution_vector,
        ref int solution_nz,
        int[] solution_indices);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getReducedRow(IntPtr highs, int row, double[] row_vector, ref int row_num_nz, int[] row_indices);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_getReducedColumn(IntPtr highs, int column, double[] col_vector, ref int col_num_nz, int[] col_indices);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_clearModel(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_clearSolver(IntPtr highs);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_passColName(IntPtr highs, int column, string name);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_passRowName(IntPtr highs, int row, string name);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_writeOptions(IntPtr highs, string filename);

    [LibraryImport(HighsLibName)]
    internal static extern int Highs_writeOptionsDeviations(IntPtr highs, string filename);
    #endregion
}