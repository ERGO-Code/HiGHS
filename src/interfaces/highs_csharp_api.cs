using System;
using System.Linq;
using System.Runtime.InteropServices;

// mcs -out:highscslib.dll -t:library highs_csharp_api.cs -unsafe

// TODO: add HighsStatus enum

// TODO: add HighsBasisStatus enum

// TODO: add objective sense enum

public class HighsModel
{
   public double[] colcost;
   public double[] collower;
   public double[] colupper;
   public double[] rowlower;
   public double[] rowupper;
   public int[] astart;
   public int[] aindex;
   public double[] avalue;

   public HighsModel()
   {

   }

   public HighsModel(double[] colcost, double[] collower, double[] colupper, double[] rowlower, double[] rowupper,
   int[] astart, int[] aindex, double[] avalue)
   {
      this.colcost = colcost;
      this.collower = collower;
      this.colupper = colupper;
      this.rowlower = rowlower;
      this.rowupper = rowupper;
      this.astart = astart;
      this.aindex = aindex;
      this.avalue = avalue;
   }
}

public class HighsSolution
{
   public double[] colvalue;
   public double[] coldual;
   public double[] rowvalue;
   public double[] rowdual;

   public HighsSolution(int numcol, int numrow)
   {
      this.colvalue = new double[numcol];
      this.coldual = new double[numcol];
      this.rowvalue = new double[numrow];
      this.rowdual = new double[numrow];
   }

   public HighsSolution(double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual)
   {
      this.colvalue = colvalue;
      this.coldual = coldual;
      this.rowvalue = rowvalue;
      this.rowdual = rowdual;
   }
}

public class HighsBasis
{
   public int[] colbasisstatus;
   public int[] rowbasisstatus;

   public HighsBasis(int numcol, int numrow)
   {
      this.colbasisstatus = new int[numcol];
      this.rowbasisstatus = new int[numrow];
   }

   public HighsBasis(int[] colbasisstatus, int[] rowbasisstatus)
   {
      this.colbasisstatus = colbasisstatus;
      this.rowbasisstatus = rowbasisstatus;
   }
}

public unsafe class HighsLpSolver
{
   private void* highs;

   [DllImport("libhighs.so")]
   private static extern int Highs_call(Int32 numcol, Int32 numrow, Int32 numnz, double[] colcost,
   double[] collower, double[] colupper, double[] rowlower, double[] rowupper, int[] astart, int[] aindex, double[] avalue,
   double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual, int[] colbasisstatus, int[] rowbasisstatus);

   [DllImport("libhighs.so")]
   private static extern void* Highs_create();

   [DllImport("libhighs.so")]
   private static extern void Highs_destroy(void* highs);

   [DllImport("libhighs.so")]
   private static extern int Highs_run(void* highs);

   [DllImport("libhighs.so")]
   private static extern int Highs_readFromFile(void* highs, string filename);

   [DllImport("libhighs.so")]
   private static extern int Highs_writeToFile(void* highs, string filename);

   [DllImport("libhighs.so")]
   private static extern int Highs_loadModel(void* highs, int numcol, int numrow, int numnz, double[] colcost,
   double[] collower, double[] colupper, double[] rowlower, double[] rowupper, int[] astart, int[] aindex, double[] avalue);

   [DllImport("libhighs.so")]
   private static extern int Highs_setOptionValue(void* highs, string option, string value);

   [DllImport("libhighs.so")]
   private static extern void Highs_getSolution(void* highs, double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual);

   [DllImport("libhighs.so")]
   private static extern int Highs_getNumCols(void* highs);

   [DllImport("libhighs.so")]
   private static extern int Highs_getNumRows(void* highs);

   [DllImport("libhighs.so")]
   private static extern int Highs_getNumNz(void* highs);

   [DllImport("libhighs.so")]
   private static extern void Highs_getBasis(void* highs, int[] colstatus, int[] rowstatus);

   [DllImport("libhighs.so")]
   private static extern double Highs_getObjectiveValue(void* highs);

   [DllImport("libhighs.so")]
   private static extern int Highs_getIterationCount(void* highs);

   [DllImport("libhighs.so")]
   private static extern int Highs_addRow(void* highs, double lower, double upper, int num_new_nz, int[] indices, double[] values);

   [DllImport("libhighs.so")]
   private static extern int Highs_addRows(void* highs, int num_new_row, double[] lower, double[] upper, 
   int num_new_nz, int[] starts, int[] indices, double[] values);

   [DllImport("libhighs.so")]
   private static extern int Highs_addCol(void* highs, double cost, double lower, double upper, 
   int num_new_nz, int[] indices, double[] values);

   [DllImport("libhighs.so")]
   private static extern int Highs_addCols(void* highs, int num_new_col, double[] costs, double[] lower, double[] upper, 
   int num_new_nz, int[] starts, int[] indices, double[] values);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeObjectiveSense(void* highs, int sense);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeColCost(void* highs, int col,  double cost);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeColsCostBySet(void* highs, int num_set_entries, int[] set, double[] cost);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeColsCostByMask(void* highs, int[] mask, double[] cost);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeColBounds(void* highs, int col, double lower, double upper);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeColsBoundsByRange(void* highs, int from_col, int to_col, double[] lower, double[] upper);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeColsBoundsBySet(void* highs, int num_set_entries, int[] set, double[] lower, double[] upper);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeColsBoundsByMask(void* highs, int[] mask, double[] lower, double[] upper);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeRowBounds(void* highs, int row, double lower, double upper);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeRowsBoundsBySet(void* highs, int num_set_entries, int[] set, double[] lower, double[] upper);

   [DllImport("libhighs.so")]
   private static extern int Highs_changeRowsBoundsByMask(void* highs, int[] mask, double[] lower, double[] upper);

   [DllImport("libhighs.so")]
   private static extern int Highs_deleteColsByRange(void* highs, int from_col, int to_col);

   [DllImport("libhighs.so")]
   private static extern int Highs_deleteColsBySet(void* highs, int num_set_entries, int[] set);

   [DllImport("libhighs.so")]
   private static extern int Highs_deleteColsByMask(void* highs, int[] mask);

   [DllImport("libhighs.so")]
   private static extern int Highs_deleteRowsByRange(void* highs, int from_row, int to_row);

   [DllImport("libhighs.so")]
   private static extern int Highs_deleteRowsBySet(void* highs, int num_set_entries, int[] set);

   [DllImport("libhighs.so")]
   private static extern int Highs_deleteRowsByMask(void* highs, int[] mask);

   [DllImport("libhighs.so")]
   private static extern int Highs_getColsByRange(void *highs, int from_col, int to_col, ref int num_col, double[] costs, 
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport("libhighs.so")]
   private static extern int Highs_getColsBySet(void *highs, int num_set_entries, int[] set, ref int num_col, double[] costs, 
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport("libhighs.so")]
   private static extern int Highs_getColsByMask(void *highs, int[] mask, ref int num_col, double[] costs, 
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport("libhighs.so")]
   private static extern int Highs_getRowsByRange(void *highs, int from_row, int to_row, ref int num_row, 
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   [DllImport("libhighs.so")]
   private static extern int Highs_getRowsBySet(void *highs, int num_set_entries, int[] set, ref int num_row, 
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);
   
   [DllImport("libhighs.so")]
   private static extern int Highs_getRowsByMask(void *highs, int[] mask, ref int num_row, 
   double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   public static int call(HighsModel model, ref HighsSolution sol, ref HighsBasis bas)
   {
      int nc = model.colcost.Length;
      int nr = model.rowlower.Length;
      int nnz = model.avalue.Length;

      return HighsLpSolver.callhighs(nc, nr, nnz, model.colcost, model.collower, model.colupper,
      model.rowlower, model.rowupper, model.astart, model.aindex, model.avalue,
      sol.colvalue, sol.coldual, sol.rowvalue, sol.rowdual, bas.colbasisstatus, bas.rowbasisstatus);
   }

   public HighsLpSolver()
   {
      this.highs = HighsLpSolver.Highs_create();
   }

   ~HighsLpSolver()
   {
      HighsLpSolver.Highs_destroy(this.highs);
   }

   public int run()
   {
      return HighsLpSolver.Highs_run(this.highs);
   }

   public int readFromFile(string filename)
   {
      return HighsLpSolver.Highs_readFromFile(this.highs, filename);
   }

   public int writeToFile(string filename)
   {
      return HighsLpSolver.Highs_writeToFile(this.highs, filename);
   }

   public int loadModel(HighsModel model)
   {
      return HighsLpSolver.Highs_loadModel(this.highs, model.colcost.Length, model.rowlower.Length, model.avalue.Length,
      model.colcost, model.collower, model.colupper, model.rowlower, model.rowupper, model.astart, model.aindex, model.avalue);
   }

   public int setOptionValue(string option, string value)
   {
      return HighsLpSolver.Highs_setOptionValue(this.highs, option, value);
   }

   public int getNumCols()
   {
      return HighsLpSolver.Highs_getNumCols(this.highs);
   }

   public int getNumRows()
   {
      return HighsLpSolver.Highs_getNumRows(this.highs);
   }

   public int getNumNz()
   {
      return HighsLpSolver.Highs_getNumNz(this.highs);
   }

   public HighsSolution getSolution()
   {
      int nc = this.getNumCols();
      int nr = this.getNumRows();

      HighsSolution sol = new HighsSolution(nc, nr);
      HighsLpSolver.Highs_getSolution(this.highs, sol.colvalue, sol.coldual, sol.rowvalue, sol.rowdual);
      return sol;
   }

   public HighsBasis getBasis()
   {
      int nc = this.getNumCols();
      int nr = this.getNumRows();

      HighsBasis bas = new HighsBasis(nc, nr);
      HighsLpSolver.Highs_getBasis(this.highs, bas.colbasisstatus, bas.rowbasisstatus);
      return bas;
   }

   public double getObjectiveValue()
   {
      return HighsLpSolver.Highs_getObjectiveValue(this.highs);
   }

   public int getIterationCount()
   {
      return HighsLpSolver.Highs_getIterationCount(this.highs);
   }

   public int addRow(double lower, double upper, int[] indices, double[] values) {
      return HighsLpSolver.Highs_addRow(this.highs, lower, upper, indices.Length, indices, values);
   }

   public int addRows(double[] lower, double[] upper, int[] starts, int[] indices, double[] values) {
      return HighsLpSolver.Highs_addRows(this.highs, lower.Length, lower, upper, indices.Length, starts, indices, values);
   }

   public int addCol(double cost, double lower, double upper, int[] indices, double[] values) {
      return HighsLpSolver.Highs_addCol(this.highs, cost, lower, upper, indices.Length, indices, values);
   }

   public int addCols(double[] costs, double[] lower, double[] upper, int[] starts, int[] indices, double[] values) {
      return HighsLpSolver.Highs_addCols(this.highs, costs.Length, costs, lower, upper, indices.Length, starts, indices, values);
   }

   public int changeObjectiveSense(int sense) {
      return HighsLpSolver.Highs_changeObjectiveSense(this.highs, sense);
   }

   public int changeColCost(int col, double cost) {
      return HighsLpSolver.Highs_changeColCost(this.highs, col, cost);
   }

   public int changeColsCostBySet(int[] cols, double[] costs) {
      return HighsLpSolver.Highs_changeColsCostBySet(this.highs, cols.Length, cols, costs);
   }

   public int changeColsCostByMask(bool[] mask, double[] cost) {
      return HighsLpSolver.Highs_changeColsCostByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), cost);
   }

   public int changeColBounds(int col, double lower, double upper) {
      return HighsLpSolver.Highs_changeColBounds(this.highs, col, lower, upper);
   }

   public int changeColsBoundsByRange(int from, int to, double[] lower, double[] upper) {
      return HighsLpSolver.Highs_changeColsBoundsByRange(this.highs, from, to, lower, upper);
   }

   public int changeColsBoundsBySet(int[] cols, double[] lower, double[] upper) {
      return HighsLpSolver.Highs_changeColsBoundsBySet(this.highs, cols.Length, cols, lower, upper);
   }

   public int changeColsBoundsByMask(bool[] mask, double[] lower, double[] upper) {
      return HighsLpSolver.Highs_changeColsBoundsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), lower, upper);
   }

   public int changeRowBounds(int row, double lower, double upper) {
      return HighsLpSolver.Highs_changeRowBounds(this.highs, row, lower, upper);
   }

   public int changeRowsBoundsBySet(int[] rows, double[] lower, double[] upper) {
      return HighsLpSolver.Highs_changeRowsBoundsBySet(this.highs, rows.Length, rows, lower, upper);
   }

   public int changeRowsBoundsByMask(bool[] mask, double[] lower, double[] upper) {
      return HighsLpSolver.Highs_changeRowsBoundsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray(), lower, upper);
   }

   public int deleteColsByRange(int from, int to) {
      return HighsLpSolver.Highs_deleteColsByRange(this.highs, from, to);
   }

   public int deleteColsBySet(int[] cols) {
      return HighsLpSolver.Highs_deleteColsBySet(this.highs, cols.Length, cols);
   }

   public int deleteColsByMask(bool[] mask) {
      return HighsLpSolver.Highs_deleteColsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray());
   }

   public int deleteRowsByRange(int from, int to) {
      return HighsLpSolver.Highs_deleteRowsByRange(this.highs, from, to);
   }

   public int deleteRowsBySet(int[] rows) {
      return HighsLpSolver.Highs_deleteRowsBySet(this.highs, rows.Length, rows);
   }

   public int deleteRowsByMask(bool[] mask) {
      return HighsLpSolver.Highs_deleteRowsByMask(this.highs, mask.Select(x => x ? 1 : 0).ToArray());
   }

   // int Highs_getColsByRange(void *highs, int from_col, int to_col, ref int num_col, double[] costs, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   // [DllImport("libhighs.so")]
   // int Highs_getColsBySet(void *highs, int num_set_entries, int[] set, ref int num_col, double[] costs, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   // [DllImport("libhighs.so")]
   // int Highs_getColsByMask(void *highs, int[] mask, ref int num_col, double[] costs, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   // [DllImport("libhighs.so")]
   // int Highs_getRowsByRange(void *highs, int from_row, int to_row, ref int num_row, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);

   // [DllImport("libhighs.so")]
   // int Highs_getRowsBySet(void *highs, int num_set_entries, int[] set, ref int num_row, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);
   
   // [DllImport("libhighs.so")]
   // int Highs_getRowsByMask(void *highs, int[] mask, ref int num_row, 
   // double[] lower, double[] upper, ref int num_nz, int[] matrix_start, int[] matrix_index, double[] matrix_value);
}
