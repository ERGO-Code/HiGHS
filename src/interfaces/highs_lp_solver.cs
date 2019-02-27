using System;
using System.Runtime.InteropServices;
// mcs -out:highscslib -t:library highs_lp_solver.cs
public class HighsLpSolver {
   [DllImport("libhighs.so")]
   private static extern void callhighs(Int32 nc, Int32 nr, Int32 nnz, double[] cc, double[] cl, double[] du, double[] rl, double[] ru, int[] astart, int[] aindex, double[] avalue);

   public static void call_highs(double[] cc, double[] cl, double[] cu, double[] rl, double[] ru, int[] astart, int[] aindex, double[] avalue) {
      int nc = cc.Length;
      int nr = rl.Length;
      int nnz = aindex.Length;
      callhighs(nc, nr, nnz, cc, cl, cu, rl, ru, astart, aindex, avalue);
   }
}