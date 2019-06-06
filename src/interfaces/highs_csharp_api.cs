using System;
using System.Runtime.InteropServices;

// mcs -out:highscslib.dll -t:library highs_csharp_api.cs
public class HighsLpSolver {
   [DllImport("libhighs.so")]
   private static extern int callhighs(Int32 numcol, Int32 numrow, Int32 numnz, double[] colcost,
   double[] collower, double[] colupper, double[] rowlower, double[] rowupper, int[] astart, int[] aindex, double[] avalue,
   double[] colvalue, double[] coldual, double[] rowvalue, double[] rowdual, int[] colbasisstatus, int[] rowbasisstatus);

   public static int Highs_callhighs(double[] cc, double[] cl, double[] cu, double[] rl, double[] ru, int[] ast, int[] ai, double[] av) {
      int nc = cc.Length;
      int nr = rl.Length;
      int nnz = ai.Length;

      double[] cv = new double[nc];
      double[] cd = new double[nc];
      double[] rv = new double[nr];
      double[] rd = new double[nr];

      int[] cbs = new int[nc];
      int[] rbs = new int[nr];

      return callhighs(nc, nr, nnz, cc, cl, cu, rl, ru, ast, ai, av, cv, cd, rv, rd, cbs, rbs);
   }
}