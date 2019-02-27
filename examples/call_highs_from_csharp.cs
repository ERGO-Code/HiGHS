using System;
using System.Runtime.InteropServices;

class Program {
   [DllImport("libhighs.so")]
   private static extern void callhighs(Int32 nc, Int32 nr, Int32 nnz, double[] cc, double[] cl, double[] du, double[] rl, double[] ru, int[] astart, int[] aindex, double[] avalue);

   static void Main(string[] args) {
      double[] cc = {1, -2};
      double[] cl = {0, 0};
      double[] cu = {10, 10};
      double[] rl = {0, 0};
      double[] ru = {2, 1};
      int[] astart = {0, 2, 4};
      int[] aindex = {0, 1, 0, 1};
      double[] avalue = {1, 2, 1, 3};
      int nc = cc.Length;
      int nr = rl.Length;
      int nnz = aindex.Length;
      callhighs(nc, nr, nnz, cc, cl, cu, rl, ru, astart, aindex, avalue);
   }
}