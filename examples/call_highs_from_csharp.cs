class Program {
   static void Main(string[] args) {
      double[] cc = {1, -2};
      double[] cl = {0, 0};
      double[] cu = {10, 10};
      double[] rl = {0, 0};
      double[] ru = {2, 1};
      int[] astart = {0, 2, 4};
      int[] aindex = {0, 1, 0, 1};
      double[] avalue = {1, 2, 1, 3};
      HighsLpSolver.call_highs(cc, cl, cu, rl, ru, astart, aindex, avalue);
   }
}