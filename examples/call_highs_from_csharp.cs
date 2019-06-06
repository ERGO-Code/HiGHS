// mcs -out:cstest call_highs_from_csharp.cs -r:highscslib.dll
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
      int status = HighsLpSolver.Highs_callhighs(cc, cl, cu, rl, ru, astart, aindex, avalue);
      //Console.Writeline("Status: " + status);
   }
}