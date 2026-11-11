// Solves an LP with the HiPO solver, which needs the Highs.Native.Extras
// NuGet package next to Highs.Native. Exits with 1 if HiPO is unavailable.

using System;

using Highs;

class Program {
   static int Main(string[] args) {
      // minimize x_2 subject to x_1 + x_2 + x_3 >= 1; x >= 0
      double[] cc = {0, 1, 0};
      double[] cl = {0, 0, 0};
      double[] cu = {1.0e30, 1.0e30, 1.0e30};
      double[] rl = {1};
      double[] ru = {1.0e30};
      int[] astart = {0, 3};
      int[] aindex = {0, 1, 2};
      double[] avalue = {1, 1, 1};

      HighsModel model = new HighsModel(cc, cl, cu, rl, ru, astart, aindex, avalue, null, 0,
                                        HighsMatrixFormat.kRowwise, HighsObjectiveSense.kMinimize);

      HighsLpSolver solver = new HighsLpSolver();
      if (solver.setStringOptionValue("solver", "hipo") != HighsStatus.kOk) {
         Console.WriteLine("HiPO is unavailable: add the Highs.Native.Extras package");
         return 1;
      }
      // so that HiPO solves the LP rather than presolve
      solver.setStringOptionValue("presolve", "off");
      solver.passLp(model);
      HighsStatus status = solver.run();
      HighsModelStatus modelStatus = solver.GetModelStatus();

      Console.WriteLine("Status: " + status);
      Console.WriteLine("Modelstatus: " + modelStatus);
      return status == HighsStatus.kOk && modelStatus == HighsModelStatus.kOptimal ? 0 : 1;
   }
}
