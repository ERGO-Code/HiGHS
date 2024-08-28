// mcs -out:cstest call_highs_from_csharp.cs -r:highscslib.dll

using System;

using Highs;

class Program {
   static void Main(string[] args) {
       // Illustrate the solution of a QP, after first solving just the LP
       //
       // minimize x_2 + (1/2)(2x_1^2 - 2x_1x_3 + 0.2x_2^2 + 2x_3^2)
       //
       // subject to x_1 + x_2 + x_3 >= 1; x>=0
       double[] cc = {0, 1, 0};
      double[] cl = {0, 0, 0};
      double[] cu = {1.0e30, 1.0e30, 1.0e30};
      double[] rl = {1};
      double[] ru = {1.0e30};
      int[] astart = {0, 3};
      int[] aindex = {0, 1, 2};
      double[] avalue = {1, 1, 1};
      HighsObjectiveSense sense = HighsObjectiveSense.kMinimize;
      double offset = 0;
      HighsMatrixFormat a_format = HighsMatrixFormat.kRowwise;

      HighsModel model = new HighsModel(cc, cl, cu, rl, ru, astart, aindex, avalue, null, offset, a_format, sense);

      HighsLpSolver solver = new HighsLpSolver();

      HighsStatus status = solver.passLp(model);
      status = solver.run();
      HighsSolution sol = solver.getSolution();
      HighsBasis bas = solver.getBasis();
      HighsModelStatus modelStatus = solver.GetModelStatus();
      
      Console.WriteLine("Status: " + status);
      Console.WriteLine("Modelstatus: " + modelStatus);
   
      for (int i=0; i<sol.colvalue.Length; i++) {
         Console.WriteLine("Activity for col " + i + " = " + sol.colvalue[i]);
      }
      for (int i=0; i<sol.rowvalue.Length; i++) {
         Console.WriteLine("Activity for row " + i + " = " + sol.rowvalue[i]);
      }
      for (int i=0; i<sol.coldual.Length; i++) {
         Console.WriteLine("Reduced cost x[" + i + "] = " + sol.coldual[i]);
      }
      for (int i=0; i<sol.rowdual.Length; i++) {
         Console.WriteLine("Dual value for row " + i + " = " + sol.rowdual[i]);
      }
      for (int i=0; i<sol.colvalue.Length; i++) {
         Console.WriteLine("x" + i + " = " + sol.colvalue[i] + " is " + bas.colbasisstatus[i]);
      }
       // Add the Hessian
      int dim = 2;
      int[] qstart = {0, 2, 3};
      int[] qindex = {0, 1, 1};
      double[] qvalue = {2, -1, 2};
      HessianFormat q_format = HessianFormat.kTriangular;
      HighsHessian hessian = new HighsHessian(dim, qstart, qindex, qvalue, q_format);
      status = solver.passHessian(hessian);
      status = solver.run();
      sol = solver.getSolution();
      modelStatus = solver.GetModelStatus();
      
      Console.WriteLine("Status: " + status);
      Console.WriteLine("Modelstatus: " + modelStatus);
   
      for (int i=0; i<sol.colvalue.Length; i++) {
         Console.WriteLine("Activity for col " + i + " = " + sol.colvalue[i]);
      }
      for (int i=0; i<sol.rowvalue.Length; i++) {
         Console.WriteLine("Activity for row " + i + " = " + sol.rowvalue[i]);
      }
      for (int i=0; i<sol.coldual.Length; i++) {
         Console.WriteLine("Reduced cost x[" + i + "] = " + sol.coldual[i]);
      }
      for (int i=0; i<sol.rowdual.Length; i++) {
         Console.WriteLine("Dual value for row " + i + " = " + sol.rowdual[i]);
      }
      
   }
}
