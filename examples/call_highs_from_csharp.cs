// mcs -out:cstest call_highs_from_csharp.cs -r:highscslib.dll

using System;

using Highs;

class Program
{
   static void Main(string[] args)
   {
      string model_name = "egout";

      HighsLpSolver solver = new HighsLpSolver();

      HighsStatus status = solver.readModel(
         "C:\\Users\\galab\\code\\HiGHS\\check\\instances\\" +
         model_name + ".mps");
      Console.WriteLine("Read status: " + status);

      status = solver.run();
      HighsModelStatus modelStatus = solver.GetModelStatus();
      double objective = solver.getObjectiveValue();

      Console.WriteLine("Status: " + status);
      Console.WriteLine("Modelstatus: " + modelStatus);
      Console.WriteLine("Objective: " + modelStatus);

   }
}
