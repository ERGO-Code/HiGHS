# Highs.Native.Extras

Optional extension for [Highs.Native](https://www.nuget.org/packages/Highs.Native/) that enables the HiPO interior point solver. It contains external dependencies whose licenses differ from the MIT license of HiGHS (Apache-2.0 and BSD-3-Clause, see `THIRD_PARTY_NOTICES.md`).

Add it next to `Highs.Native`, with the same version:

```shell
dotnet add package Highs.Native
dotnet add package Highs.Native.Extras
```

Highs.Native loads the extension at run time, so no code changes are needed. Select HiPO with

```csharp
solver.setStringOptionValue("solver", "hipo");
```

Without this package, setting `solver` to `hipo` returns `HighsStatus.kError`, and the other solvers work as before.
