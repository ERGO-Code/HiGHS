REM SPDX-FileCopyrightText: 2020 Intel Corporation
REM
REM SPDX-License-Identifier: MIT

set LANGUAGE=%1
set VS_VER=%2
set SAMPLES_TAG=%3

IF "%VS_VER%"=="2017_build_tools" (
@call "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
)

IF "%VS_VER%"=="2019_build_tools" (
@call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
)

IF "%VS_VER%"=="2022_build_tools" (
@call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
)

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

git clone --depth 1 --branch %SAMPLES_TAG% https://github.com/oneapi-src/oneAPI-samples.git

if "%LANGUAGE%" == "c++" goto cpp
if "%LANGUAGE%" == "fortran" goto fortran
if "%LANGUAGE%" == "dpc++" goto dpcpp
goto exit

:cpp
cd oneAPI-samples\DirectProgramming\C++\CompilerInfrastructure\Intrinsics
icl -O2 src\intrin_dot_sample.cpp
icl -O2 src\intrin_double_sample.cpp
icl -O2 src\intrin_ftz_sample.cpp
intrin_dot_sample.exe && intrin_double_sample.exe && intrin_ftz_sample.exe
set RESULT=%ERRORLEVEL%
del intrin_dot_sample.exe intrin_double_sample.exe intrin_ftz_sample.exe
icx -O2 -msse3 src\intrin_dot_sample.cpp
icx -O2 -msse3 src\intrin_double_sample.cpp
icx -O2 -msse3 src\intrin_ftz_sample.cpp
intrin_dot_sample.exe && intrin_double_sample.exe && intrin_ftz_sample.exe
set /a RESULT=%RESULT%+%ERRORLEVEL%
goto exit

:fortran
cd oneAPI-samples\DirectProgramming\Fortran\CombinationalLogic\openmp-primes
ifort -O2 -fpp -qopenmp src\openmp_sample.f90
openmp_sample.exe
set RESULT=%ERRORLEVEL%
ifx -O2 -fpp -qopenmp src\openmp_sample.f90
openmp_sample.exe
set /a RESULT=%RESULT%+%ERRORLEVEL%
goto exit

:dpcpp
for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\tbb\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\tbb\%LATEST_VERSION%\env\vars.bat"
cd oneAPI-samples\DirectProgramming\DPC++\DenseLinearAlgebra\vector-add
nmake -f Makefile.win
nmake -f Makefile.win run
set RESULT=%ERRORLEVEL%
@REM goto exit


exit /b %RESULT%
