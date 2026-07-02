# use cmake to build the HiGHS windows .dll
try {
	if (-not (Test-Path build_windows)) {
		mkdir build_windows
	}
	cd build_windows

	# Define a list to store arguments
	$arguments = @()

	$arguments += "-A x64"
	$arguments += "-DCMAKE_SYSTEM_NAME=Windows"
	$arguments += "-DCMAKE_BUILD_TYPE=Release"
	$arguments += "-DBUILD_SHARED_LIBS=ON"

	# Compiler flags as a single string
	# $flags = "/arch:AVX512 /Ox /Ot /Oi /O2" - Intel stopped supporting AVX-512 on some CPUs, me might run into compatibility issues
	$flags = "/Ox /Ot /Oi /O2"
	$arguments += "-DCMAKE_CXX_FLAGS='$flags'"
	$arguments += "-DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded"

	# Pass the argument list to cmake
	& cmake $arguments ../..
	cmake --build . --config Release

	# # Generate a MinGW-compatible import library (libhighs.dll.a) so that
	# # MinGW consumers can link against this MSVC-built DLL without needing
	# # the MSVC .lib format. Requires dlltool from a MinGW/MSYS2 installation.
	# $dll = (Get-ChildItem -Recurse -Filter "highs.dll" | Select-Object -First 1).FullName
	# $lib = (Get-ChildItem -Recurse -Filter "highs.lib" | Select-Object -First 1).FullName
	# if ($dll -and (Get-Command dlltool -ErrorAction SilentlyContinue)) {
	# 	$libdir = Split-Path $lib
	# 	Write-Host "Generating MinGW import library: $libdir\libhighs.dll.a"
	# 	& dlltool -D (Split-Path $dll -Leaf) -l "$libdir\libhighs.dll.a" 2>&1
	# } else {
	# 	Write-Host "dlltool not found; skipping MinGW import library generation. MinGW consumers can generate it with: dlltool -D highs.dll -l libhighs.dll.a"
	# }
	}
catch {
		Write-Error "Failed to build for Windows: $($_.Exception.Message)"
}
