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
		
	# Pass the argument list to cmake
	& cmake $arguments ../.. 
	cmake --build . --config Release
	}
catch {
		Write-Error "Failed to build for Windows: $($_.Exception.Message)"
}
