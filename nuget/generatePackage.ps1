try {
    Set-Location "$PSScriptRoot\.."
    # Use regex to match and extract version numbers with any prefix
    $content = Get-Content Version.txt -Raw
    $major = ($content | Select-String -Pattern ".*_MAJOR=(\d+)").Matches.Groups[1].Value
    $minor = ($content | Select-String -Pattern ".*_MINOR=(\d+)").Matches.Groups[1].Value
    $patch = ($content | Select-String -Pattern ".*_PATCH=(\d+)").Matches.Groups[1].Value
    # check if a prereleas is set, i.e. the line with PRE_RELEASE doesn't start with # and has the value YES
	$preReleaseMatch = $content | Select-String -Pattern "[^#]PRE_RELEASE=(.+)" -AllMatches
    if ($preReleaseMatch -ne $null -and $preReleaseMatch.Matches.Count -gt 0 -and $preReleaseMatch.Matches[0].Groups[1].Success) {
        $preRelease = $preReleaseMatch.Matches[0].Groups[1].Value.Trim()
    } else {
        $preRelease = ""
    }
    # Combine and output the version
    $version = "$major.$minor.$patch"
    if ($preRelease -eq "YES") {
        # Append a pre-release suffix if needed
        $version = "$version-alpha" # Adjust this line as needed for your pre-release naming convention
    }

    Set-Location "$PSScriptRoot"
    # Now pack
    dotnet pack -c Release /p:Version=$version
}
catch {
    Write-Error $_.Exception.Message
}
