@echo off
setlocal EnableDelayedExpansion

call build.bat
if errorlevel 1 (exit /b 1)

if defined PROGRAMFILES(x86) (
    echo.Detected 64-bit Windows.
    for %%f in (
        ..\windows\generate\x64\Release\generate.exe
        ..\windows\generate\Release\generate.exe
        ..\windows\generate\x64\Debug\generate.exe
        ..\windows\generate\Debug\generate.exe
    ) do (
        if exist %%f (
            set GENEXE=%%f
            goto found_generate_exe
        )
    )
) else (
    echo.Detected 32-bit Windows.
    for %%f in (
        ..\windows\generate\Release\generate.exe
        ..\windows\generate\Debug\generate.exe
    ) do (
        if exist %%f (
            set GENEXE=%%f
            goto found_generate_exe
        )
    )
)
echo.ERROR: Cannot find generate.exe. 
echo.Use Microsoft Visual Studio to build ..\windows\generate.sln
exit /b 1
:found_generate_exe
echo.Found executable: !GENEXE!

set OUTDIR=..\source
set EPHFILE=lnxp1600p2200.405
REM set EPHURL=ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de405/!EPHFILE!
set EPHURL=https://github.com/cosinekitty/ephemeris/raw/master/!EPHFILE!

if not exist !EPHFILE! (
    echo.
    echo.Ephemeris file not found: !EPHFILE!
    echo.Trying to download for you from:
    echo.!EPHURL!
    echo.
    
    for %%x in (wget.exe) do (set wgetexe=%%~$PATH:x)
    for %%x in (curl.exe) do (set curlexe=%%~$PATH:x)
    for %%x in (md5sum.exe) do (set md5exe=%%~$PATH:x)
    
    if defined wgetexe (
        echo.Trying download using !wgetexe! ...
        !wgetexe! !EPHURL!
        if not errorlevel 1 goto verify_eph
    ) 
    
    if defined curlexe (
        echo.Trying download using !curlexe! ...
        !curlexe! -L -o !EPHFILE! !EPHURL!
        if not errorlevel 1 goto verify_eph
    )

    if exist !EPHFILE! (del !EPHFILE!)

    echo.
    echo.Could not download the ephemeris file.
    echo.Use your browser to download the above file from
    echo.the NASA ftp site into this directory.
    echo.Then run this batch file again to continue.
    exit /b 1

:verify_eph    
    if defined md5exe (
        echo.Using !md5exe! to test integrity of downloaded !EPHFILE!
        !md5exe! -c ephemeris.md5
        if errorlevel 1 (
            echo.Corrupt ephemeris file !EPHFILE! detected.
            if exist !EPHFILE! (del !EPHFILE!)
            exit /b 1
        )
    )
)

echo.
echo.Running the target code generator.
for %%d in (output temp) do (
    if not exist %%d (
        md %%d
        if not exist %%d (
            echo.ERROR: Cannot create directory %%d
            exit /b 1
        )
    )
)
if exist output\vsop_*.txt (del /q output\vsop_*.txt)
if exist temp\* (del /q temp\*)
!GENEXE! planets
if errorlevel 1 (exit /b 1)

!GENEXE! source
if errorlevel 1 (exit /b 1)

echo.
echo Running rise/set test.
node rise_set_test.js
if errorlevel 1 (exit /b 1)

echo.
echo.Validating JavaScript code.
node astro_check.js > temp/check.txt
if errorlevel 1 (exit /b 1)

!GENEXE! check temp/check.txt
if errorlevel 1 (exit /b 1)

echo.
echo.Verifying against JPL Horizons data.
call jplcheck.bat
if errorlevel 1 (exit /b 1)

echo.
echo.Running test of moon phase search.
node moon_phase_test.js
if errorlevel 1 (exit /b 1)

echo.
echo.SUCCESS.
exit /b 0
