@echo off
setlocal EnableDelayedExpansion

set FASTMODE=true
set GENEXE=bin\generate.exe
set CTESTEXE=bin\ctest.exe
if exist "!GENEXE!" (del "!GENEXE!")
call build.bat
if errorlevel 1 (exit /b 1)

if not exist "!GENEXE!" (
    echo.FATAL[run.bat]: The executable does not exist: !GENEXE!
    exit /b 1
)

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
echo.Cleaning output directories.
for %%d in (output temp) do (
    if not exist %%d (
        md %%d
        if not exist %%d (
            echo.ERROR: Cannot create directory %%d
            exit /b 1
        )
    )
)

if exist temp\* (del /q temp\*)

if !FASTMODE! == true (
    REM *** Override fast mode if any of the required files are missing.
    for %%f in (
        output\08.eph
        output\vsop_0.txt
        output\vsop_1.txt
        output\vsop_3.txt
        output\vsop_4.txt
        output\vsop_5.txt
        output\vsop_6.txt
        output\vsop_7.txt
        output\vsop_11.txt
    ) do (
        if not exist %%f (
            echo.Missing required planet model file: %%f
            set FASTMODE=false
        )
    )
)

if !FASTMODE! == false (
    if exist output\vsop_*.txt (del /q output\vsop_*.txt)
    if exist output\*.eph (del /q output\*.eph)
    echo.
    echo.Generating planet models.
    !GENEXE! planets
    if errorlevel 1 (
        echo.FATAL: !GENEXE! planets
        exit /b 1
    )
)

call makedoc.bat
if errorlevel 1 (exit /b 1)

REM -----------------------------------------------------------------------------------------
REM     makedoc.bat has generated source code as a side effect.
REM     Call build.bat AGAIN to build ctest.c (C unit tests).

echo.Re-building to get C unit test.
if exist "!CTESTEXE!" (del "!CTESTEXE!")
call build.bat
if errorlevel 1 (exit /b 1)

if not exist "!CTESTEXE!" (
    echo.FATAL[run.bat]: The executable does not exist: !CTESTEXE!
    exit /b 1
)
REM -----------------------------------------------------------------------------------------
echo.
echo.Running C# tests.
pushd dotnet\csharp_test
dotnet restore
if errorlevel 1 (exit /b 1)
dotnet run
if errorlevel 1 (exit /b 1)
popd

REM -----------------------------------------------------------------------------------------
echo.
echo.Running longitude tests.
node elong_test.js
if errorlevel 1 (exit /b 1)

for %%f in (temp\longitude_*.txt) do (
    !GENEXE! check %%f
    if errorlevel 1 (exit /b 1)
)

echo.
echo Running seasons test.
node seasons_test.js
if errorlevel 1 (exit /b 1)

echo.
echo Running rise/set test.
node rise_set_test.js
if errorlevel 1 (exit /b 1)

echo.
echo.Validating JavaScript code.
node astro_check.js > temp/js_check.txt
if errorlevel 1 (exit /b 1)

!GENEXE! check temp/js_check.txt
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
echo.Running lunar apsis tests.
node lunar_apsis_test.js
if errorlevel 1 (exit /b 1)

echo.
echo.Running visual magnitude tests.
node mag_test.js
if errorlevel 1 (exit /b 1)

REM -----------------------------------------------------------------------------------------

echo.Running C unit tests.
!CTESTEXE!
if errorlevel 1 (exit /b 1)

!GENEXE! check temp\c_check.txt
if errorlevel 1 (exit /b 1)

!CTESTEXE! diff temp\c_check.txt temp\js_check.txt
if errorlevel 1 (exit /b 1)

!CTESTEXE! seasons seasons\seasons.txt
if errorlevel 1 (exit /b 1)

!CTESTEXE! moonphase moonphase\moonphases.txt
if errorlevel 1 (exit /b 1)

!CTESTEXE! elongation
if errorlevel 1 (exit /b 1)

for %%f in (temp\c_longitude_*.txt) do (
    !GENEXE! check %%f
    if errorlevel 1 (exit /b 1)
    echo.
)

!CTESTEXE! riseset riseset\riseset.txt
if errorlevel 1 (exit /b 1)

!CTESTEXE! magnitude
if errorlevel 1 (exit /b 1)

!CTESTEXE! apsis apsides\moon.txt
if errorlevel 1 (exit /b 1)

REM -----------------------------------------------------------------------------------------

echo.Running Python tests.
test.py time
if errorlevel 1 (exit /b 1)

test.py apsis
if errorlevel 1 (exit /b 1)

test.py magnitude
if errorlevel 1 (exit /b 1)

test.py moon
if errorlevel 1 (exit /b 1)

test.py seasons seasons\seasons.txt
if errorlevel 1 (exit /b 1)

test.py moonphase moonphase\moonphases.txt
if errorlevel 1 (exit /b 1)

test.py riseset riseset\riseset.txt
if errorlevel 1 (exit /b 1)

test.py elongation
if errorlevel 1 (exit /b 1)

for %%f in (temp\py_longitude_*.txt) do (
    !GENEXE! check %%f
    if errorlevel 1 (exit /b 1)
    echo .
)

echo.Generating Python test output.
test.py astro_check > temp\py_check.txt
if errorlevel 1 (exit /b 1)

echo.Verifying Python test output.
!GENEXE! check temp\py_check.txt
if errorlevel 1 (exit /b 1)

!CTESTEXE! diff temp\py_check.txt temp\c_check.txt
if errorlevel 1 (exit /b 1)

type pass.txt
exit /b 0
