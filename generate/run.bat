@echo off
setlocal EnableDelayedExpansion

call :Download https://github.com/cosinekitty/ephemeris/raw/master/lnxp1600p2200.405 lnxp1600p2200.405 ephemeris.sha256
if errorlevel 1 (exit /b 1)

call :Download https://github.com/cosinekitty/ephemeris/raw/master/top2013/TOP2013.dat TOP2013.dat top2013.sha256
if errorlevel 1 (exit /b 1)

call :Download https://raw.githubusercontent.com/astronexus/HYG-Database/master/hygdata_v3.csv hygdata_v3.csv hygdata_v3.sha256
if errorlevel 1 (exit /b 1)

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

if exist constellation\test_input.txt (del constellation\test_input.txt)
py make_constellation_data.py
if errorlevel 1 (
    echo.Error creating constellation test data.
    exit /b 1
)

if not exist constellation\test_input.txt (
    echo.ERROR - file was not created: constellation\test_input.txt
    exit /b 1
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

echo.
echo.Validating TOP2013 code.
!GENEXE! validate_top2013
if errorlevel 1 (exit /b 1)

if !FASTMODE! == true (
    REM *** Override fast mode if any of the required files are missing.
    for %%f in (
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
    echo.
    echo.Generating planet models.
    !GENEXE! planets
    if errorlevel 1 (
        echo.FATAL: !GENEXE! planets
        exit /b 1
    )
)

if not exist output\jupiter_moons.txt (
    echo.
    echo.Optimizing Jupiter Moon models.
    !GENEXE! jmopt
    if errorlevel 1 (
        echo.FATAL: !GENEXE! jmopt
        exit /b 1
    )
)

echo.
echo.Generating apsis test data.
if exist apsides\apsis_*.txt (del apsides\apsis_*.txt)
!GENEXE! apsis
if errorlevel 1 (
    echo.FATAL: !GENEXE! apsis
    exit /b 1
)

echo.
echo.Generating eclipse data.
cd eclipse
for %%f in (
    lunar_eclipse.txt
    solar_eclipse.txt
    mercury.txt
    venus.txt
) do (
    if exist %%f (del %%f)
)
py norm.py
if errorlevel 1 (
    echo.Error normalizing eclipse test data.
    exit /b 1
)
cd ..

echo.
echo.Generating EQJ/GAL conversion test data.
!GENEXE! galeqj temp\galeqj.txt
if errorlevel 1 (
    echo.Error generating EQJ/GAL conversion test data.
    exit /b 1
)
echo.

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
dotnet run -- all
if errorlevel 1 (exit /b 1)
popd

REM -----------------------------------------------------------------------------------------
echo.
echo.Validating JavaScript code.
node test.js astro_check > temp/js_check.txt
if errorlevel 1 (exit /b 1)

!GENEXE! check temp/js_check.txt
if errorlevel 1 (exit /b 1)

echo.
echo.Verifying against JPL Horizons data.
call jplcheck.bat
if errorlevel 1 (exit /b 1)

echo.
echo.Running JavaScript unit tests.
node test.js all
if errorlevel 1 (exit /b 1)

for %%f in (temp\longitude_*.txt) do (
    !GENEXE! check %%f
    if errorlevel 1 (exit /b 1)
)

REM -----------------------------------------------------------------------------------------

echo.Running C unit tests.
!CTESTEXE! check
if errorlevel 1 (exit /b 1)

!GENEXE! check temp\c_check.txt
if errorlevel 1 (exit /b 1)

!CTESTEXE! all
if errorlevel 1 (exit /b 1)

for %%f in (temp\c_longitude_*.txt) do (
    !GENEXE! check %%f
    if errorlevel 1 (exit /b 1)
)

REM -----------------------------------------------------------------------------------------

echo.Running Python tests.
py test.py all
if errorlevel 1 (exit /b 1)

for %%f in (temp\py_longitude_*.txt) do (
    !GENEXE! check %%f
    if errorlevel 1 (exit /b 1)
)

echo.Generating Python test output.
py test.py astro_check > temp\py_check.txt
if errorlevel 1 (exit /b 1)

echo.Verifying Python test output.
!GENEXE! check temp\py_check.txt
if errorlevel 1 (exit /b 1)

REM -----------------------------------------------------------------------------------------

echo.Running Kotlin tests.
pushd ..\source\kotlin
gradlew assemble build test dokkaHtml
if errorlevel 1 (exit /b 1)
popd

REM -----------------------------------------------------------------------------------------

call diffcalc.bat
if errorlevel 1 (exit /b 1)

if exist ..\website\src\assets\documentation.json (
    REM *** documentation.json never generates the same in Windows as Linux.
    REM *** For now, hack around this by discarding any local changes.
    git checkout -- ../website/src/assets/documentation.json
)

type pass.txt
exit /b 0

REM -----------------------------------------------------------------------------------------
REM     Subroutine for downloading an external file.
REM     Some of the files needed to generate source code are large.
REM     These files are only needed by Astronomy Engine contributors,
REM     not by developers who are using the published version of Astronomy Engine.
REM     A special download process helps keep the repo size reasonable for most users.

:Download
    setlocal
    for %%x in (wget.exe) do (set wgetexe=%%~$PATH:x)
    for %%x in (curl.exe) do (set curlexe=%%~$PATH:x)
    set EPHURL=%1
    set EPHFILE=%2
    set SHAFILE=%3
    if not exist !EPHFILE! (
        echo.
        echo.Local file not found: !EPHFILE!
        echo.Trying to download for you from:
        echo.!EPHURL!
        echo.

        if defined wgetexe (
            echo.Trying download using !wgetexe! ...
            "!wgetexe!" !EPHURL!
            if not errorlevel 1 goto verify_eph
        )

        if defined curlexe (
            echo.Trying download using !curlexe! ...
            "!curlexe!" -L -o !EPHFILE! !EPHURL!
            if not errorlevel 1 goto verify_eph
        )

        if exist !EPHFILE! (del !EPHFILE!)

        echo.
        echo.Could not download the file.
        echo.Use your browser to download the above file from
        echo.the above URL into this directory.
        echo.Then run this batch file again to continue.
        exit /b 1
    )

    :verify_eph
    echo.Using checksum.bat to test integrity of downloaded !EPHFILE!
    call checksum.bat sha256 !SHAFILE!
    if errorlevel 2 (
        echo.Error verifying checksum for !EPHFILE!.
        exit /b 1
    ) else if errorlevel 1 (
        echo.Corrupt ephemeris file !EPHFILE! detected.
        if exist !EPHFILE! (del !EPHFILE!)
        exit /b 1
    )
    exit /b 0

REM -----------------------------------------------------------------------------------------
