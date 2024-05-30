@echo off
setlocal EnableDelayedExpansion

call :Download https://raw.githubusercontent.com/cosinekitty/ephemeris/master/lnxp1600p2200.405 lnxp1600p2200.405 ephemeris.sha256 || exit /b 1
call :Download https://raw.githubusercontent.com/cosinekitty/ephemeris/master/top2013/TOP2013.dat TOP2013.dat top2013.sha256 || exit /b 1
call :Download https://raw.githubusercontent.com/astronexus/HYG-Database/3964fa862d1f08f05919a35306889fa4a0afa7d6/hyg/v3/hyg_v36_1.csv hyg_v36_1.csv hyg_v36_1.sha256 || exit /b 1

cd ..
set copyright=
for /f %%f in ('git grep -l Copyright -- generate hydrogen LICENSE source/c/astronomy.h') do (
    set copyright=!copyright! %%f
)
python generate/update_copyrights.py !copyright!
cd generate

set FASTMODE=true
set GENEXE=bin\generate.exe
set CTESTEXE=bin\ctest.exe
if exist "!GENEXE!" (del "!GENEXE!")
call build.bat || exit /b 1

if not exist "!GENEXE!" (
    echo.FATAL[run.bat]: The executable does not exist: !GENEXE!
    exit /b 1
)

if exist constellation\test_input.txt (del constellation\test_input.txt)
python make_constellation_data.py || exit /b 1
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
!GENEXE! validate_top2013 || exit /b 1

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
    !GENEXE! planets|| (
        echo.FATAL: !GENEXE! planets
        exit /b 1
    )
)

if not exist output\jupiter_moons.txt (
    echo.
    echo.Optimizing Jupiter Moon models.
    !GENEXE! jmopt ||(
        echo.FATAL: !GENEXE! jmopt
        exit /b 1
    )
)

echo.
echo.Generating apsis test data.
if exist apsides\apsis_*.txt (del apsides\apsis_*.txt)
!GENEXE! apsis || exit /b 1

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
python norm.py || exit /b 1
cd ..

echo.
echo.Generating EQJ/GAL conversion test data.
!GENEXE! galeqj temp\galeqj.txt || exit /b 1
echo.

call makedoc.bat || exit /b 1

REM -----------------------------------------------------------------------------------------
REM     makedoc.bat has generated source code as a side effect.
REM     Call build.bat AGAIN to build ctest.c (C unit tests).

echo.Re-building to get C unit test.
if exist "!CTESTEXE!" (del "!CTESTEXE!")
call build.bat || exit /b 1

if not exist "!CTESTEXE!" (
    echo.FATAL[run.bat]: The executable does not exist: !CTESTEXE!
    exit /b 1
)

REM -----------------------------------------------------------------------------------------
echo.
echo.Running C# tests.
pushd dotnet\csharp_test
dotnet run -- all || exit /b 1
popd

REM -----------------------------------------------------------------------------------------
echo.
echo.Validating JavaScript code.
node test.js astro_check > temp\js_check.txt || exit /b 1

!GENEXE! check temp\js_check.txt || exit /b 1

echo.
echo.Verifying against JPL Horizons data.
call jplcheck.bat || exit /b 1

echo.
echo.Running JavaScript unit tests.
node test.js all || exit /b 1

for %%f in (temp\js_longitude_*.txt) do (
    !GENEXE! check %%f || exit /b 1
)

REM -----------------------------------------------------------------------------------------

echo.Running C unit tests.
!CTESTEXE! check || exit /b 1
!GENEXE! check temp\c_check.txt || exit /b 1
!CTESTEXE! all || exit /b 1

for %%f in (temp\c_longitude_*.txt) do (
    !GENEXE! check %%f || exit /b 1
)

REM -----------------------------------------------------------------------------------------

echo.Running Python tests.
python test.py all || exit /b 1

for %%f in (temp\py_longitude_*.txt) do (
    !GENEXE! check %%f || exit /b 1
)

echo.Generating Python test output.
python test.py astro_check > temp\py_check.txt || exit /b 1

echo.Verifying Python test output.
!GENEXE! check temp\py_check.txt || exit /b 1

REM -----------------------------------------------------------------------------------------

echo.Running Kotlin tests.
if exist temp\k_check.txt (
    del temp\k_check.txt || exit /b 1
)

pushd ..\source\kotlin
if exist build (
    rd /s/q build || exit /b 1
    if exist build (
        echo.run.bat: FATAL: source/kotlin/build directory still exists after attempted deletion.
        exit /b 1
    )
)
call gradlew.bat assemble build test dokkaGfm jar || exit /b 1
popd

cd kotlindoc
python format_kotlin_doc.py || exit /b 1
cd ..

!GENEXE! check temp\k_check.txt || exit /b 1

for %%f in (temp\k_longitude_*.txt) do (
    !GENEXE! check %%f || exit /b 1
)

REM -----------------------------------------------------------------------------------------

call diffcalc.bat || exit /b 1

REM -----------------------------------------------------------------------------------------

echo.Validating Python demos.
pushd ..\demo\python
call demotest.bat || exit /b 1
popd

REM -----------------------------------------------------------------------------------------

echo.Validating Java demos.
pushd ..\demo\java
call demotest.bat || exit /b 1
popd

REM -----------------------------------------------------------------------------------------

echo.Validating Kotlin demos.
pushd ..\demo\kotlin
call demotest.bat || exit /b 1
popd

REM -----------------------------------------------------------------------------------------

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
    ) else (
        py checksum.py sha256 !SHAFILE! && exit /b 0
        REM     Sometimes we have to upgrade to a new version an external file.
        REM     When that happens, we just change the checksum file.
        REM     Then we notice here that the checksum doesn't match,
        REM     delete the file, and re-download it.
        echo.
        echo.The local version of !EPHFILE! is obsolete or corrupt.
        del !EPHFILE!
    )

    echo.Attempting download from:
    echo.!EPHURL!
    echo.

    if defined wgetexe (
        echo.Trying download using !wgetexe! ...
        "!wgetexe!" !EPHURL! && goto verify_eph
    )

    if defined curlexe (
        echo.Trying download using !curlexe! ...
        "!curlexe!" -L -o !EPHFILE! !EPHURL! && goto verify_eph
    )

    if exist !EPHFILE! (del !EPHFILE!)

    echo.
    echo.Could not download the file.
    echo.Use your browser to download the above file from
    echo.the above URL into this directory.
    echo.Then run this batch file again to continue.
    exit /b 1

    :verify_eph
    echo.Verifying integrity of file: !EPHFILE!
    py checksum.py sha256 !SHAFILE! || (
        echo.Corrupt file !EPHFILE! detected.
        if exist !EPHFILE! (del !EPHFILE!)
        exit /b 1
    )
    exit /b 0

REM -----------------------------------------------------------------------------------------
