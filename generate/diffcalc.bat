@echo off
setlocal EnableDelayedExpansion
set CTESTEXE=bin\ctest.exe

if not exist !CTESTEXE! (
    echo.FATAL[diffcalc]: executable does not exist: !CTESTEXE!
    exit /b 1
)

echo.Diffing calculations.

set /a FAILCOUNT = 0

!CTESTEXE! diff 5.3e-15 temp\c_check.txt dotnet\csharp_test\csharp_check.txt
if errorlevel 1 (set /a FAILCOUNT += 1)

!CTESTEXE! diff 6.3e-15 temp\c_check.txt temp\k_check.txt
if errorlevel 1 (set /a FAILCOUNT += 1)

!CTESTEXE! diff 5.7e-15 temp\c_check.txt temp\js_check.txt
if errorlevel 1 (set /a FAILCOUNT += 1)

!CTESTEXE! diff 7.11e-16 temp\c_check.txt temp\py_check.txt
if errorlevel 1 (set /a FAILCOUNT += 1)

!CTESTEXE! diff 5.7e-15 temp\js_check.txt temp\py_check.txt
if errorlevel 1 (set /a FAILCOUNT += 1)

if !FAILCOUNT! NEQ 0 (
    echo.diffcalc: *** FAILED !FAILCOUNT! TESTS ***
    exit /b 1
)

echo.diffcalc: PASS
echo.
exit /b 0
