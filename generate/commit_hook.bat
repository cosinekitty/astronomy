@echo off
setlocal EnableDelayedExpansion
REM --------------------------------------------------------------------------------
REM     Astronomy Engine - GitHub Actions steps for Windows.
REM     This batch file is executed on every push to GitHub.
REM --------------------------------------------------------------------------------

cd %~dp0
set ASTRONOMY_ENGINE_PYTHON_PROBLEM=1
py test.py precision || exit /b 1
set ASTRONOMY_ENGINE_PYTHON_PROBLEM=

REM Change to project/repo root directory.
cd %~dp0\..
echo.commit_hook: Repo root = %cd%

set DOXYGENURL=https://raw.githubusercontent.com/cosinekitty/ephemeris/master/doxygen-1.9.5.windows.x64.bin.zip
md bin
cd bin
echo.commit_hook: Downloading: !DOXYGENURL!
if exist doxygen.zip ( del doxygen.zip )
curl -L -o doxygen.zip !DOXYGENURL! || exit /b 1
echo.commit_hook: Installing Doxygen.
7z x doxygen.zip || exit /b 1
del doxygen.zip

REM change to 'generate' directory, which is where this batch file is located.
cd %~dp0
echo.commit_hook: Running Astronomy Engine tests.
del output\vsop*.txt output\*.eph output\jupiter_moons.txt
call run.bat || exit /b 1
call verify_clean.bat || exit /b 1
echo.commit_hook: SUCCESS
exit /b 0
