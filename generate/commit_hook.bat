@echo off
setlocal EnableDelayedExpansion
REM --------------------------------------------------------------------------------
REM     Astronomy Engine - GitHub Actions steps for Windows.
REM     This batch file is executed on every push to GitHub.
REM --------------------------------------------------------------------------------

REM Change to project/repo root directory.
cd %~dp0\..
echo.commit_hook: Repo root = %cd%

echo.commit_hook: Scraping Doxygen download URL.
set DOXYGENURL=
for /f %%x in ('py generate/doxygen_download_link.py') do (
    set DOXYGENURL=%%x
)
if not defined DOXYGENURL (
    echo.commit_hook: FAIL: could not determine URL to download Doxygen binaries for Windows.
    exit /b 1
)
md bin
cd bin
echo.commit_hook: Downloading: !DOXYGENURL!
curl -o doxygen.zip !DOXYGENURL! || exit /b 1
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
