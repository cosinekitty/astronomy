@echo off
REM *** Helps me quickly fix syntax errors while writing new Kotlin code in Windows. ***
pushd ..\..\generate
bin\generate.exe source
if errorlevel 1 (exit /b 1)
popd
call gradlew.bat assemble build test dokkaHtml
if errorlevel 1 (exit /b 1)
exit /b 0