@echo off
setlocal EnableDelayedExpansion

call verify_clean.bat
if errorlevel 1 (exit /b 1)

call run.bat
if errorlevel 1 (exit /b 1)

call verify_clean.bat
if errorlevel 1 (exit /b 1)

git push
if errorlevel 1 (exit /b 1)

exit /b 0
