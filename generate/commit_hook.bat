@echo off
del output\vsop*.txt output\*.eph output\jupiter_moons.txt
call run.bat
if errorlevel 1 (exit /b 1)
call verify_clean.bat
if errorlevel 1 (exit /b 1)
exit /b 0
