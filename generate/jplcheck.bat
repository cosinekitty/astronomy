@echo off
setlocal EnableDelayedExpansion

for %%f in (horizons\*.txt) do (
    node jpl_horizons_check.js %%f
    if errorlevel 1 (exit /b 1)
)

echo.jplcheck: Finished
exit /b 0
