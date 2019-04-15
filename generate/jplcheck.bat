@echo off
setlocal EnableDelayedExpansion

if exist jpl_summary.txt (del jpl_summary.txt)
for %%f in (horizons\*.txt) do (
    node jpl_horizons_check.js %%f > summary.txt
    if errorlevel 1 (exit /b 1)
    type summary.txt
    type summary.txt >> jpl_summary.txt
    del summary.txt
)

node jpl_horizons_check.js tally jpl_summary.txt
if errorlevel 1 (exit /b 1)
del jpl_summary.txt
echo.jplcheck: Finished
exit /b 0
