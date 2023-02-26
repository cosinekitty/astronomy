@echo off
setlocal EnableDelayedExpansion

for /L %%y in (1800,10,2100) do (
    if not exist %%y.json (
        echo.Downloading moon phase data for year %%y
        curl -o %%y.json "https://api.usno.navy.mil/moon/phase?year=%%y" || (
            echo.Something bad happened.
            exit /b 1
        )
    )
)

node parse_moon_phases.js || exit /b 1
exit /b 0
