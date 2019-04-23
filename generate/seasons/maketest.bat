@echo off
setlocal EnableDelayedExpansion

for /L %%y in (1800,1,2100) do (
    if not exist %%y.json (
        echo.Downloading seasons data for year %%y
        curl -o %%y.json "https://api.usno.navy.mil/seasons?year=%%y"
        if errorlevel 1 (
            echo.Something bad happened.
            exit /b 1
        )
    )
)

node parse_seasons.js
if errorlevel 1 (exit /b 1)
exit /b 0
