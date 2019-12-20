@echo off
setlocal EnableDelayedExpansion
dotnet build --output !CD!\exe
if errorlevel 1 (
    echo.Error building csdown
    exit /b 1
)
