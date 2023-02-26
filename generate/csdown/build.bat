@echo off
setlocal EnableDelayedExpansion
REM https://github.com/dotnet/sdk/issues/30624#issuecomment-1432118204
dotnet build --property:OutputPath=!CD!\exe || (
    echo.Error building csdown
    exit /b 1
)
