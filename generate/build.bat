@echo off
setlocal EnableDelayedExpansion
set msbuild=%ProgramFiles(x86)%\MSBuild\14.0\Bin\amd64\msbuild.exe
if not exist "!msbuild!" (
    echo.ERROR: cannot find msbuild.exe
    exit /b 1
)

echo.Building C code

set OUTDIR=%cd%\bin\
"!msbuild!" ..\windows\generate\generate.sln /p:Configuration=Release /v:quiet /nologo /p:clp=Summary
if errorlevel 1 (
    echo.BUILD FAILED.
    exit /b 1
)

echo.Build succeeded.
exit /b 0