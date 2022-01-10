@echo off
setlocal EnableDelayedExpansion

for %%x in (msbuild.exe) do (set msbuild=%%~$PATH:x)
if defined msbuild (
    echo Found pre-defined msbuild = !msbuild!
) else (
    set msbuild=%ProgramFiles(x86)%\MSBuild\14.0\Bin\amd64\msbuild.exe
    echo Defined msbuild = !msbuild!
)

if not exist "!msbuild!" (
    echo.ERROR: cannot find msbuild.exe
    exit /b 1
)

echo.Building C code

set OUTDIR=%cd%\bin\
"!msbuild!" windows\generate\generate.sln /p:Configuration=Release /v:quiet /nologo /p:clp=Summary
if errorlevel 1 (
    echo.BUILD FAILED.
    exit /b 1
)

echo.Build succeeded.
exit /b 0