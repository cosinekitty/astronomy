@echo off
setlocal EnableDelayedExpansion

for %%x in (msbuild.exe) do (set msbuild=%%~$PATH:x)
if not defined msbuild set msbuild=c:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\MSBuild\Current\Bin\MSBuild.exe
if not exist "!msbuild!" (
    echo.ERROR: cannot find msbuild.exe
    exit /b 1
)

echo.Building C code

set OUTDIR=%cd%\bin\
"!msbuild!" windows\generate\generate.sln /p:Configuration=Release /v:quiet /nologo /p:clp=Summary || (
    echo.BUILD FAILED.
    exit /b 1
)

echo.Build succeeded.
exit /b 0