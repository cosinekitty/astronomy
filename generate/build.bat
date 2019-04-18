@echo off
setlocal EnableDelayedExpansion
set msbuild=%ProgramFiles(x86)%\MSBuild\14.0\Bin\amd64\msbuild.exe
if not exist "!msbuild!" (
    echo.ERROR: cannot find msbuild.exe
    exit /b 1
)

"!msbuild!" ..\windows\generate\generate.sln /t:Rebuild /p:Configuration=Release /v:quiet /nologo /p:clp=Summary
