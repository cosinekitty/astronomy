@echo off
if defined PROGRAMFILES(x86) (
    echo.Detected 64-bit Windows.
    for %%f in (
        ..\windows\generate\x64\Release\generate.exe
        ..\windows\generate\Release\generate.exe
        ..\windows\generate\x64\Debug\generate.exe
        ..\windows\generate\Debug\generate.exe
    ) do (
        if exist %%f (
            set GENEXE=%%f
            goto found_generate_exe
        )
    )
) else (
    echo.Detected 32-bit Windows.
    for %%f in (
        ..\windows\generate\Release\generate.exe
        ..\windows\generate\Debug\generate.exe
    ) do (
        if exist %%f (
            set GENEXE=%%f
            goto found_generate_exe
        )
    )
)
echo.ERROR: Cannot find generate.exe. 
echo.Use Microsoft Visual Studio to build ..\windows\generate.sln
exit /b 1
:found_generate_exe
echo.Found executable: %GENEXE%
exit /b 0