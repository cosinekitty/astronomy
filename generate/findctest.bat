@echo off
if defined PROGRAMFILES(x86) (
    echo.Detected 64-bit Windows.
    for %%f in (
        ..\windows\generate\x64\Release\ctest.exe
        ..\windows\generate\Release\ctest.exe
        ..\windows\generate\x64\Debug\ctest.exe
        ..\windows\generate\Debug\ctest.exe
    ) do (
        if exist %%f (
            set CTESTEXE=%%f
            goto found_ctest_exe
        )
    )
) else (
    echo.Detected 32-bit Windows.
    for %%f in (
        ..\windows\generate\Release\ctest.exe
        ..\windows\generate\Debug\ctest.exe
    ) do (
        if exist %%f (
            set CTESTEXE=%%f
            goto found_ctest_exe
        )
    )
)
echo.ERROR: Cannot find ctest.exe. 
echo.Use Microsoft Visual Studio to build ..\windows\generate.sln
exit /b 1
:found_ctest_exe
echo.Found executable: %CTESTEXE%
exit /b 0