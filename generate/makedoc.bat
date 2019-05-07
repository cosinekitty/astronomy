@echo off
setlocal EnableDelayedExpansion

echo.
if not defined GENEXE (
    call findgenexe.bat
    if errorlevel 1 (exit /b 1)
)

if not exist "!GENEXE!" (
    echo.FATAL[makedoc.bat]: The executable does not exist: !GENEXE!
    exit /b 1
)

echo.Generating target code.
!GENEXE! source
if errorlevel 1 (
    echo.FATAL:  !GENEXE! source
    exit /b 1
)

echo.
echo.Making documentation files in Markdown format.
call jsdoc2md --separators --template ../jsdoc2md/js.hbs --files ..\source\js\astronomy.js > ..\source\js\README.md
if errorlevel 1 (
    echo.FATAL: error in jsdoc2md
    exit /b 1
)

node eol_hack.js ..\source\js\README.md
if errorlevel 1 (
    echo.FATAL: error in eol_hack.js
    exit /b 1
)

echo.Making documentation in HTML format for local viewing.
if exist html (
    rd /s/q html
)
call jsdoc ../source/js/astronomy.js --destination html
if errorlevel 1 (
    echo.FATAL: error in jsdoc
)

exit /b 0