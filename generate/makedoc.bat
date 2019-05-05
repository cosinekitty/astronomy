@echo off
echo.
echo.Making documentation files.
call jsdoc2md --template ../jsdoc2md/js.hbs --files ..\source\js\astronomy.js > ..\source\js\README.md
if errorlevel 1 (
    echo.FATAL: error in jsdoc2md
    exit /b 1
)

node eol_hack.js ..\source\js\README.md
if errorlevel 1 (
    echo.FATAL: error in eol_hack.js
    exit /b 1
)

exit /b 0