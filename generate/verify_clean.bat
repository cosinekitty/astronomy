@echo off
setlocal EnableDelayedExpansion

set clean=true
for /f %%x in ('git status --porcelain') do (
    set clean=false
)

if !clean! == false (
    echo.verify_clean: There are local file changes - build is not clean:
    git status
    git diff
    exit /b 1
)

echo.verify_clean: OK
exit /b 0