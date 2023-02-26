@echo off
setlocal EnableDelayedExpansion
call verify_clean.bat || exit /b 1
call run.bat || exit /b 1
call verify_clean.bat || exit /b 1
git push || exit /b 1
exit /b 0
