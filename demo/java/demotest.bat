@echo off
setlocal EnableDelayedExpansion

call gradlew.bat jar
if errorlevel 1 (
    echo Cannot build jar file.
    exit /b 1
)

echo.
echo.Java demos: starting...
echo.

java -jar build/libs/astronomy-demo-0.0.1.jar now
if errorlevel 1 (
    echo Error running unit test: now
    exit /b 1
)

echo.
echo.Java demos: PASS
echo.
exit /b 0
