@echo off
setlocal EnableDelayedExpansion
set RunDemo=java -jar build/libs/astronomy-demo-1.0.0.jar

echo.
echo.Java demos: starting...
echo.

if exist build ( rd /s/q build )
if exist test ( rd /s/q test )
md test

call gradlew.bat jar test
if errorlevel 1 (
    echo Cannot build/test jar file.
    type build\test-results\test\TEST-io.github.cosinekitty.astronomy.demo.MainTests.xml
    exit /b 1
)

REM ----------------------------------------------------------------------------------
REM Now

!RunDemo! now
if errorlevel 1 (
    echo Error running Java demo: now
    exit /b 1
)

REM ----------------------------------------------------------------------------------
REM MoonPhase

!RunDemo! moonphase 2019-06-15T09:15:32.987Z > test\moonphase.txt
if errorlevel 1 (
    echo Error running Java demo: moonphase
    exit /b 1
)
fc correct\moonphase.txt test\moonphase.txt
if errorlevel 1 (
    echo Incorrect output for Java demo: moonphase
    exit /b 1
)

REM ----------------------------------------------------------------------------------
REM Seasons

!RunDemo! seasons 2019 > test\seasons.txt
if errorlevel 1 (
    echo Error running Java demo: seasons
    exit /b 1
)
fc correct\seasons.txt test\seasons.txt
if errorlevel 1 (
    echo Incorrect output for Java demo: seasons
    exit /b 1
)

REM ----------------------------------------------------------------------------------

echo.
echo.Java demos: PASS
echo.
exit /b 0
