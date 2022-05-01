@echo off
setlocal EnableDelayedExpansion
set RunDemo=java -jar build/libs/AstronomyDemo-1.0.0.jar

echo.
echo.Kotlin demos: starting...
echo.

if exist build ( rd /s/q build )
if exist test ( rd /s/q test )
md test

call gradlew.bat jar
if errorlevel 1 (
    echo Cannot build Kotlin demo application.
    exit /b 1
)

REM ----------------------------------------------------------------------------------
REM MoonPhase

!RunDemo! moonphase 2019-06-15T09:15:32.987Z > test/moonphase.txt
if errorlevel 1 (
    echo Error running Kotlin demo: moonphase
    exit /b 1
)
fc correct\moonphase.txt test\moonphase.txt
if errorlevel 1 (
    echo Incorrect output for Kotlin demo: moonphase
    exit /b 1
)

REM ----------------------------------------------------------------------------------
REM Seasons

!RunDemo! seasons 2019 > test/seasons.txt
if errorlevel 1 (
    echo Error running Kotlin demo: seasons
    exit /b 1
)
fc correct\seasons.txt test\seasons.txt
if errorlevel 1 (
    echo Incorrect output for Kotlin demo: seasons
    exit /b 1
)

REM ----------------------------------------------------------------------------------

echo.
echo.Kotlin demos: PASS
echo.
exit /b 0
