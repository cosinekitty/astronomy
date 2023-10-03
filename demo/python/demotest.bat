@echo off
setlocal EnableDelayedExpansion

if not exist test (md test)
if exist test\*.txt (del test\*.txt)

call :TestDemo solar_time +38.88 -77.03 2023-02-12T17:00:00Z || exit /b 1
call :TestDemo jupiter_moons 2021-04-16T00:26:18Z || exit /b 1
call :TestDemo constellation 2021-06-01T00:00:00Z || exit /b 1
call :TestDemo moonphase 2019-06-15T09:15:32.987Z || exit /b 1
call :TestDemo riseset +45.6 -90.7 2018-11-30T17:55:07.234Z || exit /b 1
call :TestDemo positions +45.6 -90.7 2018-11-30T17:55:07.234Z || exit /b 1
call :TestDemo seasons 2019 || exit /b 1
call :TestDemo culminate +30 -90 2015-02-28T00:00:00Z || exit /b 1
call :TestDemo horizon +25.5 -85.3 2016-12-25T12:30:45Z || exit /b 1
call :TestDemo lunar_eclipse 1988-01-01 || exit /b 1
call :TestDemo lunar_angles 2021-05-15 || exit /b 1
call :TestDemo galactic 38.92056 -77.0658 22.793498 197.070510 2025-04-06T00:00:00Z || exit /b 1
call :TestDemo triangulate 48.16042 24.49986 2019 18 7 48.27305 24.36401 662 83 12 || exit /b 1
call :TestDemo stars_near_moon 30 -80 2021-11-08T23:00:00Z || exit /b 1
call :TestDemo ecliptic_of_date 53.6375 9.9981 2023-10-28T20:24:33.489Z || exit /b 1
call :TestDemo solar_eclipse 28.5 -82.5 2023-01-01 || exit /b 1

call :CameraTest _nz -41.17 175.5 2023-03-24T06:30:00Z || exit /b 1
call :CameraTest _fl +29 -81 2023-03-25T00:00:00Z || exit /b 1
call :CameraTest _ca 51.05 -114.07 2023-03-24T02:07:31.000Z || exit /b 1

echo.Testing Python demo: gravity
for /L %%x in (0, 1, 90) do (
    py gravity.py %%x 0 >> test/gravity.txt || Fail "Error running gravity.py."
)
fc correct\gravity.txt test\gravity.txt || exit /b 1
echo.PASS: Python demos
exit /b 0


:TestDemo
    set name=%1
    set args=%2 %3 %4 %5 %6 %7 %8 %9
    REM *** Hack to handle more than 8 arguments needed by triangulate.py
    shift
    set args=!args! %9
    shift
    set args=!args! %9
    echo.Testing Python demo: !name! !args!
    mypy --strict !name!.py || exit /b 1
    echo.Running !name!.py !args!
    py !name!.py !args! > test\!name!.txt || (
        echo.FATAL[demotest.bat] - error returned by !name!.py
        exit /b 1
    )
    fc correct\!name!.txt test\!name!.txt || exit /b 1
    exit /b 0

:CameraTest
    set name=camera
    set suffix=%1
    set args=%2 %3 %4 %5 %6 %7 %8 %9
    echo.Testing Python demo: !name! !args!
    mypy --strict !name!.py || exit /b 1
    echo.Running !name!.py !args!
    py !name!.py !args! > test\!name!!suffix!.txt || (
        echo.FATAL[demotest.bat] - error returned by !name!.py
        exit /b 1
    )
    fc correct\!name!!suffix!.txt test\!name!!suffix!.txt || exit /b 1
    exit /b 0
