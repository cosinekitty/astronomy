@echo off
setlocal EnableDelayedExpansion
REM -------------------------------------------------------------------------
REM    checksum.bat - by Don Cross / cosinekitty@gmail.com
REM
REM    A wrapper around the standard Windows certutil.exe utility
REM    that allows verifying checkum listings produced by the
REM    Linux utilities sha256sum, md5sum, ....
REM -------------------------------------------------------------------------

if "%1" == "" goto usage
if "%2" == "" goto usage
if not "%3" == "" goto usage
set hashfunc=%1
set textfile=%2

set certutil=
for %%x in (certutil.exe) do (set certutil=%%~$PATH:x)
if not defined certutil (
    echo ERROR: could not find certutil.exe
    exit /b 2
)

if not exist !textfile! (
    echo.ERROR - textfile "!textfile!" does not exist.
    exit /b 2
)

set /a failcount = 0
for /f "tokens=1,2" %%a in (!textfile!) do (
    if not exist %%b (
        echo.ERROR - File does not exist: %%b
        set /a failcount += 1
    ) else (
        certutil.exe -hashfile %%b !hashfunc! > checksum.tmp
        if errorlevel 1 (
            type checksum.tmp
            echo.ERROR - unexpected error !ErrorLevel! returned by certutil.exe.
            exit /b 2
        )

        set /a lnum = 0
        set sum=
        for /f %%x in (checksum.tmp) do (
            set /a lnum += 1
            if !lnum! == 2 (
                set sum=%%x
            )
        )
        if not defined sum (
            echo.ERROR - missing checksum in certutil.exe output for file %%b
            exit /b 2
        )
        if !sum! == %%a (
            echo.%%b : OK
        ) else (
            type checksum.tmp
            echo.expected   checksum = %%a
            echo.calculated checksum = !sum!
            echo.%%b : FAILURE
            set /a failcount += 1
        )
    )
)
if exist checksum.tmp (del checksum.tmp)

if !failcount! gtr 0 (
    echo.ERROR - !failcount! checksum failures
    exit /b 1
)

echo.SUCCESS
exit /b 0

REM -------------------------------------------------------------------------
:usage
echo.checksum.bat by Don Cross - https://github.com/cosinekitty/astronomy
echo.
echo.USAGE: checksum.bat hashfunc textfile
echo.
echo.where
echo.
echo.hashfunc = sha256, md5, ...
echo.textfile = text file produced by Linux sha256sum, md5sum, ...
echo.
echo.Processes the lines of text inside textfile, and verifies
echo.the checksums inside each.
echo.
echo.Return code:
echo.0 = all checksums verified correctly
echo.1 = at least one checksum failure
echo.2 = some other error (see printed output)
echo.
exit /b 2
