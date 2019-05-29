@echo off
setlocal EnableDelayedExpansion

echo.
if not defined GENEXE (
    set GENEXE=%cd%\bin\generate.exe
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

if exist generate_c_docs (
    echo.Generating C documentation.
    cd ..\source\c
    for %%d in (html xml latex) do (
        if exist %%d (
            echo.Deleting directory: %%d
            rd /s/q %%d
        )
    )
    doxygen Doxyfile > doxygen.log
    if errorlevel 1 (
        echo.Error in doxygen
        exit /b 1
    )
    if not exist xml (
        echo.Missing directory xml. Was supposed to be created by doxygen.
        exit /b 1
    )
    cd xml
    echo.Merging doxygen xml files.
    xsltproc combine.xslt index.xml > all.xml
    if errorlevel 1 (
        echo.Error in xsltproc. That means we could not merge doxygen xml files.
        exit /b 1
    )
    cd ..
    echo.Translating doxygen XML to Markdown.
    node ..\..\hydrogen\hydrogen.js ..\..\hydrogen\c_prefix.md .\xml\all.xml
    if errorlevel 1 (
        echo.Error in hydrogen.js.
        exit /b 1
    )
    cd ..\..\generate
    node eol_hack.js ..\source\c\README.md
    if errorlevel 1 (
        echo.Error in eol_hack.js
        exit /b 1
    )
) else (
    echo.Skipping generation of C documentation because file generate_c_docs does not exist.
)

exit /b 0