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

echo.Trimming trailing whitespace in source code.
for %%f in (template\astronomy.c ..\source\c\astronomy.h) do (
    node trimspace.js %%f
    if errorlevel 1 (exit /b 1)
)

echo.Generating target code.
!GENEXE! source
if errorlevel 1 (
    echo.FATAL:  !GENEXE! source
    exit /b 1
)

echo.Trimming trailing whitespace in target code.
for %%f in (..\source\c\astronomy.c ..\source\js\astronomy.js ..\source\python\astronomy.py ..\source\csharp\astronomy.cs) do (
    node trimspace.js %%f
    if errorlevel 1 (exit /b 1)
)

REM   C# is a special case. We have to compile the code to get the documentation.
echo.Building C# code to get documentation.
cd dotnet\csharp_test
dotnet build --output !CD!\exe
if errorlevel 1 (exit /b 1)
cd ..\..\csdown
call build.bat
if errorlevel 1 (exit /b 1)
dotnet exe\csdown.dll csharp_prefix.md ..\dotnet\csharp_test\exe\astronomy.dll ..\dotnet\csharp_test\exe\astronomy.xml ..\..\source\csharp\README.md
if errorlevel 1 (
    echo.Error generating C# documentation.
    exit /b 1
)
cd ..
node eol_hack.js ..\source\csharp\README.md
if errorlevel 1 (
    echo.FATAL: error in eol_hack.js [C#]
    exit /b 1
)

echo.Install NodeJS dev dependencies
if not exist node_modules (
    call npm ci
    if errorlevel 1 (
        echo.Error installing NodeJS dev dependencies.
        exit /b 1
    )
)

echo.Compiling TypeScript to JavaScript.
call npm run build
if errorlevel 1 (
    echo.Error in typescript compiler.
    exit /b 1
)

echo.Bundling JavaScript code for Browser.
call npm run build:browser
if errorlevel 1 (
    echo.Error building browser bundle
    exit /b 1
)

echo.Minifying JavaScript code.
call npm run minify
if errorlevel 1 (
    echo.Error minifying astronomy.js
    exit /b 1
)

call npm run minify:browser
if errorlevel 1 (
    echo.Error minifying astronomy.browser.js
    exit /b 1
)

for %%f in (
    ..\source\js\astronomy.js
    ..\source\js\astronomy.min.js
    ..\source\js\astronomy.browser.js
    ..\source\js\astronomy.browser.min.js
) do (
    node eol_hack.js %%f
    if errorlevel 1 (
        echo.ERROR cleaning newlines in file: %%f
        exit /b 1
    )
)

if exist ..\website (
    echo.Generating JS documentation in JSON format.
    call npm run docs:json
    if errorlevel 1 (
        echo.Error generating JSON documentation.
        exit /b 1
    )
    node jsdoc_strip_path.js
    if errorlevel 1 (
        echo.Error stripping absolute paths.
        exit /b 1
    )
)

echo.Generating JS documentation in Markdown format.
call npm run docs:md
if errorlevel 1 (
    echo.Error generating JS documentation.
    exit /b 1
)

node eol_hack.js ..\source\js\README.md
if errorlevel 1 (exit /b 1)

check_internal_links.py ..\source\js\README.md
if errorlevel 1 (exit /b 1)

echo.Making documentation in HTML format for local viewing.
if exist html (
    rd /s/q html
)
call npm run docs:html
if errorlevel 1 (
    echo.FATAL: error in jsdoc
)

if exist disable_generate_c_docs (
    echo.Generation of C documentation is disabled.
) else (
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
    if not exist ..\..\generate\hydrogen\node_modules\xml2js (
        pushd ..\..\generate\hydrogen
        npm ci
        if errorlevel 1 (
            echo.Error installing dependencies for hydrogen.
            exit /b 1
        )
        popd
    )
    echo.Translating doxygen XML to Markdown.
    node ..\..\generate\hydrogen\hydrogen.js ..\..\generate\hydrogen\c_prefix.md .\xml\all.xml
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
)

echo.Generating Python documentation.
pydown\pydown.py pydown\py_prefix.md ..\source\python\astronomy.py ..\source\python\README.md
if errorlevel 1 (exit /b 1)

echo.Making redundant copies of source in demo folders.

copy ..\source\js\astronomy.browser.js ..\demo\browser\
if errorlevel 1 (exit /b 1)

copy ..\source\js\astronomy.js ..\demo\nodejs\
if errorlevel 1 (exit /b 1)

copy ..\source\python\astronomy.py ..\demo\python\
if errorlevel 1 (exit /b 1)

exit /b 0
