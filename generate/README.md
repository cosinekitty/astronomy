# Source Generator

This directory contains code and data for generating the various
language implementations of Astronomy Engine.

It is only needed by contributors who want to make enhancements
or fix bugs.  People who just want to use Astronomy Engine
for a given programming langauge can safely ignore this directory
and use the [source code that has already been generated for that language](../source).

---

# Linux and macOS

## Tool setup

The following tools are required for developers:
- gcc
- Node.js
- npm (needed to install other tools like TypeScript, jsdoc2md, ...)
- Python 3.7+
- Microsoft .NET 6.0 SDK
- doxygen
- xsltproc
- coreutils

Change into the directory `hydrogen` and execute:  `npm init`.

## Build process

Once you have all the tools installed and configured, you are ready to proceed.

Change into the `generate` directory (this directory) and run the bash script
`./run` to rebuild all code, generate all documentation, and run all the unit tests.

---

# Windows

## Tool setup

The following tools are required for developers:
- Microsoft .NET 6.0 SDK
- Node.js
- npm (needed to install other tools like TypeScript, jsdoc2md, ...)
- Python 3.7+
- doxygen
- xsltproc  (Follow instructions at https://www.zlatkovic.com/libxml.en.html)

    Hint for 64-bit Windows: You will need to download the following archives from [here](https://www.zlatkovic.com/pub/libxml/64bit/): 
    - iconv-1.14-win32-x86_64.7z
    - libtool-2.4.6-win32-x86_64.7z
    - libxml2-2.9.3-win32-x86_64.7z
    - libxslt-1.1.28-win32-x86_64.7z
    - mingwrt-5.2.0-win32-x86_64.7z
    
    Use [7-Zip](https://www.7-zip.org/) to expand these archives into a newly-created
    empty directory. Just unzip them all to the same place. Verify that the `bin`
    directory beneath your new directory contains `xsltproc.exe` and a bunch of DLLs.
    Either copy the contents of `bin` to somewhere in your `PATH`, or add this `bin`
    directory to your `PATH`. Verify that you can run `xsltproc.exe` from the command
    prompt without any popups about missing DLLs.

## Build process

Once you have all the tools installed and configured, you are ready to proceed.

Change into the `generate` directory (this directory) and run the
batch file `run.bat` to rebuild all code, generate all documentation,
and run all the unit tests.
