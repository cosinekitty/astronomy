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
- gcc and g++
- Node.js
- npm (needed to install other tools like TypeScript, jsdoc2md, ...)
- Python 3.7+ and the following modules:
    - [mypy](https://pypi.org/project/mypy/)
- Microsoft .NET 7 SDK
- doxygen
- xsltproc
- coreutils
- cppcheck
- Java Developer Kit (JDK), for Kotlin.
    - Hint for quick start: install [Android Developer Studio](https://developer.android.com/studio)
    or the Community version of [IntelliJ IDEA](https://www.jetbrains.com/idea/).
    - Set the environment variable `JAVA_HOME` to point to the JDK that comes bundled with either.
    - For example, if you install Android Developer Studio to `/home/yourname/android_studio`,
    you can add the following to your `.bash_aliases` file: `export JAVA_HOME=/home/yourname/android_studio/jre`.

## Build process

Once you have all the tools installed and configured, you are ready to proceed.

Change into the `generate` directory (this directory) and run the bash script
`./run` to rebuild all code, generate all documentation, and run all the unit tests.

---

# Windows

## Tool setup

The following tools are required for developers:
- Microsoft .NET 7 SDK
- Node.js
- npm (needed to install other tools like TypeScript, jsdoc2md, ...)
- Python 3.7+ and the following modules:
    - [mypy](https://pypi.org/project/mypy/)
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

- Java Developer Kit (JDK), for Kotlin.
    - Hint for quick start: install [Android Developer Studio](https://developer.android.com/studio)
    or the Community version of [IntelliJ IDEA](https://www.jetbrains.com/idea/).
    - Set the environment variable `JAVA_HOME` to point to the JDK that comes bundled with either.
    - For example, if you install IntelliJ IDEA to its default location, you can set your user
      environment variable `JAVA_HOME` to the following (adjusting for your actual version):

      `C:\Program Files\JetBrains\IntelliJ IDEA Community Edition 2021.3.3\jbr`


## Build process

Once you have all the tools installed and configured, you are ready to proceed.

Change into the `generate` directory (this directory) and run the
batch file `run.bat` to rebuild all code, generate all documentation,
and run all the unit tests.
