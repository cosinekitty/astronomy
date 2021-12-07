# Source Generator

This directory contains code and data for generating the various
language implementations of Astronomy Engine.

It is only needed by contributors who want to make enhancements
or fix bugs.  People who just want to use Astronomy Engine
for a given programming langauge can safely ignore this directory
and use the [source code that has already been generated for that language](../source).

---

# Linux

## Tool setup

The following tools are required for developers:
- gcc
- Node.js
- npm
- Python 3.7+
- Microsoft .NET 6.0 SDK
- jsdoc2md
- doxygen
- xsltproc
- graphviz

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
- npm
- Python 3.7+
- jsdoc2md
- doxygen
- xsltproc  (Follow instructions at https://www.zlatkovic.com/libxml.en.html)
- graphviz

Change into the directory `hydrogen` and execute:  `npm init`.

## Build process

Once you have all the tools installed and configured, you are ready to proceed.

Change into the `generate` directory (this directory) and run the
batch file `run.bat` to rebuild all code, generate all documentation,
and run all the unit tests.

---

# Mac

I could use some help getting this to work on the Mac. If you are looking
for an open source project to help with, and you have a Mac, here is your chance!
This should be similar to the Linux steps, but may require some tweaks.
See issue #142 for more information.

