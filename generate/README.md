# Source Generator

This directory contains code and data for generating the various
language implementations of Astronomy Engine.

It is only needed by contributors who want to make enhancements
or fix bugs.  People who just want to use the astronomy calculator
for a given programming langauge can safely ignore this directory
and use the source code that has already been generated for that language.

---

# Linux

## Tool setup

The following tools are required for developers:
- gcc
- Node.js
- npm

## Optional documentation generation

The following steps are optional, if you want to generate documentation.

Install these packages:

- jsdoc2md
- doxygen
- xsltproc
- graphviz

In the directory `generate`, create a file `generate_c_docs`; it doesn't matter what it contains.
This enables creating documentation for the C code.

Change into the directory `hydrogen` and execute:  `npm init`.

## Build process

Once you have all the tools installed and configured, you are ready to proceed.

Change into the `generate` directory (this directory) and run the bash script
`./run` to rebuild all code, generate all documentation, and run all the unit tests.

---

# Windows

## Tool setup

The following tools are required for developers:
- Microsoft Visual Studio (2015 or later Community Edition is freely available and works fine).
- Node.js
- npm

## Optional documentation generation

The following steps are optional, if you want to generate documentation.

Install these packages:

- jsdoc2md
- doxygen
- xsltproc  (Follow instructions at https://www.zlatkovic.com/libxml.en.html)
- graphviz

In the directory `generate`, create a file `generate_c_docs`; it doesn't matter what it contains.
This enables creating documentation for the C code.

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
Creating issues is welcome. Creating pull requests is even *more* welcome!
This should be similar to the Linux steps, but may require some tweaks.
