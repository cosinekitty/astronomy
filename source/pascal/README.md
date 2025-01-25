# Astronomy Engine for FreePascal language
This is the programming reference for the [FreePascal](https://www.freepascal.org/) version of Astronomy Engine. It can be used directly from FreePascal (and probably also Delphi) programs. Other programming languages are supported. See the [home page](https://github.com/cosinekitty/astronomy) for more info.

---

## Quick Start
To include Astronomy Engine in your own FreePascal program, all you need is the
file `Astronomy.pas` from this directory.

To get started quickly, here are some [examples](../../demo/pascal/).

## Pascal developer note
The `Astronomy.pas` Pascal unit (module) is a [wrapper library](https://en.wikipedia.org/wiki/Wrapper_library) that bridges the C and Pascal programming languages so that the Astronomy Engine C/C++ library can be used in the Pascal language. The C library binary (`Astronomy.dll` on Windows) is also required to use this Pascal module. This means one needs to put the binary file in the same directory where the Pascal program executable is located.

To build Astronomy Engine C/C++ library (for Windows) check the [documentation](../c/README.Windows.md).
