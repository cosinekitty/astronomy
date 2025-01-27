# Astronomy Engine (C/C++) for Windows platform

We include the project file (.vcxproj) and solution file (.sln) so that you can easily build the Astronomy Engine on the Windows platform using Visual Studio.

## Build

Open and build the `Astronomy.sln` file in Visual Studio (2022 and newer).
Please note that the `Desktop development with C++` workload is required to build the Astronomy Engine in Visual Studio.
A dynamic link library file `Astronomy.dll` is created on a successful build as an output.

## Language bindings possibility

The Astronomy Engine (C/C++) project is built in a way that allows to create [language binding](https://en.wikipedia.org/wiki/Language_binding) for other programming languages that support _external functions_. This way one can call C/C++ functions from Pascal language for example.
