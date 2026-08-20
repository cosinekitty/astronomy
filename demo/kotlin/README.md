# Astronomy Engine examples in Kotlin Multiplatform

This is a Kotlin Multiplatform console application that contains a collection of
examples of how to use the Kotlin version of Astronomy Engine.

The demo targets four platforms:
- **JVM** — standard Java Virtual Machine
- **WASM/WASI** — WebAssembly with WASI, runs via Node.js
- **Linux native** — compiled to a native ELF binary
- **Windows native** — compiled to a native PE executable

The demo program is split into separate source files, one for each topic, as listed below.
Each topic is chosen by a command line parameter.

## Building

Build all targets:

```
./gradlew assemble
```

This requires the library at `../../source/kotlin` to be available (it is referenced
via Gradle composite build).

## Running

### JVM

```
./gradlew jvmRun --args="seasons 2026"
./gradlew jvmRun --args="moonphase"
./gradlew jvmRun --args="positions 29 -81"
```

### WASM/WASI (via Node.js)

```
./gradlew wasmWasiNodeProductionRun
```

### Linux native

```
./build/bin/linuxX64/releaseExecutable/AstronomyDemo.kexe seasons 2026
```

### Windows native

```
build\bin\mingwX64\releaseExecutable\AstronomyDemo.exe seasons 2026
```

## Usage

Run without arguments (or with `--args=""` on JVM) to see usage text listing
all available demo commands and their parameters.

---

### [Constellation.kt](src/commonMain/kotlin/Constellation.kt)
This demo finds what constellation the Moon is in at a given time.
It also shows how to do a binary search to find the moment in time
when the Moon moves across the border between constellations.

### [JupiterMoons.kt](src/commonMain/kotlin/JupiterMoons.kt)
Calculates the coordinates of Jupiter and its four major moons
(Io, Europa, Ganymede, and Callisto) as seen from the Earth
at a given date and time. This demo illustrates how to correct
for the delay caused by the time it takes for light to reach
the Earth from the Jupiter system.

### [LunarEclipse.kt](src/commonMain/kotlin/LunarEclipse.kt)
Searches for the first 10 lunar eclipses (partial or total)
that occur after the specified time. Penumbral lunar eclipses
are ignored, as these are difficult to observe in practice.

### [MoonPhase.kt](src/commonMain/kotlin/MoonPhase.kt)
This example shows how to determine the Moon's current phase,
and how to predict when the next 10 quarter phases will occur.

### [Positions.kt](src/commonMain/kotlin/Positions.kt)
Given an observer's geographic latitude and longitude,
and an optional date and time, this demo displays the
equatorial and horizontal coordinates of the Sun, Moon, and planets.

### [RiseSetCulm.kt](src/commonMain/kotlin/RiseSetCulm.kt)
Calculates rise, set, and culmination times of the Sun and Moon.
Culmination is when a body reaches the highest
point in an observer's sky as it crosses the meridian.

### [Seasons.kt](src/commonMain/kotlin/Seasons.kt)
Calculates the equinoxes and solstices for a given calendar year.

### [SolarTime.kt](src/commonMain/kotlin/SolarTime.kt)
Displays the true solar time for a given observer location.
