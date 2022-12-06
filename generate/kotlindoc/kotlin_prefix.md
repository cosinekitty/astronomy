# Astronomy Engine (Kotlin)

[![JitPack build status](https://jitpack.io/v/cosinekitty/astronomy.svg)](https://jitpack.io/#cosinekitty/astronomy)

---

## Quick Start

Here are some [Kotlin examples](../../demo/kotlin/) and [Java examples](../../demo/java/) to get you started.

Feel free to [start a new discussion topic](https://github.com/cosinekitty/astronomy/discussions)
if you need some help in your astronomy-related project, either with this code or astronomy
computation concepts in general.

To include Astronomy Engine in your project, add this in your root `build.gradle.kts` at the end of repositories section:
```kotlin
allprojects {
    repositories {
        ...
        maven("https://jitpack.io")
    }
}
```

Now add the dependency:
```kotlin
dependencies {
    implementation("io.github.cosinekitty:astronomy:[ASTRONOMY_ENGINE_VERSION]")
}
```

For other build tools support have a look at [this](https://jitpack.io/#cosinekitty/astronomy).

---

## Contents

- [Coordinate Transforms](#coords)
- [Gravity Simulator](#gravsim)
- [Types](#types)
- [Functions](#functions)
- [Properties](#properties)

---

<a name="coords"></a>
## Coordinate Transforms

The following orientation systems are supported.
Astronomy Engine can convert a vector from any of these orientations to any of the others.
It also allows converting from a vector to spherical (angular) coordinates and back,
within a given orientation. Note the 3-letter codes for each of the orientation systems;
these are used in function and type names.

- **EQJ = Equatorial J2000**: Uses the Earth's equator on January 1, 2000, at noon UTC.
- **EQD = Equator of Date**: Uses the Earth's equator on a given date and time, adjusted for precession and nutation.
- **ECT = True Ecliptic of Date**: Uses the true orbital plane and equator of the Earth on the given date.
- **ECL = Mean J2000 Ecliptic**: Uses the plane of the Earth's orbit around the Sun in the year 2000. The x-axis is referenced against the J2000 mean equinox.
- **HOR = Horizontal**: Uses the viewpoint of an observer at a specific location on the Earth at a given date and time.
- **GAL = Galactic**: Based on the IAU 1958 definition of galactic coordinates.

<a name="gravsim"></a>
## Gravity Simulator

Astronomy Engine provides a [GravitySimulator](doc/-gravity-simulator/index.md) class
that allows you to model the trajectories of one or more small bodies like asteroids,
comets, or coasting spacecraft. If you know an initial position vector
and velocity vector for a small body, the gravity simulator can incrementally
simulate the pull of gravity on it from the Sun and planets, to calculate its
movement through the Solar System.

---
