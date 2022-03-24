# Astronomy Engine (Kotlin)

[![JitPack build status](https://jitpack.io/v/cosinekitty/astronomy.svg)](https://jitpack.io/#cosinekitty/astronomy)

---

## Quick Start

Add this in your root `build.gradle.kts` at the end of repositories section:
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
    implementation("io.github.cosinekitty:astronomy:0.0.1")
}
```

For other build tools support have a look at [this](https://jitpack.io/#cosinekitty/astronomy).

---

## Contents

- [Coordinate Transforms](#coords)
- [Reference](#reference)

---

<a name="coords"></a>
## Coordinate Transforms

The following five orientation systems are supported.
Astronomy Engine can convert a vector from any of these orientations to any of the others.
It also allows converting from a vector to spherical (angular) coordinates and back,
within a given orientation. Note the 3-letter codes for each of the orientation systems;
these are used in function and type names.

- **EQJ = Equatorial J2000**: Uses the Earth's equator on January 1, 2000, at noon UTC.
- **EQD = Equator of-date**: Uses the Earth's equator on a given date and time, adjusted for precession and nutation.
- **ECL = Ecliptic**: Uses the mean plane of the Earth's orbit around the Sun. The x-axis is referenced against the J2000 equinox.
- **HOR = Horizontal**: Uses the viewpoint of an observer at a specific location on the Earth at a given date and time.
- **GAL = Galactic**: Based on the IAU 1958 definition of galactic coordinates.

---

<a name="reference"></a>
## Reference

(More content coming here soon.)

