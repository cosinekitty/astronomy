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


## Types

| Name | Summary |
|---|---|
| [Aberration](doc/-aberration/index.md) | [jvm]<br>enum [Aberration](doc/-aberration/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Aberration](doc/-aberration/index.md)&gt; <br>Aberration calculation options. |
| [Astronomy](doc/-astronomy/index.md) | [jvm]<br>object [Astronomy](doc/-astronomy/index.md) |
| [AstroTime](doc/-astro-time/index.md) | [jvm]<br>class [AstroTime](doc/-astro-time/index.md)<br>A date and time used for astronomical calculations. |
| [AstroVector](doc/-astro-vector/index.md) | [jvm]<br>data class [AstroVector](doc/-astro-vector/index.md)(x: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), y: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), z: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), t: [AstroTime](doc/-astro-time/index.md)) |
| [Body](doc/-body/index.md) | [jvm]<br>enum [Body](doc/-body/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Body](doc/-body/index.md)&gt; <br>The enumeration of celestial bodies supported by Astronomy Engine. |
| [Direction](doc/-direction/index.md) | [jvm]<br>enum [Direction](doc/-direction/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Direction](doc/-direction/index.md)&gt; <br>Selects whether to search for a rising event or a setting event for a celestial body. |
| [EquatorEpoch](doc/-equator-epoch/index.md) | [jvm]<br>enum [EquatorEpoch](doc/-equator-epoch/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[EquatorEpoch](doc/-equator-epoch/index.md)&gt; <br>Selects the date for which the Earth's equator is to be used for representing equatorial coordinates. |
| [Equatorial](doc/-equatorial/index.md) | [jvm]<br>class [Equatorial](doc/-equatorial/index.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vec: [AstroVector](doc/-astro-vector/index.md))<br>Equatorial angular and cartesian coordinates. |
| [JupiterMoonsInfo](doc/-jupiter-moons-info/index.md) | [jvm]<br>class [JupiterMoonsInfo](doc/-jupiter-moons-info/index.md)(moon: [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[StateVector](doc/-state-vector/index.md)&gt;)<br>Holds the positions and velocities of Jupiter's major 4 moons. |
| [Observer](doc/-observer/index.md) | [jvm]<br>data class [Observer](doc/-observer/index.md)(latitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), longitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), height: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>The location of an observer on (or near) the surface of the Earth. |
| [Refraction](doc/-refraction/index.md) | [jvm]<br>enum [Refraction](doc/-refraction/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Refraction](doc/-refraction/index.md)&gt; <br>Selects whether to correct for atmospheric refraction, and if so, how. |
| [RotationMatrix](doc/-rotation-matrix/index.md) | [jvm]<br>class [RotationMatrix](doc/-rotation-matrix/index.md)(rot: [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[DoubleArray](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double-array/index.html)&gt;)<br>A rotation matrix that can be used to transform one coordinate system to another. |
| [Spherical](doc/-spherical/index.md) | [jvm]<br>data class [Spherical](doc/-spherical/index.md)(lat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), lon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Spherical coordinates: latitude, longitude, distance. |
| [StateVector](doc/-state-vector/index.md) | [jvm]<br>data class [StateVector](doc/-state-vector/index.md)(x: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), y: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), z: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vx: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vy: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vz: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), t: [AstroTime](doc/-astro-time/index.md)) |
| [Visibility](doc/-visibility/index.md) | [jvm]<br>enum [Visibility](doc/-visibility/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Visibility](doc/-visibility/index.md)&gt; <br>Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening. |

## Functions

| Name | Summary |
|---|---|
| [degreesToRadians](doc/degrees-to-radians.md) | [jvm]<br>fun [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html).[degreesToRadians](doc/degrees-to-radians.md)(): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Convert an angle expressed in degrees to an angle expressed in radians. |
| [radiansToDegrees](doc/radians-to-degrees.md) | [jvm]<br>fun [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html).[radiansToDegrees](doc/radians-to-degrees.md)(): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Convert an angle expressed in radians to an angle expressed in degrees. |
| [times](doc/times.md) | [jvm]<br>operator fun [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html).[times](doc/times.md)(vec: [AstroVector](doc/-astro-vector/index.md)): [AstroVector](doc/-astro-vector/index.md) |

## Properties

| Name | Summary |
|---|---|
| [DEG2RAD](doc/-d-e-g2-r-a-d.md) | [jvm]<br>const val [DEG2RAD](doc/-d-e-g2-r-a-d.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 0.017453292519943295<br>The factor to convert degrees to radians = pi/180. |
| [RAD2DEG](doc/-r-a-d2-d-e-g.md) | [jvm]<br>const val [RAD2DEG](doc/-r-a-d2-d-e-g.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 57.29577951308232<br>The factor to convert radians to degrees = 180/pi. |
