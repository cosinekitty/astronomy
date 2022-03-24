//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)

# Package io.github.cosinekitty.astronomy

## Types

| Name | Summary |
|---|---|
| [Aberration](-aberration/index.md) | [jvm]<br>enum [Aberration](-aberration/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Aberration](-aberration/index.md)&gt; <br>Aberration calculation options. |
| [Astronomy](-astronomy/index.md) | [jvm]<br>object [Astronomy](-astronomy/index.md) |
| [AstroTime](-astro-time/index.md) | [jvm]<br>class [AstroTime](-astro-time/index.md)<br>A date and time used for astronomical calculations. |
| [AstroVector](-astro-vector/index.md) | [jvm]<br>data class [AstroVector](-astro-vector/index.md)(x: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), y: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), z: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), t: [AstroTime](-astro-time/index.md)) |
| [Body](-body/index.md) | [jvm]<br>enum [Body](-body/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Body](-body/index.md)&gt; <br>The enumeration of celestial bodies supported by Astronomy Engine. |
| [Direction](-direction/index.md) | [jvm]<br>enum [Direction](-direction/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Direction](-direction/index.md)&gt; <br>Selects whether to search for a rising event or a setting event for a celestial body. |
| [EquatorEpoch](-equator-epoch/index.md) | [jvm]<br>enum [EquatorEpoch](-equator-epoch/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[EquatorEpoch](-equator-epoch/index.md)&gt; <br>Selects the date for which the Earth's equator is to be used for representing equatorial coordinates. |
| [Equatorial](-equatorial/index.md) | [jvm]<br>class [Equatorial](-equatorial/index.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vec: [AstroVector](-astro-vector/index.md))<br>Equatorial angular and cartesian coordinates. |
| [JupiterMoonsInfo](-jupiter-moons-info/index.md) | [jvm]<br>class [JupiterMoonsInfo](-jupiter-moons-info/index.md)(moon: [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[StateVector](-state-vector/index.md)&gt;)<br>Holds the positions and velocities of Jupiter's major 4 moons. |
| [Observer](-observer/index.md) | [jvm]<br>data class [Observer](-observer/index.md)(latitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), longitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), height: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>The location of an observer on (or near) the surface of the Earth. |
| [Refraction](-refraction/index.md) | [jvm]<br>enum [Refraction](-refraction/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Refraction](-refraction/index.md)&gt; <br>Selects whether to correct for atmospheric refraction, and if so, how. |
| [RotationMatrix](-rotation-matrix/index.md) | [jvm]<br>class [RotationMatrix](-rotation-matrix/index.md)(rot: [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[DoubleArray](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double-array/index.html)&gt;)<br>A rotation matrix that can be used to transform one coordinate system to another. |
| [Spherical](-spherical/index.md) | [jvm]<br>data class [Spherical](-spherical/index.md)(lat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), lon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Spherical coordinates: latitude, longitude, distance. |
| [StateVector](-state-vector/index.md) | [jvm]<br>data class [StateVector](-state-vector/index.md)(x: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), y: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), z: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vx: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vy: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vz: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), t: [AstroTime](-astro-time/index.md)) |
| [Visibility](-visibility/index.md) | [jvm]<br>enum [Visibility](-visibility/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Visibility](-visibility/index.md)&gt; <br>Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening. |

## Functions

| Name | Summary |
|---|---|
| [times](times.md) | [jvm]<br>operator fun [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html).[times](times.md)(vec: [AstroVector](-astro-vector/index.md)): [AstroVector](-astro-vector/index.md) |

## Properties

| Name | Summary |
|---|---|
| [DEG2RAD](-d-e-g2-r-a-d.md) | [jvm]<br>const val [DEG2RAD](-d-e-g2-r-a-d.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 0.017453292519943295<br>The factor to convert degrees to radians = pi/180. |
| [RAD2DEG](-r-a-d2-d-e-g.md) | [jvm]<br>const val [RAD2DEG](-r-a-d2-d-e-g.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 57.29577951308232<br>The factor to convert radians to degrees = 180/pi. |
