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
| [ApsisInfo](doc/-apsis-info/index.md) | [jvm]<br>class [ApsisInfo](doc/-apsis-info/index.md)(time: [AstroTime](doc/-astro-time/index.md), kind: [ApsisKind](doc/-apsis-kind/index.md), distAu: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), distKm: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>An apsis event: pericenter (closest approach) or apocenter (farthest distance). |
| [ApsisKind](doc/-apsis-kind/index.md) | [jvm]<br>enum [ApsisKind](doc/-apsis-kind/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[ApsisKind](doc/-apsis-kind/index.md)&gt; <br>The type of apsis: pericenter (closest approach) or apocenter (farthest distance). |
| [Astronomy](doc/-astronomy/index.md) | [jvm]<br>object [Astronomy](doc/-astronomy/index.md)<br>The main container of astronomy calculation functions. |
| [AstroTime](doc/-astro-time/index.md) | [jvm]<br>class [AstroTime](doc/-astro-time/index.md)<br>A date and time used for astronomical calculations. |
| [AstroVector](doc/-astro-vector/index.md) | [jvm]<br>data class [AstroVector](doc/-astro-vector/index.md)(x: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), y: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), z: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), t: [AstroTime](doc/-astro-time/index.md)) |
| [AxisInfo](doc/-axis-info/index.md) | [jvm]<br>class [AxisInfo](doc/-axis-info/index.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), spin: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), north: [AstroVector](doc/-astro-vector/index.md))<br>Information about a body's rotation axis at a given time. |
| [Body](doc/-body/index.md) | [jvm]<br>enum [Body](doc/-body/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Body](doc/-body/index.md)&gt; <br>The enumeration of celestial bodies supported by Astronomy Engine. |
| [ConstellationInfo](doc/-constellation-info/index.md) | [jvm]<br>class [ConstellationInfo](doc/-constellation-info/index.md)(symbol: [String](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-string/index.html), name: [String](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-string/index.html), ra1875: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec1875: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Reports the constellation that a given celestial point lies within. |
| [Direction](doc/-direction/index.md) | [jvm]<br>enum [Direction](doc/-direction/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Direction](doc/-direction/index.md)&gt; <br>Selects whether to search for a rising event or a setting event for a celestial body. |
| [EclipseEvent](doc/-eclipse-event/index.md) | [jvm]<br>class [EclipseEvent](doc/-eclipse-event/index.md)(time: [AstroTime](doc/-astro-time/index.md), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Holds a time and the observed altitude of the Sun at that time. |
| [EclipseKind](doc/-eclipse-kind/index.md) | [jvm]<br>enum [EclipseKind](doc/-eclipse-kind/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[EclipseKind](doc/-eclipse-kind/index.md)&gt; <br>The different kinds of lunar/solar eclipses. |
| [Ecliptic](doc/-ecliptic/index.md) | [jvm]<br>data class [Ecliptic](doc/-ecliptic/index.md)(vec: [AstroVector](doc/-astro-vector/index.md), elat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), elon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Ecliptic angular and Cartesian coordinates. |
| [ElongationInfo](doc/-elongation-info/index.md) | [jvm]<br>class [ElongationInfo](doc/-elongation-info/index.md)(time: [AstroTime](doc/-astro-time/index.md), visibility: [Visibility](doc/-visibility/index.md), elongation: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), eclipticSeparation: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Contains information about the visibility of a celestial body at a given date and time. |
| [EquatorEpoch](doc/-equator-epoch/index.md) | [jvm]<br>enum [EquatorEpoch](doc/-equator-epoch/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[EquatorEpoch](doc/-equator-epoch/index.md)&gt; <br>Selects the date for which the Earth's equator is to be used for representing equatorial coordinates. |
| [Equatorial](doc/-equatorial/index.md) | [jvm]<br>class [Equatorial](doc/-equatorial/index.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vec: [AstroVector](doc/-astro-vector/index.md))<br>Equatorial angular and cartesian coordinates. |
| [GlobalSolarEclipseInfo](doc/-global-solar-eclipse-info/index.md) | [jvm]<br>class [GlobalSolarEclipseInfo](doc/-global-solar-eclipse-info/index.md)(kind: [EclipseKind](doc/-eclipse-kind/index.md), peak: [AstroTime](doc/-astro-time/index.md), distance: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), latitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), longitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Reports the time and geographic location of the peak of a solar eclipse. |
| [HourAngleInfo](doc/-hour-angle-info/index.md) | [jvm]<br>class [HourAngleInfo](doc/-hour-angle-info/index.md)(time: [AstroTime](doc/-astro-time/index.md), hor: [Topocentric](doc/-topocentric/index.md))<br>Information about a celestial body crossing a specific hour angle. |
| [IllumInfo](doc/-illum-info/index.md) | [jvm]<br>class [IllumInfo](doc/-illum-info/index.md)(time: [AstroTime](doc/-astro-time/index.md), mag: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), phaseAngle: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), phaseFraction: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), helioDist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), ringTilt: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Information about the brightness and illuminated shape of a celestial body. |
| [InvalidBodyException](doc/-invalid-body-exception/index.md) | [jvm]<br>class [InvalidBodyException](doc/-invalid-body-exception/index.md)(body: [Body](doc/-body/index.md)) : [Exception](https://docs.oracle.com/javase/8/docs/api/java/lang/Exception.html)<br>An invalid body was specified for the given function. |
| [JupiterMoonsInfo](doc/-jupiter-moons-info/index.md) | [jvm]<br>class [JupiterMoonsInfo](doc/-jupiter-moons-info/index.md)(moon: [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[StateVector](doc/-state-vector/index.md)&gt;)<br>Holds the positions and velocities of Jupiter's major 4 moons. |
| [LibrationInfo](doc/-libration-info/index.md) | [jvm]<br>data class [LibrationInfo](doc/-libration-info/index.md)(elat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), elon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), mlat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), mlon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), distKm: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), diamDeg: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Lunar libration angles, returned by #Astronomy.libration. |
| [LocalSolarEclipseInfo](doc/-local-solar-eclipse-info/index.md) | [jvm]<br>class [LocalSolarEclipseInfo](doc/-local-solar-eclipse-info/index.md)(kind: [EclipseKind](doc/-eclipse-kind/index.md), partialBegin: [EclipseEvent](doc/-eclipse-event/index.md), totalBegin: [EclipseEvent](doc/-eclipse-event/index.md), peak: [EclipseEvent](doc/-eclipse-event/index.md), totalEnd: [EclipseEvent](doc/-eclipse-event/index.md), partialEnd: [EclipseEvent](doc/-eclipse-event/index.md))<br>Information about a solar eclipse as seen by an observer at a given time and geographic location. |
| [LunarEclipseInfo](doc/-lunar-eclipse-info/index.md) | [jvm]<br>class [LunarEclipseInfo](doc/-lunar-eclipse-info/index.md)(kind: [EclipseKind](doc/-eclipse-kind/index.md), peak: [AstroTime](doc/-astro-time/index.md), sdPenum: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), sdPartial: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), sdTotal: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Information about a lunar eclipse. |
| [MoonQuarterInfo](doc/-moon-quarter-info/index.md) | [jvm]<br>class [MoonQuarterInfo](doc/-moon-quarter-info/index.md)(quarter: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), time: [AstroTime](doc/-astro-time/index.md))<br>A lunar quarter event (new moon, first quarter, full moon, or third quarter) along with its date and time. |
| [NodeEventInfo](doc/-node-event-info/index.md) | [jvm]<br>class [NodeEventInfo](doc/-node-event-info/index.md)(time: [AstroTime](doc/-astro-time/index.md), kind: [NodeEventKind](doc/-node-event-kind/index.md))<br>Information about an ascending or descending node of a body. |
| [NodeEventKind](doc/-node-event-kind/index.md) | [jvm]<br>enum [NodeEventKind](doc/-node-event-kind/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[NodeEventKind](doc/-node-event-kind/index.md)&gt; <br>Indicates whether a crossing through the ecliptic plane is ascending or descending. |
| [Observer](doc/-observer/index.md) | [jvm]<br>data class [Observer](doc/-observer/index.md)(latitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), longitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), height: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>The location of an observer on (or near) the surface of the Earth. |
| [Refraction](doc/-refraction/index.md) | [jvm]<br>enum [Refraction](doc/-refraction/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Refraction](doc/-refraction/index.md)&gt; <br>Selects whether to correct for atmospheric refraction, and if so, how. |
| [RotationMatrix](doc/-rotation-matrix/index.md) | [jvm]<br>class [RotationMatrix](doc/-rotation-matrix/index.md)(rot: [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[DoubleArray](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double-array/index.html)&gt;)<br>A rotation matrix that can be used to transform one coordinate system to another. |
| [SearchContext](doc/-search-context/index.md) | [jvm]<br>interface [SearchContext](doc/-search-context/index.md) |
| [SeasonsInfo](doc/-seasons-info/index.md) | [jvm]<br>class [SeasonsInfo](doc/-seasons-info/index.md)(marEquinox: [AstroTime](doc/-astro-time/index.md), junSolstice: [AstroTime](doc/-astro-time/index.md), sepEquinox: [AstroTime](doc/-astro-time/index.md), decSolstice: [AstroTime](doc/-astro-time/index.md))<br>The dates and times of changes of season for a given calendar year. |
| [Spherical](doc/-spherical/index.md) | [jvm]<br>data class [Spherical](doc/-spherical/index.md)(lat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), lon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Spherical coordinates: latitude, longitude, distance. |
| [StateVector](doc/-state-vector/index.md) | [jvm]<br>data class [StateVector](doc/-state-vector/index.md)(x: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), y: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), z: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vx: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vy: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vz: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), t: [AstroTime](doc/-astro-time/index.md)) |
| [Topocentric](doc/-topocentric/index.md) | [jvm]<br>data class [Topocentric](doc/-topocentric/index.md)(azimuth: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Coordinates of a celestial body as seen by a topocentric observer. |
| [TransitInfo](doc/-transit-info/index.md) | [jvm]<br>class [TransitInfo](doc/-transit-info/index.md)(start: [AstroTime](doc/-astro-time/index.md), peak: [AstroTime](doc/-astro-time/index.md), finish: [AstroTime](doc/-astro-time/index.md), separation: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Information about a transit of Mercury or Venus, as seen from the Earth. |
| [Visibility](doc/-visibility/index.md) | [jvm]<br>enum [Visibility](doc/-visibility/index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Visibility](doc/-visibility/index.md)&gt; <br>Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening. |

## Functions

| Name | Summary |
|---|---|
| [degreesToRadians](doc/degrees-to-radians.md) | [jvm]<br>fun [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html).[degreesToRadians](doc/degrees-to-radians.md)(): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Convert an angle expressed in degrees to an angle expressed in radians. |
| [radiansToDegrees](doc/radians-to-degrees.md) | [jvm]<br>fun [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html).[radiansToDegrees](doc/radians-to-degrees.md)(): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Convert an angle expressed in radians to an angle expressed in degrees. |
| [times](doc/times.md) | [jvm]<br>operator fun [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html).[times](doc/times.md)(vec: [AstroVector](doc/-astro-vector/index.md)): [AstroVector](doc/-astro-vector/index.md)<br>Multiply a scalar by a vector, yielding another vector. |

## Properties

| Name | Summary |
|---|---|
| [C_AUDAY](doc/-c_-a-u-d-a-y.md) | [jvm]<br>const val [C_AUDAY](doc/-c_-a-u-d-a-y.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 173.1446326846693<br>The speed of light in AU/day. |
| [CALLISTO_RADIUS_KM](doc/-c-a-l-l-i-s-t-o_-r-a-d-i-u-s_-k-m.md) | [jvm]<br>const val [CALLISTO_RADIUS_KM](doc/-c-a-l-l-i-s-t-o_-r-a-d-i-u-s_-k-m.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 2410.3<br>The mean radius of Jupiter's moon Callisto, expressed in kilometers. |
| [DEG2RAD](doc/-d-e-g2-r-a-d.md) | [jvm]<br>const val [DEG2RAD](doc/-d-e-g2-r-a-d.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 0.017453292519943295<br>The factor to convert degrees to radians = pi/180. |
| [EUROPA_RADIUS_KM](doc/-e-u-r-o-p-a_-r-a-d-i-u-s_-k-m.md) | [jvm]<br>const val [EUROPA_RADIUS_KM](doc/-e-u-r-o-p-a_-r-a-d-i-u-s_-k-m.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 1560.8<br>The mean radius of Jupiter's moon Europa, expressed in kilometers. |
| [GANYMEDE_RADIUS_KM](doc/-g-a-n-y-m-e-d-e_-r-a-d-i-u-s_-k-m.md) | [jvm]<br>const val [GANYMEDE_RADIUS_KM](doc/-g-a-n-y-m-e-d-e_-r-a-d-i-u-s_-k-m.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 2631.2<br>The mean radius of Jupiter's moon Ganymede, expressed in kilometers. |
| [HOUR2RAD](doc/-h-o-u-r2-r-a-d.md) | [jvm]<br>const val [HOUR2RAD](doc/-h-o-u-r2-r-a-d.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 0.26179938779914946<br>The factor to convert sidereal hours to radians = pi/12. |
| [IO_RADIUS_KM](doc/-i-o_-r-a-d-i-u-s_-k-m.md) | [jvm]<br>const val [IO_RADIUS_KM](doc/-i-o_-r-a-d-i-u-s_-k-m.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 1821.6<br>The mean radius of Jupiter's moon Io, expressed in kilometers. |
| [JUPITER_EQUATORIAL_RADIUS_KM](doc/-j-u-p-i-t-e-r_-e-q-u-a-t-o-r-i-a-l_-r-a-d-i-u-s_-k-m.md) | [jvm]<br>const val [JUPITER_EQUATORIAL_RADIUS_KM](doc/-j-u-p-i-t-e-r_-e-q-u-a-t-o-r-i-a-l_-r-a-d-i-u-s_-k-m.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 71492.0<br>The equatorial radius of Jupiter, expressed in kilometers. |
| [JUPITER_MEAN_RADIUS_KM](doc/-j-u-p-i-t-e-r_-m-e-a-n_-r-a-d-i-u-s_-k-m.md) | [jvm]<br>const val [JUPITER_MEAN_RADIUS_KM](doc/-j-u-p-i-t-e-r_-m-e-a-n_-r-a-d-i-u-s_-k-m.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 69911.0<br>The volumetric mean radius of Jupiter, expressed in kilometers. |
| [JUPITER_POLAR_RADIUS_KM](doc/-j-u-p-i-t-e-r_-p-o-l-a-r_-r-a-d-i-u-s_-k-m.md) | [jvm]<br>const val [JUPITER_POLAR_RADIUS_KM](doc/-j-u-p-i-t-e-r_-p-o-l-a-r_-r-a-d-i-u-s_-k-m.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 66854.0<br>The polar radius of Jupiter, expressed in kilometers. |
| [KM_PER_AU](doc/-k-m_-p-e-r_-a-u.md) | [jvm]<br>const val [KM_PER_AU](doc/-k-m_-p-e-r_-a-u.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 1.4959787069098932E8<br>The number of kilometers in one astronomical unit (AU). |
| [RAD2DEG](doc/-r-a-d2-d-e-g.md) | [jvm]<br>const val [RAD2DEG](doc/-r-a-d2-d-e-g.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 57.29577951308232<br>The factor to convert radians to degrees = 180/pi. |
| [RAD2HOUR](doc/-r-a-d2-h-o-u-r.md) | [jvm]<br>const val [RAD2HOUR](doc/-r-a-d2-h-o-u-r.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) = 3.819718634205488<br>The factor to convert radians to sidereal hours = 12/pi. |
