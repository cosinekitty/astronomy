//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[AxisInfo](index.md)

# AxisInfo

[jvm]\
class [AxisInfo](index.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), spin: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), north: [AstroVector](../-astro-vector/index.md))

Information about a body's rotation axis at a given time.

This structure is returned by #Astronomy.rotationAxis to report the orientation of a body's rotation axis at a given moment in time. The axis is specified by the direction in space that the body's north pole points, using angular equatorial coordinates in the J2000 system (EQJ).

Thus ra is the right ascension, and dec is the declination, of the body's north pole vector at the given moment in time. The north pole of a body is defined as the pole that lies on the north side of the [Solar System's invariable plane](https://en.wikipedia.org/wiki/Invariable_plane), regardless of the body's direction of rotation.

The spin field indicates the angular position of a prime meridian arbitrarily recommended for the body by the International Astronomical Union (IAU).

The fields ra, dec, and spin correspond to the variables α0, δ0, and W, respectively, from [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).

The field north is a unit vector pointing in the direction of the body's north pole. It is expressed in the equatorial J2000 system (EQJ).

## Constructors

| | |
|---|---|
| [AxisInfo](-axis-info.md) | [jvm]<br>fun [AxisInfo](-axis-info.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), spin: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), north: [AstroVector](../-astro-vector/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [dec](dec.md) | [jvm]<br>val [dec](dec.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The J2000 declination of the body's north pole direction, in degrees. |
| [north](north.md) | [jvm]<br>val [north](north.md): [AstroVector](../-astro-vector/index.md)<br>A J2000 dimensionless unit vector pointing in the direction of the body's north pole. |
| [ra](ra.md) | [jvm]<br>val [ra](ra.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The J2000 right ascension of the body's north pole direction, in sidereal hours. |
| [spin](spin.md) | [jvm]<br>val [spin](spin.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Rotation angle of the body's prime meridian, in degrees. |
