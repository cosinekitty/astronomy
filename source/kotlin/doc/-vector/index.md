//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Vector](index.md)

# Vector

data class [Vector](index.md)(x: Double, y: Double, z: Double, t: [Time](../-time/index.md))

A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).

## Constructors

| | |
|---|---|
| [Vector](-vector.md)<br>fun [Vector](-vector.md)(x: Double, y: Double, z: Double, t: [Time](../-time/index.md)) |

## Functions

| Name | Summary |
|---|---|
| [angleWith](angle-with.md)<br>fun [angleWith](angle-with.md)(other: [Vector](index.md)): Double<br>Calculates the angle in degrees (0..180) between two vectors. |
| [div](div.md)<br>operator fun [div](div.md)(denom: Double): [Vector](index.md)<br>Divides a vector by a scalar. |
| [dot](dot.md)<br>infix fun [dot](dot.md)(other: [Vector](index.md)): Double<br>Takes the dot product of two vectors. |
| [length](length.md)<br>fun [length](length.md)(): Double<br>The total distance in AU represented by this vector. |
| [minus](minus.md)<br>operator fun [minus](minus.md)(other: [Vector](index.md)): [Vector](index.md)<br>Subtracts one vector from another. Both operands must have identical times. |
| [plus](plus.md)<br>operator fun [plus](plus.md)(other: [Vector](index.md)): [Vector](index.md)<br>Adds two vectors. Both operands must have identical times. |
| [toEquatorial](to-equatorial.md)<br>fun [toEquatorial](to-equatorial.md)(): [Equatorial](../-equatorial/index.md)<br>Given an equatorial vector, calculates equatorial angular coordinates. |
| [toHorizontal](to-horizontal.md)<br>fun [toHorizontal](to-horizontal.md)(refraction: [Refraction](../-refraction/index.md)): [Spherical](../-spherical/index.md)<br>Converts Cartesian coordinates to horizontal coordinates. |
| [toObserver](to-observer.md)<br>fun [toObserver](to-observer.md)(equator: [EquatorEpoch](../-equator-epoch/index.md)): [Observer](../-observer/index.md)<br>Calculates the geographic location corresponding to a geocentric equatorial vector. |
| [toSpherical](to-spherical.md)<br>fun [toSpherical](to-spherical.md)(): [Spherical](../-spherical/index.md)<br>Converts Cartesian coordinates to spherical coordinates. |
| [unaryMinus](unary-minus.md)<br>operator fun [unaryMinus](unary-minus.md)(): [Vector](index.md)<br>Negates a vector; the same as multiplying the vector by the scalar -1. |
| [withTime](with-time.md)<br>fun [withTime](with-time.md)(time: [Time](../-time/index.md)): [Vector](index.md)<br>Creates a new vector with the same coordinates but a different time. |

## Properties

| Name | Summary |
|---|---|
| [t](t.md)<br>val [t](t.md): [Time](../-time/index.md)<br>The date and time at which this vector is valid. |
| [x](x.md)<br>val [x](x.md): Double<br>A Cartesian x-coordinate expressed in AU. |
| [y](y.md)<br>val [y](y.md): Double<br>A Cartesian y-coordinate expressed in AU. |
| [z](z.md)<br>val [z](z.md): Double<br>A Cartesian z-coordinate expressed in AU. |
