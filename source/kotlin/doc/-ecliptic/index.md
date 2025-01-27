//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Ecliptic](index.md)

# Ecliptic

data class [Ecliptic](index.md)(vec: [Vector](../-vector/index.md), elat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), elon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html))

Ecliptic angular and Cartesian coordinates.

Coordinates of a celestial body as seen from the center of the Sun (heliocentric), oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).

## Constructors

| | |
|---|---|
| [Ecliptic](-ecliptic.md)<br>fun [Ecliptic](-ecliptic.md)(vec: [Vector](../-vector/index.md), elat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), elon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [elat](elat.md)<br>val [elat](elat.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Latitude in degrees north (positive) or south (negative) of the ecliptic plane. |
| [elon](elon.md)<br>val [elon](elon.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Longitude in degrees around the ecliptic plane prograde from the equinox. |
| [vec](vec.md)<br>val [vec](vec.md): [Vector](../-vector/index.md)<br>Cartesian ecliptic vector, with components as follows: x: the direction of the equinox along the ecliptic plane. y: in the ecliptic plane 90 degrees prograde from the equinox. z: perpendicular to the ecliptic plane. Positive is north. |
