//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Spherical](index.md)

# Spherical

data class [Spherical](index.md)(lat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), lon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html))

Spherical coordinates: latitude, longitude, distance.

## Constructors

| | |
|---|---|
| [Spherical](-spherical.md)<br>fun [Spherical](-spherical.md)(lat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), lon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)) |

## Functions

| Name | Summary |
|---|---|
| [toVector](to-vector.md)<br>fun [toVector](to-vector.md)(time: [Time](../-time/index.md)): [Vector](../-vector/index.md)<br>Converts spherical coordinates to Cartesian coordinates. |
| [toVectorFromHorizon](to-vector-from-horizon.md)<br>fun [toVectorFromHorizon](to-vector-from-horizon.md)(time: [Time](../-time/index.md), refraction: [Refraction](../-refraction/index.md)): [Vector](../-vector/index.md)<br>Given apparent angular horizontal coordinates, calculate the unrefracted horizontal vector. |

## Properties

| Name | Summary |
|---|---|
| [dist](dist.md)<br>val [dist](dist.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Distance in AU. |
| [lat](lat.md)<br>val [lat](lat.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>The latitude angle: -90..+90 degrees. |
| [lon](lon.md)<br>val [lon](lon.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>The longitude angle: 0..360 degrees. |
