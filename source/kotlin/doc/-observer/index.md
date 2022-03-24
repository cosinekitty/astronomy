//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Observer](index.md)

# Observer

[jvm]\
data class [Observer](index.md)(latitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), longitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), height: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

The location of an observer on (or near) the surface of the Earth.

This object is passed to functions that calculate phenomena as observed from a particular place on the Earth.

## Constructors

| | |
|---|---|
| [Observer](-observer.md) | [jvm]<br>fun [Observer](-observer.md)(latitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), longitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), height: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [height](height.md) | [jvm]<br>val [height](height.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The height above (positive) or below (negative) sea level, expressed in meters. |
| [latitude](latitude.md) | [jvm]<br>val [latitude](latitude.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Geographic latitude in degrees north (positive) or south (negative) of the equator. |
| [longitude](longitude.md) | [jvm]<br>val [longitude](longitude.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England. |
