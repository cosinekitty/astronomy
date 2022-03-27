//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[AstroVector](index.md)

# AstroVector

[jvm]\
data class [AstroVector](index.md)(x: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), y: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), z: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), t: [AstroTime](../-astro-time/index.md))

## Functions

| Name | Summary |
|---|---|
| [angleWith](angle-with.md) | [jvm]<br>fun [angleWith](angle-with.md)(other: [AstroVector](index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) |
| [div](div.md) | [jvm]<br>operator fun [div](div.md)(denom: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [AstroVector](index.md) |
| [dot](dot.md) | [jvm]<br>infix fun [dot](dot.md)(other: [AstroVector](index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) |
| [length](length.md) | [jvm]<br>fun [length](length.md)(): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The total distance in AU represented by this vector. |
| [minus](minus.md) | [jvm]<br>operator fun [minus](minus.md)(other: [AstroVector](index.md)): [AstroVector](index.md) |
| [plus](plus.md) | [jvm]<br>operator fun [plus](plus.md)(other: [AstroVector](index.md)): [AstroVector](index.md) |
| [toEquatorial](to-equatorial.md) | [jvm]<br>fun [toEquatorial](to-equatorial.md)(): [Equatorial](../-equatorial/index.md)<br>Given an equatorial vector, calculates equatorial angular coordinates. |
| [toHorizontal](to-horizontal.md) | [jvm]<br>fun [toHorizontal](to-horizontal.md)(refraction: [Refraction](../-refraction/index.md)): [Spherical](../-spherical/index.md)<br>Converts Cartesian coordinates to horizontal coordinates. |
| [toSpherical](to-spherical.md) | [jvm]<br>fun [toSpherical](to-spherical.md)(): [Spherical](../-spherical/index.md)<br>Converts Cartesian coordinates to spherical coordinates. |
| [unaryMinus](unary-minus.md) | [jvm]<br>operator fun [unaryMinus](unary-minus.md)(): [AstroVector](index.md) |

## Properties

| Name | Summary |
|---|---|
| [t](t.md) | [jvm]<br>val [t](t.md): [AstroTime](../-astro-time/index.md)<br>The date and time at which this vector is valid. |
| [x](x.md) | [jvm]<br>val [x](x.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A Cartesian x-coordinate expressed in AU. |
| [y](y.md) | [jvm]<br>val [y](y.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A Cartesian y-coordinate expressed in AU. |
| [z](z.md) | [jvm]<br>val [z](z.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A Cartesian z-coordinate expressed in AU. |
