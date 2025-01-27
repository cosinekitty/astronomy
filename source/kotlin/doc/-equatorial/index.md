//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Equatorial](index.md)

# Equatorial

class [Equatorial](index.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), vec: [Vector](../-vector/index.md))

Equatorial angular and cartesian coordinates.

Coordinates of a celestial body as seen from the Earth (geocentric or topocentric, depending on context), oriented with respect to the projection of the Earth's equator onto the sky.

## Constructors

| | |
|---|---|
| [Equatorial](-equatorial.md)<br>fun [Equatorial](-equatorial.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), vec: [Vector](../-vector/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [dec](dec.md)<br>val [dec](dec.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Declination in degrees. |
| [dist](dist.md)<br>val [dist](dist.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Distance to the celestial body in AU. |
| [ra](ra.md)<br>val [ra](ra.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Right ascension in sidereal hours. |
| [vec](vec.md)<br>val [vec](vec.md): [Vector](../-vector/index.md)<br>Equatorial coordinates in cartesian vector form: x = March equinox, y = June solstice, z = north. |
