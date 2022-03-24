//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Equatorial](index.md)

# Equatorial

[jvm]\
class [Equatorial](index.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vec: [AstroVector](../-astro-vector/index.md))

Equatorial angular and cartesian coordinates.

Coordinates of a celestial body as seen from the Earth (geocentric or topocentric, depending on context), oriented with respect to the projection of the Earth's equator onto the sky.

## Constructors

| | |
|---|---|
| [Equatorial](-equatorial.md) | [jvm]<br>fun [Equatorial](-equatorial.md)(ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), vec: [AstroVector](../-astro-vector/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [dec](dec.md) | [jvm]<br>val [dec](dec.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Declination in degrees. |
| [dist](dist.md) | [jvm]<br>val [dist](dist.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Distance to the celestial body in AU. |
| [ra](ra.md) | [jvm]<br>val [ra](ra.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Right ascension in sidereal hours. |
| [vec](vec.md) | [jvm]<br>val [vec](vec.md): [AstroVector](../-astro-vector/index.md)<br>Equatorial coordinates in cartesian vector form: x = March equinox, y = June solstice, z = north. |
