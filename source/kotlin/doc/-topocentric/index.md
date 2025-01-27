//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Topocentric](index.md)

# Topocentric

data class [Topocentric](index.md)(azimuth: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html))

Coordinates of a celestial body as seen by a topocentric observer.

Horizontal and equatorial coordinates seen by an observer on or near the surface of the Earth (a topocentric observer). Optionally corrected for atmospheric refraction.

## Constructors

| | |
|---|---|
| [Topocentric](-topocentric.md)<br>fun [Topocentric](-topocentric.md)(azimuth: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [altitude](altitude.md)<br>val [altitude](altitude.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Angle in degrees above (positive) or below (negative) the observer's horizon. |
| [azimuth](azimuth.md)<br>val [azimuth](azimuth.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West. |
| [dec](dec.md)<br>val [dec](dec.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Declination in degrees. |
| [ra](ra.md)<br>val [ra](ra.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Right ascension in sidereal hours. |
