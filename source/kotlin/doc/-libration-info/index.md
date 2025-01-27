//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[LibrationInfo](index.md)

# LibrationInfo

data class [LibrationInfo](index.md)(elat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), elon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), mlat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), mlon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), distanceKm: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), diamDeg: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html))

Lunar libration angles, returned by [libration](../libration.md).

## Constructors

| | |
|---|---|
| [LibrationInfo](-libration-info.md)<br>fun [LibrationInfo](-libration-info.md)(elat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), elon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), mlat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), mlon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), distanceKm: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), diamDeg: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [diamDeg](diam-deg.md)<br>val [diamDeg](diam-deg.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth. |
| [distanceKm](distance-km.md)<br>val [distanceKm](distance-km.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Distance between the centers of the Earth and Moon in kilometers. |
| [elat](elat.md)<br>val [elat](elat.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Sub-Earth libration ecliptic latitude angle, in degrees. |
| [elon](elon.md)<br>val [elon](elon.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Sub-Earth libration ecliptic longitude angle, in degrees. |
| [mlat](mlat.md)<br>val [mlat](mlat.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Moon's geocentric ecliptic latitude. |
| [mlon](mlon.md)<br>val [mlon](mlon.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>Moon's geocentric ecliptic longitude. |
