//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[LibrationInfo](index.md)

# LibrationInfo

[jvm]\
data class [LibrationInfo](index.md)(elat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), elon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), mlat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), mlon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), distKm: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), diamDeg: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

Lunar libration angles, returned by Astronomy.libration.

## Constructors

| | |
|---|---|
| [LibrationInfo](-libration-info.md) | [jvm]<br>fun [LibrationInfo](-libration-info.md)(elat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), elon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), mlat: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), mlon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), distKm: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), diamDeg: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [diamDeg](diam-deg.md) | [jvm]<br>val [diamDeg](diam-deg.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth. |
| [distKm](dist-km.md) | [jvm]<br>val [distKm](dist-km.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Distance between the centers of the Earth and Moon in kilometers. |
| [elat](elat.md) | [jvm]<br>val [elat](elat.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Sub-Earth libration ecliptic latitude angle, in degrees. |
| [elon](elon.md) | [jvm]<br>val [elon](elon.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Sub-Earth libration ecliptic longitude angle, in degrees. |
| [mlat](mlat.md) | [jvm]<br>val [mlat](mlat.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Moon's geocentric ecliptic latitude. |
| [mlon](mlon.md) | [jvm]<br>val [mlon](mlon.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Moon's geocentric ecliptic longitude. |
