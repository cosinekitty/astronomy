//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[IllumInfo](index.md)

# IllumInfo

[jvm]\
class [IllumInfo](index.md)(time: [AstroTime](../-astro-time/index.md), mag: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), phaseAngle: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), phaseFraction: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), helioDist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), ringTilt: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

Information about the brightness and illuminated shape of a celestial body.

Returned by the functions Astronomy.illumination and Astronomy.searchPeakMagnitude to report the visual magnitude and illuminated fraction of a celestial body at a given date and time.

## Constructors

| | |
|---|---|
| [IllumInfo](-illum-info.md) | [jvm]<br>fun [IllumInfo](-illum-info.md)(time: [AstroTime](../-astro-time/index.md), mag: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), phaseAngle: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), phaseFraction: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), helioDist: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), ringTilt: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [helioDist](helio-dist.md) | [jvm]<br>val [helioDist](helio-dist.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The distance between the Sun and the body at the observation time. |
| [mag](mag.md) | [jvm]<br>val [mag](mag.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The visual magnitude of the body. Smaller values are brighter. |
| [phaseAngle](phase-angle.md) | [jvm]<br>val [phaseAngle](phase-angle.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth. |
| [phaseFraction](phase-fraction.md) | [jvm]<br>val [phaseFraction](phase-fraction.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>A value in the range 0.0, 1.0 indicating what fraction of the body's apparent disc is illuminated, as seen from the Earth. |
| [ringTilt](ring-tilt.md) | [jvm]<br>val [ringTilt](ring-tilt.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>For Saturn, the tilt angle in degrees of its rings as seen from Earth. For all other bodies, 0.0. |
| [time](time.md) | [jvm]<br>val [time](time.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the observation. |
