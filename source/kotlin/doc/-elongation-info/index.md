//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[ElongationInfo](index.md)

# ElongationInfo

class [ElongationInfo](index.md)(time: [Time](../-time/index.md), visibility: [Visibility](../-visibility/index.md), elongation: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), eclipticSeparation: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

Contains information about the visibility of a celestial body at a given date and time.

See [elongation](elongation.md) for more detailed information about the members of this class. See also [searchMaxElongation](../search-max-elongation.md) for how to search for maximum elongation events.

## Constructors

| | |
|---|---|
| [ElongationInfo](-elongation-info.md)<br>fun [ElongationInfo](-elongation-info.md)(time: [Time](../-time/index.md), visibility: [Visibility](../-visibility/index.md), elongation: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), eclipticSeparation: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [eclipticSeparation](ecliptic-separation.md)<br>val [eclipticSeparation](ecliptic-separation.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth. |
| [elongation](elongation.md)<br>val [elongation](elongation.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The angle in degrees between the body and the Sun, as seen from the Earth. |
| [time](time.md)<br>val [time](time.md): [Time](../-time/index.md)<br>The date and time of the observation. |
| [visibility](visibility.md)<br>val [visibility](visibility.md): [Visibility](../-visibility/index.md)<br>Whether the body is best seen in the morning or the evening. |
