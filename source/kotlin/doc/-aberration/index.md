//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Aberration](index.md)

# Aberration

[jvm]\
enum [Aberration](index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Aberration](index.md)&gt; 

Aberration calculation options.

[Aberration](https://en.wikipedia.org/wiki/Aberration_of_light) is an effect causing the apparent direction of an observed body to be shifted due to transverse movement of the Earth with respect to the rays of light coming from that body. This angular correction can be anywhere from 0 to about 20 arcseconds, depending on the position of the observed body relative to the instantaneous velocity vector of the Earth.

Some Astronomy Engine functions allow optional correction for aberration by passing in a value of this enumerated type.

Aberration correction is useful to improve accuracy of coordinates of apparent locations of bodies seen from the Earth. However, because aberration affects not only the observed body (such as a planet) but the surrounding stars, aberration may be unhelpful (for example) for determining exactly when a planet crosses from one constellation to another.

## Entries

| | |
|---|---|
| [None](-none/index.md) | [jvm]<br>[None](-none/index.md)()<br>Do not correct for aberration. |
| [Corrected](-corrected/index.md) | [jvm]<br>[Corrected](-corrected/index.md)()<br>Request correction for aberration. |

## Properties

| Name | Summary |
|---|---|
| [name](../-node-event-kind/-ascending/index.md#-372974862%2FProperties%2F-1216412040) | [jvm]<br>val [name](../-node-event-kind/-ascending/index.md#-372974862%2FProperties%2F-1216412040): [String](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-string/index.html) |
| [ordinal](../-node-event-kind/-ascending/index.md#-739389684%2FProperties%2F-1216412040) | [jvm]<br>val [ordinal](../-node-event-kind/-ascending/index.md#-739389684%2FProperties%2F-1216412040): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html) |
