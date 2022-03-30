//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[ApsisInfo](index.md)

# ApsisInfo

[jvm]\
class [ApsisInfo](index.md)(time: [AstroTime](../-astro-time/index.md), kind: [ApsisKind](../-apsis-kind/index.md), distAu: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), distKm: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

An apsis event: pericenter (closest approach) or apocenter (farthest distance).

For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an event where the orbiting body reaches its closest or farthest point from the primary body. The closest approach is called *pericenter* and the farthest point is *apocenter*.

More specific terminology is common for particular orbiting bodies. The Moon's closest approach to the Earth is called *perigee* and its farthest point is called *apogee*. The closest approach of a planet to the Sun is called *perihelion* and the furthest point is called *aphelion*.

This data structure is returned by Astronomy.searchLunarApsis and Astronomy.nextLunarApsis to iterate through consecutive alternating perigees and apogees.

## Constructors

| | |
|---|---|
| [ApsisInfo](-apsis-info.md) | [jvm]<br>fun [ApsisInfo](-apsis-info.md)(time: [AstroTime](../-astro-time/index.md), kind: [ApsisKind](../-apsis-kind/index.md), distAu: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), distKm: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [distAu](dist-au.md) | [jvm]<br>val [distAu](dist-au.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The distance between the centers of the bodies in astronomical units. |
| [distKm](dist-km.md) | [jvm]<br>val [distKm](dist-km.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The distance between the centers of the bodies in kilometers. |
| [kind](kind.md) | [jvm]<br>val [kind](kind.md): [ApsisKind](../-apsis-kind/index.md)<br>Whether this is a pericenter or apocenter event. |
| [time](time.md) | [jvm]<br>val [time](time.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the apsis. |
