//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Body](index.md)

# Body

[jvm]\
enum [Body](index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Body](index.md)&gt; 

The enumeration of celestial bodies supported by Astronomy Engine.

## Entries

| | |
|---|---|
| [SSB](-s-s-b/index.md) | [jvm]<br>[SSB](-s-s-b/index.md)(null, null, null)<br>The Solar System Barycenter. |
| [EMB](-e-m-b/index.md) | [jvm]<br>[EMB](-e-m-b/index.md)(EARTH_GM + MOON_GM, MEAN_SYNODIC_MONTH, null)<br>The Earth/Moon Barycenter. |
| [Moon](-moon/index.md) | [jvm]<br>[Moon](-moon/index.md)(MOON_GM, MEAN_SYNODIC_MONTH, null)<br>The Earth's natural satellite, the Moon. |
| [Sun](-sun/index.md) | [jvm]<br>[Sun](-sun/index.md)(SUN_GM, null, null)<br>The Sun. |
| [Pluto](-pluto/index.md) | [jvm]<br>[Pluto](-pluto/index.md)(PLUTO_GM, 90560.0, null)<br>The planet Pluto. |
| [Neptune](-neptune/index.md) | [jvm]<br>[Neptune](-neptune/index.md)(NEPTUNE_GM, NEPTUNE_ORBITAL_PERIOD, VsopModel(vsopLonNeptune, vsopLatNeptune, vsopRadNeptune))<br>The planet Neptune. |
| [Uranus](-uranus/index.md) | [jvm]<br>[Uranus](-uranus/index.md)(URANUS_GM, 30685.4, VsopModel(vsopLonUranus, vsopLatUranus, vsopRadUranus))<br>The planet Uranus. |
| [Saturn](-saturn/index.md) | [jvm]<br>[Saturn](-saturn/index.md)(SATURN_GM, 10759.22, VsopModel(vsopLonSaturn, vsopLatSaturn, vsopRadSaturn))<br>The planet Saturn. |
| [Jupiter](-jupiter/index.md) | [jvm]<br>[Jupiter](-jupiter/index.md)(JUPITER_GM, 4332.589, VsopModel(vsopLonJupiter, vsopLatJupiter, vsopRadJupiter))<br>The planet Jupiter. |
| [Mars](-mars/index.md) | [jvm]<br>[Mars](-mars/index.md)(MARS_GM, 686.98, VsopModel(vsopLonMars, vsopLatMars, vsopRadMars))<br>The planet Mars. |
| [Earth](-earth/index.md) | [jvm]<br>[Earth](-earth/index.md)(EARTH_GM, EARTH_ORBITAL_PERIOD, VsopModel(vsopLonEarth, vsopLatEarth, vsopRadEarth))<br>The planet Earth. Some functions that accept a Body parameter will fail if passed this value because they assume that an observation is being made from the Earth, and therefore the Earth is not a target of observation. |
| [Venus](-venus/index.md) | [jvm]<br>[Venus](-venus/index.md)(VENUS_GM, 224.701, VsopModel(vsopLonVenus, vsopLatVenus, vsopRadVenus))<br>The planet Venus. |
| [Mercury](-mercury/index.md) | [jvm]<br>[Mercury](-mercury/index.md)(MERCURY_GM, 87.969, VsopModel(vsopLonMercury, vsopLatMercury, vsopRadMercury))<br>The planet Mercury. |

## Properties

| Name | Summary |
|---|---|
| [name](../-node-event-kind/-invalid/index.md#-372974862%2FProperties%2F-1216412040) | [jvm]<br>val [name](../-node-event-kind/-invalid/index.md#-372974862%2FProperties%2F-1216412040): [String](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-string/index.html) |
| [ordinal](../-node-event-kind/-invalid/index.md#-739389684%2FProperties%2F-1216412040) | [jvm]<br>val [ordinal](../-node-event-kind/-invalid/index.md#-739389684%2FProperties%2F-1216412040): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html) |
