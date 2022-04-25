//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[JupiterMoonsInfo](index.md)

# JupiterMoonsInfo

class [JupiterMoonsInfo](index.md)(moon: [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[StateVector](../-state-vector/index.md)&gt;)

Holds the positions and velocities of Jupiter's major 4 moons.

The [jupiterMoons](../jupiter-moons.md) function returns an object of this type to report position and velocity vectors for Jupiter's largest 4 moons Io, Europa, Ganymede, and Callisto. Each position vector is relative to the center of Jupiter. Both position and velocity are oriented in the EQJ system (that is, using Earth's equator at the J2000 epoch). The positions are expressed in astronomical units (AU), and the velocities in AU/day.

## Constructors

| | |
|---|---|
| [JupiterMoonsInfo](-jupiter-moons-info.md)<br>fun [JupiterMoonsInfo](-jupiter-moons-info.md)(moon: [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[StateVector](../-state-vector/index.md)&gt;) |

## Properties

| Name | Summary |
|---|---|
| [moon](moon.md)<br>val [moon](moon.md): [Array](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-array/index.html)&lt;[StateVector](../-state-vector/index.md)&gt;<br>An array of state vectors for each of the 4 moons, in the following order: 0 = Io, 1 = Europa, 2 = Ganymede, 3 = Callisto. |
