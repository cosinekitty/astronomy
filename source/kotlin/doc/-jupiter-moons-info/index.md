//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[JupiterMoonsInfo](index.md)

# JupiterMoonsInfo

class [JupiterMoonsInfo](index.md)(io: [StateVector](../-state-vector/index.md), europa: [StateVector](../-state-vector/index.md), ganymede: [StateVector](../-state-vector/index.md), callisto: [StateVector](../-state-vector/index.md))

Holds the positions and velocities of Jupiter's major 4 moons.

The [jupiterMoons](../jupiter-moons.md) function returns an object of this type to report position and velocity vectors for Jupiter's largest 4 moons Io, Europa, Ganymede, and Callisto. Each position vector is relative to the center of Jupiter. Both position and velocity are oriented in the EQJ system (that is, using Earth's equator at the J2000 epoch). The positions are expressed in astronomical units (AU), and the velocities in AU/day.

## Constructors

| | |
|---|---|
| [JupiterMoonsInfo](-jupiter-moons-info.md)<br>fun [JupiterMoonsInfo](-jupiter-moons-info.md)(io: [StateVector](../-state-vector/index.md), europa: [StateVector](../-state-vector/index.md), ganymede: [StateVector](../-state-vector/index.md), callisto: [StateVector](../-state-vector/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [callisto](callisto.md)<br>val [callisto](callisto.md): [StateVector](../-state-vector/index.md)<br>The position and velocity of Jupiter's moon Callisto. |
| [europa](europa.md)<br>val [europa](europa.md): [StateVector](../-state-vector/index.md)<br>The position and velocity of Jupiter's moon Europa. |
| [ganymede](ganymede.md)<br>val [ganymede](ganymede.md): [StateVector](../-state-vector/index.md)<br>The position and velocity of Jupiter's moon Ganymede. |
| [io](io.md)<br>val [io](io.md): [StateVector](../-state-vector/index.md)<br>The position and velocity of Jupiter's moon Io. |
