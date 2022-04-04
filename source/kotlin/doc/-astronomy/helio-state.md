//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[helioState](helio-state.md)

# helioState

[jvm]\
fun [helioState](helio-state.md)(body: [Body](../-body/index.md), time: [AstroTime](../-astro-time/index.md)): [StateVector](../-state-vector/index.md)

Calculates heliocentric position and velocity vectors for the given body.

Given a body and a time, calculates the position and velocity vectors for the center of that body at that time, relative to the center of the Sun. The vectors are expressed in equatorial J2000 coordinates (EQJ). If you need the position vector only, it is more efficient to call [Astronomy.helioVector](helio-vector.md). The Sun's center is a non-inertial frame of reference. In other words, the Sun experiences acceleration due to gravitational forces, mostly from the larger planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum, kinetic energy, or other quantities that require a non-accelerating frame of reference, consider using Astronomy.baryState instead.

#### Return

    A state vector that contains heliocentric position and velocity vectors.
    The positions are expressed in AU.
    The velocities are expressed in AU/day.

## Parameters

jvm

| | |
|---|---|
| body | The celestial body whose heliocentric state vector is to be calculated. Supported values are Body.Sun, Body.Moon, Body.EMB, Body.SSB, and all planets: Body.Mercury, Body.Venus, Body.Earth, Body.Mars, Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto. |
| time | The date and time for which to calculate position and velocity. |
