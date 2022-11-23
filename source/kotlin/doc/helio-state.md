//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[helioState](helio-state.md)

# helioState

fun [helioState](helio-state.md)(body: [Body](-body/index.md), time: [Time](-time/index.md)): [StateVector](-state-vector/index.md)

Calculates heliocentric position and velocity vectors for the given body.

Given a body and a time, calculates the position and velocity vectors for the center of that body at that time, relative to the center of the Sun. The vectors are expressed in equatorial J2000 coordinates (EQJ). If you need the position vector only, it is more efficient to call [helioVector](helio-vector.md). The Sun's center is a non-inertial frame of reference. In other words, the Sun experiences acceleration due to gravitational forces, mostly from the larger planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum, kinetic energy, or other quantities that require a non-accelerating frame of reference, consider using [baryState](bary-state.md) instead.

#### Return

A state vector that contains heliocentric position and velocity vectors. The positions are expressed in AU. The velocities are expressed in AU/day.

## Parameters

| | |
|---|---|
| body | The celestial body whose heliocentric state vector is to be calculated. Supported values are [Body.Sun](-body/-sun/index.md), [Body.Moon](-body/-moon/index.md), [Body.EMB](-body/-e-m-b/index.md), [Body.SSB](-body/-s-s-b/index.md), and all planets: [Body.Mercury](-body/-mercury/index.md), [Body.Venus](-body/-venus/index.md), [Body.Earth](-body/-earth/index.md), [Body.Mars](-body/-mars/index.md), [Body.Jupiter](-body/-jupiter/index.md), [Body.Saturn](-body/-saturn/index.md), [Body.Uranus](-body/-uranus/index.md), [Body.Neptune](-body/-neptune/index.md), [Body.Pluto](-body/-pluto/index.md). Also allowed to be a user-defined star created by [defineStar](define-star.md). |
| time | The date and time for which to calculate position and velocity. |
