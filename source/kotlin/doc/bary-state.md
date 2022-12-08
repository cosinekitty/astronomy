//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[baryState](bary-state.md)

# baryState

fun [baryState](bary-state.md)(body: [Body](-body/index.md), time: [Time](-time/index.md)): [StateVector](-state-vector/index.md)

Calculates barycentric position and velocity vectors for the given body.

Given a body and a time, calculates the barycentric position and velocity vectors for the center of that body at that time. The vectors are expressed in J2000 mean equator coordinates (EQJ).

#### Return

The barycentric position and velocity vectors of the body.

## Parameters

| | |
|---|---|
| body | The celestial body whose barycentric state vector is to be calculated. Supported values are [Body.Sun](-body/-sun/index.md), [Body.Moon](-body/-moon/index.md), [Body.EMB](-body/-e-m-b/index.md), [Body.SSB](-body/-s-s-b/index.md), and all planets: [Body.Mercury](-body/-mercury/index.md), [Body.Venus](-body/-venus/index.md), [Body.Earth](-body/-earth/index.md), [Body.Mars](-body/-mars/index.md), [Body.Jupiter](-body/-jupiter/index.md), [Body.Saturn](-body/-saturn/index.md), [Body.Uranus](-body/-uranus/index.md), [Body.Neptune](-body/-neptune/index.md), [Body.Pluto](-body/-pluto/index.md). |
| time | The date and time for which to calculate position and velocity. |
