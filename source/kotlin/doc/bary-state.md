//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[baryState](bary-state.md)

# baryState

fun [baryState](bary-state.md)(body: [Body](-body/index.md), time: [Time](-time/index.md)): [StateVector](-state-vector/index.md)

Calculates barycentric position and velocity vectors for the given body.

Given a body and a time, calculates the barycentric position and velocity vectors for the center of that body at that time. The vectors are expressed in equatorial J2000 coordinates (EQJ).

#### Return

The barycentric position and velocity vectors of the body.

## Parameters

| | |
|---|---|
| body | The celestial body whose barycentric state vector is to be calculated.     Supported values are `Body.Sun`, `Body.Moon`, `Body.EMB`, `Body.SSB`, and all planets:     `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`,     `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`. |
| time | The date and time for which to calculate position and velocity. |
