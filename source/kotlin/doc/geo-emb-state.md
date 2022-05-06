//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[geoEmbState](geo-emb-state.md)

# geoEmbState

fun [geoEmbState](geo-emb-state.md)(time: [Time](-time/index.md)): [StateVector](-state-vector/index.md)

Calculates the geocentric position and velocity of the Earth/Moon barycenter.

Given a time of observation, calculates the geocentric position and velocity vectors of the Earth/Moon barycenter (EMB). The position (x, y, z) components are expressed in AU (astronomical units). The velocity (vx, vy, vz) components are expressed in AU/day.

#### Return

The EMB's position and velocity vectors in geocentric J2000 equatorial coordinates.

## Parameters

| | |
|---|---|
| time | The date and time for which to calculate the EMB vectors. |
