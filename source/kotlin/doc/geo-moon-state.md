//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[geoMoonState](geo-moon-state.md)

# geoMoonState

[jvm]\
fun [geoMoonState](geo-moon-state.md)(time: [AstroTime](-astro-time/index.md)): [StateVector](-state-vector/index.md)

Calculates equatorial geocentric position and velocity of the Moon at a given time.

Given a time of observation, calculates the Moon's position and velocity vectors. The position and velocity are of the Moon's center relative to the Earth's center. The position (x, y, z) components are expressed in AU (astronomical units). The velocity (vx, vy, vz) components are expressed in AU/day. The coordinates are oriented with respect to the Earth's equator at the J2000 epoch. In Astronomy Engine, this orientation is called EQJ. If you need the Moon's position only, and not its velocity, it is much more efficient to use [geoMoon](geo-moon.md) instead.

#### Return

The Moon's position and velocity vectors in J2000 equatorial coordinates (EQJ).

## Parameters

jvm

| | |
|---|---|
| time | The date and time for which to calculate the Moon's position and velocity. |
