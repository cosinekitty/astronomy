//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[geoMoon](geo-moon.md)

# geoMoon

[jvm]\
fun [geoMoon](geo-moon.md)(time: [AstroTime](../-astro-time/index.md)): [AstroVector](../-astro-vector/index.md)

Calculates equatorial geocentric position of the Moon at a given time.

Given a time of observation, calculates the Moon's position vector. The vector indicates the Moon's center relative to the Earth's center. The vector components are expressed in AU (astronomical units). The coordinates are oriented with respect to the Earth's equator at the J2000 epoch. In Astronomy Engine, this orientation is called EQJ.

## Parameters

jvm

| | |
|---|---|
| time | The date and time for which to calculate the Moon's position. |
