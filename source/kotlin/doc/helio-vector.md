//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[helioVector](helio-vector.md)

# helioVector

[jvm]\
fun [helioVector](helio-vector.md)(body: [Body](-body/index.md), time: [AstroTime](-astro-time/index.md)): [AstroVector](-astro-vector/index.md)

Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.

This function calculates the position of the given celestial body as a vector, using the center of the Sun as the origin.  The result is expressed as a Cartesian vector in the J2000 equatorial system: the coordinates are based on the mean equator of the Earth at noon UTC on 1 January 2000.

The position is not corrected for light travel time or aberration. This is different from the behavior of [geoVector](geo-vector.md).

If given an invalid value for body, this function will throw an [InvalidBodyException](-invalid-body-exception/index.md).

#### Return

The heliocentric position vector of the center of the given body.

## Parameters

jvm

| | |
|---|---|
| body | A body for which to calculate a heliocentric position:     the Sun, Moon, EMB, SSB, or any of the planets. |
| time | The date and time for which to calculate the position. |
