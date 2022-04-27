//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[geoVector](geo-vector.md)

# geoVector

fun [geoVector](geo-vector.md)(body: [Body](-body/index.md), time: [Time](-time/index.md), aberration: [Aberration](-aberration/index.md)): [Vector](-vector/index.md)

Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.

This function calculates the position of the given celestial body as a vector, using the center of the Earth as the origin.  The result is expressed as a Cartesian vector in the J2000 equatorial system: the coordinates are based on the mean equator of the Earth at noon UTC on 1 January 2000.

If given an invalid value for body, this function will throw an exception.

Unlike [helioVector](helio-vector.md), this function always corrects for light travel time. This means the position of the body is "back-dated" by the amount of time it takes light to travel from that body to an observer on the Earth.

Also, the position can optionally be corrected for [aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect causing the apparent direction of the body to be shifted due to transverse movement of the Earth with respect to the rays of light coming from that body.

#### Return

A geocentric position vector of the center of the given body.

## Parameters

| | |
|---|---|
| body | A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets. |
| time | The date and time for which to calculate the position. |
| aberration | [Aberration.Corrected](-aberration/-corrected/index.md) to correct for aberration, or [Aberration.None](-aberration/-none/index.md) to leave uncorrected. |
