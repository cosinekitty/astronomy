//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[sunPosition](sun-position.md)

# sunPosition

[jvm]\
fun [sunPosition](sun-position.md)(time: [Time](-time/index.md)): [Ecliptic](-ecliptic/index.md)

Calculates geocentric ecliptic coordinates for the Sun.

This function calculates the position of the Sun as seen from the Earth. The returned value includes both Cartesian and spherical coordinates. The x-coordinate and longitude values in the returned object are based on the *true equinox of date*: one of two points in the sky where the instantaneous plane of the Earth's equator at the given date and time (the *equatorial plane*) intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*). By convention, the apparent location of the Sun at the March equinox is chosen as the longitude origin and x-axis direction, instead of the one for September.

sunPosition corrects for precession and nutation of the Earth's axis in order to obtain the exact equatorial plane at the given time.

This function can be used for calculating changes of seasons: equinoxes and solstices. In fact, the function [seasons](seasons.md) does use this function for that purpose.

#### Return

The ecliptic coordinates of the Sun using the Earth's true equator of date.

## Parameters

jvm

| | |
|---|---|
| time | The date and time for which to calculate the Sun's position. |
