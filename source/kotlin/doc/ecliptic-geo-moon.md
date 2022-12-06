//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[eclipticGeoMoon](ecliptic-geo-moon.md)

# eclipticGeoMoon

fun [eclipticGeoMoon](ecliptic-geo-moon.md)(time: [Time](-time/index.md)): [Spherical](-spherical/index.md)

Calculates spherical ecliptic geocentric position of the Moon.

Given a time of observation, calculates the Moon's geocentric position in ecliptic spherical coordinates. Provides the ecliptic latitude and longitude in degrees, and the geocentric distance in astronomical units (AU).

The ecliptic angles are measured in "ECT": relative to the true ecliptic plane and equatorial plane at the specified time. This means the Earth's equator is corrected for precession and nutation, and the plane of the Earth's orbit is corrected for gradual obliquity drift.

This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954, which in turn derives from E. W. Brown's lunar theories from the early twentieth century. It is adapted from Turbo Pascal code from the book [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210) by Montenbruck and Pfleger.

To calculate an equatorial J2000 vector instead, use [geoMoon](geo-moon.md).

#### Return

The Moon's distance, ecliptic latitude, and ecliptic longitude, expressed in true equinox of date.

## Parameters

| | |
|---|---|
| time | The date and time for which to calculate the Moon's position. |
