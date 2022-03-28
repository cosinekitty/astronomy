//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[eclipticGeoMoon](ecliptic-geo-moon.md)

# eclipticGeoMoon

[jvm]\
fun [eclipticGeoMoon](ecliptic-geo-moon.md)(time: [AstroTime](../-astro-time/index.md)): [Spherical](../-spherical/index.md)

Calculates spherical ecliptic geocentric position of the Moon.

Given a time of observation, calculates the Moon's geocentric position in ecliptic spherical coordinates. Provides the ecliptic latitude and longitude in degrees, and the geocentric distance in astronomical units (AU). The ecliptic longitude is measured relative to the equinox of date.

This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954, which in turn derives from E. W. Brown's lunar theories from the early twentieth century. It is adapted from Turbo Pascal code from the book [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210) by Montenbruck and Pfleger.

To calculate an equatorial J2000 vector instead, use #Astronomy.geoMoon.
