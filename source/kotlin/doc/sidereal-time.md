//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[siderealTime](sidereal-time.md)

# siderealTime

fun [siderealTime](sidereal-time.md)(time: [Time](-time/index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)

Calculates Greenwich Apparent Sidereal Time (GAST).

Given a date and time, this function calculates the rotation of the Earth, represented by the equatorial angle of the Greenwich prime meridian with respect to distant stars (not the Sun, which moves relative to background stars by almost one degree per day). This angle is called Greenwich Apparent Sidereal Time (GAST). GAST is measured in sidereal hours in the half-open range [0, 24). When GAST = 0, it means the prime meridian is aligned with the of-date equinox, corrected at that time for precession and nutation of the Earth's axis. In this context, the *equinox* is the direction in space where the Earth's orbital plane (the ecliptic) intersects with the plane of the Earth's equator, at the location on the Earth's orbit of the (seasonal) March equinox. As the Earth rotates, GAST increases from 0 up to 24 sidereal hours, then starts over at 0. To convert to degrees, multiply the return value by 15.

## Parameters

| | |
|---|---|
| time | The date and time for which to find GAST. As an optimization, this function caches the sideral time value in time, unless it has already been cached, in which case the cached value is reused. |
