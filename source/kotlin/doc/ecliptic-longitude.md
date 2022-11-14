//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[eclipticLongitude](ecliptic-longitude.md)

# eclipticLongitude

fun [eclipticLongitude](ecliptic-longitude.md)(body: [Body](-body/index.md), time: [Time](-time/index.md)): Double

Calculates heliocentric ecliptic longitude of a body based on the J2000 equinox.

This function calculates the angle around the plane of the Earth's orbit of a celestial body, as seen from the center of the Sun. The angle is measured prograde (in the direction of the Earth's orbit around the Sun) in degrees from the J2000 equinox. The ecliptic longitude is always in the range [0, 360).

#### Return

The ecliptic longitude in degrees of the given body at the given time.

## Parameters

| | |
|---|---|
| body | A body other than the Sun. |
| time | The date and time at which the body's ecliptic longitude is to be calculated. |
