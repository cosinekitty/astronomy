//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[pairLongitude](pair-longitude.md)

# pairLongitude

fun [pairLongitude](pair-longitude.md)(body1: [Body](-body/index.md), body2: [Body](-body/index.md), time: [Time](-time/index.md)): Double

Returns one body's ecliptic longitude with respect to another, as seen from the Earth.

This function determines where one body appears around the ecliptic plane (the plane of the Earth's orbit around the Sun) as seen from the Earth, relative to the another body's apparent position. The function returns an angle in the half-open range [0, 360) degrees. The value is the ecliptic longitude of body1 relative to the ecliptic longitude of body2.

The angle is 0 when the two bodies are at the same ecliptic longitude as seen from the Earth. The angle increases in the prograde direction (the direction that the planets orbit the Sun and the Moon orbits the Earth).

When the angle is 180 degrees, it means the two bodies appear on opposite sides of the sky for an Earthly observer.

Neither body1 nor body2 is allowed to be [Body.Earth](-body/-earth/index.md). If this happens, the function throws an exception.

#### Return

An angle in the range [0, 360), expressed in degrees.

## Parameters

| | |
|---|---|
| body1 | The first body, whose longitude is to be found relative to the second body. |
| body2 | The second body, relative to which the longitude of the first body is to be found. |
| time | The date and time of the observation. |
