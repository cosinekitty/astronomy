//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[illumination](illumination.md)

# illumination

fun [illumination](illumination.md)(body: [Body](-body/index.md), time: [Time](-time/index.md)): [IlluminationInfo](-illumination-info/index.md)

Finds visual magnitude, phase angle, and other illumination information about a celestial body.

This function calculates information about how bright a celestial body appears from the Earth, reported as visual magnitude, which is a smaller (or even negative) number for brighter objects and a larger number for dimmer objects.

For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction of the body appears illuminated as seen from the Earth. For example, when the phase angle is near zero, it means the body appears "full" as seen from the Earth.  A phase angle approaching 180 degrees means the body appears as a thin crescent as seen from the Earth.  A phase angle of 90 degrees means the body appears "half full". For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it, so it doesn't have a phase angle.

When the body is Saturn, the returned structure contains a field ringTilt that holds the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means the rings appear edge-on, and are thus nearly invisible from the Earth. The ringTilt holds 0 for all bodies other than Saturn.

## Parameters

| | |
|---|---|
| body | The Sun, Moon, or any planet other than the Earth. |
| time | The date and time of the observation. |
