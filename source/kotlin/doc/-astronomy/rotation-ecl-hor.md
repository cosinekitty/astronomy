//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[rotationEclHor](rotation-ecl-hor.md)

# rotationEclHor

[jvm]\
fun [rotationEclHor](rotation-ecl-hor.md)(time: [AstroTime](../-astro-time/index.md), observer: [Observer](../-observer/index.md)): [RotationMatrix](../-rotation-matrix/index.md)

Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECL = ecliptic system, using equator at J2000 epoch. Target: HOR = horizontal system.

## Parameters

jvm

| | |
|---|---|
| time | The date and time of the desired horizontal orientation. |
| observer | A location near the Earth's mean sea level that defines the observer's horizon. |
