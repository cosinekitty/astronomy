//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[rotationEqjHor](rotation-eqj-hor.md)

# rotationEqjHor

[jvm]\
fun [rotationEqjHor](rotation-eqj-hor.md)(time: [AstroTime](../-astro-time/index.md), observer: [Observer](../-observer/index.md)): [RotationMatrix](../-rotation-matrix/index.md)

Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using the equator at the J2000 epoch. Target: HOR = horizontal system.

## Parameters

jvm

| | |
|---|---|
| time | The date and time of the observation. |
| observer | A location near the Earth's mean sea level that defines the observer's horizon. |
