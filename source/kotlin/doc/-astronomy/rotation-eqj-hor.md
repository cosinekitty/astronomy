//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[rotationEqjHor](rotation-eqj-hor.md)

# rotationEqjHor

[jvm]\
fun [rotationEqjHor](rotation-eqj-hor.md)(time: [AstroTime](../-astro-time/index.md), observer: [Observer](../-observer/index.md)): [RotationMatrix](../-rotation-matrix/index.md)

Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using the equator at the J2000 epoch. Target: HOR = horizontal system.

#### Return

    A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
    The components of the horizontal vector are:
    x = north, y = west, z = zenith (straight up from the observer).
    These components are chosen so that the "right-hand rule" works for the vector
    and so that north represents the direction where azimuth = 0.

## Parameters

jvm

| | |
|---|---|
| time | The date and time of the observation. |
| observer | A location near the Earth's mean sea level that defines the observer's horizon. |
