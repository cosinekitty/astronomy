//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[rotationEclHor](rotation-ecl-hor.md)

# rotationEclHor

fun [rotationEclHor](rotation-ecl-hor.md)(time: [Time](-time/index.md), observer: [Observer](-observer/index.md)): [RotationMatrix](-rotation-matrix/index.md)

Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECL = ecliptic system, using equator at J2000 epoch. Target: HOR = horizontal system.

#### Return

A rotation matrix that converts ECL to HOR at time and for observer. The components of the horizontal vector are: x = north, y = west, z = zenith (straight up from the observer). These components are chosen so that the "right-hand rule" works for the vector and so that north represents the direction where azimuth = 0.

## Parameters

| | |
|---|---|
| time | The date and time of the desired horizontal orientation. |
| observer | A location near the Earth's mean sea level that defines the observer's horizon. |
