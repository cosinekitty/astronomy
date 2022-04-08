//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[rotationEclEqd](rotation-ecl-eqd.md)

# rotationEclEqd

[jvm]\
fun [rotationEclEqd](rotation-ecl-eqd.md)(time: [AstroTime](-astro-time/index.md)): [RotationMatrix](-rotation-matrix/index.md)

Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECL = ecliptic system, using equator at J2000 epoch. Target: EQD = equatorial system, using equator of date.

#### Return

A rotation matrix that converts ECL to EQD.

## Parameters

jvm

| | |
|---|---|
| time | The date and time of the desired equator. |
