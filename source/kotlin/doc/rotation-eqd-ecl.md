//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[rotationEqdEcl](rotation-eqd-ecl.md)

# rotationEqdEcl

fun [rotationEqdEcl](rotation-eqd-ecl.md)(time: [Time](-time/index.md)): [RotationMatrix](-rotation-matrix/index.md)

Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean ecliptic (ECL).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQD = equatorial system, using equator of date. Target: ECL = ecliptic system, using equator at J2000 epoch.

#### Return

A rotation matrix that converts EQD to ECL.

## Parameters

| | |
|---|---|
| time | The date and time of the source equator. |
