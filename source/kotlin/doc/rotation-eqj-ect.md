//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[rotationEqjEct](rotation-eqj-ect.md)

# rotationEqjEct

fun [rotationEqjEct](rotation-eqj-ect.md)(time: [Time](-time/index.md)): [RotationMatrix](-rotation-matrix/index.md)

Calculates a rotation matrix from J2000 mean equator (EQJ) to true ecliptic of date (ECT).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using equator at J2000 epoch. Target: ECT = ecliptic system, using true equinox of the specified date/time.

## Parameters

| | |
|---|---|
| time | The date and time at which the Earth's equator defines the target orientation. |
