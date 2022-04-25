//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[rotationEqjEqd](rotation-eqj-eqd.md)

# rotationEqjEqd

fun [rotationEqjEqd](rotation-eqj-eqd.md)(time: [Time](-time/index.md)): [RotationMatrix](-rotation-matrix/index.md)

Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQJ = equatorial system, using equator at J2000 epoch. Target: EQD = equatorial system, using equator of the specified date/time.

## Parameters

jvm

| | |
|---|---|
| time | The date and time at which the Earth's equator defines the target orientation. |
