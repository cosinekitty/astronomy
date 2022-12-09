//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[rotationEctEqj](rotation-ect-eqj.md)

# rotationEctEqj

fun [rotationEctEqj](rotation-ect-eqj.md)(time: [Time](-time/index.md)): [RotationMatrix](-rotation-matrix/index.md)

Calculates a rotation matrix from true ecliptic of date (ECT) to J2000 mean equator (EQJ).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: ECT = ecliptic system, using true equinox of the specified date/time. Target: EQJ = equatorial system, using equator at J2000 epoch.

## Parameters

| | |
|---|---|
| time | The date and time at which the Earth's equator defines the target orientation. |
