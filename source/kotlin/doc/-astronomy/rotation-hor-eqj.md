//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[rotationHorEqj](rotation-hor-eqj.md)

# rotationHorEqj

[jvm]\
fun [rotationHorEqj](rotation-hor-eqj.md)(time: [AstroTime](../-astro-time/index.md), observer: [Observer](../-observer/index.md)): [RotationMatrix](../-rotation-matrix/index.md)

Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ). This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: HOR = horizontal system (x=North, y=West, z=Zenith). Target: EQJ = equatorial system, using equator at the J2000 epoch.

## Parameters

jvm

| | |
|---|---|
| time | The date and time of the observation. |
| observer | A location near the Earth's mean sea level that defines the observer's horizon. |
