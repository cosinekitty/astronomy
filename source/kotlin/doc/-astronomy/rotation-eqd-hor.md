//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[rotationEqdHor](rotation-eqd-hor.md)

# rotationEqdHor

[jvm]\
fun [rotationEqdHor](rotation-eqd-hor.md)(time: [AstroTime](../-astro-time/index.md), observer: [Observer](../-observer/index.md)): [RotationMatrix](../-rotation-matrix/index.md)

Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: EQD = equatorial system, using equator of the specified date/time. Target: HOR = horizontal system.

## Parameters

jvm

| | |
|---|---|
| time | The date and time at which the Earth's equator applies. |
| observer | A location near the Earth's mean sea level that defines the observer's horizon. |
