//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[rotationHorEqd](rotation-hor-eqd.md)

# rotationHorEqd

[jvm]\
fun [rotationHorEqd](rotation-hor-eqd.md)(time: [AstroTime](../-astro-time/index.md), observer: [Observer](../-observer/index.md)): [RotationMatrix](../-rotation-matrix/index.md)

Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).

This is one of the family of functions that returns a rotation matrix for converting from one orientation to another. Source: HOR = horizontal system (x=North, y=West, z=Zenith). Target: EQD = equatorial system, using equator of the specified date/time.

#### Return

    A rotation matrix that converts HOR to EQD at `time` and for `observer`.

## Parameters

jvm

| | |
|---|---|
| time | The date and time at which the Earth's equator applies. |
| observer | A location near the Earth's mean sea level that defines the observer's horizon. |
