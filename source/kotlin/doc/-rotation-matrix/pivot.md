//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[RotationMatrix](index.md)/[pivot](pivot.md)

# pivot

fun [pivot](pivot.md)(axis: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-int/index.html), angle: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)): [RotationMatrix](index.md)

Re-orients the rotation matrix by pivoting it by an angle around one of its axes.

Given this rotation matrix, a selected coordinate axis, and an angle in degrees, this function pivots the rotation matrix by that angle around that coordinate axis. The function returns a new rotation matrix; it does not mutate this matrix.

For example, if you have rotation matrix that converts ecliptic coordinates (ECL) to horizontal coordinates (HOR), but you really want to convert ECL to the orientation of a telescope camera pointed at a given body, you can use RotationMatrix.pivot twice: (1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the western axis by the body's altitude angle. The resulting rotation matrix will then reorient ECL coordinates to the orientation of your telescope camera.

## Parameters

| | |
|---|---|
| axis | An integer that selects which coordinate axis to rotate around: 0 = x, 1 = y, 2 = z. Any other value will cause an exception. |
| angle | An angle in degrees indicating the amount of rotation around the specified axis. Positive angles indicate rotation counterclockwise as seen from the positive direction along that axis, looking towards the origin point of the orientation system. Any finite number of degrees is allowed, but best precision will result from keeping angle in the range -360, +360. |
