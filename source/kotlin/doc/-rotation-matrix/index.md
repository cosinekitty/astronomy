//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[RotationMatrix](index.md)

# RotationMatrix

class [RotationMatrix](index.md)(rot: Array&lt;DoubleArray&gt;)

A rotation matrix that can be used to transform one coordinate system to another.

## Constructors

| | |
|---|---|
| [RotationMatrix](-rotation-matrix.md)<br>fun [RotationMatrix](-rotation-matrix.md)(a00: Double, a01: Double, a02: Double, a10: Double, a11: Double, a12: Double, a20: Double, a21: Double, a22: Double) |
| [RotationMatrix](-rotation-matrix.md)<br>fun [RotationMatrix](-rotation-matrix.md)(rot: Array&lt;DoubleArray&gt;) |

## Types

| Name | Summary |
|---|---|
| [Companion](-companion/index.md)<br>object [Companion](-companion/index.md) |

## Functions

| Name | Summary |
|---|---|
| [combine](combine.md)<br>infix fun [combine](combine.md)(other: [RotationMatrix](index.md)): [RotationMatrix](index.md)<br>Creates a rotation based on applying one rotation followed by another. |
| [inverse](inverse.md)<br>fun [inverse](inverse.md)(): [RotationMatrix](index.md)<br>Calculates the inverse of a rotation matrix. |
| [pivot](pivot.md)<br>fun [pivot](pivot.md)(axis: Int, angle: Double): [RotationMatrix](index.md)<br>Re-orients the rotation matrix by pivoting it by an angle around one of its axes. |
| [rotate](rotate.md)<br>fun [rotate](rotate.md)(state: [StateVector](../-state-vector/index.md)): [StateVector](../-state-vector/index.md)<br>Applies a rotation to a state vector, yielding a rotated state vector.<br>[jvm]<br>fun [rotate](rotate.md)(vec: [Vector](../-vector/index.md)): [Vector](../-vector/index.md)<br>Applies a rotation to a vector, yielding a rotated vector. |

## Properties

| Name | Summary |
|---|---|
| [rot](rot.md)<br>val [rot](rot.md): Array&lt;DoubleArray&gt;<br>A 3x3 array of numbers to initialize the rotation matrix. |
