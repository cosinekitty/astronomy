//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[RotationMatrix](index.md)/[rotate](rotate.md)

# rotate

fun [rotate](rotate.md)(vec: [Vector](../-vector/index.md)): [Vector](../-vector/index.md)

Applies a rotation to a vector, yielding a rotated vector.

This function transforms a vector in one orientation to a vector in another orientation.

## Parameters

| | |
|---|---|
| vec | The vector whose orientation is to be changed. |

fun [rotate](rotate.md)(state: [StateVector](../-state-vector/index.md)): [StateVector](../-state-vector/index.md)

Applies a rotation to a state vector, yielding a rotated state vector.

This function transforms a state vector in one orientation to a state vector in another orientation. The resulting state vector has both position and velocity reoriented.

## Parameters

| | |
|---|---|
| state | The state vector whose orientation is to be changed.     The value of `state` is not changed; the return value is a new state vector object. |
