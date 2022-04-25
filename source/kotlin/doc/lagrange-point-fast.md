//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[lagrangePointFast](lagrange-point-fast.md)

# lagrangePointFast

fun [lagrangePointFast](lagrange-point-fast.md)(point: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), majorState: [StateVector](-state-vector/index.md), majorMass: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), minorState: [StateVector](-state-vector/index.md), minorMass: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [StateVector](-state-vector/index.md)

Calculates one of the 5 Lagrange points from body masses and state vectors.

Given a more massive "major" body and a much less massive "minor" body, calculates one of the five Lagrange points in relation to the minor body's orbit around the major body. The parameter point is an integer that selects the Lagrange point as follows:

1 = the Lagrange point between the major body and minor body. 2 = the Lagrange point on the far side of the minor body. 3 = the Lagrange point on the far side of the major body. 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position. 5 = the Lagrange point 60 degrees behind the minor body's orbital position.

The caller passes in the state vector and mass for both bodies. The state vectors can be in any orientation and frame of reference. The body masses are expressed as GM products, where G = the universal gravitation constant and M = the body's mass. Thus the units for major_mass and minor_mass must be au^3/day^2. Use [massProduct](mass-product.md) to obtain GM values for various solar system bodies.

The function returns the state vector for the selected Lagrange point using the same orientation as the state vector parameters majorState and minorState, and the position and velocity components are with respect to the major body's center.

Consider calling [lagrangePoint](lagrange-point.md), instead of this function, for simpler usage in most cases.

#### Return

The position and velocity of the selected Lagrange point with respect to the major body's center.

## Parameters

jvm

| | |
|---|---|
| point | An integer 1..5 that selects which of the Lagrange points to calculate. |
| majorState | The state vector of the major (more massive) of the pair of bodies. |
| majorMass | The mass product GM of the major body. |
| minorState | The state vector of the minor (less massive) of the pair of bodies. |
| minorMass | The mass product GM of the minor body. |
