//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[lagrangePoint](lagrange-point.md)

# lagrangePoint

fun [lagrangePoint](lagrange-point.md)(point: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), time: [Time](-time/index.md), majorBody: [Body](-body/index.md), minorBody: [Body](-body/index.md)): [StateVector](-state-vector/index.md)

Calculates one of the 5 Lagrange points for a pair of co-orbiting bodies.

Given a more massive "major" body and a much less massive "minor" body, calculates one of the five Lagrange points in relation to the minor body's orbit around the major body. The parameter point is an integer that selects the Lagrange point as follows:

1 = the Lagrange point between the major body and minor body. 2 = the Lagrange point on the far side of the minor body. 3 = the Lagrange point on the far side of the major body. 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position. 5 = the Lagrange point 60 degrees behind the minor body's orbital position.

The function returns the state vector for the selected Lagrange point in equatorial J2000 coordinates (EQJ), with respect to the center of the major body.

To calculate Sun/Earth Lagrange points, pass in [Body.Sun](-body/-sun/index.md) for majorBody and [Body.EMB](-body/-e-m-b/index.md) (Earth/Moon barycenter) for minorBody. For Lagrange points of the Sun and any other planet, pass in just that planet (e.g. [Body.Jupiter](-body/-jupiter/index.md)) for minorBody. To calculate Earth/Moon Lagrange points, pass in [Body.Earth](-body/-earth/index.md) and [Body.Moon](-body/-moon/index.md) for the major and minor bodies respectively.

In some cases, it may be more efficient to call [lagrangePointFast](lagrange-point-fast.md), especially when the state vectors have already been calculated, or are needed for some other purpose.

#### Return

The position and velocity of the selected Lagrange point with respect to the major body's center.

## Parameters

| | |
|---|---|
| point | An integer 1..5 that selects which of the Lagrange points to calculate. |
| time | The time for which the Lagrange point is to be calculated. |
| majorBody | The more massive of the co-orbiting bodies: [Body.Sun](-body/-sun/index.md) or [Body.Earth](-body/-earth/index.md). |
| minorBody | The less massive of the co-orbiting bodies. See main remarks. |
