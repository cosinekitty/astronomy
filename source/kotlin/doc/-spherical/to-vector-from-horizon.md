//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Spherical](index.md)/[toVectorFromHorizon](to-vector-from-horizon.md)

# toVectorFromHorizon

[jvm]\
fun [toVectorFromHorizon](to-vector-from-horizon.md)(time: [Time](../-time/index.md), refraction: [Refraction](../-refraction/index.md)): [Vector](../-vector/index.md)

Given apparent angular horizontal coordinates, calculate the unrefracted horizontal vector.

Assumes this contains apparent horizontal coordinates: lat holds the refracted azimuth angle, lon holds the azimuth in degrees clockwise from north, and dist holds the distance from the observer to the object in AU.

#### Return

A vector in the horizontal system: x = north, y = west, and z = zenith (up).

## Parameters

jvm

| | |
|---|---|
| time | The date and time of the observation. This is needed because the returned     [Vector] requires a valid time value when passed to certain other functions. |
| refraction | The refraction option used to model atmospheric lensing. See [refractionAngle].     This specifies how refraction is to be removed from the altitude stored in `this.lat`. |
