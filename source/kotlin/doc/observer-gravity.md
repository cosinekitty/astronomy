//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[observerGravity](observer-gravity.md)

# observerGravity

fun [observerGravity](observer-gravity.md)(latitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), height: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)

Calculates the gravitational acceleration experienced by an observer on the Earth.

This function implements the WGS 84 Ellipsoidal Gravity Formula. The result is a combination of inward gravitational acceleration with outward centrifugal acceleration, as experienced by an observer in the Earth's rotating frame of reference. The resulting value increases toward the Earth's poles and decreases toward the equator, consistent with changes of the weight measured by a spring scale of a fixed mass moved to different latitudes and heights on the Earth.

#### Return

The effective gravitational acceleration expressed in meters per second squared.

## Parameters

jvm

| | |
|---|---|
| latitude | The latitude of the observer in degrees north or south of the equator.     By formula symmetry, positive latitudes give the same answer as negative     latitudes, so the sign does not matter. |
| height | The height above the sea level geoid in meters.     No range checking is done; however, accuracy is only valid in the     range 0 to 100000 meters. |
