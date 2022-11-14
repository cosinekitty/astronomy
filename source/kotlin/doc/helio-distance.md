//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[helioDistance](helio-distance.md)

# helioDistance

fun [helioDistance](helio-distance.md)(body: [Body](-body/index.md), time: [Time](-time/index.md)): Double

Calculates the distance between a body and the Sun at a given time.

Given a date and time, this function calculates the distance between the center of body and the center of the Sun, expressed in AU. For the planets Mercury through Neptune, this function is significantly more efficient than calling [helioVector](helio-vector.md) followed by taking the length of the resulting vector.

#### Return

The heliocentric distance in AU.

## Parameters

| | |
|---|---|
| body | A body for which to calculate a heliocentric distance: the Sun, Moon, EMB, SSB, or any of the planets. |
| time | The date and time for which to calculate the distance. |
