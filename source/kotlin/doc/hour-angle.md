//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[hourAngle](hour-angle.md)

# hourAngle

fun [hourAngle](hour-angle.md)(body: [Body](-body/index.md), time: [Time](-time/index.md), observer: [Observer](-observer/index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)

Finds the hour angle of a body for a given observer and time.

The *hour angle* of a celestial body indicates its position in the sky with respect to the Earth's rotation. The hour angle depends on the location of the observer on the Earth. The hour angle is 0 when the body's center reaches its highest angle above the horizon in a given day. The hour angle increases by 1 unit for every sidereal hour that passes after that point, up to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates the number of hours that have passed since the most recent time that the body has culminated, or reached its highest point.

#### Return

The real-valued hour angle of the body in the half-open range [0, 24).

## Parameters

| | |
|---|---|
| body | The body whose observed hour angle is to be found. |
| time | The time of the observation. |
| observer | The geographic location where the observation takes place. |
