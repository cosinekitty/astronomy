//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchHourAngle](search-hour-angle.md)

# searchHourAngle

fun [searchHourAngle](search-hour-angle.md)(body: [Body](-body/index.md), observer: [Observer](-observer/index.md), hourAngle: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), startTime: [Time](-time/index.md), direction: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html) = +1): [HourAngleInfo](-hour-angle-info/index.md)

Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.

The *hour angle* of a celestial body indicates its position in the sky with respect to the Earth's rotation. The hour angle depends on the location of the observer on the Earth. The hour angle is 0 when the body reaches its highest angle above the horizon in a given day. The hour angle increases by 1 unit for every sidereal hour that passes after that point, up to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates the number of hours that have passed since the most recent time that the body has culminated, or reached its highest point.

This function searches for the next time a celestial body reaches the given hour angle after the date and time specified by startTime. To find when a body culminates, pass 0 for hourAngle. To find when a body reaches its lowest point in the sky, pass 12 for hourAngle.

Note that, especially close to the Earth's poles, a body as seen on a given day may always be above the horizon or always below the horizon, so the caller cannot assume that a culminating object is visible nor that an object is below the horizon at its minimum altitude.

On success, the function reports the date and time, along with the horizontal coordinates of the body at that time, as seen by the given observer.

#### Return

The time when the body reaches the hour angle, and the horizontal coordinates of the body at that time.

## Parameters

| | |
|---|---|
| body | The Sun, Moon, any planet other than the Earth, or a user-defined star that was created by a call to [defineStar](define-star.md). |
| observer | A location on or near the surface of the Earth where the observer is located. |
| hourAngle | An hour angle value in the range [0, 24) indicating the number of sidereal hours after the body's most recent culmination. |
| startTime | The date and time at which to start the search. |
| direction | The direction in time to perform the search: a positive value searches forward in time, a negative value searches backward in time. The function throws an exception if direction is zero. |
