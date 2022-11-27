//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchAltitude](search-altitude.md)

# searchAltitude

fun [searchAltitude](search-altitude.md)(body: [Body](-body/index.md), observer: [Observer](-observer/index.md), direction: [Direction](-direction/index.md), startTime: [Time](-time/index.md), limitDays: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Time](-time/index.md)?

Finds the next time the center of a body reaches a given altitude.

Finds when the center of the given body ascends or descends through a given altitude angle, as seen by an observer at the specified location on the Earth. By using the appropriate combination of direction and altitude parameters, this function can be used to find when civil, nautical, or astronomical twilight begins (dawn) or ends (dusk).

Civil dawn begins before sunrise when the Sun ascends through 6 degrees below the horizon. To find civil dawn, pass Direction.Rise for direction and -6 for altitude.

Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon. To find civil dusk, pass Direction.Set for direction and -6 for altitude.

Nautical twilight is similar to civil twilight, only the altitude value should be -12 degrees.

Astronomical twilight uses -18 degrees as the altitude value.

By convention for twilight time calculations, the altitude is not corrected for atmospheric refraction. This is because the target altitudes are below the horizon, and refraction is not directly observable.

searchAltitude is not intended to find rise/set times of a body for two reasons: (1) Rise/set times of the Sun or Moon are defined by their topmost visible portion, not their centers. (2) Rise/set times are affected significantly by atmospheric refraction. Therefore, it is better to use [searchRiseSet](search-rise-set.md) to find rise/set times, which corrects for both of these considerations.

searchAltitude will not work reliably for altitudes at or near the body's maximum or minimum altitudes. To find the time a body reaches minimum or maximum altitude angles, use [searchHourAngle](search-hour-angle.md).

#### Return

The date and time of the altitude event, or null if no such event occurs within the specified time window.

## Parameters

| | |
|---|---|
| body | The Sun, Moon, any planet other than the Earth, or a user-defined star that was created by a call to [defineStar](define-star.md). |
| observer | The location where observation takes place. |
| direction | Either [Direction.Rise](-direction/-rise/index.md) to find an ascending altitude event or [Direction.Set](-direction/-set/index.md) to find a descending altitude event. |
| startTime | The date and time at which to start the search. |
| limitDays | Limits how many days to search for the body reaching the altitude angle, and defines the direction in time to search. When limitDays is positive, the search is performed into the future, after startTime. When negative, the search is performed into the past, before startTime. To limit the search to the same day, you can use a value of 1 day. In cases where you want to find the altitude event no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. |
| altitude | The desired altitude angle of the body's center above (positive) or below (negative) the observer's local horizon, expressed in degrees. Must be in the range -90, +90. |
