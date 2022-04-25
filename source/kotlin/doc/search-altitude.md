//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchAltitude](search-altitude.md)

# searchAltitude

fun [searchAltitude](search-altitude.md)(body: [Body](-body/index.md), observer: [Observer](-observer/index.md), direction: [Direction](-direction/index.md), startTime: [Time](-time/index.md), limitDays: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Time](-time/index.md)?

Finds the next time a body reaches a given altitude.

Finds when the given body ascends or descends through a given altitude angle, as seen by an observer at the specified location on the Earth. By using the appropriate combination of direction and altitude parameters, this function can be used to find when civil, nautical, or astronomical twilight begins (dawn) or ends (dusk).

Civil dawn begins before sunrise when the Sun ascends through 6 degrees below the horizon. To find civil dawn, pass Direction.Rise for direction and -6 for altitude.

Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon. To find civil dusk, pass Direction.Set for direction and -6 for altitude.

Nautical twilight is similar to civil twilight, only the altitude value should be -12 degrees.

Astronomical twilight uses -18 degrees as the altitude value.

#### Return

The date and time of the altitude event, or null if no such event occurs within the specified time window.

## Parameters

jvm

| | |
|---|---|
| body | The Sun, Moon, or any planet other than the Earth. |
| observer | The location where observation takes place. |
| direction | Either Direction.Rise to find an ascending altitude event or Direction.Set to find a descending altitude event. |
| startTime | The date and time at which to start the search. |
| limitDays | The fractional number of days after dateStart that limits when the altitude event is to be found. Must be a positive number. |
| altitude | The desired altitude angle of the body's center above (positive) or below (negative) the observer's local horizon, expressed in degrees. Must be in the range -90, +90. |
