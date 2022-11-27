//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchRiseSet](search-rise-set.md)

# searchRiseSet

fun [searchRiseSet](search-rise-set.md)(body: [Body](-body/index.md), observer: [Observer](-observer/index.md), direction: [Direction](-direction/index.md), startTime: [Time](-time/index.md), limitDays: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Time](-time/index.md)?

Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.

This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth. Rise time is when the body first starts to be visible above the horizon. For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon. Set time is the moment when the body appears to vanish below the horizon. Therefore, this function adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).

This function corrects for a typical value of atmospheric refraction, which causes celestial bodies to appear higher above the horizon than they would if the Earth had no atmosphere. Astronomy Engine uses a correction of 34 arcminutes. Real-world refraction varies based on air temperature, pressure, and humidity; such weather-based conditions are outside the scope of Astronomy Engine.

Note that rise or set may not occur in every 24 hour period. For example, near the Earth's poles, there are long periods of time where the Sun stays below the horizon, never rising. Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day. This is because the Moon sets nearly an hour later each day due to orbiting the Earth a significant amount during each rotation of the Earth. Therefore callers must not assume that the function will always succeed.

#### Return

On success, returns the date and time of the rise or set time as requested. If the function returns null, it means the rise or set event does not occur within limitDays days of startTime. This is a normal condition, not an error.

## Parameters

| | |
|---|---|
| body | The Sun, Moon, any planet other than the Earth, or a user-defined star that was created by a call to [defineStar](define-star.md). |
| observer | The location where observation takes place. |
| direction | Either [Direction.Rise](-direction/-rise/index.md) to find a rise time or [Direction.Set](-direction/-set/index.md) to find a set time. |
| startTime | The date and time at which to start the search. |
| limitDays | Limits how many days to search for a rise or set time, and defines the direction in time to search. When limitDays is positive, the search is performed into the future, after startTime. When negative, the search is performed into the past, before startTime. To limit a rise or set time to the same day, you can use a value of 1 day. In cases where you want to find the next rise or set time no matter how far in the future (for example, for an observer near the south pole), you can pass in a larger value like 365. |
