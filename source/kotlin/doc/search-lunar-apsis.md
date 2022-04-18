//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchLunarApsis](search-lunar-apsis.md)

# searchLunarApsis

[jvm]\
fun [searchLunarApsis](search-lunar-apsis.md)(startTime: [Time](-time/index.md)): [ApsisInfo](-apsis-info/index.md)

Finds the date and time of the Moon's perigee or apogee.

Given a date and time to start the search in startTime, this function finds the next date and time that the center of the Moon reaches the closest or farthest point in its orbit with respect to the center of the Earth, whichever comes first after startTime.

The closest point is called *perigee* and the farthest point is called *apogee*. The word *apsis* refers to either event.

To iterate through consecutive alternating perigee and apogee events, call searchLunarApsis once, then use the return value to call [nextLunarApsis](next-lunar-apsis.md). After that, keep feeding the previous return value from Astronomy.NextLunarApsis into another call of Astronomy.NextLunarApsis as many times as desired.

## Parameters

jvm

| | |
|---|---|
| startTime | The date and time at which to start searching for the next perigee or apogee. |
