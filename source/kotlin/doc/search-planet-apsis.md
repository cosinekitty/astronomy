//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchPlanetApsis](search-planet-apsis.md)

# searchPlanetApsis

fun [searchPlanetApsis](search-planet-apsis.md)(body: [Body](-body/index.md), startTime: [Time](-time/index.md)): [ApsisInfo](-apsis-info/index.md)

Finds the first aphelion or perihelion for a planet after a given time.

Given a date and time to start the search in startTime, this function finds the next date and time that the center of the specified planet reaches the closest or farthest point in its orbit with respect to the center of the Sun, whichever comes first, after startTime.

The closest point is called *perihelion* and the farthest point is called *aphelion*. The word *apsis* refers to either event.

To iterate through consecutive alternating perihelion and aphelion events, call searchPlanetApsis once, then use the return value to call [nextPlanetApsis](next-planet-apsis.md). After that, keep feeding the previous return value from nextPlanetApsis into another call of nextPlanetApsis as many times as desired.

## Parameters

| | |
|---|---|
| body | The planet for which to find the next perihelion/aphelion event.     Not allowed to be [Body.Sun] or [Body.Moon]. |
| startTime | The date and time at which to start searching for the next perihelion or aphelion. |
