//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchSunLongitude](search-sun-longitude.md)

# searchSunLongitude

fun [searchSunLongitude](search-sun-longitude.md)(targetLon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), startTime: [Time](-time/index.md), limitDays: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)): [Time](-time/index.md)?

Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.

This function finds the moment in time, if any exists in the given time window, that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.

This function can be used to determine equinoxes and solstices. However, it is usually more convenient and efficient to call [seasons](seasons.md) to calculate all equinoxes and solstices for a given calendar year.

The function searches the window of time specified by startTime and startTime+limitDays. The search will return null if the Sun never reaches the longitude targetLon or if the window is so large that the longitude ranges more than 180 degrees within it. It is recommended to keep the window smaller than 10 days when possible.
