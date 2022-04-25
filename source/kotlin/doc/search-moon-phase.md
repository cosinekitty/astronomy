//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchMoonPhase](search-moon-phase.md)

# searchMoonPhase

fun [searchMoonPhase](search-moon-phase.md)(targetLon: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), startTime: [Time](-time/index.md), limitDays: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Time](-time/index.md)?

Searches for the time that the Moon reaches a specified phase.

Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic longitude with respect to the Sun's geocentric ecliptic longitude. When the Moon and the Sun have the same longitude, that is defined as a new moon. When their longitudes are 180 degrees apart, that is defined as a full moon.

This function searches for any value of the lunar phase expressed as an angle in degrees in the range [0, 360).

If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter) it is much easier to call the functions [searchMoonQuarter](search-moon-quarter.md) and [nextMoonQuarter](next-moon-quarter.md). This function is useful for finding general phase angles outside those four quarters.

#### Return

If successful, returns the date and time the moon reaches the phase specified by targetlon. This function will return null if the phase does not occur within limitDays of startTime; that is, if the search window is too small.

## Parameters

jvm

| | |
|---|---|
| targetLon | The difference in geocentric longitude between the Sun and Moon     that specifies the lunar phase being sought. This can be any value     in the range [0, 360).  Certain values have conventional names:     0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter. |
| startTime | The beginning of the time window in which to search for the Moon reaching the specified phase. |
| limitDays | The number of days after `startTime` that limits the time window for the search. |
