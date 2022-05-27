//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[localSolarEclipsesAfter](local-solar-eclipses-after.md)

# localSolarEclipsesAfter

fun [localSolarEclipsesAfter](local-solar-eclipses-after.md)(startTime: [Time](-time/index.md), observer: [Observer](-observer/index.md)): [Sequence](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin.sequences/-sequence/index.html)&lt;[LocalSolarEclipseInfo](-local-solar-eclipse-info/index.md)&gt;

Enumerates a series of consecutive local solar eclipses that occur after a given time.

This function enables iteration through an unlimited number of consecutive local solar eclipses starting at a given time. This is a convenience wrapper around [searchLocalSolarEclipse](search-local-solar-eclipse.md) and [nextLocalSolarEclipse](next-local-solar-eclipse.md).

## Parameters

| | |
|---|---|
| startTime | The date and time for starting the search for a series of local solar eclipses. |
| observer | The geographic location of the observer. |
