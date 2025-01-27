//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[globalSolarEclipsesAfter](global-solar-eclipses-after.md)

# globalSolarEclipsesAfter

fun [globalSolarEclipsesAfter](global-solar-eclipses-after.md)(startTime: [Time](-time/index.md)): [Sequence](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin.sequences/-sequence/index.html)&lt;[GlobalSolarEclipseInfo](-global-solar-eclipse-info/index.md)&gt;

Enumerates a series of consecutive global solar eclipses that occur after a given time.

This function enables iteration through an unlimited number of consecutive global solar eclipses starting at a given time. This is a convenience wrapper around [searchGlobalSolarEclipse](search-global-solar-eclipse.md) and [nextGlobalSolarEclipse](next-global-solar-eclipse.md).

## Parameters

| | |
|---|---|
| startTime | The date and time for starting the search for a series of global solar eclipses. |
