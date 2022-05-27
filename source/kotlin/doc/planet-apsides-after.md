//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[planetApsidesAfter](planet-apsides-after.md)

# planetApsidesAfter

fun [planetApsidesAfter](planet-apsides-after.md)(body: [Body](-body/index.md), startTime: [Time](-time/index.md)): [Sequence](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin.sequences/-sequence/index.html)&lt;[ApsisInfo](-apsis-info/index.md)&gt;

Enumerates a series of consecutive planetary perihelia/aphelia events.

This function enables iteration through an unlimited number of consecutive planetary apsides. This is a convenience wrapper around [searchPlanetApsis](search-planet-apsis.md) and [nextPlanetApsis](next-planet-apsis.md).

## Parameters

| | |
|---|---|
| body | The planet for which to find a series of perihelia/aphelia events. Not allowed to be [Body.Sun](-body/-sun/index.md) or [Body.Moon](-body/-moon/index.md). |
| startTime | The date and time for starting the search for a series of apsides. |
