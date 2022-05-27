//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[lunarApsidesAfter](lunar-apsides-after.md)

# lunarApsidesAfter

fun [lunarApsidesAfter](lunar-apsides-after.md)(startTime: [Time](-time/index.md)): [Sequence](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin.sequences/-sequence/index.html)&lt;[ApsisInfo](-apsis-info/index.md)&gt;

Enumerates a series of consecutive lunar apsides that occur after a given time.

This function enables iteration through an unlimited number of consecutive lunar perigees/apogees starting at a given time. This is a convenience wrapper around [searchLunarApsis](search-lunar-apsis.md) and [nextLunarApsis](next-lunar-apsis.md).

## Parameters

| | |
|---|---|
| startTime | The date and time for starting the search for a series of lunar apsides. |
