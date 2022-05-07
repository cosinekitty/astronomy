//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[moonQuartersAfter](moon-quarters-after.md)

# moonQuartersAfter

fun [moonQuartersAfter](moon-quarters-after.md)(startTime: [Time](-time/index.md)): [Iterator](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin.collections/-iterator/index.html)&lt;[MoonQuarterInfo](-moon-quarter-info/index.md)&gt;

Enumerates a series of consecutive moon quarter phase events.

This function enables iteration through an unlimited number of consecutive lunar quarter phases starting at a given time. This is a convenience wrapper around [searchMoonQuarter](search-moon-quarter.md) and [nextMoonQuarter](next-moon-quarter.md).

## Parameters

| | |
|---|---|
| startTime | The date and time for starting the search for a series of quarter phases. |
