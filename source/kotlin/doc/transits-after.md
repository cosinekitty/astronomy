//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[transitsAfter](transits-after.md)

# transitsAfter

fun [transitsAfter](transits-after.md)(body: [Body](-body/index.md), startTime: [Time](-time/index.md)): [Sequence](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin.sequences/-sequence/index.html)&lt;[TransitInfo](-transit-info/index.md)&gt;

Enumerates a series of consecutive transits of Mercury or Venus.

This function enables iteration through a series of consecutive transits of Mercury or Venus that occur after a specified time.

This is a convenience wrapper around [searchTransit](search-transit.md) and [nextTransit](next-transit.md).

## Parameters

| | |
|---|---|
| body | The planet for which to enumerate transits. Must be [Body.Mercury](-body/-mercury/index.md) or [Body.Venus](-body/-venus/index.md). |
| startTime | The date and time for starting the search for a series of transits. |
