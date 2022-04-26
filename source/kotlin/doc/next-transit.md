//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[nextTransit](next-transit.md)

# nextTransit

fun [nextTransit](next-transit.md)(body: [Body](-body/index.md), prevTransitTime: [Time](-time/index.md)): [TransitInfo](-transit-info/index.md)

Searches for another transit of Mercury or Venus.

After calling [searchTransit](search-transit.md) to find a transit of Mercury or Venus, this function finds the next transit after that. Keep calling this function as many times as you want to keep finding more transits.

## Parameters

| | |
|---|---|
| body | The planet whose transit is to be found. Must be [Body.Mercury] or [Body.Venus]. |
| prevTransitTime | A date and time near the previous transit. |
