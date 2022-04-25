//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchTransit](search-transit.md)

# searchTransit

fun [searchTransit](search-transit.md)(body: [Body](-body/index.md), startTime: [Time](-time/index.md)): [TransitInfo](-transit-info/index.md)

Searches for the first transit of Mercury or Venus after a given date.

Finds the first transit of Mercury or Venus after a specified date. A transit is when an inferior planet passes between the Sun and the Earth so that the silhouette of the planet is visible against the Sun in the background. To continue the search, pass the finish time in the returned object to [nextTransit](next-transit.md).

## Parameters

jvm

| | |
|---|---|
| body | The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`. |
| startTime | The date and time for starting the search for a transit. |
