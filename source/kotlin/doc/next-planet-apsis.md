//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[nextPlanetApsis](next-planet-apsis.md)

# nextPlanetApsis

[jvm]\
fun [nextPlanetApsis](next-planet-apsis.md)(body: [Body](-body/index.md), apsis: [ApsisInfo](-apsis-info/index.md)): [ApsisInfo](-apsis-info/index.md)

Finds the next planetary perihelion or aphelion event in a series.

This function requires an [ApsisInfo](-apsis-info/index.md) value obtained from a call to [searchPlanetApsis](search-planet-apsis.md) or nextPlanetApsis. Given an aphelion event, this function finds the next perihelion event, and vice versa. See [searchPlanetApsis](search-planet-apsis.md) for more details.

## Parameters

jvm

| | |
|---|---|
| body | The planet for which to find the next perihelion/aphelion event.     Not allowed to be `Body.Sun` or `Body.Moon`.     Must match the body passed into the call that produced the `apsis` parameter. |
| apsis | An apsis event obtained from a call to [searchPlanetApsis] or `nextPlanetApsis`. |
