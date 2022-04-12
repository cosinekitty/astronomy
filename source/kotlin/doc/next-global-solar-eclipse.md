//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[nextGlobalSolarEclipse](next-global-solar-eclipse.md)

# nextGlobalSolarEclipse

[jvm]\
fun [nextGlobalSolarEclipse](next-global-solar-eclipse.md)(prevEclipseTime: [Time](-time/index.md)): [GlobalSolarEclipseInfo](-global-solar-eclipse-info/index.md)

Searches for the next global solar eclipse in a series.

After using [searchGlobalSolarEclipse](search-global-solar-eclipse.md) to find the first solar eclipse in a series, you can call this function to find the next consecutive solar eclipse. Pass in the peak value from the [GlobalSolarEclipseInfo](-global-solar-eclipse-info/index.md) returned by the previous call to searchGlobalSolarEclipse or nextGlobalSolarEclipse to find the next solar eclipse.

#### Return

Information about the next consecutive solar eclipse.

## Parameters

jvm

| | |
|---|---|
| prevEclipseTime | A date and time near a new moon. Solar eclipse search will start at the next new moon. |
