//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchGlobalSolarEclipse](search-global-solar-eclipse.md)

# searchGlobalSolarEclipse

fun [searchGlobalSolarEclipse](search-global-solar-eclipse.md)(startTime: [Time](-time/index.md)): [GlobalSolarEclipseInfo](-global-solar-eclipse-info/index.md)

Searches for a solar eclipse visible anywhere on the Earth's surface.

This function finds the first solar eclipse that occurs after startTime. A solar eclipse may be partial, annular, or total. See [GlobalSolarEclipseInfo](-global-solar-eclipse-info/index.md) for more information. To find a series of solar eclipses, call this function once, then keep calling [nextGlobalSolarEclipse](next-global-solar-eclipse.md) as many times as desired, passing in the peak value returned from the previous call.

See [globalSolarEclipsesAfter](global-solar-eclipses-after.md) for convenient iteration of consecutive eclipses.

#### Return

Information about the first solar eclipse after startTime.

## Parameters

| | |
|---|---|
| startTime | The date and time for starting the search for a solar eclipse. |
