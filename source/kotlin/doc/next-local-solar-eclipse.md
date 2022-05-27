//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[nextLocalSolarEclipse](next-local-solar-eclipse.md)

# nextLocalSolarEclipse

fun [nextLocalSolarEclipse](next-local-solar-eclipse.md)(prevEclipseTime: [Time](-time/index.md), observer: [Observer](-observer/index.md)): [LocalSolarEclipseInfo](-local-solar-eclipse-info/index.md)

Searches for the next local solar eclipse in a series.

After using [searchLocalSolarEclipse](search-local-solar-eclipse.md) to find the first solar eclipse in a series, you can call this function to find the next consecutive solar eclipse. Pass in the peak value from the [LocalSolarEclipseInfo](-local-solar-eclipse-info/index.md) returned by the previous call to searchLocalSolarEclipse or nextLocalSolarEclipse to find the next solar eclipse.

See [localSolarEclipsesAfter](local-solar-eclipses-after.md) for convenient iteration of consecutive eclipses.

#### Return

Information about the next solar eclipse visible at the specified observer location.

## Parameters

| | |
|---|---|
| prevEclipseTime | A date and time near a new moon. Solar eclipse search will start at the next new moon. |
| observer | The geographic location of the observer. |
