//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchLocalSolarEclipse](search-local-solar-eclipse.md)

# searchLocalSolarEclipse

fun [searchLocalSolarEclipse](search-local-solar-eclipse.md)(startTime: [Time](-time/index.md), observer: [Observer](-observer/index.md)): [LocalSolarEclipseInfo](-local-solar-eclipse-info/index.md)

Searches for a solar eclipse visible at a specific location on the Earth's surface.

This function finds the first solar eclipse that occurs after startTime. A solar eclipse may be partial, annular, or total. See [LocalSolarEclipseInfo](-local-solar-eclipse-info/index.md) for more information.

To find a series of solar eclipses, call this function once, then keep calling #Astronomy.NextLocalSolarEclipse as many times as desired, passing in the peak value returned from the previous call.

IMPORTANT: An eclipse reported by this function might be partly or completely invisible to the observer due to the time of day. See [LocalSolarEclipseInfo](-local-solar-eclipse-info/index.md) for more information about this topic.

#### Return

Information about the first solar eclipse visible at the specified observer location.

## Parameters

| | |
|---|---|
| startTime | The date and time for starting the search for a solar eclipse. |
| observer | The geographic location of the observer. |
