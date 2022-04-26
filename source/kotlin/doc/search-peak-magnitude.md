//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchPeakMagnitude](search-peak-magnitude.md)

# searchPeakMagnitude

fun [searchPeakMagnitude](search-peak-magnitude.md)(body: [Body](-body/index.md), startTime: [Time](-time/index.md)): [IlluminationInfo](-illumination-info/index.md)

Searches for the date and time Venus will next appear brightest as seen from the Earth.

This function searches for the date and time Venus appears brightest as seen from the Earth. Currently only Venus is supported for the body parameter, though this could change in the future. Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from the Earth, so peak magnitude events have little practical value for that planet. Planets other than Venus and Mercury reach peak magnitude at opposition, which can be found using #Astronomy.SearchRelativeLongitude. The Moon reaches peak magnitude at full moon, which can be found using [searchMoonQuarter](search-moon-quarter.md) or [searchMoonPhase](search-moon-phase.md). The Sun reaches peak magnitude at perihelion, which occurs each year in January. However, the difference is minor and has little practical value.

## Parameters

| | |
|---|---|
| body | Currently only [Body.Venus] is allowed. Any other value causes an exception. |
| startTime | The date and time to start searching for the next peak magnitude event. |
