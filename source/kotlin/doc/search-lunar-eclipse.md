//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchLunarEclipse](search-lunar-eclipse.md)

# searchLunarEclipse

fun [searchLunarEclipse](search-lunar-eclipse.md)(startTime: [Time](-time/index.md)): [LunarEclipseInfo](-lunar-eclipse-info/index.md)

Searches for a lunar eclipse.

This function finds the first lunar eclipse that occurs after startTime. A lunar eclipse may be penumbral, partial, or total. See [LunarEclipseInfo](-lunar-eclipse-info/index.md) for more information. To find a series of lunar eclipses, call this function once, then keep calling [nextLunarEclipse](next-lunar-eclipse.md) as many times as desired, passing in the center value returned from the previous call.

#### Return

Information about the first lunar eclipse that occurs after startTime.

## Parameters

| | |
|---|---|
| startTime | The date and time for starting the search for a lunar eclipse. |
