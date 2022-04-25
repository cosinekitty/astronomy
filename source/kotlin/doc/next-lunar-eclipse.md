//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[nextLunarEclipse](next-lunar-eclipse.md)

# nextLunarEclipse

fun [nextLunarEclipse](next-lunar-eclipse.md)(prevEclipseTime: [Time](-time/index.md)): [LunarEclipseInfo](-lunar-eclipse-info/index.md)

Searches for the next lunar eclipse in a series.

After using [searchLunarEclipse](search-lunar-eclipse.md) to find the first lunar eclipse in a series, you can call this function to find the next consecutive lunar eclipse. Pass in the center value from the [LunarEclipseInfo](-lunar-eclipse-info/index.md) returned by the previous call to searchLunarEclipse or nextLunarEclipse to find the next lunar eclipse.

#### Return

Information about the next lunar eclipse in a series.

## Parameters

| | |
|---|---|
| prevEclipseTime | A time near a full moon. Lunar eclipse search will start at the next full moon. |
