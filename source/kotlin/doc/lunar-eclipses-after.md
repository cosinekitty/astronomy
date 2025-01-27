//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[lunarEclipsesAfter](lunar-eclipses-after.md)

# lunarEclipsesAfter

fun [lunarEclipsesAfter](lunar-eclipses-after.md)(startTime: [Time](-time/index.md)): [Sequence](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin.sequences/-sequence/index.html)&lt;[LunarEclipseInfo](-lunar-eclipse-info/index.md)&gt;

Enumerates a series of consecutive lunar eclipses that occur after a given time.

This function enables iteration through an unlimited number of consecutive lunar eclipses starting at a given time.

This is a convenience wrapper around [searchLunarEclipse](search-lunar-eclipse.md) and [nextLunarEclipse](next-lunar-eclipse.md).

## Parameters

| | |
|---|---|
| startTime | The date and time for starting the search for a series of lunar eclipses. |
