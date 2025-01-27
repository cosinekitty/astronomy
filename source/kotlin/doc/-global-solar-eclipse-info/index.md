//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[GlobalSolarEclipseInfo](index.md)

# GlobalSolarEclipseInfo

class [GlobalSolarEclipseInfo](index.md)(kind: [EclipseKind](../-eclipse-kind/index.md), obscuration: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), peak: [Time](../-time/index.md), distance: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), latitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), longitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html))

Reports the time and geographic location of the peak of a solar eclipse.

Returned by [searchGlobalSolarEclipse](../search-global-solar-eclipse.md) or [nextGlobalSolarEclipse](../next-global-solar-eclipse.md) to report information about a solar eclipse event.

The eclipse is classified as partial, annular, or total, depending on the maximum amount of the Sun's disc obscured, as seen at the peak location on the surface of the Earth.

The kind field thus holds EclipseKind.Partial, EclipseKind.Annular, or EclipseKind.Total. A total eclipse is when the peak observer sees the Sun completely blocked by the Moon. An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface to completely block the Sun; instead, the Sun takes on a ring-shaped appearance. A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth observes either a total or annular eclipse.

If kind is EclipseKind.Total or EclipseKind.Annular, the latitude and longitude fields give the geographic coordinates of the center of the Moon's shadow projected onto the daytime side of the Earth at the instant of the eclipse's peak. If kind has any other value, latitude and longitude are undefined and should not be used.

For total or annular eclipses, the obscuration field holds the fraction (0, 1] of the Sun's apparent disc area that is blocked from view by the Moon's silhouette, as seen by an observer located at the geographic coordinates latitude, longitude at the darkest time peak. The value will always be 1 for total eclipses, and less than 1 for annular eclipses. For partial eclipses, obscuration is undefined and should not be used. This is because there is little practical use for an obscuration value of a partial eclipse without supplying a particular observation location. Developers who wish to find an obscuration value for partial solar eclipses should therefore use [searchLocalSolarEclipse](../search-local-solar-eclipse.md) and provide the geographic coordinates of an observer.

## Constructors

| | |
|---|---|
| [GlobalSolarEclipseInfo](-global-solar-eclipse-info.md)<br>fun [GlobalSolarEclipseInfo](-global-solar-eclipse-info.md)(kind: [EclipseKind](../-eclipse-kind/index.md), obscuration: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), peak: [Time](../-time/index.md), distance: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), latitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html), longitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [distance](distance.md)<br>val [distance](distance.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers. |
| [kind](kind.md)<br>val [kind](kind.md): [EclipseKind](../-eclipse-kind/index.md)<br>The type of solar eclipse: EclipseKind.Partial, EclipseKind.Annular, or EclipseKind.Total. |
| [latitude](latitude.md)<br>val [latitude](latitude.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>The geographic latitude at the center of the peak eclipse shadow. |
| [longitude](longitude.md)<br>val [longitude](longitude.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>The geographic longitude at the center of the peak eclipse shadow. |
| [obscuration](obscuration.md)<br>val [obscuration](obscuration.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)<br>The peak fraction of the Sun's apparent disc area obscured by the Moon (total and annular eclipses only). |
| [peak](peak.md)<br>val [peak](peak.md): [Time](../-time/index.md)<br>The date and time when the solar eclipse is darkest. This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center. |
