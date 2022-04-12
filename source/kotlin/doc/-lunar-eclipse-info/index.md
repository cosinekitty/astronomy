//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[LunarEclipseInfo](index.md)

# LunarEclipseInfo

[jvm]\
class [LunarEclipseInfo](index.md)(kind: [EclipseKind](../-eclipse-kind/index.md), peak: [Time](../-time/index.md), sdPenum: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), sdPartial: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), sdTotal: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

Information about a lunar eclipse.

Returned by searchLunarEclipse or nextLunarEclipse to report information about a lunar eclipse event. When a lunar eclipse is found, it is classified as penumbral, partial, or total. Penumbral eclipses are difficult to observe, because the Moon is only slightly dimmed by the Earth's penumbra; no part of the Moon touches the Earth's umbra. Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra. Total eclipses occur when the entire Moon passes into the Earth's umbra.

The kind field thus holds EclipseKind.Penumbral, EclipseKind.Partial, or EclipseKind.Total, depending on the kind of lunar eclipse found.

Field peak holds the date and time of the center of the eclipse, when it is at its peak.

Fields sdPenum, sdPartial, and sdTotal hold the semi-duration of each phase of the eclipse, which is half of the amount of time the eclipse spends in each phase (expressed in minutes), or 0.0 if the eclipse never reaches that phase. By converting from minutes to days, and subtracting/adding with peak, the caller may determine the date and time of the beginning/end of each eclipse phase.

## Constructors

| | |
|---|---|
| [LunarEclipseInfo](-lunar-eclipse-info.md) | [jvm]<br>fun [LunarEclipseInfo](-lunar-eclipse-info.md)(kind: [EclipseKind](../-eclipse-kind/index.md), peak: [Time](../-time/index.md), sdPenum: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), sdPartial: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), sdTotal: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [kind](kind.md) | [jvm]<br>val [kind](kind.md): [EclipseKind](../-eclipse-kind/index.md)<br>The type of lunar eclipse found. |
| [peak](peak.md) | [jvm]<br>val [peak](peak.md): [Time](../-time/index.md)<br>The time of the eclipse at its peak. |
| [sdPartial](sd-partial.md) | [jvm]<br>val [sdPartial](sd-partial.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The semi-duration of the partial phase in minutes, or 0.0 if none. |
| [sdPenum](sd-penum.md) | [jvm]<br>val [sdPenum](sd-penum.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The semi-duration of the penumbral phase in minutes. |
| [sdTotal](sd-total.md) | [jvm]<br>val [sdTotal](sd-total.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The semi-duration of the total phase in minutes, or 0.0 if none. |
