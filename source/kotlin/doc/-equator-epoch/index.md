//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[EquatorEpoch](index.md)

# EquatorEpoch

[jvm]\
enum [EquatorEpoch](index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[EquatorEpoch](index.md)&gt; 

Selects the date for which the Earth's equator is to be used for representing equatorial coordinates.

The Earth's equator is not always in the same plane due to precession and nutation.

Sometimes it is useful to have a fixed plane of reference for equatorial coordinates across different calendar dates.  In these cases, a fixed *epoch*, or reference time, is helpful. Astronomy Engine provides the J2000 epoch for such cases.  This refers to the plane of the Earth's orbit as it was on noon UTC on 1 January 2000.

For some other purposes, it is more helpful to represent coordinates using the Earth's equator exactly as it is on that date. For example, when calculating rise/set times or horizontal coordinates, it is most accurate to use the orientation of the Earth's equator at that same date and time. For these uses, Astronomy Engine allows *of-date* calculations.

## Entries

| | |
|---|---|
| [OfDate](-of-date/index.md) | [jvm]<br>[OfDate](-of-date/index.md)()<br>Represent equatorial coordinates using the Earth's equator at the given date and time. |
| [J2000](-j2000/index.md) | [jvm]<br>[J2000](-j2000/index.md)()<br>Represent equatorial coordinates in the J2000 epoch. |

## Properties

| Name | Summary |
|---|---|
| [name](../-eclipse-kind/-none/index.md#-372974862%2FProperties%2F-1216412040) | [jvm]<br>val [name](../-eclipse-kind/-none/index.md#-372974862%2FProperties%2F-1216412040): [String](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-string/index.html) |
| [ordinal](../-eclipse-kind/-none/index.md#-739389684%2FProperties%2F-1216412040) | [jvm]<br>val [ordinal](../-eclipse-kind/-none/index.md#-739389684%2FProperties%2F-1216412040): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html) |
