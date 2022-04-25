//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[LocalSolarEclipseInfo](index.md)

# LocalSolarEclipseInfo

class [LocalSolarEclipseInfo](index.md)(kind: [EclipseKind](../-eclipse-kind/index.md), partialBegin: [EclipseEvent](../-eclipse-event/index.md), totalBegin: [EclipseEvent](../-eclipse-event/index.md)?, peak: [EclipseEvent](../-eclipse-event/index.md), totalEnd: [EclipseEvent](../-eclipse-event/index.md)?, partialEnd: [EclipseEvent](../-eclipse-event/index.md))

Information about a solar eclipse as seen by an observer at a given time and geographic location.

Returned by [searchLocalSolarEclipse](../search-local-solar-eclipse.md) or [nextLocalSolarEclipse](../next-local-solar-eclipse.md) to report information about a solar eclipse as seen at a given geographic location.

When a solar eclipse is found, it is classified as partial, annular, or total. The kind field thus holds EclipseKind.Partial, EclipseKind.Annular, or EclipseKind.Total. A partial solar eclipse is when the Moon does not line up directly enough with the Sun to completely block the Sun's light from reaching the observer. An annular eclipse occurs when the Moon's disc is completely visible against the Sun but the Moon is too far away to completely block the Sun's light; this leaves the Sun with a ring-like appearance. A total eclipse occurs when the Moon is close enough to the Earth and aligned with the Sun just right to completely block all sunlight from reaching the observer.

There are 5 "event" fields, each of which contains a time and a solar altitude. Field peak holds the date and time of the center of the eclipse, when it is at its peak. The fields partialBegin and partialEnd are always set, and indicate when the eclipse begins/ends. If the eclipse reaches totality or becomes annular, totalBegin and totalEnd indicate when the total/annular phase begins/ends. When an event field is valid, the caller must also check its altitude field to see whether the Sun is above the horizon at the time indicated by the time field. </remarks>

## Constructors

| | |
|---|---|
| [LocalSolarEclipseInfo](-local-solar-eclipse-info.md)<br>fun [LocalSolarEclipseInfo](-local-solar-eclipse-info.md)(kind: [EclipseKind](../-eclipse-kind/index.md), partialBegin: [EclipseEvent](../-eclipse-event/index.md), totalBegin: [EclipseEvent](../-eclipse-event/index.md)?, peak: [EclipseEvent](../-eclipse-event/index.md), totalEnd: [EclipseEvent](../-eclipse-event/index.md)?, partialEnd: [EclipseEvent](../-eclipse-event/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [kind](kind.md)<br>val [kind](kind.md): [EclipseKind](../-eclipse-kind/index.md)<br>The type of solar eclipse: EclipseKind.Partial, EclipseKind.Annular, or EclipseKind.Total. |
| [partialBegin](partial-begin.md)<br>val [partialBegin](partial-begin.md): [EclipseEvent](../-eclipse-event/index.md)<br>The time and Sun altitude at the beginning of the eclipse. |
| [partialEnd](partial-end.md)<br>val [partialEnd](partial-end.md): [EclipseEvent](../-eclipse-event/index.md)<br>The time and Sun altitude at the end of the eclipse. |
| [peak](peak.md)<br>val [peak](peak.md): [EclipseEvent](../-eclipse-event/index.md)<br>The time and Sun altitude when the eclipse reaches its peak. |
| [totalBegin](total-begin.md)<br>val [totalBegin](total-begin.md): [EclipseEvent](../-eclipse-event/index.md)?<br>If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise null. |
| [totalEnd](total-end.md)<br>val [totalEnd](total-end.md): [EclipseEvent](../-eclipse-event/index.md)?<br>If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise null. |
