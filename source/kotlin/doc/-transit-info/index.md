//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[TransitInfo](index.md)

# TransitInfo

class [TransitInfo](index.md)(start: [Time](../-time/index.md), peak: [Time](../-time/index.md), finish: [Time](../-time/index.md), separation: Double)

Information about a transit of Mercury or Venus, as seen from the Earth.

Returned by [searchTransit](../search-transit.md) or [nextTransit](../next-transit.md) to report information about a transit of Mercury or Venus. A transit is when Mercury or Venus passes between the Sun and Earth so that the other planet is seen in silhouette against the Sun.

The start field reports the moment in time when the planet first becomes visible against the Sun in its background. The peak field reports when the planet is most aligned with the Sun, as seen from the Earth. The finish field reports the last moment when the planet is visible against the Sun in its background.

The calculations are performed from the point of view of a geocentric observer.

## Constructors

| | |
|---|---|
| [TransitInfo](-transit-info.md)<br>fun [TransitInfo](-transit-info.md)(start: [Time](../-time/index.md), peak: [Time](../-time/index.md), finish: [Time](../-time/index.md), separation: Double) |

## Properties

| Name | Summary |
|---|---|
| [finish](finish.md)<br>val [finish](finish.md): [Time](../-time/index.md)<br>Date and time at the end of the transit. |
| [peak](peak.md)<br>val [peak](peak.md): [Time](../-time/index.md)<br>Date and time of the peak of the transit. |
| [separation](separation.md)<br>val [separation](separation.md): Double<br>Angular separation in arcminutes between the centers of the Sun and the planet at time peak. |
| [start](start.md)<br>val [start](start.md): [Time](../-time/index.md)<br>Date and time at the beginning of the transit. |
