//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[HourAngleInfo](index.md)

# HourAngleInfo

[jvm]\
class [HourAngleInfo](index.md)(time: [AstroTime](../-astro-time/index.md), hor: [Topocentric](../-topocentric/index.md))

Information about a celestial body crossing a specific hour angle.

Returned by the function [searchHourAngle](../search-hour-angle.md) to report information about a celestial body crossing a certain hour angle as seen by a specified topocentric observer.

## Constructors

| | |
|---|---|
| [HourAngleInfo](-hour-angle-info.md) | [jvm]<br>fun [HourAngleInfo](-hour-angle-info.md)(time: [AstroTime](../-astro-time/index.md), hor: [Topocentric](../-topocentric/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [hor](hor.md) | [jvm]<br>val [hor](hor.md): [Topocentric](../-topocentric/index.md)<br>Apparent coordinates of the body at the time it crosses the specified hour angle. |
| [time](time.md) | [jvm]<br>val [time](time.md): [AstroTime](../-astro-time/index.md)<br>The date and time when the body crosses the specified hour angle. |
