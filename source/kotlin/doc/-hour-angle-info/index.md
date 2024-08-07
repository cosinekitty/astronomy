//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[HourAngleInfo](index.md)

# HourAngleInfo

class [HourAngleInfo](index.md)(time: [Time](../-time/index.md), hor: [Topocentric](../-topocentric/index.md))

Information about a celestial body crossing a specific hour angle.

Returned by the function [searchHourAngle](../search-hour-angle.md) to report information about a celestial body crossing a certain hour angle as seen by a specified topocentric observer.

## Constructors

| | |
|---|---|
| [HourAngleInfo](-hour-angle-info.md)<br>fun [HourAngleInfo](-hour-angle-info.md)(time: [Time](../-time/index.md), hor: [Topocentric](../-topocentric/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [hor](hor.md)<br>val [hor](hor.md): [Topocentric](../-topocentric/index.md)<br>Apparent coordinates of the body at the time it crosses the specified hour angle. |
| [time](time.md)<br>val [time](time.md): [Time](../-time/index.md)<br>The date and time when the body crosses the specified hour angle. |
