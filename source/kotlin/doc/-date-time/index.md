//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[DateTime](index.md)

# DateTime

class [DateTime](index.md)(year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Double)

A universal time resolved into UTC calendar date and time fields.

## Constructors

| | |
|---|---|
| [DateTime](-date-time.md)<br>fun [DateTime](-date-time.md)(year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Double) |

## Functions

| Name | Summary |
|---|---|
| [toDays](to-days.md)<br>fun [toDays](to-days.md)(): Double<br>Convert this date and time to the floating point number of days since the J2000 epoch. |
| [toString](to-string.md)<br>open override fun [toString](to-string.md)(): String<br>Converts this DateTime to ISO 8601 format, expressed in UTC with millisecond resolution. |
| [toTime](to-time.md)<br>fun [toTime](to-time.md)(): [Time](../-time/index.md)<br>Convert this date and time to a [Time](../-time/index.md) value that can be used for astronomy calculations. |

## Properties

| Name | Summary |
|---|---|
| [day](day.md)<br>val [day](day.md): Int<br>The day of the month, 1..31. |
| [hour](hour.md)<br>val [hour](hour.md): Int<br>The hour of the day, 0..23. |
| [minute](minute.md)<br>val [minute](minute.md): Int<br>The minute value 0..59. |
| [month](month.md)<br>val [month](month.md): Int<br>The month value 1=January, ..., 12=December. |
| [second](second.md)<br>val [second](second.md): Double<br>The floating point second value in the half-open range [0, 60). |
| [year](year.md)<br>val [year](year.md): Int<br>The integer year value. |
