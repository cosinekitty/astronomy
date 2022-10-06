//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[DateTime](index.md)

# DateTime

class [DateTime](index.md)(year: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), month: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), day: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), hour: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), minute: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), second: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

A universal time resolved into UTC calendar date and time fields.

## Constructors

| | |
|---|---|
| [DateTime](-date-time.md)<br>fun [DateTime](-date-time.md)(year: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), month: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), day: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), hour: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), minute: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), second: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |

## Functions

| Name | Summary |
|---|---|
| [toDays](to-days.md)<br>fun [toDays](to-days.md)(): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Convert this date and time to the floating point number of days since the J2000 epoch. |
| [toString](to-string.md)<br>open override fun [toString](to-string.md)(): [String](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-string/index.html)<br>Converts this DateTime to ISO 8601 format, expressed in UTC with millisecond resolution. |
| [toTime](to-time.md)<br>fun [toTime](to-time.md)(): [Time](../-time/index.md)<br>Convert this date and time to a [Time](../-time/index.md) value that can be used for astronomy calculations. |

## Properties

| Name | Summary |
|---|---|
| [day](day.md)<br>val [day](day.md): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html)<br>The day of the month, 1..31. |
| [hour](hour.md)<br>val [hour](hour.md): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html)<br>The hour of the day, 0..23. |
| [minute](minute.md)<br>val [minute](minute.md): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html)<br>The minute value 0..59. |
| [month](month.md)<br>val [month](month.md): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html)<br>The month value 1=January, ..., 12=December. |
| [second](second.md)<br>val [second](second.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The floating point second value in the half-open range [0, 60). |
| [year](year.md)<br>val [year](year.md): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html)<br>The integer year value. |
