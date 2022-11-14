//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Time](index.md)

# Time

class [Time](index.md) : Comparable&lt;[Time](index.md)&gt; 

A date and time used for astronomical calculations.

## Constructors

| | |
|---|---|
| [Time](-time.md)<br>fun [Time](-time.md)(ut: Double) |
| [Time](-time.md)<br>fun [Time](-time.md)(year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Double)<br>Creates a Time object from a UTC year, month, day, hour, minute and second. |

## Types

| Name | Summary |
|---|---|
| [Companion](-companion/index.md)<br>object [Companion](-companion/index.md) |

## Functions

| Name | Summary |
|---|---|
| [addDays](add-days.md)<br>fun [addDays](add-days.md)(days: Double): [Time](index.md)<br>Calculates the sum or difference of an [Time](index.md) with a specified floating point number of days. |
| [compareTo](compare-to.md)<br>open operator override fun [compareTo](compare-to.md)(other: [Time](index.md)): Int<br>Compares the chronological order of two Time values. |
| [toDateTime](to-date-time.md)<br>fun [toDateTime](to-date-time.md)(): [DateTime](../-date-time/index.md)<br>Resolves this Time into year, month, day, hour, minute, second. |
| [toMillisecondsSince1970](to-milliseconds-since1970.md)<br>fun [toMillisecondsSince1970](to-milliseconds-since1970.md)(): Long<br>Converts this Time to the integer number of millseconds since 1970. |
| [toString](to-string.md)<br>open override fun [toString](to-string.md)(): String<br>Converts this Time to ISO 8601 format, expressed in UTC with millisecond resolution. |

## Properties

| Name | Summary |
|---|---|
| [tt](tt.md)<br>val [tt](tt.md): Double<br>Terrestrial Time days since noon on January 1, 2000. |
| [ut](ut.md)<br>val [ut](ut.md): Double<br>UT1/UTC number of days since noon on January 1, 2000. |
