//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Time](index.md)

# Time

class [Time](index.md) : [Comparable](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-comparable/index.html)&lt;[Time](index.md)&gt; 

A date and time used for astronomical calculations.

## Constructors

| | |
|---|---|
| [Time](-time.md)<br>fun [Time](-time.md)(ut: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |
| [Time](-time.md)<br>fun [Time](-time.md)(year: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), month: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), day: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), hour: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), minute: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), second: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Creates a Time object from a UTC year, month, day, hour, minute and second. |

## Types

| Name | Summary |
|---|---|
| [Companion](-companion/index.md)<br>object [Companion](-companion/index.md) |

## Functions

| Name | Summary |
|---|---|
| [addDays](add-days.md)<br>fun [addDays](add-days.md)(days: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Time](index.md)<br>Calculates the sum or difference of an [Time](index.md) with a specified floating point number of days. |
| [compareTo](compare-to.md)<br>open operator override fun [compareTo](compare-to.md)(other: [Time](index.md)): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html)<br>Compares the chronological order of two Time values. |
| [nutationEps](nutation-eps.md)<br>fun [nutationEps](nutation-eps.md)(): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) |
| [nutationPsi](nutation-psi.md)<br>fun [nutationPsi](nutation-psi.md)(): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html) |
| [toDateTime](to-date-time.md)<br>fun [toDateTime](to-date-time.md)(): [DateTime](../-date-time/index.md)<br>Resolves this Time into year, month, day, hour, minute, second. |
| [toMillisecondsSince1970](to-milliseconds-since1970.md)<br>fun [toMillisecondsSince1970](to-milliseconds-since1970.md)(): [Long](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-long/index.html)<br>Converts this Time to the integer number of millseconds since 1970. |
| [toString](to-string.md)<br>open override fun [toString](to-string.md)(): [String](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-string/index.html)<br>Converts this Time to ISO 8601 format, expressed in UTC with millisecond resolution. |

## Properties

| Name | Summary |
|---|---|
| [tt](tt.md)<br>val [tt](tt.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Terrestrial Time days since noon on January 1, 2000. |
| [ut](ut.md)<br>val [ut](ut.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>UT1/UTC number of days since noon on January 1, 2000. |
