//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[AstroTime](index.md)

# AstroTime

[jvm]\
class [AstroTime](index.md)

A date and time used for astronomical calculations.

## Constructors

| | |
|---|---|
| [AstroTime](-astro-time.md) | [jvm]<br>fun [AstroTime](-astro-time.md)(ut: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |
| [AstroTime](-astro-time.md) | [jvm]<br>fun [AstroTime](-astro-time.md)(d: [Date](https://docs.oracle.com/javase/8/docs/api/java/util/Date.html))<br>Creates an AstroTime object from a Date object. |
| [AstroTime](-astro-time.md) | [jvm]<br>fun [AstroTime](-astro-time.md)(d: [Calendar](https://docs.oracle.com/javase/8/docs/api/java/util/Calendar.html))<br>Creates an AstroTime object from a Calendar object. |
| [AstroTime](-astro-time.md) | [jvm]<br>fun [AstroTime](-astro-time.md)(year: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), month: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), day: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), hour: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), minute: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), second: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))<br>Creates an AstroTime object from a UTC year, month, day, hour, minute and second. |

## Types

| Name | Summary |
|---|---|
| [Companion](-companion/index.md) | [jvm]<br>object [Companion](-companion/index.md) |

## Functions

| Name | Summary |
|---|---|
| [addDays](add-days.md) | [jvm]<br>fun [addDays](add-days.md)(days: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [AstroTime](index.md)<br>Calculates the sum or difference of an [AstroTime](index.md) with a specified floating point number of days. |
| [toDate](to-date.md) | [jvm]<br>fun [toDate](to-date.md)(): [Date](https://docs.oracle.com/javase/8/docs/api/java/util/Date.html)<br>Converts this object to a native Date equivalent. |
| [toString](to-string.md) | [jvm]<br>open override fun [toString](to-string.md)(): [String](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-string/index.html)<br>Converts this AstroTime to ISO 8601 format, expressed in UTC with millisecond resolution. |

## Properties

| Name | Summary |
|---|---|
| [tt](tt.md) | [jvm]<br>val [tt](tt.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Terrestrial Time days since noon on January 1, 2000. |
| [ut](ut.md) | [jvm]<br>val [ut](ut.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>UT1/UTC number of days since noon on January 1, 2000. |
