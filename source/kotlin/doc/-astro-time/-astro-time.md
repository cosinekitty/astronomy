//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[AstroTime](index.md)/[AstroTime](-astro-time.md)

# AstroTime

[jvm]\
fun [AstroTime](-astro-time.md)(ut: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

[jvm]\
fun [AstroTime](-astro-time.md)(d: [Date](https://docs.oracle.com/javase/8/docs/api/java/util/Date.html))

Creates an AstroTime object from a Date object.

## Parameters

jvm

| | |
|---|---|
| d | The date and time to be converted to AstroTime format. |

[jvm]\
fun [AstroTime](-astro-time.md)(d: [Calendar](https://docs.oracle.com/javase/8/docs/api/java/util/Calendar.html))

Creates an AstroTime object from a Calendar object.

## Parameters

jvm

| | |
|---|---|
| d | The date and time to be converted to AstroTime format. |

[jvm]\
fun [AstroTime](-astro-time.md)(year: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), month: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), day: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), hour: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), minute: [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html), second: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

Creates an AstroTime object from a UTC year, month, day, hour, minute and second.

## Parameters

jvm

| | |
|---|---|
| year | The UTC year value. |
| month | The UTC month value 1..12. |
| day | The UTC day of the month 1..31. |
| hour | The UTC hour value 0..23. |
| minute | The UTC minute value 0..59. |
| second | The UTC second in the half-open range [0, 60). |
