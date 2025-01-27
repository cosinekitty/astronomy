//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Time](index.md)/[addDays](add-days.md)

# addDays

fun [addDays](add-days.md)(days: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-double/index.html)): [Time](index.md)

Calculates the sum or difference of an [Time](index.md) with a specified floating point number of days.

Sometimes we need to adjust a given [Time](index.md) value by a certain amount of time. This function adds the given real number of days in days to the date and time in this object.

More precisely, the result's Universal Time field ut is exactly adjusted by days and the Terrestrial Time field tt is adjusted for the resulting UTC date and time, using a best-fit piecewise polynomial model devised by [Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).

#### Return

A date and time that is conceptually equal to time + days.

## Parameters

| | |
|---|---|
| days | A floating point number of days by which to adjust time. May be negative, 0, or positive. |
