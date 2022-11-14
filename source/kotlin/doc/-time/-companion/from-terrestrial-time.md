//[astronomy](../../../../index.md)/[io.github.cosinekitty.astronomy](../../index.md)/[Time](../index.md)/[Companion](index.md)/[fromTerrestrialTime](from-terrestrial-time.md)

# fromTerrestrialTime


@[JvmStatic](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin.jvm/-jvm-static/index.html)

fun [fromTerrestrialTime](from-terrestrial-time.md)(tt: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Time](../index.md)

Creates a Time object from a Terrestrial Time day value.

This function can be used in rare cases where a time must be based on Terrestrial Time (TT) rather than Universal Time (UT). Most developers will want to use Time(ut) with a universal time instead of this function, because usually time is based on civil time adjusted by leap seconds to match the Earth's rotation, rather than the uniformly flowing TT used to calculate solar system dynamics. In rare cases where the caller already knows TT, this function is provided to create a Time value that can be passed to Astronomy Engine functions.

## Parameters

| | |
|---|---|
| tt | The number of days after the J2000 epoch. |
