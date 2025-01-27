//[astronomy](../../../../index.md)/[io.github.cosinekitty.astronomy](../../index.md)/[Time](../index.md)/[Companion](index.md)/[fromMillisecondsSince1970](from-milliseconds-since1970.md)

# fromMillisecondsSince1970


@[JvmStatic](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin.jvm/-jvm-static/index.html)

fun [fromMillisecondsSince1970](from-milliseconds-since1970.md)(millis: [Long](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin-stdlib/kotlin/-long/index.html)): [Time](../index.md)

Creates a Time object from the number of milliseconds since the 1970 epoch.

Operating systems and runtime libraries commonly measure civil time in integer milliseconds since January 1, 1970 at 00:00 UTC (midnight). To facilitate using such values for astronomy calculations, this function converts a millsecond count into a Time object.
