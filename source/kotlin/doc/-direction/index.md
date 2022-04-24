//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Direction](index.md)

# Direction

[jvm]\
enum [Direction](index.md) : [Enum](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-enum/index.html)&lt;[Direction](index.md)&gt; 

Selects whether to search for a rising event or a setting event for a celestial body.

## Entries

| | |
|---|---|
| [Set](-set/index.md) | [jvm]<br>[Set](-set/index.md)(-1)<br>Indicates a setting event: a celestial body is observed to sink below the horizon by an observer on the Earth. |
| [Rise](-rise/index.md) | [jvm]<br>[Rise](-rise/index.md)(+1)<br>Indicates a rising event: a celestial body is observed to rise above the horizon by an observer on the Earth. |

## Properties

| Name | Summary |
|---|---|
| [name](../-node-event-kind/-ascending/index.md#-372974862%2FProperties%2F-1216412040) | [jvm]<br>val [name](../-node-event-kind/-ascending/index.md#-372974862%2FProperties%2F-1216412040): [String](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-string/index.html) |
| [ordinal](../-node-event-kind/-ascending/index.md#-739389684%2FProperties%2F-1216412040) | [jvm]<br>val [ordinal](../-node-event-kind/-ascending/index.md#-739389684%2FProperties%2F-1216412040): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html) |
| [sign](sign.md) | [jvm]<br>val [sign](sign.md): [Int](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-int/index.html)<br>A numeric value that is helpful in formulas involving rise/set. The sign is +1 for a rising event, or -1 for a setting event. |
