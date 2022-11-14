//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[moonNodesAfter](moon-nodes-after.md)

# moonNodesAfter

fun [moonNodesAfter](moon-nodes-after.md)(startTime: [Time](-time/index.md)): [Sequence](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin.sequences/-sequence/index.html)&lt;[NodeEventInfo](-node-event-info/index.md)&gt;

Enumerates a series of consecutive ascending/descending nodes of the Moon.

This function enables iteration through an unlimited number of consecutive lunar nodes starting at a given time. This is a convenience wrapper around [searchMoonNode](search-moon-node.md) and [nextMoonNode](next-moon-node.md).

## Parameters

| | |
|---|---|
| startTime | The date and time for starting the search for a series of ascending/descending nodes of the Moon. |
