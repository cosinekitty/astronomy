//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[nextMoonNode](next-moon-node.md)

# nextMoonNode

fun [nextMoonNode](next-moon-node.md)(prevNode: [NodeEventInfo](-node-event-info/index.md)): [NodeEventInfo](-node-event-info/index.md)

Searches for the next time when the Moon's center crosses through the ecliptic plane.

Call [searchMoonNode](search-moon-node.md) to find the first of a series of nodes. Then call nextMoonNode to find as many more consecutive nodes as desired.

See [moonNodesAfter](moon-nodes-after.md) for convenient iteration of consecutive nodes.

## Parameters

| | |
|---|---|
| prevNode | The previous node found from calling [searchMoonNode](search-moon-node.md) or nextMoonNode. |
