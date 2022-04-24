//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchMoonNode](search-moon-node.md)

# searchMoonNode

[jvm]\
fun [searchMoonNode](search-moon-node.md)(startTime: [Time](-time/index.md)): [NodeEventInfo](-node-event-info/index.md)

Searches for a time when the Moon's center crosses through the ecliptic plane.

Searches for the first ascending or descending node of the Moon after startTime. An ascending node is when the Moon's center passes through the ecliptic plane (the plane of the Earth's orbit around the Sun) from south to north. A descending node is when the Moon's center passes through the ecliptic plane from north to south. Nodes indicate possible times of solar or lunar eclipses, if the Moon also happens to be in the correct phase (new or full, respectively). Call searchMoonNode to find the first of a series of nodes. Then call [nextMoonNode](next-moon-node.md) to find as many more consecutive nodes as desired.

## Parameters

jvm

| | |
|---|---|
| startTime | The date and time for starting the search for an ascending or descending node of the Moon. |
