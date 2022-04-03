//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[NodeEventInfo](index.md)

# NodeEventInfo

[jvm]\
class [NodeEventInfo](index.md)(time: [AstroTime](../-astro-time/index.md), kind: [NodeEventKind](../-node-event-kind/index.md))

Information about an ascending or descending node of a body.

This object is returned by Astronomy.searchMoonNode and Astronomy.nextMoonNode to report information about the center of the Moon passing through the ecliptic plane.

## Constructors

| | |
|---|---|
| [NodeEventInfo](-node-event-info.md) | [jvm]<br>fun [NodeEventInfo](-node-event-info.md)(time: [AstroTime](../-astro-time/index.md), kind: [NodeEventKind](../-node-event-kind/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [kind](kind.md) | [jvm]<br>val [kind](kind.md): [NodeEventKind](../-node-event-kind/index.md)<br>Whether the node is ascending or descending. |
| [time](time.md) | [jvm]<br>val [time](time.md): [AstroTime](../-astro-time/index.md)<br>The time of the body's node. |
