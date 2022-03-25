//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[EclipseEvent](index.md)

# EclipseEvent

[jvm]\
class [EclipseEvent](index.md)(time: [AstroTime](../-astro-time/index.md), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

Holds a time and the observed altitude of the Sun at that time.

When reporting a solar eclipse observed at a specific location on the Earth (a "local" solar eclipse), a series of events occur. In addition to the time of each event, it is important to know the altitude of the Sun, because each event may be invisible to the observer if the Sun is below the horizon (i.e. it at night).

If altitude is negative, the event is theoretical only; it would be visible if the Earth were transparent, but the observer cannot actually see it. If altitude is positive but less than a few degrees, visibility will be impaired by atmospheric interference (sunrise or sunset conditions).

## Constructors

| | |
|---|---|
| [EclipseEvent](-eclipse-event.md) | [jvm]<br>fun [EclipseEvent](-eclipse-event.md)(time: [AstroTime](../-astro-time/index.md), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)) |

## Properties

| Name | Summary |
|---|---|
| [altitude](altitude.md) | [jvm]<br>val [altitude](altitude.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>The angular altitude of the center of the Sun above/below the horizon, at time, corrected for atmospheric refraction and expressed in degrees. |
| [time](time.md) | [jvm]<br>val [time](time.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the event. |
