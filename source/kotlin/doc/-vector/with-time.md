//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Vector](index.md)/[withTime](with-time.md)

# withTime

fun [withTime](with-time.md)(time: [Time](../-time/index.md)): [Vector](index.md)

Creates a new vector with the same coordinates but a different time.

Usually it is a mistake to add or subtract vectors corresponding to different times. The overloaded operators for adding and subtracting vectors will throw an exception if the times do not match. However, occasionally it is helpful to adjust the time associated with a vector to get around this safety check. For example, the time an event occurs may be different from the time it is observed. The geoVector function stores the observation time in the returned vector, while helioVector stores the event time.

#### Return

A cloned vector that has the same coordinates as this one, but at a different time.

## Parameters

| | |
|---|---|
| time | A time to include in a new vector with the same coordinates as this one. |
