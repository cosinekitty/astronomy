//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[search](search.md)

# search

[jvm]\
fun [search](search.md)(func: [SearchContext](../-search-context/index.md), time1: [AstroTime](../-astro-time/index.md), time2: [AstroTime](../-astro-time/index.md), toleranceSeconds: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [AstroTime](../-astro-time/index.md)?

Searches for a time at which a function's value increases through zero.

Certain astronomy calculations involve finding a time when an event occurs. Often such events can be defined as the root of a function: the time at which the function's value becomes zero.

search finds the *ascending root* of a function: the time at which the function's value becomes zero while having a positive slope. That is, as time increases, the function transitions from a negative value, through zero at a specific moment, to a positive value later. The goal of the search is to find that specific moment.

The func parameter is an instance of the interface [SearchContext](../-search-context/index.md). As an example, a caller may wish to find the moment a celestial body reaches a certain ecliptic longitude. In that case, the caller might derive a class that contains a [Body](../-body/index.md) member to specify the body and a Double to hold the target longitude. It could subtract the target longitude from the actual longitude at a given time; thus the difference would equal zero at the moment in time the planet reaches the desired longitude.

Every time it is called, func.eval returns a Double value or it throws an exception. If func.eval throws an exception, the search immediately fails and the exception is propagated to the caller. Otherwise, the search proceeds until it either finds the ascending root or fails for some reason.

The search calls func.eval repeatedly to rapidly narrow in on any ascending root within the time window specified by time1 and time2. The search never reports a solution outside this time window.

search uses a combination of bisection and quadratic interpolation to minimize the number of function calls. However, it is critical that the supplied time window be small enough that there cannot be more than one root (ascedning or descending) within it; otherwise the search can fail. Beyond that, it helps to make the time window as small as possible, ideally such that the function itself resembles a smooth parabolic curve within that window.

If an ascending root is not found, or more than one root (ascending and/or descending) exists within the window time1..time2, the search will return null.

If the search does not converge within 20 iterations, it will throw an exception.

## Parameters

jvm

| | |
|---|---|
| func | The function for which to find the time of an ascending root.     See remarks above for more details. |
| time1 | The lower time bound of the search window.     See remarks above for more details. |
| time2 | The upper time bound of the search window.     See remarks above for more details. |
| toleranceSeconds | Specifies an amount of time in seconds within which a bounded ascending root     is considered accurate enough to stop. A typical value is 1 second. |
