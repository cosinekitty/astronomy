//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Vector](index.md)/[toObserver](to-observer.md)

# toObserver

fun [toObserver](to-observer.md)(equator: [EquatorEpoch](../-equator-epoch/index.md)): [Observer](../-observer/index.md)

Calculates the geographic location corresponding to a geocentric equatorial vector.

This is the inverse function of [Observer.toVector](../-observer/to-vector.md). Given an equatorial vector from the center of the Earth to an observer on or near the Earth's surface, this function returns the geographic latitude, longitude, and elevation for that observer.

#### Return

The geographic coordinates corresponding to the vector.

## Parameters

| | |
|---|---|
| equator | Selects the date of the Earth's equator in which this vector is expressed.     The caller may select [EquatorEpoch.J2000] to use the orientation of the Earth's equator     at noon UTC on January 1, 2000, in which case this function corrects for precession     and nutation of the Earth as it was at the moment specified by the time `this.t`.     Or the caller may select [EquatorEpoch.OfDate] to use the Earth's equator at `this.t`     as the orientation. |
