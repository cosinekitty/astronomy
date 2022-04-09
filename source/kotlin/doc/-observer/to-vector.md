//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Observer](index.md)/[toVector](to-vector.md)

# toVector

[jvm]\
fun [toVector](to-vector.md)(time: [AstroTime](../-astro-time/index.md), equator: [EquatorEpoch](../-equator-epoch/index.md)): [AstroVector](../-astro-vector/index.md)

Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.

This function calculates a vector from the center of the Earth to a point on or near the surface of the Earth, expressed in equatorial coordinates. It takes into account the rotation of the Earth at the given time, along with the given latitude, longitude, and elevation of the observer.

The caller may pass a value in equator to select either [EquatorEpoch.J2000](../-equator-epoch/-j2000/index.md) for using J2000 coordinates, or [EquatorEpoch.OfDate](../-equator-epoch/-of-date/index.md) for using coordinates relative to the Earth's equator at the specified time.

The returned vector has components expressed in astronomical units (AU). To convert to kilometers, multiply the vector values by the scalar value [KM_PER_AU](../-k-m_-p-e-r_-a-u.md).

The inverse of this function is also available: AstroVector.toObserver.

#### Return

A vector from the center of the Earth to this geographic location.

## Parameters

jvm

| | |
|---|---|
| time | The date and time for which to calculate the observer's position vector. |
| equator | Selects the date of the Earth's equator in which to express the equatorial coordinates.     The caller may select [EquatorEpoch.J2000] to use the orientation of the Earth's equator     at noon UTC on January 1, 2000, in which case this function corrects for precession     and nutation of the Earth as it was at the moment specified by the `time` parameter.     Or the caller may select [EquatorEpoch.OfDate] to use the Earth's equator at `time`     as the orientation. |
