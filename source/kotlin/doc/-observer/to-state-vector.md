//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Observer](index.md)/[toStateVector](to-state-vector.md)

# toStateVector

fun [toStateVector](to-state-vector.md)(time: [Time](../-time/index.md), equator: [EquatorEpoch](../-equator-epoch/index.md)): [StateVector](../-state-vector/index.md)

Calculates geocentric equatorial position and velocity of an observer on the surface of the Earth.

This function calculates position and velocity vectors of an observer on or near the surface of the Earth, expressed in equatorial coordinates. It takes into account the rotation of the Earth at the given time, along with the given latitude, longitude, and elevation of the observer.

The caller may pass a value in equator to select either [EquatorEpoch.J2000](../-equator-epoch/-j2000/index.md) for using J2000 coordinates, or [EquatorEpoch.OfDate](../-equator-epoch/-of-date/index.md) for using coordinates relative to the Earth's equator at the specified time.

The returned position vector has components expressed in astronomical units (AU). To convert to kilometers, multiply the vector values by the scalar value [KM_PER_AU](../-k-m_-p-e-r_-a-u.md).

The returned velocity vector is measured in AU/day.

#### Return

The position and velocity of this observer with respect to the Earth's center.

## Parameters

jvm

| | |
|---|---|
| time | The date and time for which to calculate the observer's position vector. |
| equator | Selects the date of the Earth's equator in which to express the equatorial coordinates.     The caller may select [EquatorEpoch.J2000] to use the orientation of the Earth's equator     at noon UTC on January 1, 2000, in which case this function corrects for precession     and nutation of the Earth as it was at the moment specified by the `time` parameter.     Or the caller may select [EquatorEpoch.OfDate] to use the Earth's equator at `time`     as the orientation. |
