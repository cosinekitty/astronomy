//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[equatorialToEcliptic](equatorial-to-ecliptic.md)

# equatorialToEcliptic

fun [equatorialToEcliptic](equatorial-to-ecliptic.md)(eqj: [Vector](-vector/index.md)): [Ecliptic](-ecliptic/index.md)

Converts a J2000 mean equator (EQJ) vector to a true ecliptic of date (ETC) vector and angles.

Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC on 1 January 2000), this function converts those coordinates to true ecliptic coordinates of date, which are relative to the plane of the Earth's orbit around the Sun.

#### Return

Spherical and vector coordinates expressed in true ecliptic coordinates of date (ECT)..

## Parameters

| | |
|---|---|
| eqj | Equatorial coordinates in the J2000 frame of reference. You can call [geoVector](geo-vector.md) to obtain suitable equatorial coordinates. |
