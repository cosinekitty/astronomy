//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[equatorialToEcliptic](equatorial-to-ecliptic.md)

# equatorialToEcliptic

fun [equatorialToEcliptic](equatorial-to-ecliptic.md)(equ: [Vector](-vector/index.md)): [Ecliptic](-ecliptic/index.md)

Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.

Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates, which are relative to the plane of the Earth's orbit around the Sun.

#### Return

Ecliptic coordinates in the J2000 frame of reference (ECL).

## Parameters

jvm

| | |
|---|---|
| equ | Equatorial coordinates in the J2000 frame of reference.     You can call [geoVector] to obtain suitable equatorial coordinates. |
