//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Vector](index.md)/[toHorizontal](to-horizontal.md)

# toHorizontal

[jvm]\
fun [toHorizontal](to-horizontal.md)(refraction: [Refraction](../-refraction/index.md)): [Spherical](../-spherical/index.md)

Converts Cartesian coordinates to horizontal coordinates.

Given a horizontal Cartesian vector, returns horizontal azimuth and altitude. *IMPORTANT:* This function differs from [Vector.toSpherical](to-spherical.md) in two ways:

- 
   toSpherical returns a lon value that represents azimuth defined counterclockwise from north (e.g., west = +90), but this function represents a clockwise rotation (e.g., east = +90). The difference is because toSpherical is intended to preserve the vector "right-hand rule", while this function defines azimuth in a more traditional way as used in navigation and cartography.
- 
   This function optionally corrects for atmospheric refraction, while toSpherical does not.

The returned object contains the azimuth in lon. It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.

The altitude is stored in lat.

The distance to the observed object is stored in dist, and is expressed in astronomical units (AU).

## Parameters

jvm

| | |
|---|---|
| refraction | `Refraction.None`: no atmospheric refraction correction is performed.     `Refraction.Normal`: correct altitude for atmospheric refraction.     `Refraction.JplHor`: for JPL Horizons compatibility testing only; not recommended for normal use. |
