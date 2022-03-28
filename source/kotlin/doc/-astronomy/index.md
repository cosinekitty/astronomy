//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)

# Astronomy

[jvm]\
object [Astronomy](index.md)

The main container of astronomy calculation functions.

## Functions

| Name | Summary |
|---|---|
| [eclipticGeoMoon](ecliptic-geo-moon.md) | [jvm]<br>fun [eclipticGeoMoon](ecliptic-geo-moon.md)(time: [AstroTime](../-astro-time/index.md)): [Spherical](../-spherical/index.md)<br>Calculates spherical ecliptic geocentric position of the Moon. |
| [equatorFromVector](equator-from-vector.md) | [jvm]<br>fun [equatorFromVector](equator-from-vector.md)(vector: [AstroVector](../-astro-vector/index.md)): [Equatorial](../-equatorial/index.md)<br>Given an equatorial vector, calculates equatorial angular coordinates. |
| [geoMoon](geo-moon.md) | [jvm]<br>fun [geoMoon](geo-moon.md)(time: [AstroTime](../-astro-time/index.md)): [AstroVector](../-astro-vector/index.md)<br>Calculates equatorial geocentric position of the Moon at a given time. |
| [helioDistance](helio-distance.md) | [jvm]<br>fun [helioDistance](helio-distance.md)(body: [Body](../-body/index.md), time: [AstroTime](../-astro-time/index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Calculates the distance between a body and the Sun at a given time. |
| [helioVector](helio-vector.md) | [jvm]<br>fun [helioVector](helio-vector.md)(body: [Body](../-body/index.md), time: [AstroTime](../-astro-time/index.md)): [AstroVector](../-astro-vector/index.md)<br>Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system. |
| [inverseRefractionAngle](inverse-refraction-angle.md) | [jvm]<br>fun [inverseRefractionAngle](inverse-refraction-angle.md)(refraction: [Refraction](../-refraction/index.md), bentAltitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Calculates the inverse of an atmospheric refraction angle. |
| [massProduct](mass-product.md) | [jvm]<br>fun [massProduct](mass-product.md)(body: [Body](../-body/index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Returns the product of mass and universal gravitational constant of a Solar System body. |
| [refractionAngle](refraction-angle.md) | [jvm]<br>fun [refractionAngle](refraction-angle.md)(refraction: [Refraction](../-refraction/index.md), altitude: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction. |
| [rotationAxis](rotation-axis.md) | [jvm]<br>fun [rotationAxis](rotation-axis.md)(body: [Body](../-body/index.md), time: [AstroTime](../-astro-time/index.md)): [AxisInfo](../-axis-info/index.md)<br>Calculates information about a body's rotation axis at a given time. |
| [siderealTime](sidereal-time.md) | [jvm]<br>fun [siderealTime](sidereal-time.md)(time: [AstroTime](../-astro-time/index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Calculates Greenwich Apparent Sidereal Time (GAST). |
