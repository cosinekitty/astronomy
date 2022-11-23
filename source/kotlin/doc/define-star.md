//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[defineStar](define-star.md)

# defineStar

fun [defineStar](define-star.md)(body: [Body](-body/index.md), ra: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), dec: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html), distanceLightYears: [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html))

Assign equatorial coordinates to a user-defined star.

Some Astronomy Engine functions allow their body parameter to be a user-defined fixed point in the sky, loosely called a "star". This function assigns a right ascension, declination, and distance to one of the eight user-defined stars Body.Star1..Body.Star8.

Stars are not valid until defined. Once defined, they retain their definition until re-defined by another call to defineStar.

## Parameters

| | |
|---|---|
| body | One of the eight user-defined star identifiers: Body.Star1, Body.Star2, ..., Body.Star8. |
| ra | The right ascension to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ). The value is in units of sidereal hours, and must be within the half-open range [0, 24). |
| dec | The right ascension to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ). The value is in units of degrees north (positive) or south (negative) of the J2000 equator, and must be within the closed range -90, +90. |
| distanceLightYears | The distance between the star and the Sun, expressed in light-years. This value is used to calculate the tiny parallax shift as seen by an observer on Earth. If you don't know the distance to the star, using a large value like 1000 will generally work well. The minimum allowed distance is 1 light-year, which is required to provide certain internal optimizations. |
