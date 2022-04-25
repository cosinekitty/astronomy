//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Time](index.md)/[tt](tt.md)

# tt

val [tt](tt.md): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)

Terrestrial Time days since noon on January 1, 2000.

Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000. In this system, days are not based on Earth rotations, but instead by the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html) divided by 86400. Unlike ut, tt increases uniformly without adjustments for changes in the Earth's rotation.

The value in tt is used for calculations of movements not involving the Earth's rotation, such as the orbits of planets around the Sun, or the Moon around the Earth.

Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).
