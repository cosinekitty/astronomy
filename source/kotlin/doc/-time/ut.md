//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Time](index.md)/[ut](ut.md)

# ut

val [ut](ut.md): Double

UT1/UTC number of days since noon on January 1, 2000.

The floating point number of days of Universal Time since noon UTC January 1, 2000. Astronomy Engine approximates UTC and UT1 as being the same thing, although they are not exactly equivalent; UTC and UT1 can disagree by up to plus or minus 0.9 seconds. This approximation is sufficient for the accuracy requirements of Astronomy Engine.

Universal Time Coordinate (UTC) is the international standard for legal and civil timekeeping and replaces the older Greenwich Mean Time (GMT) standard. UTC is kept in sync with unpredictable observed changes in the Earth's rotation by occasionally adding leap seconds as needed.

UT1 is an idealized time scale based on observed rotation of the Earth, which gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun, large scale weather events like hurricanes, and internal seismic and convection effects. Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC is adjusted by a scheduled whole number of leap seconds as needed.

The value in ut is appropriate for any calculation involving the Earth's rotation, such as calculating rise/set times, culumination, and anything involving apparent sidereal time.

Before the era of atomic timekeeping, days based on the Earth's rotation were often known as *mean solar days*.
