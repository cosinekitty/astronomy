//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[Astronomy](index.md)/[jupiterMoons](jupiter-moons.md)

# jupiterMoons

[jvm]\
fun [jupiterMoons](jupiter-moons.md)(time: [AstroTime](../-astro-time/index.md)): [JupiterMoonsInfo](../-jupiter-moons-info/index.md)

Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.

Calculates position and velocity vectors for Jupiter's moons Io, Europa, Ganymede, and Callisto, at the given date and time. The vectors are jovicentric (relative to the center of Jupiter). Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ). The position components are expressed in astronomical units (AU), and the velocity components are in AU/day.

To convert to heliocentric position vectors, call [Astronomy.helioVector](helio-vector.md) with Body.Jupiter to get Jupiter's heliocentric position, then add the jovicentric positions. Likewise, you can call [Astronomy.geoVector](geo-vector.md) to convert to geocentric positions; however, you will have to manually correct for light travel time from the Jupiter system to Earth to figure out what time to pass to jupiterMoons to get an accurate picture of how Jupiter and its moons look from Earth.
