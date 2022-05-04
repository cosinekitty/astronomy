import io.github.cosinekitty.astronomy.*

/**
 * Print positions of Jupiter's moons as seen from Earth.
 *
 * Calculates the coordinates of Jupiter and its four major moons
 * (Io, Europa, Ganymede, and Callisto) as seen from the Earth
 * at a given date and time. This function illustrates how to correct
 * for the delay caused by the time it takes for light to reach
 * the Earth from the Jupiter system.
 *
 * Without light travel time correction, observations with a telescope
 * will differ from the calculations by plus or minus 8 minutes over
 * the span of a year, because of the amount of time it takes for light
 * to travel across the diameter of the Earth's orbit.
 *
 * @param time
 * The date and time of the observation.
 */
internal fun `Jupiter moons demo`(time: Time): Int {
    // Call geoVector to calculate the geocentric position of Jupiter.
    // geoVector corrects for light travel time.
    // That means it returns a vector to where Jupiter appears to be
    // in the sky, when the light left Jupiter to travel toward the
    // Earth to arrive here at the specified time. This is different from
    // where Jupiter is at that time.

    println("Calculations for: $time")

    val jv = geoVector(Body.Jupiter, time, Aberration.Corrected)

    // Calculate the amount of time it took light to reach the Earth from Jupiter.
    // The distance to Jupiter (AU) divided by the speed of light (AU/day) = time in days.
    val lightTravelDays = jv.length() / C_AUDAY
    println()
    println("It took light %.2f minutes to reach the Earth from Jupiter.".format(lightTravelDays * MINUTES_PER_DAY))
    println()

    // The jupiterMoons function calculates positions of Jupiter's moons without
    // correcting for light travel time. Correct for light travel by backdating
    // by the given amount of light travel time.
    val backdate = time.addDays(-lightTravelDays)

    val jm = jupiterMoons(backdate)

    // Tricky: the `+` operator for adding `Vector` will throw an exception
    // if the vectors do not have matching times. We work around this
    // by using `withTime` to clone each moon's position vector to have
    // a different time. This is a manual override to work around a safety check.

    printBody("Jupiter",  jv)
    printBody("Io",       jv + jm.io.position().withTime(jv.t))
    printBody("Europa",   jv + jm.europa.position().withTime(jv.t))
    printBody("Ganymede", jv + jm.ganymede.position().withTime(jv.t))
    printBody("Callisto", jv + jm.callisto.position().withTime(jv.t))
    println()

    return 0
}


private fun printBody(name: String, geovec: Vector) {
    // Convert the geocentric vector into equatorial coordinates.
    val equ = geovec.toEquatorial()
    println("%-8s   RA %10.6f   DEC %10.6f  %10.6f AU".format(name, equ.ra, equ.dec, equ.dist))
}

