import io.github.cosinekitty.astronomy.*

/**
 * Find current constellation of the Moon, and its constellation changes over the next 30 days.
 *
 * @param startTime
 * The date and time to start the constellation search.
 */
internal fun `Constellations demo`(startTime: Time): Int {
    val stopTime = startTime.addDays(30.0)

    // There are 12 zodiac constellations, and the moon takes
    // about 27.3 days in its sidereal period. Therefore, there
    // are roughly 2.2 days per constellation. We will sample
    // the Moon's constellation once every 0.1 days to reduce
    // the chance of missing a brief transition through a small
    // part of a constellation.
    val dayIncrement = 0.1

    var t1 = startTime
    var c1 = bodyConstellation(Body.Moon, t1)
    println("$t1 : The Moon starts in ${c1.name}.")

    while (t1 < stopTime) {     // note that Time objects can be compared chronologically
        val t2 = t1.addDays(dayIncrement)
        val c2 = bodyConstellation(Body.Moon, t2)
        if (c1.symbol == c2.symbol) {
            // No constellation change in this time step. Try again on the next time step.
            c1 = c2
            t1 = t2
        } else {
            // The body moved from one constellation to another during this time step.
            // Narrow in on the exact moment of the transition by doing a binary search.
            val change = findConstellationChange(Body.Moon, c1, t1, t2)
            println("${change.time} : The Moon enters ${change.con.name}.")
            c1 = change.con
            t1 = change.time
        }
    }

    return 0
}


private fun bodyConstellation(body: Body, time: Time): ConstellationInfo {
    // Find a vector from the center of the Earth to the center of the body.
    val vec: Vector = geoVector(body, time, Aberration.Corrected)

    // Convert cartesian vector to spherical angular coordinates.
    val equ: Equatorial = vec.toEquatorial()

    // Use the right ascension and declination to find the constellation.
    return constellation(equ.ra, equ.dec)
}


private class ConstellationEvent(
    val time: Time,
    val con: ConstellationInfo
)


private fun findConstellationChange(
    body: Body,
    c1: ConstellationInfo,
    startTime: Time,
    endTime: Time
): ConstellationEvent {
    val tolerance = 0.1 / SECONDS_PER_DAY   // one tenth of a second, expressed in days
    var t1 = startTime
    var t2 = endTime

    // Do a binary search for when the constellation changes.
    while (true) {
        // Calculate the width of the search window.
        val dt = t2.ut - t1.ut

        // Let tx = the time halfway between t1 and t2.
        val tx = t1.addDays(dt / 2.0)

        // What constellation is the body in at tx?
        val cx = bodyConstellation(body, tx)

        // Is it in the same constellation as at time t1, or a different one?
        if (cx.symbol == c1.symbol) {
            // Still in the same constellation.
            // Narrow the search window to [tx, t2].
            t1 = tx
        } else {
            if (dt < tolerance) {
                // We found the time of transition to within 0.1 seconds.
                // Always end the search inside the new constellation.
                return ConstellationEvent(tx, cx)
            }
            // The constellation changed some time in the range [t1, tx].
            // Keep searching.
            t2 = tx
        }
    }
}
