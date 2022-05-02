import io.github.cosinekitty.astronomy.*

private val bodyList = arrayOf(
    Body.Sun, Body.Moon, Body.Mercury, Body.Venus, Body.Mars,
    Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto
)

/**
 * Print a table of solar system body positions in equatorial and horizontal coordinates.
 *
 * @param observer
 * The geographic location for which to calculate apparent celestial body positions.
 *
 * @param time
 * The date and time for which to calculate positions.
 */
internal fun `Celestial body positions demo`(observer: Observer, time: Time): Int {
    println("UTC date = $time")
    println()
    println("BODY           RA      DEC       AZ      ALT")
    for (body in bodyList) {
        val equ_2000: Equatorial = equator(body, time, observer, EquatorEpoch.J2000, Aberration.Corrected)
        val equ_ofdate: Equatorial = equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected)
        val hor: Topocentric = horizon(time, observer, equ_ofdate.ra, equ_ofdate.dec, Refraction.Normal)
        println("%-8s %8.2f %8.2f %8.2f %8.2f".format(body, equ_2000.ra, equ_2000.dec, hor.azimuth, hor.altitude))
    }
    return 0
}
