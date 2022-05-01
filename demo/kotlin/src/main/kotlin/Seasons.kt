import io.github.cosinekitty.astronomy.*

/**
 * Demonstration of calculating the equinoxes and solstices for calendar year.
 *
 * @param year
 * The calendar year value for which to find equinoxes and solstices.
 */
internal fun `Seasons demo`(year: Int): Int {
    val s = seasons(year)
    println("March equinox     : ${s.marchEquinox}")
    println("June solstice     : ${s.juneSolstice}")
    println("September equinox : ${s.septemberEquinox}")
    println("December solstice : ${s.decemberSolstice}")
    return 0
}
