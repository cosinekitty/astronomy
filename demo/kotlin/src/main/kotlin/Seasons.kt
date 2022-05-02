import io.github.cosinekitty.astronomy.*

/**
 * Print the times of the equinoxes and solstices for a calendar year.
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
