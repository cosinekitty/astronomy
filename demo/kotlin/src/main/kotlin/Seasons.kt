import io.github.cosinekitty.astronomy.*

internal fun demoSeasons(year: Int): Int {
    val s = seasons(year)
    println("March equinox     : ${s.marchEquinox}")
    println("June solstice     : ${s.juneSolstice}")
    println("September equinox : ${s.septemberEquinox}")
    println("December solstice : ${s.decemberSolstice}")
    return 0
}
