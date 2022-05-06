import io.github.cosinekitty.astronomy.*

/**
 * Searches for the next 10 partial/total lunar eclipses.
 *
 * @param startTime
 * The date and time after which to start searching for lunar eclipses.
 */
internal fun `Lunar eclipse demo`(startTime: Time): Int {
    generateSequence(searchLunarEclipse(startTime)) { nextLunarEclipse(it.peak) }
        .filter { it.kind != EclipseKind.Penumbral }
        .take(10)
        .forEach(::printEclipse)
    return 0
}

private fun printEclipse(e: LunarEclipseInfo) {
    // Calculate beginning/ending of different phases
    // of an eclipse by subtracting/adding the peak time
    // with the number of minutes indicated by the "semi-duration"
    // fields sdPartial and sdTotal.

    val p1 = e.peak.addDays(-e.sdPartial / MINUTES_PER_DAY)
    println("$p1  Partial eclipse begins.")
    if (e.sdTotal > 0.0) {
        val t1 = e.peak.addDays(-e.sdTotal / MINUTES_PER_DAY)
        println("$t1  Total eclipse begins.")
    }
    println("${e.peak}  Peak of ${e.kind.toString().lowercase()} eclipse.")
    if (e.sdTotal > 0.0) {
        val t2 = e.peak.addDays(+e.sdTotal / MINUTES_PER_DAY)
        println("$t2  Total eclipse ends.")
    }
    val p2 = e.peak.addDays(+e.sdPartial / MINUTES_PER_DAY)
    println("$p2  Partial eclipse ends.")
    println()
}
