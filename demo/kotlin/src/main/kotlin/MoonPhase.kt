import io.github.cosinekitty.astronomy.*

/**
 * Display Moon phase information and the next 10 quarter phases.
 *
 * @param time
 * The date and time for which to display Moon phase information,
 * and when to initiate a search for the subsequent quarter phases.
 */
internal fun `Moon Phase demo`(time: Time): Int {
    // Calculate the Moon's ecliptic phase angle,
    // which ranges from 0 to 360 degrees.
    //   0 degrees = new mooon,
    //  90 degrees = first quarter,
    // 180 degrees = full moon,
    // 270 degrees = third quarter
    val phase = moonPhase(time)
    println("$time : Moon's ecliptic phase angle = %.3f degrees.".format(phase))

    // Calculate the fraction of the Moon's disc that appears
    // illuminated, as seen from the Earth.
    val illum = illumination(Body.Moon, time)
    println("$time : Moon's illuminated fraction = %.2f%%.".format(100.0 * illum.phaseFraction))

    // Predict when the next 10 lunar quarter phases will happen.
    println()
    println("The next 10 lunar quarters are:")
    var mq = searchMoonQuarter(time)
    for (i in 0..9) {
        if (i > 0) {
            mq = nextMoonQuarter(mq)
        }
        println("${mq.time} : ${quarterName(mq.quarter)}")
    }

    return 0
}

private fun quarterName(quarter: Int): String =
    when (quarter) {
        0 -> "New Moon"
        1 -> "First Quarter"
        2 -> "Full Moon"
        3 -> "Third Quarter"
        else -> "INVALID QUARTER"
    }
