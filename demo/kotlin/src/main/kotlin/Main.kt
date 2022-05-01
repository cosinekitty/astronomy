import java.time.Instant
import java.time.format.DateTimeParseException
import io.github.cosinekitty.astronomy.*

private val usageText = """
Command line arguments:

    moonphase [yyyy-mm-ddThh:mm:ssZ]
        Calculates the Moon's ecliptic phase and illumination percentage
        for a given date and time, or for the computer's current date and
        time if none is given on the command line.
        Also finds the dates and times of the subsequent 10 quarter phases.

    seasons year
        Given an integer year number, display the solstices and equinoxes for that year.

"""

internal fun printUsage(): Int {
    println(usageText)
    return 1
}

internal class Demo(
    val name: String,
    val minArgs: Int,
    val maxArgs: Int,
    val func: (Array<String>) -> Int
)

fun main(args: Array<String>) {
    System.exit(runDemo(args))
}

internal fun runDemo(args: Array<String>): Int {
    if (args.size > 0) {
        val verb = args[0];
        for (demo in demoList) {
            if (demo.name == verb) {
                if (args.size < demo.minArgs || args.size > demo.maxArgs) {
                    println("ERROR: Incorrect number of command-line arguments.")
                    return 1
                }
                try {
                    return demo.func(args)
                } catch (e: DateTimeParseException) {
                    println("ERROR: Invalid date/time format on command line.")
                    return 1
                }
            }
        }
    }
    return printUsage()
}

internal fun parseTime(args: Array<String>, index: Int): Time {
    val millis = (
        if (index >= 0 && index < args.size)
            Instant.parse(args[index]).toEpochMilli()
        else
            System.currentTimeMillis()
    )
    return Time.fromMillisecondsSince1970(millis)
}

internal val demoList = arrayOf(
    Demo("moonphase", 1, 2) { args -> demoMoonPhase(parseTime(args, 1)) },
    Demo("seasons", 2, 2)  { args -> demoSeasons(args[1].toInt()) }
)
