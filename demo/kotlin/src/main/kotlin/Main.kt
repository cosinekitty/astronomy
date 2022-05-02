import java.time.Instant
import java.time.format.DateTimeParseException
import io.github.cosinekitty.astronomy.*

private val usageText = """
Command line arguments:

    jupiter_moons [yyyy-mm-ddThh:mm:ssZ]
        Calculates the coordinates of Jupiter and its four major moons
        (Io, Europa, Ganymede, and Callisto) as seen from the Earth
        at a given date and time. This demo illustrates how to correct
        for the delay caused by the time it takes for light to reach
        the Earth from the Jupiter system.

    moonphase [yyyy-mm-ddThh:mm:ssZ]
        Calculates the Moon's ecliptic phase and illumination percentage
        for a given date and time, or for the computer's current date and
        time if none is given on the command line.
        Also finds the dates and times of the subsequent 10 quarter phases.

    positions latitude longitude [yyyy-mm-ddThh:mm:ssZ]
        Displays the equatorial and horizontal coordinates of
        the Sun, Moon, and planets, as seen from a given
        geographic location. Uses the date and time specified on
        the command line, if present. Otherwise, uses the computer's
        current date and time.

    seasons year
        Given an integer year number, displays the solstices and equinoxes for that year.
        The year must be in the range 0000..9999.

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

class DemoException(message:String): Exception(message)

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
                } catch (_: DateTimeParseException) {
                    println("ERROR: Invalid date/time format on command line.")
                    return 1
                } catch (e: DemoException) {
                    println("ERROR: ${e.message}")
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

internal fun parseNumber(name: String, text: String, minValue: Double, maxValue: Double): Double {
    try {
        val value = text.toDouble()
        if (!value.isFinite() || value < minValue || value > maxValue) {
            throw DemoException("Value for $name is out of range.")
        }
        return value
    } catch (_: NumberFormatException) {
        throw DemoException("Invalid numeric format '$text' for $name.")
    }
}

internal fun parseYear(text: String): Int {
    try {
        val year = text.toInt()
        if (year < 0 || year > 9999) {
            throw DemoException("Year must be in the range 0000..9999.")
        }
        return year
    } catch (_: NumberFormatException) {
        throw DemoException("Invalid numeric format '$text' for year.")
    }
}

internal fun parseObserver(args: Array<String>, index: Int): Observer {
    val latitude = parseNumber("latitude", args[index], -90.0, +90.0)
    val longitude = parseNumber("longitude", args[index+1], -180.0, +180.0)
    return Observer(latitude, longitude, 0.0)
}

internal val demoList = arrayOf(
    Demo("jupiter_moons", 1, 2) { args ->
        `Jupiter moons demo`(
            parseTime(args, 1)
        )
    },
    Demo("moonphase", 1, 2) { args ->
        `Moon Phase demo`(
            parseTime(args, 1)
        )
    },
    Demo("positions", 3, 4) { args ->
        `Celestial body positions demo`(
            parseObserver(args, 1),
            parseTime(args, 3)
        )
    },
    Demo("seasons", 2, 2)  { args ->
        `Seasons demo`(
            parseYear(args[1])
        )
    }
)
