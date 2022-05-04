import java.time.Instant
import java.time.format.DateTimeParseException
import io.github.cosinekitty.astronomy.*
import kotlin.system.exitProcess

private const val usageText = """
In all demos that include [yyyy-mm-ddThh:mm:ssZ]
on the command line, that pattern indicates an optional date/time.
If the date/time is specified, it is used for the calculation.
If absent, the computer's current date and time is used.

Command line arguments:

    constellation [yyyy-mm-ddThh:mm:ssZ]
        Finds what constellation the Moon is in at a given time.
        Then it finds the Moon's constellation changes over the
        subsequent 30 days.

    jupiter_moons [yyyy-mm-ddThh:mm:ssZ]
        Calculates the coordinates of Jupiter and its four major moons
        (Io, Europa, Ganymede, and Callisto) as seen from the Earth
        at a given date and time. This demo illustrates how to correct
        for the delay caused by the time it takes for light to reach
        the Earth from the Jupiter system.

    moonphase [yyyy-mm-ddThh:mm:ssZ]
        Calculates the Moon's ecliptic phase and illumination percentage
        for a given date and time. Also finds the dates and times of
        the subsequent 10 quarter phases of the Moon.

    positions latitude longitude [yyyy-mm-ddThh:mm:ssZ]
        Displays the equatorial and horizontal coordinates of
        the Sun, Moon, and planets, as seen from a given
        geographic location.

    riseset latitude longitude [yyyy-mm-ddThh:mm:ssZ]
        Displays the times of the following 6 events:
        sunrise, sun culmination, sunset,
        moonrise, moon culmination, moonset.
        Culmination is when a body reaches the highest
        point in an observer's sky as it crosses the meridian.
        The specified time is the starting point of the search.
        The displayed events are those that occur first after that time.
        The events are displayed in chronological order.

    seasons year
        Given an integer year number, displays the solstices and
        equinoxes for that year. The year must be in the range 0000..9999.

"""

private fun printUsage(): Int {
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
    exitProcess(runDemo(args))
}

class DemoException(message:String): Exception(message)

private fun runDemo(args: Array<String>): Int {
    if (args.isEmpty()) return printUsage()
    val verb = args[0]
    val demo = demoList.firstOrNull { it.name == verb }
    if (demo == null) {
        println("ERROR: Unknown command '$verb'.")
        return 1
    }
    if (args.size < demo.minArgs || args.size > demo.maxArgs) {
        println("ERROR: Incorrect number of command-line arguments.")
        return 1
    }
    return try {
        demo.func(args)
    } catch (_: DateTimeParseException) {
        println("ERROR: Invalid date/time format on command line.")
        1
    } catch (e: DemoException) {
        println("ERROR: ${e.message}")
        1
    }
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

internal val demoList = listOf(
    Demo("constellation", 1, 2) { args ->
        `Constellations demo`(
            parseTime(args, 1)
        )
    },
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
    Demo("riseset", 3, 4) { args ->
        `Rise Set Culmination demo`(
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
