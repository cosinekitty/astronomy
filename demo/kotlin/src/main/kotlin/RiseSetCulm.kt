import kotlin.collections.sorted
import io.github.cosinekitty.astronomy.*

/**
 * Display rise/set/culmination times for the Sun and Moon.
 *
 * @param observer
 * The geographic location of the observer.
 *
 * @param startTime
 * The date and time after which to search for the first
 * occurrences of rise/set/culmination times of the Sun and Moon.
 */
internal fun `Rise Set Culmination demo`(observer: Observer, startTime: Time): Int {
    // At extreme latitudes (close to the poles), there can be long
    // periods of time where the Sun/Moon never rise/set.
    // Allow the search to proceed for 2 days, which will usually find events
    // for most populated latitudes.
    // In other applications (e.g. a daily calendar),
    // the developer may wish to use a 1-day search limit.
    // Note that the Moon rises/sets about 1 hour later each day, so it is
    // common for the Moon to not have a rise or set event in a given day.
    val dayLimit = 2.0
    val eventList = listOfNotNull(
        // Rise/set times may or may not occur within any finite search window.
        maybe("sunrise",  searchRiseSet(Body.Sun,  observer, Direction.Rise, startTime, dayLimit)),
        maybe("sunset",   searchRiseSet(Body.Sun,  observer, Direction.Set,  startTime, dayLimit)),
        maybe("moonrise", searchRiseSet(Body.Moon, observer, Direction.Rise, startTime, dayLimit)),
        maybe("moonset",  searchRiseSet(Body.Moon, observer, Direction.Set,  startTime, dayLimit)),

        // Culmination times can always be found regardless of latitude.
        AstroEvent("sunculm",  searchHourAngle(Body.Sun,  observer, 0.0, startTime)),
        AstroEvent("moonculm", searchHourAngle(Body.Moon, observer, 0.0, startTime)),
    )

    // Print the list in sorted order.
    eventList.sorted().forEach(AstroEvent::display)

    return 0
}


private fun maybe(name: String, time: Time?): AstroEvent? = time?.let { AstroEvent(name, it) }


private class AstroEvent(
    val name: String,
    val time: Time,
    val altitude: Double? = null
): Comparable<AstroEvent> {
    // Allow a list of AstroEvent to be sorted chronologically.
    override operator fun compareTo(other: AstroEvent): Int = this.time.compareTo(other.time)

    constructor(name: String, hourAngleInfo: HourAngleInfo)
        : this(name, hourAngleInfo.time, hourAngleInfo.hor.altitude)

    fun display() {
        print("%-8s : %s".format(name, time))
        if (altitude != null) {
            print("   altitude = %5.2f".format(altitude))
        }
        println()
    }
}
