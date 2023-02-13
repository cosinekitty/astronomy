import io.github.cosinekitty.astronomy.*
import kotlin.math.roundToInt

/**
 * Calculate true solar time for a given observer and UTC time.
 *
 * @param observer
 * The geographic location for which to calculate apparent celestial body positions.
 *
 * @param time
 * The date and time for which to calculate positions.
 */
internal fun `Solar true time`(observer: Observer, time: Time): Int {
    val ha = hourAngle(Body.Sun, time, observer)
    val solarTimeHours = (ha + 12.0) % 24.0
    var milli = (solarTimeHours * 3.6e+6).roundToInt()
    var second = milli / 1000
    milli %= 1000
    var minute = second / 60
    second %= 60
    var hour = minute / 60
    minute %= 60
    hour %= 24

    println("True solar time = %7.4f hours (%02d:%02d:%02d.%03d)".format(
        solarTimeHours, hour, minute, second, milli))

    return 0
}
