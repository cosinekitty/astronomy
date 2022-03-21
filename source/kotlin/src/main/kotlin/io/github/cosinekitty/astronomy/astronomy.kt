package io.github.cosinekitty.astronomy

import java.text.SimpleDateFormat
import java.util.*
import kotlin.math.absoluteValue
import kotlin.math.roundToLong
import kotlin.math.sqrt

/*
    Astronomy Engine for Kotlin / JVM.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019-2022 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/


private const val DAYS_PER_TROPICAL_YEAR = 365.24217


/**
 * The enumeration of celestial bodies supported by Astronomy Engine.
 */
enum class Body {
    /**
     * A placeholder value representing an invalid or unknown celestial body.
     */
    Invalid,

    /**
     * The planet Mercury.
     */
    Mercury,

    /**
     * The planet Venus.
     */
    Venus,

    /**
     * The planet Earth.
     * Some functions that accept a `Body` parameter will fail if passed this value
     * because they assume that an observation is being made from the Earth,
     * and therefore the Earth is not a target of observation.
     */
    Earth,

    /**
     * The planet Mars.
     */
    Mars,

    /**
     * The planet Jupiter.
     */
    Jupiter,

    /**
     * The planet Saturn.
     */
    Saturn,

    /**
     * The planet Uranus.
     */
    Uranus,

    /**
     * The planet Neptune.
     */
    Neptune,

    /**
     * The planet Pluto.
     */
    Pluto,

    /**
     * The Sun.
     */
    Sun,

    /**
     * The Earth's natural satellite, the Moon.
     */
    Moon,

    /**
     * The Earth/Moon Barycenter.
     */
    EMB,

    /**
     * The Solar System Barycenter.
     */
    SSB,
}


/**
 * A date and time used for astronomical calculations.
 */
class AstroTime private constructor(
    /**
     * UT1/UTC number of days since noon on January 1, 2000.
     *
     * The floating point number of days of Universal Time since noon UTC January 1, 2000.
     * Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
     * not exactly equivalent; UTC and UT1 can disagree by up to plus or minus 0.9 seconds.
     * This approximation is sufficient for the accuracy requirements of Astronomy Engine.
     *
     * Universal Time Coordinate (UTC) is the international standard for legal and civil
     * timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
     * UTC is kept in sync with unpredictable observed changes in the Earth's rotation
     * by occasionally adding leap seconds as needed.
     *
     * UT1 is an idealized time scale based on observed rotation of the Earth, which
     * gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
     * large scale weather events like hurricanes, and internal seismic and convection effects.
     * Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
     * is adjusted by a scheduled whole number of leap seconds as needed.
     *
     * The value in `ut` is appropriate for any calculation involving the Earth's rotation,
     * such as calculating rise/set times, culumination, and anything involving apparent
     * sidereal time.
     *
     * Before the era of atomic timekeeping, days based on the Earth's rotation
     * were often known as *mean solar days*.
     */
    val ut: Double,

    /**
     * Terrestrial Time days since noon on January 1, 2000.
     *
     * Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
     * In this system, days are not based on Earth rotations, but instead by
     * the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html)
     * divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments
     * for changes in the Earth's rotation.
     *
     * The value in `tt` is used for calculations of movements not involving the Earth's rotation,
     * such as the orbits of planets around the Sun, or the Moon around the Earth.
     *
     * Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).
     */
    val tt: Double
) {
    /*
     * For internal use only. Used to optimize Earth tilt calculations.
     */
    internal var psi = Double.NaN

    /*
     * For internal use only. Used to optimize Earth tilt calculations.
     */
    internal var eps = Double.NaN

    /*
     * For internal use only. Lazy-caches sidereal time (Earth rotation).
     */
    internal var st = Double.NaN

    constructor(ut: Double) : this(ut, Astronomy.terrestrialTime(ut))

    /**
     * Creates an `AstroTime` object from a `Date` object.
     *
     * @param d The date and time to be converted to AstroTime format.
     */
    constructor(d: Date) : this((d.time - origin.time).toDouble() / MILLIS_PER_DAY)

    /**
     * Creates an `AstroTime` object from a `Calendar` object.
     *
     * @param d The date and time to be converted to AstroTime format.
     */
    constructor(d: Calendar) : this(d.time)

    /**
     * Creates an `AstroTime` object from a UTC year, month, day, hour, minute and second.
     *
     * @param year The UTC year value.
     * @param month The UTC month value 1..12.
     * @param day The UTC day of the month 1..31.
     * @param hour The UTC hour value 0..23.
     * @param minute The UTC minute value 0..59.
     * @param second The UTC second value 0..59.
     */
    constructor(year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Int) : this(
        GregorianCalendar(TimeZone.getTimeZone("UTC")).also {
            it.set(year, month - 1, day, hour, minute, second)
            it.set(Calendar.MILLISECOND, 0)
        }
    )

    /**
     * Converts this object to Java `Date` object.
     *
     * @return a UTC `Date` object for this `AstroTime` value.
     */
    fun toDate(): Date = Date(origin.time + (ut * MILLIS_PER_DAY).roundToLong())

    /**
     * Converts this `AstroTime` to ISO 8601 format, expressed in UTC with millisecond resolution.
     *
     * @return Example: "2019-08-30T17:45:22.763".
     */
    override fun toString(): String = dateFormat.format(toDate())

    /**
     * Calculates the sum or difference of an #AstroTime with a specified floating point number of days.
     *
     * Sometimes we need to adjust a given #AstroTime value by a certain amount of time.
     * This function adds the given real number of days in `days` to the date and time in this object.
     *
     * More precisely, the result's Universal Time field `ut` is exactly adjusted by `days` and
     * the Terrestrial Time field `tt` is adjusted for the resulting UTC date and time,
     * using a best-fit piecewise polynomial model devised by
     * [Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).
     *
     * @param days A floating point number of days by which to adjust `time`. May be negative, 0, or positive.
     * @return A date and time that is conceptually equal to `time + days`.
     */
    fun addDays(days: Double): AstroTime = AstroTime(ut + days)

    companion object {
        private val origin = GregorianCalendar(TimeZone.getTimeZone("UTC")).also {
            it.set(2000, 0, 1, 12, 0, 0)
            it.set(Calendar.MILLISECOND, 0)
        }.time

        private const val MILLIS_PER_DAY = 24 * 3600 * 1000

        private val dateFormat = SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.sss'Z'").also {
            it.timeZone = TimeZone.getTimeZone("UTC")
        }

        /**
         * Creates an `AstroTime` object from a Terrestrial Time day value.
         *
         * This function can be used in rare cases where a time must be based
         * on Terrestrial Time (TT) rather than Universal Time (UT).
         * Most developers will want to invoke `new AstroTime(ut)` with a universal time
         * instead of this function, because usually time is based on civil time adjusted
         * by leap seconds to match the Earth's rotation, rather than the uniformly
         * flowing TT used to calculate solar system dynamics. In rare cases
         * where the caller already knows TT, this function is provided to create
         * an `AstroTime` value that can be passed to Astronomy Engine functions.
         *
         * @param tt The number of days after the J2000 epoch.
         */
        fun fromTerrestrialTime(tt: Double): AstroTime = AstroTime(Astronomy.universalTime(tt), tt)
    }
}


object Astronomy {
    private const val SECONDS_PER_DAY = 24 * 3600
    internal fun terrestrialTime(ut: Double): Double = ut + deltaT(ut) / SECONDS_PER_DAY

    internal fun deltaT(ut: Double): Double {
        /*
            Fred Espenak writes about Delta-T generically here:
            https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
            https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html

            He provides polynomial approximations for distant years here:
            https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html

            They start with a year value 'y' such that y=2000 corresponds
            to the UTC Date 15-January-2000. Convert difference in days
            to mean tropical years.
        */
        val u: Double
        val u2: Double
        val u3: Double
        val u4: Double
        val u5: Double
        val u6: Double
        val u7: Double
        val y = 2000 + (ut - 14) / DAYS_PER_TROPICAL_YEAR
        return when {
            y < -500 -> {
                u = (y - 1820) / 100
                -20.0 + (32.0 * u * u)
            }
            y < 500.0 -> {
                u = y / 100
                u2 = u * u; u3 = u * u2; u4 = u2 * u2; u5 = u2 * u3; u6 = u3 * u3
                10583.6 - 1014.41 * u + 33.78311 * u2 - 5.952053 * u3 - 0.1798452 * u4 + 0.022174192 * u5 + 0.0090316521 * u6
            }
            y < 1600.0 -> {
                u = (y - 1000) / 100
                u2 = u * u; u3 = u * u2; u4 = u2 * u2; u5 = u2 * u3; u6 = u3 * u3
                1574.2 - 556.01 * u + 71.23472 * u2 + 0.319781 * u3 - 0.8503463 * u4 - 0.005050998 * u5 + 0.0083572073 * u6
            }
            y < 1700.0 -> {
                u = y - 1600
                u2 = u * u; u3 = u * u2
                120.0 - 0.9808 * u - 0.01532 * u2 + u3 / 7129
            }
            y < 1800.0 -> {
                u = y - 1700
                u2 = u * u; u3 = u * u2; u4 = u2 * u2
                8.83 + 0.1603 * u - 0.0059285 * u2 + 0.00013336 * u3 - u4 / 1174000
            }
            y < 1860.0 -> {
                u = y - 1800
                u2 = u * u; u3 = u * u2; u4 = u2 * u2; u5 = u2 * u3; u6 = u3 * u3; u7 = u3 * u4
                13.72 - 0.332447 * u + 0.0068612 * u2 + 0.0041116 * u3 - 0.00037436 * u4 + 0.0000121272 * u5 - 0.0000001699 * u6 + 0.000000000875 * u7
            }
            y < 1900.0 -> {
                u = y - 1860
                u2 = u * u; u3 = u * u2; u4 = u2 * u2; u5 = u2 * u3
                7.62 + 0.5737 * u - 0.251754 * u2 + 0.01680668 * u3 - 0.0004473624 * u4 + u5 / 233174
            }
            y < 1920.0 -> {
                u = y - 1900
                u2 = u * u; u3 = u * u2; u4 = u2 * u2
                -2.79 + 1.494119 * u - 0.0598939 * u2 + 0.0061966 * u3 - 0.000197 * u4
            }
            y < 1941.0 -> {
                u = y - 1920
                u2 = u * u; u3 = u * u2
                21.20 + 0.84493 * u - 0.076100 * u2 + 0.0020936 * u3
            }
            y < 1961 -> {
                u = y - 1950
                u2 = u * u; u3 = u * u2
                29.07 + 0.407 * u - u2 / 233 + u3 / 2547
            }
            y < 1986.0 -> {
                u = y - 1975
                u2 = u * u; u3 = u * u2
                45.45 + 1.067 * u - u2 / 260 - u3 / 718
            }
            y < 2005 -> {
                u = y - 2000
                u2 = u * u; u3 = u * u2; u4 = u2 * u2; u5 = u2 * u3
                63.86 + 0.3345 * u - 0.060374 * u2 + 0.0017275 * u3 + 0.000651814 * u4 + 0.00002373599 * u5
            }
            y < 2050 -> {
                u = y - 2000
                62.92 + 0.32217 * u + 0.005589 * u * u
            }
            y < 2150 -> {
                u = (y - 1820) / 100
                -20 + 32 * u * u - 0.5628 * (2150 - y)
            }
            /* all years after 2150 */
            else -> {
                u = (y - 1820) / 100
                -20 + (32 * u * u)
            }
        }
    }

    internal fun universalTime(tt: Double): Double {
        // This is the inverse function of TerrestrialTime.
        // This is an iterative numerical solver, but because
        // the relationship between UT and TT is almost perfectly linear,
        // it converges extremely fast (never more than 3 iterations).

        // dt = tt - ut
        var dt = terrestrialTime(tt) - tt
        while (true) {
            val ut = tt - dt
            val ttCheck = terrestrialTime(ut)
            val err = ttCheck - tt
            if (err.absoluteValue < 1.0e-12) return ut
            dt += err
        }
    }
}

internal data class TerseVector(val x: Double, val y: Double, val z: Double) {
    // fun toAstroVector(time: AstroTime): AstroVector = AstroVector(x, y, z, time)

    operator fun plus(other: TerseVector): TerseVector =
        TerseVector(x + other.x, y + other.y, z + other.z)

    operator fun rem(other: TerseVector): TerseVector =
        TerseVector(x - other.x, y - other.y, z - other.z)

    operator fun times(other: Double): TerseVector =
        TerseVector(x * other, y * other, z * other)

    operator fun div(other: Double): TerseVector =
        TerseVector(x / other, y / other, z / other)

    val quadrature get() = x * x + y * y + z * z
    val magnitude get() = sqrt(quadrature)

    companion object {
        val zero = TerseVector(0.0, 0.0, 0.0)
    }
}
