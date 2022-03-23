package io.github.cosinekitty.astronomy

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

import java.text.SimpleDateFormat
import java.util.*
import kotlin.math.absoluteValue
import kotlin.math.roundToLong
import kotlin.math.sqrt


private val TimeZoneUtc = TimeZone.getTimeZone("UTC")
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
     * @param second The UTC second in the half-open range [0, 60).
     */
    constructor(year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Double) : this(
            GregorianCalendar(TimeZoneUtc).also {
                val totalMillis: Int = Math.floor(second * 1000.0).toInt()
                it.set(year, month - 1, day, hour, minute, totalMillis / 1000)
                it.set(Calendar.MILLISECOND, totalMillis % 1000)
            }
    )

    /**
     * Converts this object to a native `Date` equivalent.
     *
     * @return a UTC `Date` object for this `AstroTime` value.
     */
    fun toDate(): Date = Date(origin.time + (ut * MILLIS_PER_DAY).roundToLong())

    /**
     * Converts this `AstroTime` to ISO 8601 format, expressed in UTC with millisecond resolution.
     *
     * @return Example: "2019-08-30T17:45:22.763Z".
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
    fun addDays(days: Double) = AstroTime(ut + days)

    companion object {
        private val origin = GregorianCalendar(TimeZoneUtc).also {
            it.set(2000, 0, 1, 12, 0, 0)
            it.set(Calendar.MILLISECOND, 0)
        }.time

        private const val MILLIS_PER_DAY = 24 * 3600 * 1000

        private val dateFormat = SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS'Z'").also {
            it.timeZone = TimeZoneUtc
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


internal data class TerseVector(val x: Double, val y: Double, val z: Double) {
    fun toAstroVector(time: AstroTime) =
            AstroVector(x, y, z, time)

    operator fun plus(other: TerseVector) =
            TerseVector(x + other.x, y + other.y, z + other.z)

    operator fun minus(other: TerseVector) =
            TerseVector(x - other.x, y - other.y, z - other.z)

    operator fun times(other: Double) =
            TerseVector(x * other, y * other, z * other)

    operator fun div(other: Double) =
            TerseVector(x / other, y / other, z / other)

    val quadrature get() = x * x + y * y + z * z
    val magnitude get() = sqrt(quadrature)

    companion object {
        val zero = TerseVector(0.0, 0.0, 0.0)
    }
}

data class AstroVector(
    /**
     * A Cartesian x-coordinate expressed in AU.
     */
    val x: Double,

    /**
     * A Cartesian y-coordinate expressed in AU.
     */
    val y: Double,

    /**
     * A Cartesian z-coordinate expressed in AU.
     */
    val z: Double,

    /**
     * The date and time at which this vector is valid.
     */
    val t: AstroTime
) {
    /**
     * The total distance in AU represented by this vector.
     */
    fun length() =
        sqrt(x*x + y*y + z*z)

    private fun verifyIdenticalTimes(otherTime: AstroTime): AstroTime {
        if (t.tt != otherTime.tt)
            throw IllegalArgumentException("Attempt to operate on two vectors from different times.")
        return t
    }

    operator fun plus(other: AstroVector) =
            AstroVector(x + other.x, y + other.y, z + other.z, verifyIdenticalTimes(other.t))

    operator fun minus(other: AstroVector) =
            AstroVector(x - other.x, y - other.y, z - other.z, verifyIdenticalTimes(other.t))

    operator fun unaryMinus() =
            AstroVector(-x, -y, -z, t)

    operator fun times(other: AstroVector): Double {    // scalar dot product
        verifyIdenticalTimes(other.t)
        return x*other.x + y*other.y + z*other.z
    }

    operator fun div(denom: Double) =
            AstroVector(x/denom, y/denom, z/denom, t)
}

operator fun Double.times(vec: AstroVector) =
        AstroVector(this*vec.x, this*vec.y, this*vec.z, vec.t)


data class StateVector(
    /**
     * A Cartesian position x-coordinate expressed in AU.
     */
    val x: Double,

    /**
     * A Cartesian position y-coordinate expressed in AU.
     */
    val y: Double,

    /**
     * A Cartesian position z-coordinate expressed in AU.
     */
    val z: Double,

    /**
     * A Cartesian velocity x-component expressed in AU/day.
     */
    val vx: Double,

    /**
     * A Cartesian velocity y-component expressed in AU/day.
     */
    val vy: Double,

    /**
     * A Cartesian velocity z-component expressed in AU/day.
     */
    val vz: Double,

    /**
     * The date and time at which this vector is valid.
     */
    val t: AstroTime
) {

    /**
     * Combines a position vector and a velocity vector into a single state vector.
     *
     * @param pos   A position vector.
     * @param vel   A velocity vector.
     * @param time  The common time that represents the given position and velocity.
     */
    constructor(pos: AstroVector, vel: AstroVector, time: AstroTime)
        : this(pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, time)

    /**
     * Returns the position vector associated with this state vector.
     */
    val position get() = AstroVector(x, y, z, t)

    /**
     * Returns the velocity vector associated with this state vector.
     */
    val velocity get() = AstroVector(vx, vy, vz, t)
}


/**
 * Holds the positions and velocities of Jupiter's major 4 moons.
 *
 * The #Astronomy.JupiterMoons function returns an object of this type
 * to report position and velocity vectors for Jupiter's largest 4 moons
 * Io, Europa, Ganymede, and Callisto. Each position vector is relative
 * to the center of Jupiter. Both position and velocity are oriented in
 * the EQJ system (that is, using Earth's equator at the J2000 epoch).
 * The positions are expressed in astronomical units (AU),
 * and the velocities in AU/day.
 */
class JupiterMoonsInfo(
    /**
     * An array of state vectors for each of the 4 moons, in the following order:
     * 0 = Io, 1 = Europa, 2 = Ganymede, 3 = Callisto.
     */
    val moon: Array<StateVector>
)


/**
 * A rotation matrix that can be used to transform one coordinate system to another.
 */
class RotationMatrix(
    /**
     * A 3x3 array of numbers to initialize the rotation matrix.
     */
    val rot: Array<DoubleArray>
) {
    init {
        if (rot.size != 3 || rot[0].size != 3 || rot[1].size != 3 || rot[2].size != 3)
            throw IllegalArgumentException("Rotation matrix must be a 3x3 array.")
    }
}


/**
 * Spherical coordinates: latitude, longitude, distance.
 */
data class Spherical(
    /**
     * The latitude angle: -90..+90 degrees.
     */
    val lat: Double,

    /**
     * The longitude angle: 0..360 degrees.
     */
    val lon: Double,

    /**
     * Distance in AU.
     */
    val dist: Double
)

/**
 * The location of an observer on (or near) the surface of the Earth.
 *
 * This object is passed to functions that calculate phenomena as observed
 * from a particular place on the Earth.
 */
data class Observer(
    /**
     * Geographic latitude in degrees north (positive) or south (negative) of the equator.
     */
    val latitude: Double,

    /**
     * Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.
     */
    val longitude: Double,

    /**
     * The height above (positive) or below (negative) sea level, expressed in meters.
     */
    val height: Double
)


/**
 * Selects the date for which the Earth's equator is to be used for representing equatorial coordinates.
 *
 * The Earth's equator is not always in the same plane due to precession and nutation.
 *
 * Sometimes it is useful to have a fixed plane of reference for equatorial coordinates
 * across different calendar dates.  In these cases, a fixed *epoch*, or reference time,
 * is helpful. Astronomy Engine provides the J2000 epoch for such cases.  This refers
 * to the plane of the Earth's orbit as it was on noon UTC on 1 January 2000.
 *
 * For some other purposes, it is more helpful to represent coordinates using the Earth's
 * equator exactly as it is on that date. For example, when calculating rise/set times
 * or horizontal coordinates, it is most accurate to use the orientation of the Earth's
 * equator at that same date and time. For these uses, Astronomy Engine allows *of-date*
 * calculations.
 */
enum class EquatorEpoch {
    /**
     * Represent equatorial coordinates in the J2000 epoch.
     */
    J2000,

    /**
     * Represent equatorial coordinates using the Earth's equator at the given date and time.
     */
    OfDate,
}


/**
 * Aberration calculation options.
 *
 * [Aberration](https://en.wikipedia.org/wiki/Aberration_of_light) is an effect
 * causing the apparent direction of an observed body to be shifted due to transverse
 * movement of the Earth with respect to the rays of light coming from that body.
 * This angular correction can be anywhere from 0 to about 20 arcseconds,
 * depending on the position of the observed body relative to the instantaneous
 * velocity vector of the Earth.
 *
 * Some Astronomy Engine functions allow optional correction for aberration by
 * passing in a value of this enumerated type.
 *
 * Aberration correction is useful to improve accuracy of coordinates of
 * apparent locations of bodies seen from the Earth.
 * However, because aberration affects not only the observed body (such as a planet)
 * but the surrounding stars, aberration may be unhelpful (for example)
 * for determining exactly when a planet crosses from one constellation to another.
 */
enum class Aberration {
    /**
     * Request correction for aberration.
     */
    Corrected,

    /**
     * Do not correct for aberration.
     */
    None,
}


/**
 * Selects whether to correct for atmospheric refraction, and if so, how.
 */
enum class Refraction {
    /**
     * No atmospheric refraction correction (airless).
     */
    None,

    /**
     * Recommended correction for standard atmospheric refraction.
     */
    Normal,

    /**
     * Used only for compatibility testing with JPL Horizons online tool.
     */
    JplHor,
}


/**
 * Selects whether to search for a rising event or a setting event for a celestial body.
 */
enum class Direction {
    /**
     * Indicates a rising event: a celestial body is observed to rise above the horizon by an observer on the Earth.
     */
    Rise,

    /**
     * Indicates a setting event: a celestial body is observed to sink below the horizon by an observer on the Earth.
     */
    Set,
}


/**
 * Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.
 */
enum class Visibility {
    /**
     * The body is best visible in the morning, before sunrise.
     */
    Morning,

    /**
     * The body is best visible in the evening, after sunset.
     */
    Evening,
}


/**
 * Equatorial angular and cartesian coordinates.
 *
 * Coordinates of a celestial body as seen from the Earth
 * (geocentric or topocentric, depending on context),
 * oriented with respect to the projection of the Earth's equator onto the sky.
 */
class Equatorial(
    /**
     * Right ascension in sidereal hours.
     */
    val ra: Double,

    /**
     * Declination in degrees.
     */
    val dec: Double,

    /**
     * Distance to the celestial body in AU.
     */
    val dist: Double,

    /**
     * Equatorial coordinates in cartesian vector form: x = March equinox, y = June solstice, z = north.
     */
    val vec: AstroVector
)


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
