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
@file:JvmName("Astronomy")

package io.github.cosinekitty.astronomy

import java.text.SimpleDateFormat
import java.util.*

import kotlin.math.absoluteValue
import kotlin.math.abs
import kotlin.math.acos
import kotlin.math.atan2
import kotlin.math.ceil
import kotlin.math.cos
import kotlin.math.floor
import kotlin.math.hypot
import kotlin.math.min
import kotlin.math.PI
import kotlin.math.roundToLong
import kotlin.math.sin
import kotlin.math.sqrt
import kotlin.math.tan

/**
 * An invalid body was specified for the given function.
 */
class InvalidBodyException(body: Body) : Exception("Invalid body: $body")

/**
 * The Earth is not allowed as the body parameter.
 */
class EarthNotAllowedException() : Exception("The Earth is not allowed as the body parameter.")

/**
 * An unexpected internal error occurred in Astronomy Engine
 */
class InternalError(message: String) : Exception(message)

/**
 * The factor to convert degrees to radians = pi/180.
 */
const val DEG2RAD = 0.017453292519943296

/**
 * Convert an angle expressed in degrees to an angle expressed in radians.
 */
fun Double.degreesToRadians() = this * DEG2RAD

internal fun dsin(degrees: Double) = sin(degrees.degreesToRadians())
internal fun dcos(degrees: Double) = cos(degrees.degreesToRadians())
internal fun dtan(degrees: Double) = tan(degrees.degreesToRadians())

/**
 * The factor to convert radians to degrees = 180/pi.
 */
const val RAD2DEG = 57.295779513082321

/**
 * The factor to convert radians to sidereal hours = 12/pi.
 */
const val RAD2HOUR = 3.819718634205488

/**
 * The factor to convert sidereal hours to radians = pi/12.
 */
const val HOUR2RAD = 0.2617993877991494365


/**
 * The equatorial radius of Jupiter, expressed in kilometers.
 */
const val JUPITER_EQUATORIAL_RADIUS_KM = 71492.0

/**
 * The polar radius of Jupiter, expressed in kilometers.
 */
const val JUPITER_POLAR_RADIUS_KM = 66854.0

/**
 * The volumetric mean radius of Jupiter, expressed in kilometers.
 */
const val JUPITER_MEAN_RADIUS_KM = 69911.0

/**
 * The mean radius of Jupiter's moon Io, expressed in kilometers.
 */
const val IO_RADIUS_KM = 1821.6

/**
 * The mean radius of Jupiter's moon Europa, expressed in kilometers.
 */
const val EUROPA_RADIUS_KM = 1560.8

/**
 * The mean radius of Jupiter's moon Ganymede, expressed in kilometers.
 */
const val GANYMEDE_RADIUS_KM = 2631.2

/**
 * The mean radius of Jupiter's moon Callisto, expressed in kilometers.
 */
const val CALLISTO_RADIUS_KM = 2410.3

/**
 * The speed of light in AU/day.
 */
const val C_AUDAY = 173.1446326846693


/**
 * Convert an angle expressed in radians to an angle expressed in degrees.
 */
fun Double.radiansToDegrees() = this * RAD2DEG

/**
 * The number of kilometers in one astronomical unit (AU).
 */
const val KM_PER_AU = 1.4959787069098932e+8


private val TimeZoneUtc = TimeZone.getTimeZone("UTC")
private const val DAYS_PER_TROPICAL_YEAR = 365.24217
private const val DAYS_PER_MILLENNIUM = 365250.0

private const val ASEC360 = 1296000.0
private const val ASEC2RAD = 4.848136811095359935899141e-6
private const val PI2 = 2.0 * PI
private const val ARC = 3600.0 * 180.0 / PI       // arcseconds per radian
private const val SUN_RADIUS_KM  = 695700.0
private const val SUN_RADIUS_AU  = SUN_RADIUS_KM / KM_PER_AU
private const val EARTH_FLATTENING = 0.996647180302104
private const val EARTH_FLATTENING_SQUARED = EARTH_FLATTENING * EARTH_FLATTENING
private const val EARTH_EQUATORIAL_RADIUS_KM = 6378.1366
private const val EARTH_EQUATORIAL_RADIUS_AU = EARTH_EQUATORIAL_RADIUS_KM / KM_PER_AU
private const val EARTH_POLAR_RADIUS_KM = EARTH_EQUATORIAL_RADIUS_KM * EARTH_FLATTENING
private const val EARTH_MEAN_RADIUS_KM = 6371.0    // mean radius of the Earth's geoid, without atmosphere
private const val EARTH_ATMOSPHERE_KM = 88.0       // effective atmosphere thickness for lunar eclipses
private const val EARTH_ECLIPSE_RADIUS_KM = EARTH_MEAN_RADIUS_KM + EARTH_ATMOSPHERE_KM
private const val MOON_EQUATORIAL_RADIUS_KM = 1738.1
private const val MOON_MEAN_RADIUS_KM       = 1737.4
private const val MOON_POLAR_RADIUS_KM      = 1736.0
private const val MOON_EQUATORIAL_RADIUS_AU = (MOON_EQUATORIAL_RADIUS_KM / KM_PER_AU)
private const val ANGVEL = 7.2921150e-5
private const val SECONDS_PER_DAY = 24.0 * 3600.0
private const val SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592
private const val MEAN_SYNODIC_MONTH = 29.530588     // average number of days for Moon to return to the same phase
private const val EARTH_ORBITAL_PERIOD = 365.256
private const val NEPTUNE_ORBITAL_PERIOD = 60189.0
private const val REFRACTION_NEAR_HORIZON = 34.0 / 60.0  // degrees of refractive "lift" seen for objects near horizon
private const val ASEC180 = 180.0 * 60.0 * 60.0          // arcseconds per 180 degrees (or pi radians)
private const val AU_PER_PARSEC = (ASEC180 / PI)         // exact definition of how many AU = one parsec
private const val EARTH_MOON_MASS_RATIO = 81.30056

/*
    Masses of the Sun and outer planets, used for:
    (1) Calculating the Solar System Barycenter
    (2) Integrating the movement of Pluto

    https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf

    Page 10 in the above document describes the constants used in the DE405 ephemeris.
    The following are G*M values (gravity constant * mass) in [au^3 / day^2].
    This side-steps issues of not knowing the exact values of G and masses M[i];
    the products GM[i] are known extremely accurately.
*/
private const val SUN_GM     = 0.2959122082855911e-03
private const val MERCURY_GM = 0.4912547451450812e-10
private const val VENUS_GM   = 0.7243452486162703e-09
private const val EARTH_GM   = 0.8887692390113509e-09
private const val MARS_GM    = 0.9549535105779258e-10
private const val JUPITER_GM = 0.2825345909524226e-06
private const val SATURN_GM  = 0.8459715185680659e-07
private const val URANUS_GM  = 0.1292024916781969e-07
private const val NEPTUNE_GM = 0.1524358900784276e-07
private const val PLUTO_GM   = 0.2188699765425970e-11
private const val MOON_GM = EARTH_GM / EARTH_MOON_MASS_RATIO


private fun Double.withMinDegreeValue(min: Double): Double {
    var deg = this
    while (deg < min)
        deg += 360.0
    while (deg >= min + 360.0)
        deg -= 360.0
    return deg
}

private fun Double.withMaxDegreeValue(max: Double): Double {
    var deg = this
    while (deg <= max - 360.0)
        deg += 360.0
    while (deg > max)
        deg -= 360.0
    return deg
}


private fun toggleAzimuthDirection(az: Double) = (360.0 - az).withMinDegreeValue(0.0)
private fun longitudeOffset(diff: Double) = diff.withMaxDegreeValue(+180.0)
private fun normalizeLongitude(lon: Double) = lon.withMinDegreeValue(0.0)

/**
 * The enumeration of celestial bodies supported by Astronomy Engine.
 */
enum class Body(
    internal val massProduct: Double?,
    internal val vsopModel: VsopModel?
) {
    /**
     * The planet Mercury.
     */
    Mercury(
        MERCURY_GM,
        VsopModel(vsopLonMercury, vsopLatMercury, vsopRadMercury)
    ),

    /**
     * The planet Venus.
     */
    Venus(
        VENUS_GM,
        VsopModel(vsopLonVenus, vsopLatVenus, vsopRadVenus)
    ),

    /**
     * The planet Earth.
     * Some functions that accept a `Body` parameter will fail if passed this value
     * because they assume that an observation is being made from the Earth,
     * and therefore the Earth is not a target of observation.
     */
    Earth(
        EARTH_GM,
        VsopModel(vsopLonEarth, vsopLatEarth, vsopRadEarth)
    ),

    /**
     * The planet Mars.
     */
    Mars(
        MARS_GM,
        VsopModel(vsopLonMars, vsopLatMars, vsopRadMars)
    ),

    /**
     * The planet Jupiter.
     */
    Jupiter(
        JUPITER_GM,
        VsopModel(vsopLonJupiter, vsopLatJupiter, vsopRadJupiter)
    ),

    /**
     * The planet Saturn.
     */
    Saturn(
        SATURN_GM,
        VsopModel(vsopLonSaturn, vsopLatSaturn, vsopRadSaturn)
    ),

    /**
     * The planet Uranus.
     */
    Uranus(
        URANUS_GM,
        VsopModel(vsopLonUranus, vsopLatUranus, vsopRadUranus)
    ),

    /**
     * The planet Neptune.
     */
    Neptune(
        NEPTUNE_GM,
        VsopModel(vsopLonNeptune, vsopLatNeptune, vsopRadNeptune)
    ),

    /**
     * The planet Pluto.
     */
    Pluto(
        PLUTO_GM,
        null,
    ),

    /**
     * The Sun.
     */
    Sun(
        SUN_GM,
        null
    ),

    /**
     * The Earth's natural satellite, the Moon.
     */
    Moon(
        MOON_GM,
        null
    ),

    /**
     * The Earth/Moon Barycenter.
     */
    EMB(
        EARTH_GM + MOON_GM,
        null
    ),

    /**
     * The Solar System Barycenter.
     */
    SSB(
        null,
        null
    ),
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

    constructor(ut: Double) : this(ut, terrestrialTime(ut))

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
     * @return A UTC `Date` object for this `AstroTime` value.
     */
    fun toDate(): Date = Date(origin.time + (ut * MILLIS_PER_DAY).roundToLong())

    /**
     * Converts this `AstroTime` to ISO 8601 format, expressed in UTC with millisecond resolution.
     *
     * @return Example: "2019-08-30T17:45:22.763Z".
     */
    override fun toString(): String = dateFormat.format(toDate())

    /**
     * Calculates the sum or difference of an [AstroTime] with a specified floating point number of days.
     *
     * Sometimes we need to adjust a given [AstroTime] value by a certain amount of time.
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

    internal fun julianCenturies() = tt / 36525.0
    internal fun julianMillennia() = tt / DAYS_PER_MILLENNIUM

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
        fun fromTerrestrialTime(tt: Double): AstroTime = AstroTime(universalTime(tt), tt)
    }
}


internal data class TerseVector(var x: Double, var y: Double, var z: Double) {
    fun toAstroVector(time: AstroTime) =
        AstroVector(x, y, z, time)

    operator fun plus(other: TerseVector) =
        TerseVector(x + other.x, y + other.y, z + other.z)

    operator fun minus(other: TerseVector) =
        TerseVector(x - other.x, y - other.y, z - other.z)

    operator fun unaryMinus() =
        TerseVector(-x, -y, -z)

    operator fun times(other: Double) =
        TerseVector(x * other, y * other, z * other)

    operator fun div(denom: Double) =
        TerseVector(x / denom, y / denom, z / denom)

    fun mean(other: TerseVector) =
        TerseVector((x + other.x) / 2.0, (y + other.y) / 2.0, (z + other.z) / 2.0)

    fun quadrature() = (x * x) + (y * y) + (z * z)
    fun magnitude() = sqrt(quadrature())

    fun decrement(other: TerseVector) {
        x -= other.x
        y -= other.y
        z -= other.z
    }

    fun increment(other: TerseVector) {
        x += other.x
        y += other.y
        z += other.z
    }

    fun mix(ramp: Double, other: TerseVector) {
        x = (1.0 - ramp)*x + ramp*other.x
        y = (1.0 - ramp)*y + ramp*other.y
        z = (1.0 - ramp)*z + ramp*other.z
    }

    companion object {
        fun zero() = TerseVector(0.0, 0.0, 0.0)
    }
}

internal operator fun Double.times(vec: TerseVector) =
    TerseVector(this * vec.x, this * vec.y, this * vec.z)


internal fun verifyIdenticalTimes(t1: AstroTime, t2: AstroTime): AstroTime {
    if (t1.tt != t2.tt)
        throw IllegalArgumentException("Attempt to operate on two vectors from different times.")
    return t1
}


/**
 * A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).
 */
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
    fun length() = sqrt((x * x) + (y * y) + (z * z))

    /**
     * Adds two vectors. Both operands must have identical times.
     */
    operator fun plus(other: AstroVector) =
        AstroVector(x + other.x, y + other.y, z + other.z, verifyIdenticalTimes(t, other.t))

    /**
     * Subtracts one vector from another. Both operands must have identical times.
     */
    operator fun minus(other: AstroVector) =
        AstroVector(x - other.x, y - other.y, z - other.z, verifyIdenticalTimes(t, other.t))

    /**
     * Negates a vector; the same as multiplying the vector by the scalar -1.
     */
    operator fun unaryMinus() =
        AstroVector(-x, -y, -z, t)

    /**
     * Takes the dot product of two vectors.
     */
    infix fun dot(other: AstroVector): Double {
        verifyIdenticalTimes(t, other.t)
        return x*other.x + y*other.y + z*other.z
    }

    /**
     * Calculates the angle in degrees (0..180) between two vectors.
     */
    fun angleWith(other: AstroVector): Double {
        val d = (this dot other) / (length() * other.length())
        return when {
            d <= -1.0 -> 180.0
            d >= +1.0 -> 0.0
            else -> acos(d).radiansToDegrees()
        }
    }

    /**
     * Divides a vector by a scalar.
     */
    operator fun div(denom: Double) =
        AstroVector(x/denom, y/denom, z/denom, t)

    /**
     * Converts Cartesian coordinates to spherical coordinates.
     *
     * Given a Cartesian vector, returns latitude, longitude, and distance.
     */
    fun toSpherical(): Spherical {
        val xyproj = x*x + y*y
        val dist = sqrt(xyproj + z*z)
        var lat: Double
        var lon: Double
        if (xyproj == 0.0) {
            if (z == 0.0) {
                // Indeterminate coordinates; pos vector has zero length.
                throw IllegalArgumentException("Cannot convert zero-length vector to spherical coordinates.")
            }
            lon = 0.0
            lat = if (z < 0.0) -90.0 else +90.0
        } else {
            lon = atan2(y, x).radiansToDegrees().withMinDegreeValue(0.0)
            lat = atan2(z, sqrt(xyproj)).radiansToDegrees()
        }
        return Spherical(lat, lon, dist)
    }

    /**
     * Given an equatorial vector, calculates equatorial angular coordinates.
     */
    fun toEquatorial(): Equatorial {
        val sphere = toSpherical()
        return Equatorial(sphere.lon / 15.0, sphere.lat, sphere.dist, this)
    }

    /**
     * Converts Cartesian coordinates to horizontal coordinates.
     *
     * Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.
     * *IMPORTANT:* This function differs from [AstroVector.toSpherical] in two ways:
     * - `toSpherical` returns a `lon` value that represents azimuth defined counterclockwise
     *   from north (e.g., west = +90), but this function represents a clockwise rotation
     *   (e.g., east = +90). The difference is because `toSpherical` is intended
     *   to preserve the vector "right-hand rule", while this function defines azimuth in a more
     *   traditional way as used in navigation and cartography.
     * - This function optionally corrects for atmospheric refraction, while `toSpherical`
     *   does not.
     *
     * The returned object contains the azimuth in `lon`.
     * It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.
     *
     * The altitude is stored in `lat`.
     *
     * The distance to the observed object is stored in `dist`,
     * and is expressed in astronomical units (AU).
     *
     * @param refraction
     *      `Refraction.None`: no atmospheric refraction correction is performed.
     *      `Refraction.Normal`: correct altitude for atmospheric refraction.
     *      `Refraction.JplHor`: for JPL Horizons compatibility testing only; not recommended for normal use.
     */
    fun toHorizontal(refraction: Refraction): Spherical {
        val sphere = toSpherical()
        return Spherical(
            sphere.lat + refractionAngle(refraction, sphere.lat),
            toggleAzimuthDirection(sphere.lon),
            sphere.dist
        )
    }
}

/**
 * Multiply a scalar by a vector, yielding another vector.
 */
operator fun Double.times(vec: AstroVector) =
        AstroVector(this*vec.x, this*vec.y, this*vec.z, vec.t)


/**
 * Represents a combined position vector and velocity vector at a given moment in time.
 */
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

    internal constructor(state: BodyState, time: AstroTime)
        : this(state.r.x, state.r.y, state.r.z, state.v.x, state.v.y, state.v.z, time)

    /**
     * Returns the position vector associated with this state vector.
     */
    fun position() = AstroVector(x, y, z, t)

    /**
     * Returns the velocity vector associated with this state vector.
     */
    fun velocity() = AstroVector(vx, vy, vz, t)

    /**
     * Adds two state vetors, yielding the state vector sum.
     */
    operator fun plus(other: StateVector) =
        StateVector(
            x + other.x,
            y + other.y,
            z + other.z,
            vx + other.vx,
            vy + other.vy,
            vz + other.vz,
            verifyIdenticalTimes(t, other.t)
        )

    /**
     * Divides a state vector by a scalar.
     */
    operator fun div(denom: Double) =
        StateVector(
            x / denom,
            y / denom,
            z / denom,
            vx / denom,
            vy / denom,
            vz / denom,
            t
        )
}


/**
 * Holds the positions and velocities of Jupiter's major 4 moons.
 *
 * The [jupiterMoons] function returns an object of this type
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
        if (rot.size != 3 || rot.any { it.size != 3 })
            throw IllegalArgumentException("Rotation matrix must be a 3x3 array.")
    }

    constructor(
        a00: Double, a01: Double, a02: Double,
        a10: Double, a11: Double, a12: Double,
        a20: Double, a21: Double, a22: Double
    ) : this(
        arrayOf(
            doubleArrayOf(a00, a01, a02),
            doubleArrayOf(a10, a11, a12),
            doubleArrayOf(a20, a21, a22)
        )
    )

    /**
     * Calculates the inverse of a rotation matrix.
     *
     * Returns a rotation matrix that performs the reverse transformation
     * that this rotation matrix performs.
     */
    fun inverse() = RotationMatrix(
        rot[0][0], rot[1][0], rot[2][0],
        rot[0][1], rot[1][1], rot[2][1],
        rot[0][2], rot[1][2], rot[2][2]
    )

    /**
     * Applies a rotation to a vector, yielding a rotated vector.
     *
     * This function transforms a vector in one orientation to a vector
     * in another orientation.
     *
     * @param vec
     *      The vector whose orientation is to be changed.
     */
    fun rotate(vec: AstroVector) = AstroVector(
        rot[0][0]*vec.x + rot[1][0]*vec.y + rot[2][0]*vec.z,
        rot[0][1]*vec.x + rot[1][1]*vec.y + rot[2][1]*vec.z,
        rot[0][2]*vec.x + rot[1][2]*vec.y + rot[2][2]*vec.z,
        vec.t
    )

    /**
     * Applies a rotation to a state vector, yielding a rotated state vector.
     *
     * This function transforms a state vector in one orientation to a state vector in another orientation.
     * The resulting state vector has both position and velocity reoriented.
     *
     * @param state
     *      The state vector whose orientation is to be changed.
     *      The value of `state` is not changed; the return value is a new state vector object.
     */
    fun rotate(state: StateVector) = StateVector(
        rotate(state.position()),
        rotate(state.velocity()),
        state.t
    )

    /**
     * Creates a rotation based on applying one rotation followed by another.
     *
     * Given two rotation matrices, returns a combined rotation matrix that is
     * equivalent to rotating based on this matrix, followed by the matrix `other`.
     */
    infix fun combine(other: RotationMatrix) = RotationMatrix (
        other.rot[0][0]*rot[0][0] + other.rot[1][0]*rot[0][1] + other.rot[2][0]*rot[0][2],
        other.rot[0][1]*rot[0][0] + other.rot[1][1]*rot[0][1] + other.rot[2][1]*rot[0][2],
        other.rot[0][2]*rot[0][0] + other.rot[1][2]*rot[0][1] + other.rot[2][2]*rot[0][2],
        other.rot[0][0]*rot[1][0] + other.rot[1][0]*rot[1][1] + other.rot[2][0]*rot[1][2],
        other.rot[0][1]*rot[1][0] + other.rot[1][1]*rot[1][1] + other.rot[2][1]*rot[1][2],
        other.rot[0][2]*rot[1][0] + other.rot[1][2]*rot[1][1] + other.rot[2][2]*rot[1][2],
        other.rot[0][0]*rot[2][0] + other.rot[1][0]*rot[2][1] + other.rot[2][0]*rot[2][2],
        other.rot[0][1]*rot[2][0] + other.rot[1][1]*rot[2][1] + other.rot[2][1]*rot[2][2],
        other.rot[0][2]*rot[2][0] + other.rot[1][2]*rot[2][1] + other.rot[2][2]*rot[2][2]
    )

    /**
     * Re-orients the rotation matrix by pivoting it by an angle around one of its axes.
     *
     * Given this rotation matrix, a selected coordinate axis, and an angle in degrees,
     * this function pivots the rotation matrix by that angle around that coordinate axis.
     * The function returns a new rotation matrix; it does not mutate this matrix.
     *
     * For example, if you have rotation matrix that converts ecliptic coordinates (ECL)
     * to horizontal coordinates (HOR), but you really want to convert ECL to the orientation
     * of a telescope camera pointed at a given body, you can use `RotationMatrix.pivot` twice:
     * (1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the
     * western axis by the body's altitude angle. The resulting rotation matrix will then
     * reorient ECL coordinates to the orientation of your telescope camera.
     *
     * @param axis
     *      An integer that selects which coordinate axis to rotate around:
     *      0 = x, 1 = y, 2 = z. Any other value will cause an exception.
     *
     * @param angle
     *      An angle in degrees indicating the amount of rotation around the specified axis.
     *      Positive angles indicate rotation counterclockwise as seen from the positive
     *      direction along that axis, looking towards the origin point of the orientation system.
     *      Any finite number of degrees is allowed, but best precision will result from keeping
     *      `angle` in the range [-360, +360].
     */
    fun pivot(axis: Int, angle: Double): RotationMatrix {
        if (axis < 0 || axis > 2)
            throw IllegalArgumentException("Invalid coordinate axis $axis. Must be 0..2.")

        if (!angle.isFinite())
            throw IllegalArgumentException("Angle must be a finite number.")

        val radians = angle.degreesToRadians()
        val c = cos(radians)
        val s = sin(radians)

        // We need to maintain the "right-hand" rule, no matter which
        // axis was selected. That means we pick (i, j, k) axis order
        // such that the following vector cross product is satisfied:
        // i x j = k
        val i = (axis + 1) % 3
        val j = (axis + 2) % 3
        val k = axis

        val piv = arrayOf(DoubleArray(3), DoubleArray(3), DoubleArray(3))

        piv[i][i] = c*rot[i][i] - s*rot[i][j]
        piv[i][j] = s*rot[i][i] + c*rot[i][j]
        piv[i][k] = rot[i][k]
        piv[j][i] = c*rot[j][i] - s*rot[j][j]
        piv[j][j] = s*rot[j][i] + c*rot[j][j]
        piv[j][k] = rot[j][k]
        piv[k][i] = c*rot[k][i] - s*rot[k][j]
        piv[k][j] = s*rot[k][i] + c*rot[k][j]
        piv[k][k] = rot[k][k]

        return RotationMatrix(piv)
    }

    companion object {
        /**
         * The identity rotation matrix.
         *
         * A matrix that has no effect on orientation.
         * This matrix can be the starting point for other operations,
         * such as calling a series of [RotationMatrix.combine] or [RotationMatrix.pivot].
         */
        val identity = RotationMatrix (
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        )
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
) {
    /**
     * Converts spherical coordinates to Cartesian coordinates.
     *
     * Given spherical coordinates and a time at which they are valid,
     * returns a vector of Cartesian coordinates. The returned value
     * includes the time, as required by the type [AstroVector].
     *
     * @param time
     *      The time that should be included in the return value.
     */
    fun toVector(time: AstroTime): AstroVector {
        val radlat = lat.degreesToRadians()
        val radlon = lon.degreesToRadians()
        val rcoslat = dist * cos(radlat)
        return AstroVector(
            rcoslat * cos(radlon),
            rcoslat * sin(radlon),
            dist * sin(radlat),
            time
        )
    }

    /**
     * Given apparent angular horizontal coordinates, calculate the unrefracted horizontal vector.
     *
     * Assumes `this` contains apparent horizontal coordinates:
     * `lat` holds the refracted azimuth angle,
     * `lon` holds the azimuth in degrees clockwise from north,
     * and `dist` holds the distance from the observer to the object in AU.
     *
     * @param time
     *      The date and time of the observation. This is needed because the returned
     *      [AstroVector] requires a valid time value when passed to certain other functions.
     *
     * @param refraction
     *      The refraction option used to model atmospheric lensing. See [refractionAngle].
     *      This specifies how refraction is to be removed from the altitude stored in `this.lat`.
     *
     * @return A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).
     */
    fun toVectorFromHorizon(time: AstroTime, refraction: Refraction): AstroVector =
        Spherical(
            lat + inverseRefractionAngle(refraction, lat),
            toggleAzimuthDirection(lon),
            dist
        )
        .toVector(time)
}

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
enum class Direction(
    /**
     * A numeric value that is helpful in formulas involving rise/set.
     * The sign is +1 for a rising event, or -1 for a setting event.
     */
    val sign: Int
) {
    /**
     * Indicates a rising event: a celestial body is observed to rise above the horizon by an observer on the Earth.
     */
    Rise(+1),

    /**
     * Indicates a setting event: a celestial body is observed to sink below the horizon by an observer on the Earth.
     */
    Set(-1),
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


/**
 * Ecliptic angular and Cartesian coordinates.
 *
 * Coordinates of a celestial body as seen from the center of the Sun (heliocentric),
 * oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).
 */
data class Ecliptic(
    /**
     * Cartesian ecliptic vector, with components as follows:
     * x: the direction of the equinox along the ecliptic plane.
     * y: in the ecliptic plane 90 degrees prograde from the equinox.
     * z: perpendicular to the ecliptic plane. Positive is north.
     */
    val vec: AstroVector,

    /**
     * Latitude in degrees north (positive) or south (negative) of the ecliptic plane.
     */
    val elat: Double,

    /**
     * Longitude in degrees around the ecliptic plane prograde from the equinox.
     */
    val elon: Double
)


/**
 * Coordinates of a celestial body as seen by a topocentric observer.
 *
 * Horizontal and equatorial coordinates seen by an observer on or near
 * the surface of the Earth (a topocentric observer).
 * Optionally corrected for atmospheric refraction.
 */
data class Topocentric(
    /**
     * Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West.
     */
    val azimuth: Double,

    /**
     * Angle in degrees above (positive) or below (negative) the observer's horizon.
     */
    val altitude: Double,

    /**
     * Right ascension in sidereal hours.
     */
    val ra: Double,

    /**
     * Declination in degrees.
     */
    val dec: Double
)


/**
 * The dates and times of changes of season for a given calendar year.
 *
 * Call [seasons] to calculate this data structure for a given year.
 */
class SeasonsInfo(
    /**
     * The date and time of the March equinox for the specified year.
     */
    val marchEquinox: AstroTime,

    /**
     * The date and time of the June soltice for the specified year.
     */
    val juneSolstice: AstroTime,

    /**
     * The date and time of the September equinox for the specified year.
     */
    val septemberEquinox: AstroTime,

    /**
     * The date and time of the December solstice for the specified year.
     */
    val decemberSolstice: AstroTime
)


/**
 * A lunar quarter event (new moon, first quarter, full moon, or third quarter) along with its date and time.
 */
class MoonQuarterInfo(
    /**
    * 0=new moon, 1=first quarter, 2=full moon, 3=third quarter.
    */
    val quarter: Int,

    /**
    * The date and time of the lunar quarter.
    */
    val time: AstroTime
)


/**
 * Lunar libration angles, returned by [libration].
 */
data class LibrationInfo(
    /**
     * Sub-Earth libration ecliptic latitude angle, in degrees.
     */
    val elat: Double,

    /**
     * Sub-Earth libration ecliptic longitude angle, in degrees.
     */
    val elon: Double,

    /**
     * Moon's geocentric ecliptic latitude.
     */
    val mlat: Double,

    /**
     * Moon's geocentric ecliptic longitude.
     */
    val mlon: Double,

    /**
     * Distance between the centers of the Earth and Moon in kilometers.
     */
    val distKm: Double,

    /**
     * The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth.
     */
    val diamDeg: Double
)


/**
 * Information about a celestial body crossing a specific hour angle.
 *
 * Returned by the function [searchHourAngle] to report information about
 * a celestial body crossing a certain hour angle as seen by a specified topocentric observer.
 */
class HourAngleInfo(
    /**
     * The date and time when the body crosses the specified hour angle.
     */
    val time: AstroTime,

    /**
     * Apparent coordinates of the body at the time it crosses the specified hour angle.
     */
    val hor: Topocentric
)


/**
 * Contains information about the visibility of a celestial body at a given date and time.
 *
 * See [elongation] for more detailed information about the members of this class.
 * See also [searchMaxElongation] for how to search for maximum elongation events.
 */
class ElongationInfo(
    /**
     * The date and time of the observation.
     */
    val time: AstroTime,

    /**
     * Whether the body is best seen in the morning or the evening.
     */
    val visibility: Visibility,

    /**
     * The angle in degrees between the body and the Sun, as seen from the Earth.
     */
    val elongation: Double,

    /**
     * The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth.
     */
    val eclipticSeparation: Double
)


/**
 * The type of apsis: pericenter (closest approach) or apocenter (farthest distance).
 */
enum class ApsisKind {
    /**
     * The body is at its closest approach to the object it orbits.
     */
    Pericenter,

    /**
     * The body is at its farthest distance from the object it orbits.
     */
    Apocenter,
}


/**
 * An apsis event: pericenter (closest approach) or apocenter (farthest distance).
 *
 * For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an
 * event where the orbiting body reaches its closest or farthest point from the primary body.
 * The closest approach is called *pericenter* and the farthest point is *apocenter*.
 *
 * More specific terminology is common for particular orbiting bodies.
 * The Moon's closest approach to the Earth is called *perigee* and its farthest
 * point is called *apogee*. The closest approach of a planet to the Sun is called
 * *perihelion* and the furthest point is called *aphelion*.
 *
 * This data structure is returned by [searchLunarApsis] and [nextLunarApsis]
 * to iterate through consecutive alternating perigees and apogees.
 */
class ApsisInfo(
    /**
     * The date and time of the apsis.
     */
    val time: AstroTime,

    /**
     * Whether this is a pericenter or apocenter event.
     */
    val kind: ApsisKind,

    /**
     * The distance between the centers of the bodies in astronomical units.
     */
    val distAu: Double,

    /**
     * The distance between the centers of the bodies in kilometers.
     */
    val distKm: Double
)


/**
 * The different kinds of lunar/solar eclipses.
 */
enum class EclipseKind {
    /**
     * No eclipse found.
     */
    None,

    /**
     * A penumbral lunar eclipse. (Never used for a solar eclipse.)
     */
    Penumbral,

    /**
     * A partial lunar/solar eclipse.
     */
    Partial,

    /**
     * An annular solar eclipse. (Never used for a lunar eclipse.)
     */
    Annular,

    /**
     * A total lunar/solar eclipse.
     */
    Total,
}


/**
 * Information about a lunar eclipse.
 *
 * Returned by [searchLunarEclipse] or [nextLunarEclipse]
 * to report information about a lunar eclipse event.
 * When a lunar eclipse is found, it is classified as penumbral, partial, or total.
 * Penumbral eclipses are difficult to observe, because the Moon is only slightly dimmed
 * by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
 * Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
 * Total eclipses occur when the entire Moon passes into the Earth's umbra.
 *
 * The `kind` field thus holds `EclipseKind.Penumbral`, `EclipseKind.Partial`,
 * or `EclipseKind.Total`, depending on the kind of lunar eclipse found.
 *
 * Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
 *
 * Fields `sdPenum`, `sdPartial`, and `sdTotal` hold the semi-duration of each phase
 * of the eclipse, which is half of the amount of time the eclipse spends in each
 * phase (expressed in minutes), or 0.0 if the eclipse never reaches that phase.
 * By converting from minutes to days, and subtracting/adding with `peak`, the caller
 * may determine the date and time of the beginning/end of each eclipse phase.
 */
class LunarEclipseInfo(

    /**
     * The type of lunar eclipse found.
     */
    val kind: EclipseKind,

    /**
     * The time of the eclipse at its peak.
     */
    val peak: AstroTime,

    /**
     * The semi-duration of the penumbral phase in minutes.
     */
    val sdPenum: Double,

    /**
     * The semi-duration of the partial phase in minutes, or 0.0 if none.
     */
    val sdPartial: Double,

    /**
     * The semi-duration of the total phase in minutes, or 0.0 if none.
     */
    val sdTotal: Double
)


/**
 * Reports the time and geographic location of the peak of a solar eclipse.
 *
 * Returned by [searchGlobalSolarEclipse] or [nextGlobalSolarEclipse]
 * to report information about a solar eclipse event.
 *
 * The eclipse is classified as partial, annular, or total, depending on the
 * maximum amount of the Sun's disc obscured, as seen at the peak location
 * on the surface of the Earth.
 *
 * The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
 * A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
 * An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
 * to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
 * A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
 * observes either a total or annular eclipse.
 *
 * If `kind` is `EclipseKind.Total` or `EclipseKind.Annular`, the `latitude` and `longitude`
 * fields give the geographic coordinates of the center of the Moon's shadow projected
 * onto the daytime side of the Earth at the instant of the eclipse's peak.
 * If `kind` has any other value, `latitude` and `longitude` are undefined and should
 * not be used.
 */
class GlobalSolarEclipseInfo(
    /**
     * The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
     */
    val kind: EclipseKind,

    /**
     * The date and time when the solar eclipse is darkest.
     * This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.
     */
    val peak: AstroTime,

    /**
     * The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers.
     */
    val distance: Double,

    /**
     * The geographic latitude at the center of the peak eclipse shadow.
     */
    val latitude: Double,

    /**
     * The geographic longitude at the center of the peak eclipse shadow.
     */
    val longitude: Double
)


/**
 * Holds a time and the observed altitude of the Sun at that time.
 *
 * When reporting a solar eclipse observed at a specific location on the Earth
 * (a "local" solar eclipse), a series of events occur. In addition
 * to the time of each event, it is important to know the altitude of the Sun,
 * because each event may be invisible to the observer if the Sun is below
 * the horizon (i.e. it at night).
 *
 * If `altitude` is negative, the event is theoretical only; it would be
 * visible if the Earth were transparent, but the observer cannot actually see it.
 * If `altitude` is positive but less than a few degrees, visibility will be impaired by
 * atmospheric interference (sunrise or sunset conditions).
 */
class EclipseEvent (
    /**
     * The date and time of the event.
     */
    val time: AstroTime,

    /**
     * The angular altitude of the center of the Sun above/below the horizon, at `time`,
     * corrected for atmospheric refraction and expressed in degrees.
     */
    val altitude: Double
)


/**
 * Information about a solar eclipse as seen by an observer at a given time and geographic location.
 *
 * Returned by [searchLocalSolarEclipse] or [nextLocalSolarEclipse]
 * to report information about a solar eclipse as seen at a given geographic location.
 *
 * When a solar eclipse is found, it is classified as partial, annular, or total.
 * The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
 * A partial solar eclipse is when the Moon does not line up directly enough with the Sun
 * to completely block the Sun's light from reaching the observer.
 * An annular eclipse occurs when the Moon's disc is completely visible against the Sun
 * but the Moon is too far away to completely block the Sun's light; this leaves the
 * Sun with a ring-like appearance.
 * A total eclipse occurs when the Moon is close enough to the Earth and aligned with the
 * Sun just right to completely block all sunlight from reaching the observer.
 *
 * There are 5 "event" fields, each of which contains a time and a solar altitude.
 * Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
 * The fields `partialBegin` and `partialEnd` are always set, and indicate when
 * the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
 * `totalBegin` and `totalEnd` indicate when the total/annular phase begins/ends.
 * When an event field is valid, the caller must also check its `altitude` field to
 * see whether the Sun is above the horizon at the time indicated by the `time` field.
 * </remarks>
 */
class LocalSolarEclipseInfo (
    /**
     * The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
     */
    val kind: EclipseKind,

    /**
     * The time and Sun altitude at the beginning of the eclipse.
     */
    val partialBegin: EclipseEvent,

    /**
     * If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise invalid.
     */
    val totalBegin: EclipseEvent,

    /**
     * The time and Sun altitude when the eclipse reaches its peak.
     */
    val peak: EclipseEvent,

    /**
     * If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise invalid.
     */
    val totalEnd: EclipseEvent,

    /**
     * The time and Sun altitude at the end of the eclipse.
     */
    val partialEnd: EclipseEvent
)


/**
 * Information about a transit of Mercury or Venus, as seen from the Earth.
 *
 * Returned by [searchTransit] or [nextTransit] to report
 * information about a transit of Mercury or Venus.
 * A transit is when Mercury or Venus passes between the Sun and Earth so that
 * the other planet is seen in silhouette against the Sun.
 *
 * The `start` field reports the moment in time when the planet first becomes
 * visible against the Sun in its background.
 * The `peak` field reports when the planet is most aligned with the Sun,
 * as seen from the Earth.
 * The `finish` field reports the last moment when the planet is visible
 * against the Sun in its background.
 *
 * The calculations are performed from the point of view of a geocentric observer.
 */
class TransitInfo(
    /**
     * Date and time at the beginning of the transit.
     */
    val start: AstroTime,

    /**
     * Date and time of the peak of the transit.
     */
    val peak: AstroTime,

    /**
     * Date and time at the end of the transit.
     */
    val finish: AstroTime,

    /**
     * Angular separation in arcminutes between the centers of the Sun and the planet at time `peak`.
     */
    val separation: Double
)


internal class ShadowInfo(
    /**
     * The time of the shadow state.
     */
    val time: AstroTime,

    /**
     * Dot product of (heliocentric earth) and (geocentric moon).
     * Defines the shadow plane where the Moon is.
     */
    val u: Double,

    /**
     * km distance between center of Moon and the line passing through the centers of the Sun and Earth.
     */
    val r: Double,

    /**
     * Umbra radius in km, at the shadow plane.
     */
    val k: Double,

    /**
     * Penumbra radius in km, at the shadow plane.
     */
    val p: Double,

    /**
     * Coordinates of target body relative to shadow-casting body at `time`.
     */
    val target: AstroVector,

    /**
     * Heliocentric coordinates of shadow-casting body at `time`.
     */
    val dir: AstroVector
)


/**
 * Information about the brightness and illuminated shape of a celestial body.
 *
 * Returned by the functions [illumination] and [searchPeakMagnitude]
 * to report the visual magnitude and illuminated fraction of a celestial body at a given date and time.
 */
class IlluminationInfo(
    /**
     * The date and time of the observation.
     */
    val time: AstroTime,

    /**
     * The visual magnitude of the body. Smaller values are brighter.
     */
    val mag: Double,

    /**
     * The angle in degrees between the Sun and the Earth, as seen from the body.
     * Indicates the body's phase as seen from the Earth.
     */
    val phaseAngle: Double,

    /**
     * A value in the range [0.0, 1.0] indicating what fraction of the body's
     * apparent disc is illuminated, as seen from the Earth.
     */
    val phaseFraction: Double,

    /**
     * The distance between the Sun and the body at the observation time.
     */
    val helioDist: Double,

    /**
     * For Saturn, the tilt angle in degrees of its rings as seen from Earth.
     * For all other bodies, 0.0.
     */
    val ringTilt: Double
)


/**
 * Information about a body's rotation axis at a given time.
 *
 * This structure is returned by [rotationAxis] to report
 * the orientation of a body's rotation axis at a given moment in time.
 * The axis is specified by the direction in space that the body's north pole
 * points, using angular equatorial coordinates in the J2000 system (EQJ).
 *
 * Thus `ra` is the right ascension, and `dec` is the declination, of the
 * body's north pole vector at the given moment in time. The north pole
 * of a body is defined as the pole that lies on the north side of the
 * [Solar System's invariable plane](https://en.wikipedia.org/wiki/Invariable_plane),
 * regardless of the body's direction of rotation.
 *
 * The `spin` field indicates the angular position of a prime meridian
 * arbitrarily recommended for the body by the International Astronomical
 * Union (IAU).
 *
 * The fields `ra`, `dec`, and `spin` correspond to the variables
 * α0, δ0, and W, respectively, from
 * [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).
 *
 * The field `north` is a unit vector pointing in the direction of the body's north pole.
 * It is expressed in the equatorial J2000 system (EQJ).
 */
class AxisInfo(
    /**
     * The J2000 right ascension of the body's north pole direction, in sidereal hours.
     */
    val ra: Double,

    /**
     * The J2000 declination of the body's north pole direction, in degrees.
     */
    val dec: Double,

    /**
     * Rotation angle of the body's prime meridian, in degrees.
     */
    val spin: Double,

    /**
     * A J2000 dimensionless unit vector pointing in the direction of the body's north pole.
     */
    val north: AstroVector
)


/**
 * Indicates whether a crossing through the ecliptic plane is ascending or descending.
 */
enum class NodeEventKind {
    /**
     * Placeholder value for a missing or invalid node.
     */
    Invalid,

    /**
     * The body passes through the ecliptic plane from south to north.
     */
    Ascending,

    /**
     * The body passes through the ecliptic plane from north to south.
     */
    Descending,
}


/**
 * Information about an ascending or descending node of a body.
 *
 * This object is returned by [searchMoonNode] and [nextMoonNode]
 * to report information about the center of the Moon passing through the ecliptic plane.
 */
class NodeEventInfo(
    /**
     * The time of the body's node.
     */
    val time: AstroTime,

    /**
     * Whether the node is ascending or descending.
     */
    val kind: NodeEventKind
)


/**
 * Represents a function whose ascending root is to be found.
 *
 * This interface must be implemented for callers of [search]
 * in order to find the ascending root of a smooth function.
 * A class that implements `SearchContext` can hold state information
 * needed to evaluate the scalar function `eval`.
 */
interface SearchContext {
    /**
     * Evaluates a scalar function at a given time.
     *
     * @param time
     *      The time at which to evaluate the function.
     *
     * @return The floating point value of the scalar function at the given time.
     */
    fun eval(time: AstroTime): Double
}

//---------------------------------------------------------------------------------------

/**
 * Reports the constellation that a given celestial point lies within.
 *
 * The [constellation] function returns this object
 * to report which constellation corresponds with a given point in the sky.
 * Constellations are defined with respect to the B1875 equatorial system
 * per IAU standard. Although `constellation` requires J2000 equatorial
 * coordinates, `ConstellationInfo` contains converted B1875 coordinates for reference.
 */
class ConstellationInfo(
    /**
     * 3-character mnemonic symbol for the constellation, e.g. "Ori".
     */
    val symbol: String,

    /**
     * Full name of constellation, e.g. "Orion".
     */
    val name: String,

    /**
     * Right ascension expressed in B1875 coordinates.
     */
    val ra1875: Double,

    /**
     * Declination expressed in B1875 coordinates.
     */
    val dec1875: Double
)

internal class ConstellationText(
    val symbol: String,
    val name: String
)

internal class ConstellationBoundary(
    val index: Int,
    val raLo: Double,
    val raHi: Double,
    val decLo: Double
)

//---------------------------------------------------------------------------------------
// VSOP87: semi-analytic model of major planet positions

internal class VsopTerm(
    val amplitude: Double,
    val phase: Double,
    val frequency: Double
)

internal class VsopSeries(
    val term: Array<VsopTerm>
)

internal class VsopFormula(
    val series: Array<VsopSeries>
)

internal class VsopModel(
    val lon: VsopFormula,
    val lat: VsopFormula,
    val rad: VsopFormula
)

private fun vsopFormulaCalc(formula: VsopFormula, t: Double, clampAngle: Boolean): Double {
    var coord = 0.0
    var tpower = 1.0
    for (series in formula.series) {
        var sum = 0.0
        for (term in series.term)
            sum += term.amplitude * cos(term.phase + (t * term.frequency))
        coord +=
            if (clampAngle)
                (tpower * sum) % PI2    // improve precision: longitude angles can be hundreds of radians
            else
                tpower * sum
        tpower *= t
    }
    return coord
}

private fun vsopDistance(model: VsopModel, time: AstroTime) =
    vsopFormulaCalc(model.rad, time.julianMillennia(), false)

private fun vsopRotate(eclip: TerseVector) =
    TerseVector(
        eclip.x + 0.000000440360*eclip.y - 0.000000190919*eclip.z,
        -0.000000479966*eclip.x + 0.917482137087*eclip.y - 0.397776982902*eclip.z,
        0.397776982902*eclip.y + 0.917482137087*eclip.z
    )

private fun vsopSphereToRect(lon: Double, lat: Double, radius: Double): TerseVector {
    val rCosLat = radius * cos(lat)
    return TerseVector (
        rCosLat * cos(lon),
        rCosLat * sin(lon),
        radius * sin(lat)
    )
}

private fun calcVsop(model: VsopModel, time: AstroTime): AstroVector {
    val t = time.julianMillennia()

    // Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates.
    val lon = vsopFormulaCalc(model.lon, t, true)
    val lat = vsopFormulaCalc(model.lat, t, false)
    val rad = vsopFormulaCalc(model.rad, t, false)

    // Convert ecliptic spherical coordinates to ecliptic Cartesian coordinates.
    val eclip = vsopSphereToRect(lon, lat, rad)

    // Convert ecliptic Cartesian coordinates to equatorial Cartesian coordinates.
    // Also convert from TerseVector (coordinates only) to AstroVector (coordinates + time).
    return vsopRotate(eclip).toAstroVector(time)
}

private fun vsopDerivCalc(formula: VsopFormula, t: Double): Double {
    var tpower = 1.0        // t^s
    var dpower = 0.0        // t^(s-1)
    var deriv = 0.0
    var s: Int = 0
    for (series in formula.series) {
        var sinSum = 0.0
        var cosSum = 0.0
        for (term in series.term) {
            val angle = term.phase + (term.frequency * t)
            sinSum += term.amplitude * (term.frequency * sin(angle))
            if (s > 0)
                cosSum += term.amplitude * cos(angle)
        }
        deriv += (s * dpower * cosSum) - (tpower * sinSum)
        dpower = tpower
        tpower *= t
        ++s
    }
    return deriv
}

internal class BodyState(
    val tt: Double,
    val r: TerseVector,
    val v: TerseVector
) {
    fun decrement(other: BodyState) {
        r.decrement(other.r)
        v.decrement(other.v)
    }
}

private fun calcVsopPosVel(model: VsopModel, tt: Double): BodyState {
    val t = tt / DAYS_PER_MILLENNIUM

    // Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates.
    val lon = vsopFormulaCalc(model.lon, t, true)
    val lat = vsopFormulaCalc(model.lat, t, false)
    val rad = vsopFormulaCalc(model.rad, t, false)

    val eclipPos = vsopSphereToRect(lon, lat, rad)

    // Calculate derivatives of spherical coordinates with respect to time.
    val dlon = vsopDerivCalc(model.lon, t)      // dlon = d(lon) / dt
    val dlat = vsopDerivCalc(model.lat, t)      // dlat = d(lat) / dt
    val drad = vsopDerivCalc(model.rad, t)      // drad = d(rad) / dt

    // Use spherical coords and spherical derivatives to calculate
    // the velocity vector in rectangular coordinates.

    val coslon = cos(lon)
    val sinlon = sin(lon)
    val coslat = cos(lat)
    val sinlat = sin(lat)

    val vx = (
        + (drad * coslat * coslon)
        - (rad * sinlat * coslon * dlat)
        - (rad * coslat * sinlon * dlon)
    )

    val vy = (
        + (drad * coslat * sinlon)
        - (rad * sinlat * sinlon * dlat)
        + (rad * coslat * coslon * dlon)
    )

    val vz = (
        + (drad * sinlat)
        + (rad * coslat * dlat)
    )

    // Convert speed units from [AU/millennium] to [AU/day].
    val eclipVel = TerseVector(
        vx / DAYS_PER_MILLENNIUM,
        vy / DAYS_PER_MILLENNIUM,
        vz / DAYS_PER_MILLENNIUM
    )

    // Rotate the vectors from ecliptic to equatorial coordinates.
    val equPos = vsopRotate(eclipPos)
    val equVel = vsopRotate(eclipVel)
    return BodyState(tt, equPos, equVel)
}

private fun vsopModel(body: Body): VsopModel =
    body.vsopModel ?: throw InvalidBodyException(body)

private fun vsopHelioVector(body: Body, time: AstroTime) =
    calcVsop(vsopModel(body), time)

//---------------------------------------------------------------------------------------
// Geocentric Moon

/**
 * Emulate two-dimensional arrays from the Pascal programming language.
 *
 * The original Pascal source code used 2D arrays
 * with negative lower bounds on the indices. This is trivial
 * in Pascal, but most other languages use 0-based arrays.
 * This wrapper class adapts the original algorithm to Kotlin's
 * zero-based arrays.
 */
private class PascalArray2(
    val xmin: Int,
    val xmax: Int,
    val ymin: Int,
    val ymax: Int
) {
    private val array = Array<DoubleArray>((xmax- xmin) + 1) { _ -> DoubleArray((ymax - ymin) + 1) }
    operator fun get(x: Int, y: Int) = array[x - xmin][y - ymin]
    operator fun set(x: Int, y: Int, v: Double) { array[x - xmin][y - ymin] = v }
}

private class MoonContext(time: AstroTime) {
    // Variable names inside this class do not follow coding style on purpose.
    // They reflect the exact names from the original Pascal source code,
    // for ease of porting and consistent code generation.
    private var T: Double
    private var DGAM = 0.0
    private var DLAM: Double
    private var N = 0.0
    private var GAM1C: Double
    private var SINPI: Double
    private var L0: Double
    private var L: Double
    private var LS: Double
    private var F: Double
    private var D: Double
    private var DL0 = 0.0
    private var DL = 0.0
    private var DLS = 0.0
    private var DF = 0.0
    private var DD = 0.0
    private var DS = 0.0
    private val CO = PascalArray2(-6, 6, 1, 4)
    private val SI = PascalArray2(-6, 6, 1, 4)
    private val ARC = 3600.0 * RAD2DEG      // arcseconds per radian

    init {
        T = time.julianCenturies()
        val T2 = T*T
        DLAM = 0.0
        DS = 0.0
        GAM1C = 0.0
        SINPI = 3422.7
        longPeriodic()
        L0 = PI2 * frac(0.60643382 + (1336.85522467 * T) - (0.00000313 * T2)) + (DL0 / ARC)
        L  = PI2 * frac(0.37489701 + (1325.55240982 * T) + (0.00002565 * T2)) + (DL  / ARC)
        LS = PI2 * frac(0.99312619 + (  99.99735956 * T) - (0.00000044 * T2)) + (DLS / ARC)
        F  = PI2 * frac(0.25909118 + (1342.22782980 * T) - (0.00000892 * T2)) + (DF  / ARC)
        D  = PI2 * frac(0.82736186 + (1236.85308708 * T) - (0.00000397 * T2)) + (DD  / ARC)
        for (I in 1..4) {
            var ARG: Double
            var MAX: Int
            var FAC: Double
            when (I) {
                1    -> { ARG=L;  MAX=4; FAC=1.000002208               }
                2    -> { ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T }
                3    -> { ARG=F;  MAX=4; FAC=1.000002708+139.978*DGAM  }
                else -> { ARG=D;  MAX=6; FAC=1.0;                      }
            }
            CO[0,I] = 1.0
            CO[1,I] = cos(ARG) * FAC
            SI[0,I] = 0.0
            SI[1,I] = sin(ARG) * FAC

            for (J in 2..MAX) {
                val c1 = CO[J-1,I]
                val s1 = SI[J-1,I]
                val c2 = CO[1,I]
                val s2 = SI[1,I]
                CO[J,I] = c1*c2 - s1*s2
                SI[J,I] = s1*c2 + c1*s2
            }

            for (J in 1..MAX) {
                CO[-J,I] =  CO[J,I]
                SI[-J,I] = -SI[J,I]
            }
        }
    }

    /**
     * Sine of an angle `phi` expressed in revolutions, not radians or degrees.
     */
    private fun sine(phi: Double) = sin(PI2 * phi)

    private fun frac(x: Double) = x - floor(x)

    private fun longPeriodic() {
        val S1 = sine(0.19833 + (0.05611 * T))
        val S2 = sine(0.27869 + (0.04508 * T))
        val S3 = sine(0.16827 - (0.36903 * T))
        val S4 = sine(0.34734 - (5.37261 * T))
        val S5 = sine(0.10498 - (5.37899 * T))
        val S6 = sine(0.42681 - (0.41855 * T))
        val S7 = sine(0.14943 - (5.37511 * T))

        DL0 = ( 0.84 * S1) + (0.31 * S2) + (14.27 * S3) + (7.26 * S4) + (0.28 * S5) + (0.24 * S6)
        DL  = ( 2.94 * S1) + (0.31 * S2) + (14.27 * S3) + (9.34 * S4) + (1.12 * S5) + (0.83 * S6)
        DLS = (-6.40 * S1) - (1.89 * S6)
        DF  = ( 0.21 * S1) + (0.31 * S2) + (14.27 * S3) - (88.70*S4) - (15.30 * S5) + (0.24 * S6) - (1.86 * S7)
        DD  = DL0 - DLS
        DGAM = (
            -3332.0e-9 * sine(0.59734 - (5.37261 * T))
             -539.0e-9 * sine(0.35498 - (5.37899 * T))
              -64.0e-9 * sine(0.39943 - (5.37511 * T))
        )
    }

    // This is an ugly hack for efficiency.
    // Kotlin does not allow passing variables by reference.
    // Because this function is called *millions* of times
    // during the unit test, we can't afford to be allocating
    // tiny arrays or classes to pass a tuple (xTerm, yTerm) by reference.
    // It is much more efficient to keep (xTerm, yTerm) at object scope,
    // and allow term() to overwrite their values.
    // Callers know that each time they call term(), they must
    // access the output through (xTerm, yTerm) before calling term() again.

    var xTerm = Double.NaN
    var yTerm = Double.NaN
    private fun term(p: Int, q: Int, r: Int, s: Int) {
        xTerm = 1.0
        yTerm = 0.0
        if (p != 0) addTheta(CO[p, 1], SI[p, 1])
        if (q != 0) addTheta(CO[q, 2], SI[q, 2])
        if (r != 0) addTheta(CO[r, 3], SI[r, 3])
        if (s != 0) addTheta(CO[s, 4], SI[s, 4])
    }

    private fun addTheta(c2: Double, s2: Double) {
        val c1 = xTerm
        val s1 = yTerm
        xTerm = c1*c2 - s1*s2
        yTerm = s1*c2 + c1*s2
    }

    internal fun addSol(
        coeffl: Double,
        coeffs: Double,
        coeffg: Double,
        coeffp: Double,
        p: Int,
        q: Int,
        r: Int,
        s: Int
    ) {
        term(p, q, r, s)
        DLAM  += coeffl * yTerm
        DS    += coeffs * yTerm
        GAM1C += coeffg * xTerm
        SINPI += coeffp * xTerm
    }

    internal fun addn(coeffn: Double, p: Int, q: Int, r: Int, s: Int) {
        term(p, q, r, s)
        N += (coeffn * yTerm)
    }

    private fun solarN() {
        N = 0.0
        addn(-526.069,  0, 0, 1, -2)
        addn(  -3.352,  0, 0, 1, -4)
        addn( +44.297, +1, 0, 1, -2)
        addn(  -6.000, +1, 0, 1, -4)
        addn( +20.599, -1, 0, 1,  0)
        addn( -30.598, -1, 0, 1, -2)
        addn( -24.649, -2, 0, 1,  0)
        addn(  -2.000, -2, 0, 1, -2)
        addn( -22.571,  0,+1, 1, -2)
        addn( +10.985,  0,-1, 1, -2)
    }

    private fun planetary() {
        DLAM += (
            +0.82*sine(0.7736   -62.5512*T) + 0.31*sine(0.0466  -125.1025*T)
            +0.35*sine(0.5785   -25.1042*T) + 0.66*sine(0.4591 +1335.8075*T)
            +0.64*sine(0.3130   -91.5680*T) + 1.14*sine(0.1480 +1331.2898*T)
            +0.21*sine(0.5918 +1056.5859*T) + 0.44*sine(0.5784 +1322.8595*T)
            +0.24*sine(0.2275    -5.7374*T) + 0.28*sine(0.2965    +2.6929*T)
            +0.33*sine(0.3132    +6.3368*T)
        )
    }

    internal fun calcMoon(): Spherical {
        addSolarTerms(this)
        solarN()
        planetary()
        val S = F + DS/ARC
        val latSeconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*sin(S)-6.24*sin(3*S) + N
        return Spherical(
            latSeconds / 3600.0,
            360.0 * frac((L0+DLAM/ARC) / PI2),
            (ARC * EARTH_EQUATORIAL_RADIUS_AU) / (0.999953253 * SINPI)
        )
    }
}

//---------------------------------------------------------------------------------------
// Pluto/SSB gravitation simulator

private class BodyGravCalc(
    val tt: Double,         // J2000 terrestrial time [days]
    val r: TerseVector,     // position [au]
    val v: TerseVector,     // velocity [au/day]
    val a: TerseVector      // acceleration [au/day]
)

private class MajorBodies(
    val sun:        BodyState,
    val jupiter:    BodyState,
    val saturn:     BodyState,
    val uranus:     BodyState,
    val neptune:    BodyState
) {
    fun acceleration(smallPos: TerseVector): TerseVector = (
        accelerationIncrement(smallPos, SUN_GM,     sun.r    ) +
        accelerationIncrement(smallPos, JUPITER_GM, jupiter.r) +
        accelerationIncrement(smallPos, SATURN_GM,  saturn.r ) +
        accelerationIncrement(smallPos, URANUS_GM,  uranus.r ) +
        accelerationIncrement(smallPos, NEPTUNE_GM, neptune.r)
    )

    private fun accelerationIncrement(smallPos: TerseVector, gm: Double, majorPos: TerseVector): TerseVector {
        val delta = majorPos - smallPos
        val r2 = delta.quadrature()
        return (gm / (r2 * sqrt(r2))) * delta
    }
}

private fun updatePosition(dt: Double, r: TerseVector, v: TerseVector, a: TerseVector) =
    TerseVector(
        r.x + dt*(v.x + dt*a.x/2.0),
        r.y + dt*(v.y + dt*a.y/2.0),
        r.z + dt*(v.z + dt*a.z/2.0),
    )

private fun updateVelocity(dt: Double, v: TerseVector, a: TerseVector) =
    TerseVector(
        v.x + dt*a.x,
        v.y + dt*a.y,
        v.z + dt*a.z
    )

private fun adjustBarycenterPosVel(ssb: BodyState, tt: Double, body: Body, planetGm: Double): BodyState {
    val shift = planetGm / (planetGm + SUN_GM)
    val planet = calcVsopPosVel(vsopModel(body), tt)
    ssb.r.increment(shift * planet.r)
    ssb.v.increment(shift * planet.v)
    return planet
}

private fun majorBodyBary(tt: Double): MajorBodies {
    // Calculate the state vector of the Solar System Barycenter (SSB).
    val ssb = BodyState(tt, TerseVector.zero(), TerseVector.zero())
    val jupiter = adjustBarycenterPosVel(ssb, tt, Body.Jupiter, JUPITER_GM)
    val saturn  = adjustBarycenterPosVel(ssb, tt, Body.Saturn,  SATURN_GM )
    val uranus  = adjustBarycenterPosVel(ssb, tt, Body.Uranus,  URANUS_GM )
    val neptune = adjustBarycenterPosVel(ssb, tt, Body.Neptune, NEPTUNE_GM)

    // Convert planet state vectors from heliocentric to barycentric.
    jupiter.decrement(ssb)
    saturn.decrement(ssb)
    uranus.decrement(ssb)
    neptune.decrement(ssb)

    // Convert the heliocentric SSB to a barycentric Sun state.
    val sun = BodyState(tt, -ssb.r, -ssb.v)

    return MajorBodies(sun, jupiter, saturn, uranus, neptune)
}

private class GravSim(
    val bary: MajorBodies,
    val grav: BodyGravCalc
)

private fun simulateGravity(
    tt2: Double,            // a target time to be calculated (before or after current time tt1)
    calc1: BodyGravCalc     // [pos, vel, acc] of the simulated body at tt1
): GravSim {
    val dt = tt2 - calc1.tt

    // Calculate where the major bodies (Sun, Jupiter...Neptune) will be at tt2.
    val bary = majorBodyBary(tt2)

    // Estimate the position of the small body as if the current
    // acceleration applies across the whole time interval.
    val approxPos = updatePosition(dt, calc1.r, calc1.v, calc1.a)

    // Calculate the average acceleration of the endpoints.
    // This becomes our estimate of the mean effective acceleration
    // over the whole time interval.
    val meanAcc = bary.acceleration(approxPos).mean(calc1.a)

    // Refine the estimates of [pos, vel, acc] at tt2 using the mean acceleration.
    val pos = updatePosition(dt, calc1.r, calc1.v, meanAcc)
    val vel = updateVelocity(dt, calc1.v, meanAcc)
    val acc = bary.acceleration(pos)

    val grav = BodyGravCalc(tt2, pos, vel, acc)
    return GravSim(bary, grav)
}


private fun clampIndex(frac: Double, nsteps: Int): Int {
    val index = frac.toInt()
    return when {
        index < 0 -> 0
        index >= nsteps -> nsteps - 1
        else -> index
    }
}

private fun gravFromState(state: BodyState): GravSim {
    val bary = majorBodyBary(state.tt)
    val r = state.r + bary.sun.r
    val v = state.v + bary.sun.v
    val a = bary.acceleration(r)
    val grav = BodyGravCalc(state.tt, r, v, a)
    return GravSim(bary, grav)
}

private fun getPlutoSegment(tt: Double): List<BodyGravCalc>? {
    if (tt < plutoStateTable[0].tt || tt > plutoStateTable[PLUTO_NUM_STATES-1].tt)
        return null     // Don't bother calculating a segment. Let the caller crawl backward/forward to this time

    val segIndex = clampIndex((tt - plutoStateTable[0].tt) / PLUTO_TIME_STEP, PLUTO_NUM_STATES-1)
    synchronized (plutoCache) {
        var list = plutoCache.get(segIndex)
        if (list == null) {
            val seg = ArrayList<BodyGravCalc>()

            // The first endpoint is exact.
            var sim = gravFromState(plutoStateTable[segIndex])
            seg.add(sim.grav)

            // Simulate forwards from the lower time bound.
            var steptt = sim.grav.tt
            for (i in 1 until PLUTO_NSTEPS-1) {
                steptt += PLUTO_DT
                sim = simulateGravity(steptt, sim.grav)
                seg.add(sim.grav)
            }

            // The last endpoint is exact.
            sim = gravFromState(plutoStateTable[segIndex + 1])
            seg.add(sim.grav)

            // Simulate backwards from the upper time bound.
            val backward = ArrayList<BodyGravCalc>()
            backward.add(sim.grav)

            steptt = sim.grav.tt
            for (i in (PLUTO_NSTEPS-2) downTo 1) {
                steptt -= PLUTO_DT
                sim = simulateGravity(steptt, sim.grav)
                backward.add(sim.grav)
            }

            backward.add(seg[0])
            val reverse = backward.reversed()

            // Fade-mix the two series so that there are no discontinuities.
            for (i in (PLUTO_NSTEPS-2) downTo 1) {
                val ramp = i.toDouble() / (PLUTO_NSTEPS - 1)
                seg[i].r.mix(ramp, reverse[i].r)
                seg[i].v.mix(ramp, reverse[i].v)
                seg[i].a.mix(ramp, reverse[i].a)
            }

            list = seg.toList()
            plutoCache.set(segIndex, list)
        }
        return list
    }
}

private fun calcPlutoOneWay(
    initState: BodyState,
    targetTt: Double,
    dt: Double
) : GravSim {
    var sim = gravFromState(initState)
    val n: Int = ceil((targetTt - sim.grav.tt) / dt).toInt()
    for (i in 0 until n) {
        val tt = if (i+1 == n) targetTt else (sim.grav.tt + dt)
        sim = simulateGravity(tt, sim.grav)
    }
    return sim
}

private fun calcPluto(time: AstroTime, helio: Boolean): StateVector {
    val seg = getPlutoSegment(time.tt)
    var calc: BodyGravCalc
    var bary: MajorBodies? = null
    if (seg == null) {
        // The target time is outside the year range 0000..4000.
        // Calculate it by crawling backward from 0000 or forward from 4000.
        // FIXFIXFIX - This is super slow. Could optimize this with extra caching if needed.
        val sim =
            if (time.tt < plutoStateTable[0].tt)
                calcPlutoOneWay(plutoStateTable[0], time.tt, (-PLUTO_DT).toDouble())
            else
                calcPlutoOneWay(plutoStateTable[PLUTO_NUM_STATES-1], time.tt, (+PLUTO_DT).toDouble())

        calc = sim.grav
        bary = sim.bary
    } else {
        val left = clampIndex((time.tt - seg[0].tt) / PLUTO_DT, PLUTO_NSTEPS-1)
        val s1 = seg[left]
        val s2 = seg[left+1]

        // Find mean acceleration vector over the interval.
        val acc: TerseVector = s1.a.mean(s2.a)

        // Use Newtonian mechanics to extrapolate away from t1 in the positive time direction.
        val ra: TerseVector = updatePosition(time.tt - s1.tt, s1.r, s1.v, acc)
        val va: TerseVector = updateVelocity(time.tt - s1.tt, s1.v, acc)

        // Use Newtonian mechanics to extrapolate away from t2 in the negative time direction.
        val rb: TerseVector = updatePosition(time.tt - s2.tt, s2.r, s2.v, acc)
        val vb: TerseVector = updateVelocity(time.tt - s2.tt, s2.v, acc)

        // Use fade in/out idea to blend the two position estimates.
        val ramp = (time.tt - s1.tt)/PLUTO_DT
        calc = BodyGravCalc(
            time.tt,
            (1.0 - ramp)*ra + ramp*rb,      // ramp/mix the position vector
            (1.0 - ramp)*va + ramp*vb,      // ramp/mix the velocity vector
            TerseVector.zero()              // the acceleration isn't used
        )
    }

    if (helio) {
        // Lazy evaluate the Solar System Barycenter state if needed.
        if (bary == null)
            bary = majorBodyBary(time.tt)

        // Convert barycentric vectors to heliocentric vectors.
        calc.r.decrement(bary.sun.r)
        calc.v.decrement(bary.sun.v)
    }

    return StateVector(
        calc.r.x, calc.r.y, calc.r.z,
        calc.v.x, calc.v.y, calc.v.z,
        time
    )
}

//---------------------------------------------------------------------------------------

internal class JupiterMoon (
    val mu: Double,
    val al0: Double,
    val al1: Double,
    val a: Array<VsopTerm>,
    val l: Array<VsopTerm>,
    val z: Array<VsopTerm>,
    val zeta: Array<VsopTerm>
)


private fun jupiterMoonElemToPv(
    time: AstroTime,
    mu: Double,
    A: Double,
    AL: Double,
    K: Double,
    H: Double,
    Q: Double,
    P: Double
): StateVector {
    // Translation of FORTRAN subroutine ELEM2PV from:
    // https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

    val AN = sqrt(mu / (A*A*A))

    var CE: Double
    var SE: Double
    var DE: Double
    var EE = AL + K*sin(AL) - H*cos(AL)
    do {
        CE = cos(EE)
        SE = sin(EE)
        DE = (AL - EE + K*SE - H*CE) / (1.0 - K*CE - H*SE)
        EE += DE
    } while (DE.absoluteValue >= 1.0e-12)

    CE = cos(EE)
    SE = sin(EE)
    val DLE = H*CE - K*SE
    val RSAM1 = -K*CE - H*SE
    val ASR = 1.0/(1.0 + RSAM1)
    val PHI = sqrt(1.0 - K*K - H*H)
    val PSI = 1.0/(1.0 + PHI)
    val X1 = A*(CE - K - PSI*H*DLE)
    val Y1 = A*(SE - H + PSI*K*DLE)
    val VX1 = AN*ASR*A*(-SE - PSI*H*RSAM1)
    val VY1 = AN*ASR*A*(+CE + PSI*K*RSAM1)
    val F2 = 2.0*sqrt(1.0 - Q*Q - P*P)
    val P2 = 1.0 - 2.0*P*P
    val Q2 = 1.0 - 2.0*Q*Q
    val PQ = 2.0*P*Q

    return StateVector(
        X1*P2 + Y1*PQ,
        X1*PQ + Y1*Q2,
        (Q*Y1 - X1*P)*F2,
        VX1*P2 + VY1*PQ,
        VX1*PQ + VY1*Q2,
        (Q*VY1 - VX1*P)*F2,
        time
    )
}


private fun calcJupiterMoon(time: AstroTime, m: JupiterMoon): StateVector {
    // This is a translation of FORTRAN code by Duriez, Lainey, and Vienne:
    // https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

    val t = time.tt + 18262.5     // number of days since 1950-01-01T00:00:00Z

    // Calculate 6 orbital elements at the given time t.
    var elem0 = 0.0
    for (term in m.a)
        elem0 += term.amplitude * cos(term.phase + (t * term.frequency))

    var elem1 = m.al0 + (t * m.al1)
    for (term in m.l)
        elem1 += term.amplitude * sin(term.phase + (t * term.frequency))

    elem1 %= PI2
    if (elem1 < 0)
        elem1 += PI2

    var elem2 = 0.0
    var elem3 = 0.0
    for (term in m.z) {
        val arg = term.phase + (t * term.frequency)
        elem2 += term.amplitude * cos(arg)
        elem3 += term.amplitude * sin(arg)
    }

    var elem4 = 0.0
    var elem5 = 0.0
    for (term in m.zeta) {
        var arg = term.phase + (t * term.frequency)
        elem4 += term.amplitude * cos(arg)
        elem5 += term.amplitude * sin(arg)
    }

    // Convert the oribital elements into position vectors in the Jupiter equatorial system (JUP).
    val state = jupiterMoonElemToPv(time, m.mu, elem0, elem1, elem2, elem3, elem4, elem5)

    // Re-orient position and velocity vectors from Jupiter-equatorial (JUP) to Earth-equatorial in J2000 (EQJ).
    return rotationJupEqj.rotate(state)
}


//---------------------------------------------------------------------------------------


internal fun terrestrialTime(ut: Double): Double = ut + deltaT(ut) / SECONDS_PER_DAY

private val epoch2000 = AstroTime(0.0)

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
            10583.6 - (1014.41 * u) + (33.78311 * u2) - (5.952053 * u3) - (0.1798452 * u4) + (0.022174192 * u5) + (0.0090316521 * u6)
        }
        y < 1600.0 -> {
            u = (y - 1000) / 100
            u2 = u * u; u3 = u * u2; u4 = u2 * u2; u5 = u2 * u3; u6 = u3 * u3
            1574.2 - (556.01 * u) + (71.23472 * u2) + (0.319781 * u3) - (0.8503463 * u4) - (0.005050998 * u5) + (0.0083572073 * u6)
        }
        y < 1700.0 -> {
            u = y - 1600
            u2 = u * u; u3 = u * u2
            120.0 - (0.9808 * u) - (0.01532 * u2) + (u3 / 7129)
        }
        y < 1800.0 -> {
            u = y - 1700
            u2 = u * u; u3 = u * u2; u4 = u2 * u2
            8.83 + (0.1603 * u) - (0.0059285 * u2) + (0.00013336 * u3) - (u4 / 1174000)
        }
        y < 1860.0 -> {
            u = y - 1800
            u2 = u * u; u3 = u * u2; u4 = u2 * u2; u5 = u2 * u3; u6 = u3 * u3; u7 = u3 * u4
            13.72 - (0.332447 * u) + (0.0068612 * u2) + (0.0041116 * u3) - (0.00037436 * u4) + (0.0000121272 * u5) - (0.0000001699 * u6) + (0.000000000875 * u7)
        }
        y < 1900.0 -> {
            u = y - 1860
            u2 = u * u; u3 = u * u2; u4 = u2 * u2; u5 = u2 * u3
            7.62 + (0.5737 * u) - (0.251754 * u2) + (0.01680668 * u3) - (0.0004473624 * u4) + (u5 / 233174)
        }
        y < 1920.0 -> {
            u = y - 1900
            u2 = u * u; u3 = u * u2; u4 = u2 * u2
            -2.79 + (1.494119 * u) - (0.0598939 * u2) + (0.0061966 * u3) - (0.000197 * u4)
        }
        y < 1941.0 -> {
            u = y - 1920
            u2 = u * u; u3 = u * u2
            21.20 + (0.84493 * u) - (0.076100 * u2) + (0.0020936 * u3)
        }
        y < 1961 -> {
            u = y - 1950
            u2 = u * u; u3 = u * u2
            29.07 + (0.407 * u) - (u2 / 233) + (u3 / 2547)
        }
        y < 1986.0 -> {
            u = y - 1975
            u2 = u * u; u3 = u * u2
            45.45 + (1.067 * u) - (u2 / 260) - (u3 / 718)
        }
        y < 2005 -> {
            u = y - 2000
            u2 = u * u; u3 = u * u2; u4 = u2 * u2; u5 = u2 * u3
            63.86 + (0.3345 * u) - (0.060374 * u2) + (0.0017275 * u3) + (0.000651814 * u4) + (0.00002373599 * u5)
        }
        y < 2050 -> {
            u = y - 2000
            62.92 + (0.32217 * u) + (0.005589 * u * u)
        }
        y < 2150 -> {
            u = (y - 1820) / 100
            -20 + (32 * u * u) - (0.5628 * (2150 - y))
        }
        else -> {   // all years after 2150
            u = (y - 1820) / 100
            -20 + (32 * u * u)
        }
    }
}

internal fun universalTime(tt: Double): Double {
    // This is the inverse function of terrestrialTime.
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

/**
 * Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.
 *
 * Given an altitude angle and a refraction option, calculates
 * the amount of "lift" caused by atmospheric refraction.
 * This is the number of degrees higher in the sky an object appears
 * due to the lensing of the Earth's atmosphere.
 *
 * @param refraction
 *      The option selecting which refraction correction to use.
 *      If `Refraction.Normal`, uses a well-behaved refraction model that works well for
 *      all valid values (-90 to +90) of `altitude`.
 *      If `Refraction.JplHor`, this function returns a compatible value with the JPL Horizons tool.
 *      If any other value, including `Refraction.None`, this function returns 0.
 *
 * @param altitude
 *      An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90.
 */
fun refractionAngle(refraction: Refraction, altitude: Double): Double {
    if (altitude < -90.0 || altitude > +90.0)
        return 0.0     // no attempt to correct an invalid altitude

    var angle: Double
    if (refraction == Refraction.Normal || refraction == Refraction.JplHor) {
        // http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
        // JPL Horizons says it uses refraction algorithm from
        // Meeus "Astronomical Algorithms", 1991, p. 101-102.
        // I found the following Go implementation:
        // https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
        // This is a translation from the function "Saemundsson" there.
        // I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
        // This is important because the tangent formula below goes crazy near hd = -5.11.
        val hd = altitude.coerceAtLeast(-1.0)
        angle = (1.02 / dtan(hd + 10.3/(hd + 5.11))) / 60.0

        if (refraction == Refraction.Normal && altitude < -1.0) {
            // In "normal" mode we gradually reduce refraction toward the nadir
            // so that we never get an altitude angle less than -90 degrees.
            // When horizon angle is -1 degrees, the factor is exactly 1.
            // As altitude approaches -90 (the nadir), the fraction approaches 0 linearly.
            angle *= (altitude + 90.0) / 89.0
        }
    } else {
        // No refraction, or the refraction option is invalid.
        angle = 0.0
    }
    return angle
}

/**
 * Calculates the inverse of an atmospheric refraction angle.
 *
 * Given an observed altitude angle that includes atmospheric refraction,
 * calculates the negative angular correction to obtain the unrefracted
 * altitude. This is useful for cases where observed horizontal
 * coordinates are to be converted to another orientation system,
 * but refraction first must be removed from the observed position.
 *
 * @param refraction
 *      The option selecting which refraction correction to use.
 *
 * @param bentAltitude
 *      The apparent altitude that includes atmospheric refraction.
 *
 * @param
 *      The angular adjustment in degrees to be added to the
 *      altitude angle to remove atmospheric lensing.
 *      This will be less than or equal to zero.
 */
fun inverseRefractionAngle(refraction: Refraction, bentAltitude: Double): Double {
    if (bentAltitude < -90.0 || bentAltitude > +90.0)
        return 0.0     // no attempt to correct an invalid altitude

    // Find the pre-adjusted altitude whose refraction correction leads to 'altitude'.
    var altitude = bentAltitude - refractionAngle(refraction, bentAltitude)
    while (true) {
        // See how close we got.
        var diff = (altitude + refractionAngle(refraction, altitude)) - bentAltitude
        if (diff.absoluteValue < 1.0e-14)
            return altitude - bentAltitude
        altitude -= diff
    }
}

/**
 * Returns the product of mass and universal gravitational constant of a Solar System body.
 *
 * For problems involving the gravitational interactions of Solar System bodies,
 * it is helpful to know the product GM, where G = the universal gravitational constant
 * and M = the mass of the body. In practice, GM is known to a higher precision than
 * either G or M alone, and thus using the product results in the most accurate results.
 * This function returns the product GM in the units au^3/day^2.
 * The values come from page 10 of a
 * [JPL memorandum regarding the DE405/LE405 ephemeris](https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf).
 *
 * @param body
 *      The body for which to find the GM product.
 *      Allowed to be the Sun, Moon, EMB (Earth/Moon Barycenter), or any planet.
 *      Any other value will cause an exception to be thrown.
 *
 * @return The mass product of the given body in au^3/day^2.
 */
fun massProduct(body: Body): Double =
    body.massProduct ?: throw InvalidBodyException(body)

private enum class PrecessDirection {
    From2000,
    Into2000,
}

private fun precessionRot(time: AstroTime, dir: PrecessDirection): RotationMatrix {
    val t = time.julianCenturies()
    val eps0 = 84381.406

    val psia   = (((((-    0.0000000951  * t
                      +    0.000132851 ) * t
                      -    0.00114045  ) * t
                      -    1.0790069   ) * t
                      + 5038.481507    ) * t) * ASEC2RAD

    val omegaa = (((((+    0.0000003337  * t
                      -    0.000000467 ) * t
                      -    0.00772503  ) * t
                      +    0.0512623   ) * t
                      -    0.025754    ) * t + eps0) * ASEC2RAD

    val chia   = (((((-    0.0000000560  * t
                      +    0.000170663 ) * t
                      -    0.00121197  ) * t
                      -    2.3814292   ) * t
                      +   10.556403    ) * t) * ASEC2RAD

    val sa = sin(eps0 * ASEC2RAD)
    val ca = cos(eps0 * ASEC2RAD)
    val sb = sin(-psia)
    val cb = cos(-psia)
    val sc = sin(-omegaa)
    val cc = cos(-omegaa)
    val sd = sin(chia)
    val cd = cos(chia)

    val xx =  cd*cb - sb*sd*cc
    val yx =  cd*sb*ca + sd*cc*cb*ca - sa*sd*sc
    val zx =  cd*sb*sa + sd*cc*cb*sa + ca*sd*sc
    val xy = -sd*cb - sb*cd*cc
    val yy = -sd*sb * ca + cd*cc*cb*ca - sa*cd*sc
    val zy = -sd*sb * sa + cd*cc*cb*sa + ca*cd*sc
    val xz =  sb*sc
    val yz = -sc*cb*ca - sa*cc
    val zz = -sc*cb*sa + cc*ca

    return when (dir) {
        // Perform rotation from other epoch to J2000.0.
        PrecessDirection.Into2000 ->
            RotationMatrix(
                xx, yx, zx,
                xy, yy, zy,
                xz, yz, zz
            )

        // Perform rotation from J2000.0 to other epoch.
        PrecessDirection.From2000 ->
            RotationMatrix(
                xx, xy, xz,
                yx, yy, yz,
                zx, zy, zz
            )
    }
}

private fun precession(pos: AstroVector, dir: PrecessDirection) =
    precessionRot(pos.t, dir).rotate(pos)

private fun precessionPosVel(state: StateVector, dir: PrecessDirection) =
    precessionRot(state.t, dir).rotate(state)

private class EarthTilt(
    val tt: Double,
    val dpsi: Double,
    val deps: Double,
    val ee: Double,
    val mobl: Double,
    val tobl: Double
)

private fun iau2000b(time: AstroTime) {
    // Adapted from the NOVAS C 3.1 function of the same name.
    // We cache Earth nutation angles `psi` and `eps` inside AstroTime for efficiency.
    // If nutation has not already been calculated, these values will be NaN.
    // Lazy-evaluate both angles.

    if (time.psi.isNaN()) {
        val t = time.julianCenturies()
        val el  = ((485868.249036 + t * 1717915923.2178) % ASEC360) * ASEC2RAD
        val elp = ((1287104.79305 + t * 129596581.0481)  % ASEC360) * ASEC2RAD
        val f   = ((335779.526232 + t * 1739527262.8478) % ASEC360) * ASEC2RAD
        val d   = ((1072260.70369 + t * 1602961601.2090) % ASEC360) * ASEC2RAD
        val om  = ((450160.398036 - t * 6962890.5431)    % ASEC360) * ASEC2RAD
        var dp = 0.0
        var de = 0.0
        for (i in 76 downTo 0) {
            val arg = (
                iauRow[i].nals0*el + iauRow[i].nals1*elp +
                iauRow[i].nals2*f + iauRow[i].nals3*d +
                iauRow[i].nals4*om
            ) % PI2
            val sarg = sin(arg)
            val carg = cos(arg)
            dp += (iauRow[i].cls0 + iauRow[i].cls1*t) * sarg + iauRow[i].cls2*carg
            de += (iauRow[i].cls3 + iauRow[i].cls4*t) * carg + iauRow[i].cls5*sarg
        }
        time.psi = -0.000135 + (dp * 1.0e-7)
        time.eps = +0.000388 + (de * 1.0e-7)
    }
}

private fun meanObliquity(time: AstroTime): Double {
    val t = time.julianCenturies()
    val asec =
        ((((  -0.0000000434   * t
            -  0.000000576  ) * t
            +  0.00200340   ) * t
            -  0.0001831    ) * t
            - 46.836769     ) * t + 84381.406
    return asec / 3600
}

private fun earthTilt(time: AstroTime): EarthTilt {
    iau2000b(time)  // lazy-evaluate time.psi and time.eps
    val mobl = meanObliquity(time)
    val tobl = mobl + (time.eps / 3600)
    val ee = time.psi * dcos(mobl) / 15.0
    return EarthTilt(time.tt, time.psi, time.eps, ee, mobl, tobl)
}

private fun earthRotationAngle(time: AstroTime): Double {
    val thet1 = 0.7790572732640 + (0.00273781191135448 * time.ut)
    val thet3 = time.ut % 1.0
    val theta = 360.0 *((thet1 + thet3) % 1.0)
    return if (theta < 0.0) theta + 360.0 else theta
}

/**
 * Calculates Greenwich Apparent Sidereal Time (GAST).
 *
 * Given a date and time, this function calculates the rotation of the
 * Earth, represented by the equatorial angle of the Greenwich prime meridian
 * with respect to distant stars (not the Sun, which moves relative to background
 * stars by almost one degree per day).
 * This angle is called Greenwich Apparent Sidereal Time (GAST).
 * GAST is measured in sidereal hours in the half-open range [0, 24).
 * When GAST = 0, it means the prime meridian is aligned with the of-date equinox,
 * corrected at that time for precession and nutation of the Earth's axis.
 * In this context, the *equinox* is the direction in space where the Earth's
 * orbital plane (the ecliptic) intersects with the plane of the Earth's equator,
 * at the location on the Earth's orbit of the (seasonal) March equinox.
 * As the Earth rotates, GAST increases from 0 up to 24 sidereal hours,
 * then starts over at 0.
 * To convert to degrees, multiply the return value by 15.
 *
 * @param time
 *      The date and time for which to find GAST.
 *      As an optimization, this function caches the sideral time value in `time`,
 *      unless it has already been cached, in which case the cached value is reused.
 */
fun siderealTime(time: AstroTime): Double {
    if (time.st.isNaN()) {
        val t = time.julianCenturies()
        val eqeq = 15.0 * earthTilt(time).ee
        val theta = earthRotationAngle(time)
        val st = (eqeq + 0.014506 +
               (((( -    0.0000000368  * t
                   -    0.000029956  ) * t
                   -    0.00000044   ) * t
                   +    1.3915817    ) * t
                   + 4612.156534     ) * t)
        val gst = ((st/3600.0 + theta) % 360.0) / 15.0
        time.st = if (gst < 0.0) gst + 24.0 else gst
    }
    return time.st
}

/**
 * Calcluate the geocentric position and velocity of a topocentric observer.
 *
 * Given the geographic location of an observer and a time, calculate the position
 * and velocity of the observer, relative to the center of the Earth.
 * The state vector is expressed in the equator-of-date system (EQD).
 *
 * @param observer
 *      The latitude, longitude, and elevation of the observer.
 *
 * @param time
 *      The time of the observation.
 *
 * @return
 * An EQD state vector that holds the geocentric position and velocity
 * of the observer at the given time.
 */
private fun terra(observer: Observer, time: AstroTime): StateVector {
    val st = siderealTime(time)
    val phi = observer.latitude.degreesToRadians()
    val sinphi = sin(phi)
    val cosphi = cos(phi)
    val c = 1.0 / hypot(cosphi,  EARTH_FLATTENING * sinphi)
    val s = c * EARTH_FLATTENING_SQUARED
    val heightKm = observer.height / 1000.0
    val ach = (EARTH_EQUATORIAL_RADIUS_KM * c) + heightKm
    val ash = (EARTH_EQUATORIAL_RADIUS_KM * s) + heightKm
    val stlocl = (15.0*st + observer.longitude).degreesToRadians()
    val sinst = sin(stlocl)
    val cosst = cos(stlocl)

    return StateVector(
        ach * cosphi * cosst / KM_PER_AU,
        ach * cosphi * sinst / KM_PER_AU,
        ash * sinphi / KM_PER_AU,
        -(ANGVEL * 86400.0 / KM_PER_AU) * ach * cosphi * sinst,
        +(ANGVEL * 86400.0 / KM_PER_AU) * ach * cosphi * cosst,
        0.0,
        time
    )
}

/**
 * Calculate the geographic coordinates corresponding to a position vector.
 *
 * Given a geocentric position vector expressed in equator-of-date (EQD) coordinates,
 * this function calculates the latitude, longitude, and elevation of that location.
 * Note that the time `ovec.t` must be set correctly in order to determine the
 * Earth's rotation angle.
 *
 * This function is intended for positions known to be on or near the Earth's surface.
 *
 * @param ovec
 *      A geocentric position on or near the Earth's surface, in EQD coordinates.
 *
 * @return
 * The location on or near the Earth's surface corresponding to
 * the given position vector and time.
 */
private fun inverseTerra(ovec: AstroVector): Observer {
    var lonDeg: Double
    var latDeg: Double
    var heightKm: Double

    // Convert from AU to kilometers.
    val x = ovec.x * KM_PER_AU
    val y = ovec.y * KM_PER_AU
    val z = ovec.z * KM_PER_AU
    val p = hypot(x, y)
    if (p < 1.0e-6) {
        // Special case: within 1 millimeter of a pole!
        // Use arbitrary longitude, and latitude determined by polarity of z.
        lonDeg = 0.0
        latDeg = if (z > 0.0) +90.0 else -90.0
        // Elevation is calculated directly from z.
        heightKm = z.absoluteValue - EARTH_POLAR_RADIUS_KM
    } else {
        // Calculate exact longitude in the half-open range (-180, +180].
        val stlocl = atan2(y, x)
        lonDeg = longitudeOffset(stlocl.radiansToDegrees() - (15 * siderealTime(ovec.t)))
        // Numerically solve for exact latitude, using Newton's Method.
        val F = EARTH_FLATTENING_SQUARED
        // Start with initial latitude estimate, based on a spherical Earth.
        var lat = atan2(z, p)
        var c: Double
        var s: Double
        var denom: Double
        while (true) {
            // Calculate the error function W(lat).
            // We try to find the root of W, meaning where the error is 0.
            c = cos(lat)
            s = sin(lat)
            val factor = (F-1)*EARTH_EQUATORIAL_RADIUS_KM
            val c2 = c*c
            val s2 = s*s
            val radicand = c2 + F*s2
            denom = sqrt(radicand)
            val W = ((factor * s * c) / denom) - (z * c) + (p * s)
            if (W.absoluteValue < 1.0e-12)
                break  // The error is now negligible.
            // Error is still too large. Find the next estimate.
            // Calculate D = the derivative of W with respect to lat.
            val D = (factor * ((c2 - s2) / denom) - (s2 * c2 * (F - 1)/(factor * radicand))) + (z * s) + (p * c)
            lat -= (W / D)
        }
        // We now have a solution for the latitude in radians.
        latDeg = lat.radiansToDegrees()
        // Solve for exact height in kilometers.
        // There are two formulas I can use. Use whichever has the less risky denominator.
        val adjust = EARTH_EQUATORIAL_RADIUS_KM / denom
        heightKm =
            if (s.absoluteValue > c.absoluteValue)
                z/s - F*adjust
            else
                p/c - adjust
    }

    return Observer(latDeg, lonDeg, 1000.0 * heightKm)
}

private fun gyration(pos: AstroVector, dir: PrecessDirection) =
    when (dir) {
        PrecessDirection.Into2000 -> precession(nutation(pos, dir), dir)
        PrecessDirection.From2000 -> nutation(precession(pos, dir), dir)
    }

private fun gyrationPosVel(state: StateVector, dir: PrecessDirection) =
    when (dir) {
        PrecessDirection.Into2000 -> precessionPosVel(nutationPosVel(state, dir), dir)
        PrecessDirection.From2000 -> nutationPosVel(precessionPosVel(state, dir), dir)
    }

private fun geoPos(time: AstroTime, observer: Observer) =
    gyration(
        terra(observer, time).position(),
        PrecessDirection.Into2000
    )

private fun spin(angle: Double, pos: AstroVector): AstroVector {
    val cosang = dcos(angle)
    val sinang = dsin(angle)
    return AstroVector(
        +cosang*pos.x + sinang*pos.y,
        -sinang*pos.x + cosang*pos.y,
        pos.z,
        pos.t
    )
}

private fun nutationRot(time: AstroTime, dir: PrecessDirection): RotationMatrix {
    val tilt = earthTilt(time)
    val oblm = tilt.mobl.degreesToRadians()
    val oblt = tilt.tobl.degreesToRadians()
    val psi = tilt.dpsi * ASEC2RAD
    val cobm = cos(oblm)
    val sobm = sin(oblm)
    val cobt = cos(oblt)
    val sobt = sin(oblt)
    val cpsi = cos(psi)
    val spsi = sin(psi)

    val xx = cpsi
    val yx = -spsi * cobm
    val zx = -spsi * sobm
    val xy = spsi * cobt
    val yy = cpsi * cobm * cobt + sobm * sobt
    val zy = cpsi * sobm * cobt - cobm * sobt
    val xz = spsi * sobt
    val yz = cpsi * cobm * sobt - sobm * cobt
    val zz = cpsi * sobm * sobt + cobm * cobt

    return when (dir) {
        // Perform rotation from other epoch to J2000.0.
        PrecessDirection.Into2000 ->
            RotationMatrix(
                xx, yx, zx,
                xy, yy, zy,
                xz, yz, zz
            )

        // Perform rotation from J2000.0 to other epoch.
        PrecessDirection.From2000 ->
            RotationMatrix(
                xx, xy, xz,
                yx, yy, yz,
                zx, zy, zz
            )
    }
}

private fun nutation(pos: AstroVector, dir: PrecessDirection) =
    nutationRot(pos.t, dir).rotate(pos)

private fun nutationPosVel(state: StateVector, dir: PrecessDirection) =
    nutationRot(state.t, dir).rotate(state)

private fun eclipticToEquatorial(ecl: AstroVector): AstroVector {
    val obl = meanObliquity(ecl.t).degreesToRadians()
    val cosObl = cos(obl)
    val sinObl = sin(obl)
    return AstroVector(
        ecl.x,
        (ecl.y * cosObl) - (ecl.z * sinObl),
        (ecl.y * sinObl) + (ecl.z * cosObl),
        ecl.t
    )
}

/**
 * Given an equatorial vector, calculates equatorial angular coordinates.
 *
 * @vector
 *      A vector in an equatorial coordinate system.
 *
 * @return Angular coordinates expressed in the same equatorial system as `vector`.
 */
fun equatorFromVector(vector: AstroVector): Equatorial {
    val sphere = vector.toSpherical()
    return Equatorial(sphere.lon / 15.0, sphere.lat, sphere.dist, vector)
}


private fun earthRotationAxis(time: AstroTime): AxisInfo {
    // Unlike the other planets, we have a model of precession and nutation
    // for the Earth's axis that provides a north pole vector.
    // So calculate the vector first, then derive the (RA,DEC) angles from the vector.

    // Start with a north pole vector in equator-of-date coordinates: (0,0,1).
    val pos1 = AstroVector(0.0, 0.0, 1.0, time)

    // Convert the vector into J2000 coordinates to find the north pole direction.
    val pos2 = nutation(pos1, PrecessDirection.Into2000)
    val north = precession(pos2, PrecessDirection.Into2000)

    // Derive angular values: right ascension and declination.
    val equ = north.toEquatorial()

    // Use a modified version of the era() function that does not trim to 0..360 degrees.
    // This expression is also corrected to give the correct angle at the J2000 epoch.
    val spin = 190.41375788700253 + (360.9856122880876 * time.ut)

    return AxisInfo(equ.ra, equ.dec, spin, north)
}


/**
 * Calculates information about a body's rotation axis at a given time.
 *
 * Calculates the orientation of a body's rotation axis, along with
 * the rotation angle of its prime meridian, at a given moment in time.
 *
 * This function uses formulas standardized by the IAU Working Group
 * on Cartographics and Rotational Elements 2015 report, as described
 * in the following document:
 *
 * https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf
 *
 * See [AxisInfo] for more detailed information.
 *
 * @param body
 *      One of the following values:
 *      `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`,
 *      `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
 *
 * @param time
 *      The time at which to calculate the body's rotation axis.
 *
 * @return North pole orientation and body spin angle.
 */
fun rotationAxis(body: Body, time: AstroTime): AxisInfo {
    if (body == Body.Earth)
        return earthRotationAxis(time)

    val d = time.tt
    val T = time.julianCenturies()
    val ra: Double
    val dec: Double
    val w: Double
    when (body) {
        Body.Sun -> {
            ra = 286.13
            dec = 63.87
            w = 84.176 + (14.1844 * d)
        }

        Body.Mercury -> {
            ra = 281.0103 - (0.0328 * T)
            dec = 61.4155 - (0.0049 * T)
            w = (
                329.5988
                + (6.1385108 * d)
                + (0.01067257 * dsin((174.7910857 + 4.092335*d)))
                - (0.00112309 * dsin((349.5821714 + 8.184670*d)))
                - (0.00011040 * dsin((164.3732571 + 12.277005*d)))
                - (0.00002539 * dsin((339.1643429 + 16.369340*d)))
                - (0.00000571 * dsin((153.9554286 + 20.461675*d)))
            )
        }

        Body.Venus -> {
            ra = 272.76
            dec = 67.16
            w = 160.20 - (1.4813688 * d)
        }

        Body.Moon -> {
            // See page 8, Table 2 in:
            // https://astropedia.astrogeology.usgs.gov/alfresco/d/d/workspace/SpacesStore/28fd9e81-1964-44d6-a58b-fbbf61e64e15/WGCCRE2009reprint.pdf
            val E1  = 125.045 -  0.0529921*d
            val E2  = 250.089 -  0.1059842*d
            val E3  = 260.008 + 13.0120009*d
            val E4  = 176.625 + 13.3407154*d
            val E5  = 357.529 +  0.9856003*d
            val E6  = 311.589 + 26.4057084*d
            val E7  = 134.963 + 13.0649930*d
            val E8  = 276.617 +  0.3287146*d
            val E9  = 34.226  +  1.7484877*d
            val E10 = 15.134  -  0.1589763*d
            val E11 = 119.743 +  0.0036096*d
            val E12 = 239.961 +  0.1643573*d
            val E13 = 25.053  + 12.9590088*d

            ra = (
                269.9949 + 0.0031*T
                - 3.8787*dsin(E1)
                - 0.1204*dsin(E2)
                + 0.0700*dsin(E3)
                - 0.0172*dsin(E4)
                + 0.0072*dsin(E6)
                - 0.0052*dsin(E10)
                + 0.0043*dsin(E13)
            )

            dec = (
                66.5392 + 0.0130*T
                + 1.5419*dcos(E1)
                + 0.0239*dcos(E2)
                - 0.0278*dcos(E3)
                + 0.0068*dcos(E4)
                - 0.0029*dcos(E6)
                + 0.0009*dcos(E7)
                + 0.0008*dcos(E10)
                - 0.0009*dcos(E13)
            )

            w = (
                38.3213 + (13.17635815 - 1.4e-12*d)*d
                + 3.5610*dsin(E1)
                + 0.1208*dsin(E2)
                - 0.0642*dsin(E3)
                + 0.0158*dsin(E4)
                + 0.0252*dsin(E5)
                - 0.0066*dsin(E6)
                - 0.0047*dsin(E7)
                - 0.0046*dsin(E8)
                + 0.0028*dsin(E9)
                + 0.0052*dsin(E10)
                + 0.0040*dsin(E11)
                + 0.0019*dsin(E12)
                - 0.0044*dsin(E13)
            )
        }

        Body.Mars -> {
            ra = (
                317.269202 - 0.10927547*T
                + 0.000068 * dsin(198.991226 + 19139.4819985*T)
                + 0.000238 * dsin(226.292679 + 38280.8511281*T)
                + 0.000052 * dsin(249.663391 + 57420.7251593*T)
                + 0.000009 * dsin(266.183510 + 76560.6367950*T)
                + 0.419057 * dsin(79.398797 + 0.5042615*T)
            )

            dec = (
                54.432516 - 0.05827105*T
                + 0.000051 * dcos(122.433576 + 19139.9407476*T)
                + 0.000141 * dcos(43.058401 + 38280.8753272*T)
                + 0.000031 * dcos(57.663379 + 57420.7517205*T)
                + 0.000005 * dcos(79.476401 + 76560.6495004*T)
                + 1.591274 * dcos(166.325722 + 0.5042615*T)
            )

            w = (
                176.049863 + 350.891982443297*d
                + 0.000145 * dsin(129.071773 + 19140.0328244*T)
                + 0.000157 * dsin(36.352167 + 38281.0473591*T)
                + 0.000040 * dsin(56.668646 + 57420.9295360*T)
                + 0.000001 * dsin(67.364003 + 76560.2552215*T)
                + 0.000001 * dsin(104.792680 + 95700.4387578*T)
                + 0.584542 * dsin(95.391654 + 0.5042615*T)
            )
        }

        Body.Jupiter -> {
            val Ja = 99.360714  + 4850.4046*T
            val Jb = 175.895369 + 1191.9605*T
            val Jc = 300.323162 + 262.5475*T
            val Jd = 114.012305 + 6070.2476*T
            val Je = 49.511251  + 64.3000*T

            ra = (
                268.056595 - 0.006499*T
                + 0.000117 * dsin(Ja)
                + 0.000938 * dsin(Jb)
                + 0.001432 * dsin(Jc)
                + 0.000030 * dsin(Jd)
                + 0.002150 * dsin(Je)
            )

            dec = (
                64.495303 + 0.002413*T
                + 0.000050 * dcos(Ja)
                + 0.000404 * dcos(Jb)
                + 0.000617 * dcos(Jc)
                - 0.000013 * dcos(Jd)
                + 0.000926 * dcos(Je)
            )

            w = 284.95 + 870.536*d
        }

        Body.Saturn -> {
            ra = 40.589 - 0.036*T
            dec = 83.537 - 0.004*T
            w = 38.90 + 810.7939024*d
        }

        Body.Uranus -> {
            ra = 257.311
            dec = -15.175
            w = 203.81 - 501.1600928*d
        }

        Body.Neptune -> {
            val N = 357.85 + 52.316*T
            val sinN = dsin(N)
            ra = 299.36 + 0.70*sinN
            dec = 43.46 - 0.51*dcos(N)
            w = 249.978 + 541.1397757*d - 0.48*sinN
        }

        Body.Pluto -> {
            ra = 132.993
            dec = -6.163
            w = 302.695 + 56.3625225*d
        }

        else -> throw InvalidBodyException(body)
    }

    // Calculate the north pole vector using the given angles.
    val rcoslat = dcos(dec)
    val north = AstroVector(
        rcoslat * dcos(ra),
        rcoslat * dsin(ra),
        dsin(dec),
        time
    )

    return AxisInfo(ra / 15.0, dec, w, north)
}

/**
* Calculates spherical ecliptic geocentric position of the Moon.
*
* Given a time of observation, calculates the Moon's geocentric position
* in ecliptic spherical coordinates. Provides the ecliptic latitude and
* longitude in degrees, and the geocentric distance in astronomical units (AU).
* The ecliptic longitude is measured relative to the equinox of date.
*
* This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
* which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
* It is adapted from Turbo Pascal code from the book
* [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
* by Montenbruck and Pfleger.
*
* To calculate an equatorial J2000 vector instead, use [geoMoon].
*/
fun eclipticGeoMoon(time: AstroTime) = MoonContext(time).calcMoon()

/**
 * Calculates equatorial geocentric position of the Moon at a given time.
 *
 * Given a time of observation, calculates the Moon's position vector.
 * The vector indicates the Moon's center relative to the Earth's center.
 * The vector components are expressed in AU (astronomical units).
 * The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
 * In Astronomy Engine, this orientation is called EQJ.
 *
 * @param time
 *      The date and time for which to calculate the Moon's position.
 *
 * @return The Moon's position vector in J2000 equatorial coordinates (EQJ).
 */
fun geoMoon(time: AstroTime): AstroVector {
    val eclSphere = eclipticGeoMoon(time)
    val eclVec = eclSphere.toVector(time)
    val equVec = eclipticToEquatorial(eclVec)
    return precession(equVec, PrecessDirection.Into2000)
}

/**
 * Calculates equatorial geocentric position and velocity of the Moon at a given time.
 *
 * Given a time of observation, calculates the Moon's position and velocity vectors.
 * The position and velocity are of the Moon's center relative to the Earth's center.
 * The position (x, y, z) components are expressed in AU (astronomical units).
 * The velocity (vx, vy, vz) components are expressed in AU/day.
 * The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
 * In Astronomy Engine, this orientation is called EQJ.
 * If you need the Moon's position only, and not its velocity,
 * it is much more efficient to use [geoMoon] instead.
 *
 * @param time
 *      The date and time for which to calculate the Moon's position and velocity.
 *
 * @return The Moon's position and velocity vectors in J2000 equatorial coordinates (EQJ).
 */
fun geoMoonState(time: AstroTime): StateVector {
    // This is a hack, because trying to figure out how to derive
    // a time derivative for MoonContext.calcMoon() would be painful!
    // Calculate just before and just after the given time.
    // Average to find position, subtract to find velocity.
    val dt = 1.0e-5        // 0.864 seconds
    val t1 = time.addDays(-dt)
    val t2 = time.addDays(+dt)
    val r1 = geoMoon(t1)
    val r2 = geoMoon(t2)

    // The desired position is the average of the two calculated positions.
    // The difference of position vectors divided by the time span gives the velocity vector.
    return StateVector(
        (r1.x + r2.x) / 2.0,
        (r1.y + r2.y) / 2.0,
        (r1.z + r2.z) / 2.0,
        (r2.x - r1.x) / (2.0 * dt),
        (r2.y - r1.y) / (2.0 * dt),
        (r2.z - r1.z) / (2.0 * dt),
        time
    )
}

private fun helioEarthPos(time: AstroTime) =
    calcVsop(vsopModel(Body.Earth), time)

private fun helioEarthState(time: AstroTime) =
    StateVector(calcVsopPosVel(vsopModel(Body.Earth), time.tt), time)

private fun barycenterPosContrib(time: AstroTime, body: Body, planetGm: Double) =
    (planetGm / (planetGm + SUN_GM)) * vsopHelioVector(body, time)

private fun solarSystemBarycenterPos(time: AstroTime): AstroVector {
    val j = barycenterPosContrib(time, Body.Jupiter, JUPITER_GM)
    val s = barycenterPosContrib(time, Body.Saturn,  SATURN_GM)
    var u = barycenterPosContrib(time, Body.Uranus,  URANUS_GM)
    var n = barycenterPosContrib(time, Body.Neptune, NEPTUNE_GM)
    return AstroVector(
        j.x + s.x + u.x + n.x,
        j.y + s.y + u.y + n.y,
        j.z + s.z + u.z + n.z,
        time
    )
}

private fun barycenterStateContrib(time: AstroTime, body: Body, planetGm: Double): StateVector {
    val helioPlanet = calcVsopPosVel(vsopModel(body), time.tt)
    val factor = planetGm / (planetGm + SUN_GM)
    return StateVector(
        factor * helioPlanet.r.x,
        factor * helioPlanet.r.y,
        factor * helioPlanet.r.z,
        factor * helioPlanet.v.x,
        factor * helioPlanet.v.y,
        factor * helioPlanet.v.z,
        time
    )
}

private fun solarSystemBarycenterState(time: AstroTime): StateVector {
    val j = barycenterStateContrib(time, Body.Jupiter, JUPITER_GM)
    val s = barycenterStateContrib(time, Body.Saturn,  SATURN_GM)
    var u = barycenterStateContrib(time, Body.Uranus,  URANUS_GM)
    var n = barycenterStateContrib(time, Body.Neptune, NEPTUNE_GM)
    return StateVector(
        j.x + s.x + u.x + n.x,
        j.y + s.y + u.y + n.y,
        j.z + s.z + u.z + n.z,
        j.vx + s.vx + u.vx + n.vx,
        j.vy + s.vy + u.vy + n.vy,
        j.vz + s.vz + u.vz + n.vz,
        time
    )
}

/**
 * Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.
 *
 * This function calculates the position of the given celestial body as a vector,
 * using the center of the Sun as the origin.  The result is expressed as a Cartesian
 * vector in the J2000 equatorial system: the coordinates are based on the mean equator
 * of the Earth at noon UTC on 1 January 2000.
 *
 * The position is not corrected for light travel time or aberration.
 * This is different from the behavior of [geoVector].
 *
 * If given an invalid value for `body`, this function will throw an [InvalidBodyException].
 *
 * @param body
 *      A body for which to calculate a heliocentric position:
 *      the Sun, Moon, EMB, SSB, or any of the planets.
 *
 * @param time
 *      The date and time for which to calculate the position.
 *
 * @return The heliocentric position vector of the center of the given body.
 */
fun helioVector(body: Body, time: AstroTime): AstroVector =
    if (body.vsopModel != null)
        calcVsop(body.vsopModel, time)
    else when (body) {
        Body.Sun     -> AstroVector(0.0, 0.0, 0.0, time)
        Body.Pluto   -> calcPluto(time, true).position()
        Body.Moon    -> helioEarthPos(time) + geoMoon(time)
        Body.EMB     -> helioEarthPos(time) + (geoMoon(time) / (1.0 + EARTH_MOON_MASS_RATIO))
        Body.SSB     -> solarSystemBarycenterPos(time)
        else -> throw InvalidBodyException(body)
    }

/**
 * Calculates the distance between a body and the Sun at a given time.
 *
 * Given a date and time, this function calculates the distance between
 * the center of `body` and the center of the Sun, expressed in AU.
 * For the planets Mercury through Neptune, this function is significantly
 * more efficient than calling [helioVector] followed by taking the length
 * of the resulting vector.
 *
 * @param body
 *      A body for which to calculate a heliocentric distance:
 *      the Sun, Moon, EMB, SSB, or any of the planets.
 *
 * @param time
 *      The date and time for which to calculate the distance.
 *
 * @return The heliocentric distance in AU.
 */
fun helioDistance(body: Body, time: AstroTime): Double =
    when {
        body == Body.Sun -> 0.0
        body.vsopModel != null -> vsopDistance(body.vsopModel, time)
        else -> helioVector(body, time).length()
    }

/**
 * Calculates heliocentric position and velocity vectors for the given body.
 *
 * Given a body and a time, calculates the position and velocity
 * vectors for the center of that body at that time, relative to the center of the Sun.
 * The vectors are expressed in equatorial J2000 coordinates (EQJ).
 * If you need the position vector only, it is more efficient to call [helioVector].
 * The Sun's center is a non-inertial frame of reference. In other words, the Sun
 * experiences acceleration due to gravitational forces, mostly from the larger
 * planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum,
 * kinetic energy, or other quantities that require a non-accelerating frame
 * of reference, consider using [baryState] instead.
 *
 * @param body
 * The celestial body whose heliocentric state vector is to be calculated.
 * Supported values are `Body.Sun`, `Body.Moon`, `Body.EMB`, `Body.SSB`, and all planets:
 * `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`,
 * `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
 *
 * @param time
 *      The date and time for which to calculate position and velocity.
 *
 * @return
 * A state vector that contains heliocentric position and velocity vectors.
 * The positions are expressed in AU.
 * The velocities are expressed in AU/day.
 */
fun helioState(body: Body, time: AstroTime): StateVector =
    if (body.vsopModel != null)
        StateVector(calcVsopPosVel(body.vsopModel, time.tt), time)
    else when (body) {
        Body.Sun   -> StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time)
        Body.Pluto -> calcPluto(time, true)
        Body.Moon  -> helioEarthState(time) + geoMoonState(time)
        Body.EMB   -> helioEarthState(time) + (geoMoonState(time) / (1.0 + EARTH_MOON_MASS_RATIO))
        Body.SSB   -> solarSystemBarycenterState(time)
        else -> throw InvalidBodyException(body)
    }

/**
 * Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.
 *
 * This function calculates the position of the given celestial body as a vector,
 * using the center of the Earth as the origin.  The result is expressed as a Cartesian
 * vector in the J2000 equatorial system: the coordinates are based on the mean equator
 * of the Earth at noon UTC on 1 January 2000.
 *
 * If given an invalid value for `body`, this function will throw an exception.
 *
 * Unlike [helioVector], this function always corrects for light travel time.
 * This means the position of the body is "back-dated" by the amount of time it takes
 * light to travel from that body to an observer on the Earth.
 *
 * Also, the position can optionally be corrected for
 * [aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
 * causing the apparent direction of the body to be shifted due to transverse
 * movement of the Earth with respect to the rays of light coming from that body.
 *
 * @param body
 *      A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.
 *
 * @param time
 *      The date and time for which to calculate the position.
 *
 * @param aberration
 *      `Aberration.Corrected` to correct for aberration, or `Aberration.None` to leave uncorrected.
 *
 * @return A geocentric position vector of the center of the given body.
 */
fun geoVector(body: Body, time: AstroTime, aberration: Aberration): AstroVector {
    if (body == Body.Earth)
        return AstroVector(0.0, 0.0, 0.0, time)

    if (body == Body.Moon)
        return geoMoon(time)

    // For all other bodies, apply light travel time correction.
    // The intention is to find the apparent position of the body
    // from the Earth's point of view.

    var earth = helioEarthPos(time)

    var ltime = time
    for (iter in 0..9) {
        val helio = helioVector(body, ltime)
        if (aberration == Aberration.Corrected && iter > 0) {
            // Include aberration, so make a good first-order approximation
            // by backdating the Earth's position also.
            // This is confusing, but it works for objects within the Solar System
            // because the distance the Earth moves in that small amount of light
            // travel time (a few minutes to a few hours) is well approximated
            // by a line segment that substends the angle seen from the remote
            // body viewing Earth. That angle is pretty close to the aberration
            // angle of the moving Earth viewing the remote body.
            // In other words, both of the following approximate the aberration angle:
            //     (transverse distance Earth moves) / (distance to body)
            //     (transverse speed of Earth) / (speed of light).
            earth = helioEarthPos(ltime)
        }

        // Convert heliocentric vector to geocentric vector.
        // Tricky: we cannot use the subtraction operator because
        // it will get angry that we are using mismatching times!
        // It is intentional here that the calculation time was backdated,
        // but the observation time is not.
        var geopos = AstroVector(
            helio.x - earth.x,
            helio.y - earth.y,
            helio.z - earth.z,
            time
        )

        // Calculate the time in the past when light left the body on its way toward Earth.
        val ltime2 = time.addDays(-geopos.length() / C_AUDAY)

        // Very quickly we should converge on a solution for how far
        // in the past light must have left the planet in order to
        // reach the Earth at the given time, even though the observed
        // body was in a slightly different orbital position when
        // light left it.
        if ((ltime2.tt - ltime.tt).absoluteValue < 1.0e-9)
            return geopos

        // Otherwise we refine the estimate and try again.
        ltime = ltime2
    }

    // This should never happen. Usually the solver converges
    // after 3 iterations. We allow for 10 iterations.
    // Something is really wrong if this ever happens.
    throw InternalError("Light travel time did not converge")
}

/**
 * Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.
 *
 * Calculates topocentric equatorial coordinates in one of two different systems:
 * J2000 or true-equator-of-date, depending on the value of the `equdate` parameter.
 * Equatorial coordinates include right ascension, declination, and distance in astronomical units.
 *
 * This function corrects for light travel time: it adjusts the apparent location
 * of the observed body based on how long it takes for light to travel from the body to the Earth.
 *
 * This function corrects for *topocentric parallax*, meaning that it adjusts for the
 * angular shift depending on where the observer is located on the Earth. This is most
 * significant for the Moon, because it is so close to the Earth. However, parallax corection
 * has a small effect on the apparent positions of other bodies.
 *
 * Correction for aberration is optional, using the `aberration` parameter.
 *
 * @param body
 *      The celestial body to be observed. Not allowed to be `Body.Earth`.
 *
 * @param time
 *      The date and time at which the observation takes place.
 *
 * @param observer
 *      A location on or near the surface of the Earth.
 *
 * @param equdate
 *      Selects the date of the Earth's equator in which to express the equatorial coordinates.
 *
 * @param aberration
 *      Selects whether or not to correct for aberration.
 *
 * @return Topocentric equatorial coordinates of the celestial body.
 */
fun equator(
    body: Body,
    time: AstroTime,
    observer: Observer,
    equdate: EquatorEpoch,
    aberration: Aberration
): Equatorial {
    val gcObserver = geoPos(time, observer)
    val gc = geoVector(body, time, aberration)
    val j2000 = gc - gcObserver
    val vector = when (equdate) {
        EquatorEpoch.OfDate -> gyration(j2000, PrecessDirection.From2000)
        EquatorEpoch.J2000  -> j2000
    }
    return equatorFromVector(vector)
}

/**
 * Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.
 *
 * Given a date and time, the geographic location of an observer on the Earth, and
 * equatorial coordinates (right ascension and declination) of a celestial body,
 * this function returns horizontal coordinates (azimuth and altitude angles) for the body
 * relative to the horizon at the geographic location.
 *
 * The right ascension `ra` and declination `dec` passed in must be *equator of date*
 * coordinates, based on the Earth's true equator at the date and time of the observation.
 * Otherwise the resulting horizontal coordinates will be inaccurate.
 * Equator of date coordinates can be obtained by calling [equator], passing in
 * [EquatorEpoch.OfDate] as its `equdate` parameter. It is also recommended to enable
 * aberration correction by passing in [Aberration.Corrected] as the `aberration` parameter.
 *
 * This function optionally corrects for atmospheric refraction.
 * For most uses, it is recommended to pass [Refraction.Normal] in the `refraction` parameter to
 * correct for optical lensing of the Earth's atmosphere that causes objects
 * to appear somewhat higher above the horizon than they actually are.
 * However, callers may choose to avoid this correction by passing in [Refraction.None].
 * If refraction correction is enabled, the azimuth, altitude, right ascension, and declination
 * in the [Topocentric] object returned by this function will all be corrected for refraction.
 * If refraction is disabled, none of these four coordinates will be corrected; in that case,
 * the right ascension and declination in the returned structure will be numerically identical
 * to the respective `ra` and `dec` values passed in.
 *
 * @param time
 *      The date and time of the observation.
 *
 * @param observer
 *      The geographic location of the observer.
 *
 * @param ra
 *      The right ascension of the body in sidereal hours. See remarks above for more details.
 *
 * @param dec
 *      The declination of the body in degrees. See remarks above for more details.
 *
 * @param refraction
 *      Selects whether to correct for atmospheric refraction, and if so, which model to use.
 *      The recommended value for most uses is `Refraction.Normal`.
 *      See remarks above for more details.
 *
 * @return The body's apparent horizontal coordinates and equatorial coordinates, both optionally corrected for refraction.
 */
fun horizon(
    time: AstroTime,
    observer: Observer,
    ra: Double,
    dec: Double,
    refraction: Refraction
): Topocentric {
    val sinlat = dsin(observer.latitude)
    val coslat = dcos(observer.latitude)
    val sinlon = dsin(observer.longitude)
    val coslon = dcos(observer.longitude)
    val sindc = dsin(dec)
    val cosdc = dcos(dec)
    val sinra = dsin(ra * 15.0)
    val cosra = dcos(ra * 15.0)

    // Calculate three mutually perpendicular unit vectors
    // in equatorial coordinates: uze, une, uwe.
    //
    // uze = The direction of the observer's local zenith (straight up).
    // une = The direction toward due north on the observer's horizon.
    // uwe = The direction toward due west on the observer's horizon.
    //
    // HOWEVER, these are uncorrected for the Earth's rotation due to the time of day.
    //
    // The components of these 3 vectors are as follows:
    // x = direction from center of Earth toward 0 degrees longitude (the prime meridian) on equator.
    // y = direction from center of Earth toward 90 degrees west longitude on equator.
    // z = direction from center of Earth toward the north pole.
    val uze = AstroVector(coslat * coslon, coslat * sinlon, sinlat, time)
    val une = AstroVector(-sinlat * coslon, -sinlat * sinlon, coslat, time)
    val uwe = AstroVector(sinlon, -coslon, 0.0, time)

    // Correct the vectors uze, une, uwe for the Earth's rotation by calculating
    // sideral time. Call spin() for each uncorrected vector to rotate about
    // the Earth's axis to yield corrected unit vectors uz, un, uw.
    // Multiply sidereal hours by -15 to convert to degrees and flip eastward
    // rotation of the Earth to westward apparent movement of objects with time.
    val angle = -15.0 * siderealTime(time)
    val uz = spin(angle, uze)
    val un = spin(angle, une)
    val uw = spin(angle, uwe)

    // Convert angular equatorial coordinates (RA, DEC) to
    // cartesian equatorial coordinates in 'p', using the
    // same orientation system as uze, une, uwe.
    val p = AstroVector(cosdc * cosra, cosdc * sinra, sindc, time)

    // Use dot products of p with the zenith, north, and west
    // vectors to obtain the cartesian coordinates of the body in
    // the observer's horizontal orientation system.
    // pz = zenith component [-1, +1]
    // pn = north  component [-1, +1]
    // pw = west   component [-1, +1]
    val pz = p dot uz
    val pn = p dot un
    val pw = p dot uw

    // projHor is the "shadow" of the body vector along the observer's flat ground.
    val projHor = hypot(pn, pw)

    // Calculate az = azimuth (compass direction clockwise from East.)
    val az = (
        if (projHor > 0.0) (
            // If the body is not exactly straight up/down, it has an azimuth.
            // Invert the angle to produce degrees eastward from north.
            (-atan2(pw, pn)).radiansToDegrees().withMinDegreeValue(0.0)
        ) else (
            // The body is straight up/down, so it does not have an azimuth.
            // Report an arbitrary but reasonable value.
            0.0
        )
    )

    // zd = the angle of the body away from the observer's zenith, in degrees.
    var zd = atan2(projHor, pz).radiansToDegrees()
    var horRa = ra
    var horDec = dec

    if (refraction != Refraction.None) {
        val zd0 = zd
        val refr = refractionAngle(refraction, 90.0 - zd)
        zd -= refr

        if (refr > 0.0 && zd > 3.0e-4) {
            // Calculate refraction-corrected equatorial coordinates.
            val sinzd  = dsin(zd)
            val coszd  = dcos(zd)
            val sinzd0 = dsin(zd0)
            val coszd0 = dcos(zd0)

            val prx = ((p.x - coszd0 * uz.x) / sinzd0)*sinzd + uz.x*coszd
            val pry = ((p.y - coszd0 * uz.y) / sinzd0)*sinzd + uz.y*coszd
            val prz = ((p.z - coszd0 * uz.z) / sinzd0)*sinzd + uz.z*coszd

            val projEqu = hypot(prx, pry)

            horRa =
                if (projEqu > 0.0)
                    atan2(pry, prx).radiansToDegrees().withMinDegreeValue(0.0) / 15.0
                else
                    0.0

            horDec = atan2(prz, projEqu).radiansToDegrees()
        }
    }

    return Topocentric(az, 90.0 - zd, horRa, horDec)
}

/**
 * Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.
 *
 * Calculates position and velocity vectors for Jupiter's moons
 * Io, Europa, Ganymede, and Callisto, at the given date and time.
 * The vectors are jovicentric (relative to the center of Jupiter).
 * Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ).
 * The position components are expressed in astronomical units (AU), and the
 * velocity components are in AU/day.
 *
 * To convert to heliocentric position vectors, call [helioVector]
 * with `Body.Jupiter` to get Jupiter's heliocentric position, then
 * add the jovicentric positions. Likewise, you can call [geoVector]
 * to convert to geocentric positions; however, you will have to manually
 * correct for light travel time from the Jupiter system to Earth to
 * figure out what time to pass to `jupiterMoons` to get an accurate picture
 * of how Jupiter and its moons look from Earth.
 */
fun jupiterMoons(time: AstroTime) =
    JupiterMoonsInfo(arrayOf(
        calcJupiterMoon(time, jupiterMoonModel[0]),
        calcJupiterMoon(time, jupiterMoonModel[1]),
        calcJupiterMoon(time, jupiterMoonModel[2]),
        calcJupiterMoon(time, jupiterMoonModel[3])
    ))

/**
 * Searches for a time at which a function's value increases through zero.
 *
 * Certain astronomy calculations involve finding a time when an event occurs.
 * Often such events can be defined as the root of a function:
 * the time at which the function's value becomes zero.
 *
 * `search` finds the *ascending root* of a function: the time at which
 * the function's value becomes zero while having a positive slope. That is, as time increases,
 * the function transitions from a negative value, through zero at a specific moment,
 * to a positive value later. The goal of the search is to find that specific moment.
 *
 * The `func` parameter is an instance of the interface [SearchContext].
 * As an example, a caller may wish to find the moment a celestial body reaches a certain
 * ecliptic longitude. In that case, the caller might derive a class that contains
 * a [Body] member to specify the body and a `Double` to hold the target longitude.
 * It could subtract the target longitude from the actual longitude at a given time;
 * thus the difference would equal zero at the moment in time the planet reaches the
 * desired longitude.
 *
 * Every time it is called, `func.eval` returns a `Double` value or it throws an exception.
 * If `func.eval` throws an exception, the search immediately fails and the exception
 * is propagated to the caller. Otherwise, the search proceeds until it either finds
 * the ascending root or fails for some reason.
 *
 * The search calls `func.eval` repeatedly to rapidly narrow in on any ascending
 * root within the time window specified by `time1` and `time2`. The search never
 * reports a solution outside this time window.
 *
 * `search` uses a combination of bisection and quadratic interpolation
 * to minimize the number of function calls. However, it is critical that the
 * supplied time window be small enough that there cannot be more than one root
 * (ascedning or descending) within it; otherwise the search can fail.
 * Beyond that, it helps to make the time window as small as possible, ideally
 * such that the function itself resembles a smooth parabolic curve within that window.
 *
 * If an ascending root is not found, or more than one root
 * (ascending and/or descending) exists within the window `time1`..`time2`,
 * the search will return `null`.
 *
 * If the search does not converge within 20 iterations, it will throw an exception.
 *
 * @param func
 *      The function for which to find the time of an ascending root.
 *      See remarks above for more details.
 *
 * @param time1
 *      The lower time bound of the search window.
 *      See remarks above for more details.
 *
 * @param time2
 *      The upper time bound of the search window.
 *      See remarks above for more details.
 *
 * @param toleranceSeconds
 *      Specifies an amount of time in seconds within which a bounded ascending root
 *      is considered accurate enough to stop. A typical value is 1 second.
 *
 * @return
 * If successful, returns an [AstroTime] value indicating a date and time
 * that is within `toleranceSeconds` of an ascending root.
 * If no ascending root is found, or more than one root exists in the time
 * window `time1`..`time2`, the function returns `null`.
 */
fun search(
    func: SearchContext,
    time1: AstroTime,
    time2: AstroTime,
    toleranceSeconds: Double
): AstroTime? {
    var t1 = time1
    var t2 = time2
    val iterLimit = 20
    val toleranceDays = abs(toleranceSeconds / SECONDS_PER_DAY)
    var f1 = func.eval(t1)
    var f2 = func.eval(t2)
    var calcFmid = true
    var fmid = 0.0

    for (iter in 1..iterLimit) {
        val dt = (t2.tt - t1.tt) / 2.0
        val tmid = t1.addDays(dt)
        if (dt.absoluteValue < toleranceDays) {
            // We are close enough to the event to stop the search.
            return tmid
        }

        if (calcFmid)
            fmid = func.eval(tmid)
        else
            calcFmid = true     // we already have the correct value of fmid from the previous loop

        // Quadratic interpolation:
        // Try to find a parabola that passes through the 3 points we have sampled:
        // (t1,f1), (tmid,fmid), (t2,f2).
        val tm = tmid.ut
        val tspan = t2.ut - tmid.ut
        val q = (f2 + f1)/2.0 - fmid
        val r = (f2 - f1)/2.0
        val s = fmid
        var foundInterpolation = false
        var x = Double.NaN
        if (q == 0.0) {
            // This is a line, not a parabola.
            if (r != 0.0) {     // skip horizontal lines: they don't have a root
                x = -s / r
                foundInterpolation = (-1.0 <= x && x <= +1.0)
            }
        } else {
            // This really is a parabola. Find its roots x1, x2.
            val u = r*r - 4.0*q*s
            if (u > 0.0) {      // skip imaginary or tangent roots
                // See if there is a unique solution for x in the range [-1, +1].
                val ru = sqrt(u)
                val x1 = (-r + ru) / (2.0 * q)
                val x2 = (-r - ru) / (2.0 * q)
                val x1Valid = (-1.0 <= x1 && x1 <= +1.0)
                val x2Valid = (-1.0 <= x2 && x2 <= +1.0)
                if (x1Valid && !x2Valid) {
                    x = x1
                    foundInterpolation = true
                } else if (x2Valid && !x1Valid) {
                    x = x2
                    foundInterpolation = true
                }
            }
        }
        if (foundInterpolation) {
            val qut = tm + x*tspan
            val qslope = (2*q*x + r) / tspan
            val tq = AstroTime(qut)
            val fq = func.eval(tq)
            if (qslope != 0.0) {
                var dtGuess = abs(fq / qslope)
                if (dtGuess < toleranceDays) {
                    // The estimated time error is small enough that we can quit now.
                    return tq
                }

                // Try guessing a tighter boundary with the interpolated root at the center.
                dtGuess *= 1.2
                if (dtGuess < dt / 10.0) {
                    val tleft = tq.addDays(-dtGuess)
                    val tright = tq.addDays(+dtGuess)
                    if ((tleft.ut - t1.ut)*(tleft.ut - t2.ut) < 0.0) {
                        if ((tright.ut - t1.ut)*(tright.ut - t2.ut) < 0.0) {
                            val fleft = func.eval(tleft)
                            val fright = func.eval(tright)
                            if ((fleft < 0.0) && (fright >= 0.0)) {
                                f1 = fleft
                                f2 = fright
                                t1 = tleft
                                t2 = tright
                                fmid = fq
                                calcFmid = false    // save a little work; no need to recalculate fmid next time
                                continue
                            }
                        }
                    }
                }
            }
        }

        // The quadratic interpolation didn't work this time.
        // Use bisection: divide the region in two parts and pick
        // whichever one appears to contain a root.
        if (f1 < 0.0 && fmid >= 0.0) {
            t2 = tmid
            f2 = fmid
            continue
        }

        if (fmid < 0.0 && f2 >= 0.0) {
            t1 = tmid
            f1 = fmid
            continue
        }

        // Either there is no ascending zero-crossing in this range
        // or the search window is too wide (more than one zero-crossing).
        // Either way, the search has failed.
        return null
    }

    throw InternalError("Search did not converge within $iterLimit iterations.")
}

/**
 * Calculates geocentric ecliptic coordinates for the Sun.
 *
 * This function calculates the position of the Sun as seen from the Earth.
 * The returned value includes both Cartesian and spherical coordinates.
 * The x-coordinate and longitude values in the returned object are based
 * on the *true equinox of date*: one of two points in the sky where the instantaneous
 * plane of the Earth's equator at the given date and time (the *equatorial plane*)
 * intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*).
 * By convention, the apparent location of the Sun at the March equinox is chosen
 * as the longitude origin and x-axis direction, instead of the one for September.
 *
 * `sunPosition` corrects for precession and nutation of the Earth's axis
 * in order to obtain the exact equatorial plane at the given time.
 *
 * This function can be used for calculating changes of seasons: equinoxes and solstices.
 * In fact, the function [seasons] does use this function for that purpose.
 *
 * @param time
 *      The date and time for which to calculate the Sun's position.
 *
 * @return The ecliptic coordinates of the Sun using the Earth's true equator of date.
 */
fun sunPosition(time: AstroTime): Ecliptic {
    // Correct for light travel time from the Sun.
    // Otherwise season calculations (equinox, solstice) will all be early by about 8 minutes!
    val adjustedTime = time.addDays(-1.0 / C_AUDAY)
    val earth2000 = helioEarthPos(adjustedTime)

    // Convert heliocentric location of Earth to geocentric location of Sun.
    val sun2000 = -earth2000

    // Convert to equatorial Cartesian coordinates of date.
    val sunOfDate = gyration(sun2000, PrecessDirection.From2000)

    // Convert equatorial coordinates to ecliptic coordinates.
    val trueObliq = earthTilt(adjustedTime).tobl.degreesToRadians()
    return rotateEquatorialToEcliptic(sunOfDate, trueObliq)
}

private fun rotateEquatorialToEcliptic(pos: AstroVector, obliqRadians: Double): Ecliptic {
    val cosOb = cos(obliqRadians)
    val sinOb = sin(obliqRadians)
    val ex = +pos.x
    val ey = +pos.y*cosOb + pos.z*sinOb
    val ez = -pos.y*sinOb + pos.z*cosOb
    val xyproj = hypot(ex, ey)
    val elon =
        if (xyproj > 0.0)
            atan2(ey, ex).radiansToDegrees().withMinDegreeValue(0.0)
        else
            0.0
    val elat = atan2(ez, xyproj).radiansToDegrees()
    val vec = AstroVector(ex, ey, ez, pos.t)
    return Ecliptic(vec, elat, elon)
}

/**
 * Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.
 *
 * Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
 * on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates,
 * which are relative to the plane of the Earth's orbit around the Sun.
 *
 * @param equ
 *      Equatorial coordinates in the J2000 frame of reference.
 *      You can call [geoVector] to obtain suitable equatorial coordinates.
 *
 * @return Ecliptic coordinates in the J2000 frame of reference (ECL).
 */
fun equatorialToEcliptic(equ: AstroVector): Ecliptic =
    rotateEquatorialToEcliptic(
        equ,
        0.40909260059599012     // mean obliquity of the J2000 ecliptic in radians
    )

/**
 * Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.
 *
 * This function finds the moment in time, if any exists in the given time window,
 * that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.
 *
 * This function can be used to determine equinoxes and solstices.
 * However, it is usually more convenient and efficient to call [seasons]
 * to calculate all equinoxes and solstices for a given calendar year.
 *
 * The function searches the window of time specified by `startTime` and `startTime+limitDays`.
 * The search will return `null` if the Sun never reaches the longitude `targetLon` or
 * if the window is so large that the longitude ranges more than 180 degrees within it.
 * It is recommended to keep the window smaller than 10 days when possible.
 */
fun searchSunLongitude(targetLon: Double, startTime: AstroTime, limitDays: Double): AstroTime? {
    class Context(val targetLon: Double) : SearchContext {
        override fun eval(time: AstroTime) =
            longitudeOffset(sunPosition(time).elon - targetLon)
    }
    val context = Context(targetLon)
    val time2 = startTime.addDays(limitDays)
    return search(context, startTime, time2, 0.01)
}

/**
 * Finds both equinoxes and both solstices for a given calendar year.
 *
 * The changes of seasons are defined by solstices and equinoxes.
 * Given a calendar year number, this function calculates the
 * March and September equinoxes and the June and December solstices.
 *
 * The equinoxes are the moments twice each year when the plane of the
 * Earth's equator passes through the center of the Sun. In other words,
 * the Sun's declination is zero at both equinoxes.
 * The March equinox defines the beginning of spring in the northern hemisphere
 * and the beginning of autumn in the southern hemisphere.
 * The September equinox defines the beginning of autumn in the northern hemisphere
 * and the beginning of spring in the southern hemisphere.
 *
 * The solstices are the moments twice each year when one of the Earth's poles
 * is most tilted toward the Sun. More precisely, the Sun's declination reaches
 * its minimum value at the December solstice, which defines the beginning of
 * winter in the northern hemisphere and the beginning of summer in the southern
 * hemisphere. The Sun's declination reaches its maximum value at the June solstice,
 * which defines the beginning of summer in the northern hemisphere and the beginning
 * of winter in the southern hemisphere.
 *
 * @param year
 *      The calendar year number for which to calculate equinoxes and solstices.
 *      The value may be any integer, but only the years 1800 through 2100 have been
 *      validated for accuracy: unit testing against data from the
 *      United States Naval Observatory confirms that all equinoxes and solstices
 *      for that range of years are within 2 minutes of the correct time.
 *
 * @return
 * A [SeasonsInfo] object that contains four [AstroTime] values:
 * the March and September equinoxes and the June and December solstices.
 */
fun seasons(year: Int) =
    SeasonsInfo(
        findSeasonChange(  0.0, year,  3, 10),
        findSeasonChange( 90.0, year,  6, 10),
        findSeasonChange(180.0, year,  9, 10),
        findSeasonChange(270.0, year, 12, 10)
    )

private fun findSeasonChange(targetLon: Double, year: Int, month: Int, day: Int): AstroTime {
    var startTime = AstroTime(year, month, day, 0, 0, 0.0)
    return searchSunLongitude(targetLon, startTime, 20.0) ?:
        throw InternalError("Cannot find solution for Sun longitude $targetLon for year $year")
}

/**
 * Returns one body's ecliptic longitude with respect to another, as seen from the Earth.
 *
 * This function determines where one body appears around the ecliptic plane
 * (the plane of the Earth's orbit around the Sun) as seen from the Earth,
 * relative to the another body's apparent position.
 * The function returns an angle in the half-open range [0, 360) degrees.
 * The value is the ecliptic longitude of `body1` relative to the ecliptic
 * longitude of `body2`.
 *
 * The angle is 0 when the two bodies are at the same ecliptic longitude
 * as seen from the Earth. The angle increases in the prograde direction
 * (the direction that the planets orbit the Sun and the Moon orbits the Earth).
 *
 * When the angle is 180 degrees, it means the two bodies appear on opposite sides
 * of the sky for an Earthly observer.
 *
 * Neither `body1` nor `body2` is allowed to be `Body.Earth`.
 * If this happens, the function throws an exception.
 *
 * @param body1 The first body, whose longitude is to be found relative to the second body.
 * @param body2 The second body, relative to which the longitude of the first body is to be found.
 * @param time  The date and time of the observation.
 * @return An angle in the range [0, 360), expressed in degrees.
 */
fun pairLongitude(body1: Body, body2: Body, time: AstroTime): Double {
    if (body1 == Body.Earth || body2 == Body.Earth)
        throw EarthNotAllowedException()

    val vector1 = geoVector(body1, time, Aberration.None)
    val eclip1 = equatorialToEcliptic(vector1)

    val vector2 = geoVector(body2, time, Aberration.None)
    val eclip2 = equatorialToEcliptic(vector2)

    return normalizeLongitude(eclip1.elon - eclip2.elon)
}

/**
 * Returns the Moon's phase as an angle from 0 to 360 degrees.
 *
 * This function determines the phase of the Moon using its apparent
 * ecliptic longitude relative to the Sun, as seen from the center of the Earth.
 * Certain values of the angle have conventional definitions:
 *
 * - 0 = new moon
 * - 90 = first quarter
 * - 180 = full moon
 * - 270 = third quarter
 *
 * @param time  The date and time of the observation.
 * @return The angle as described above, a value in the range 0..360 degrees.
 */
fun moonPhase(time: AstroTime): Double =
    pairLongitude(Body.Moon, Body.Sun, time)

/**
 * Searches for the time that the Moon reaches a specified phase.
 *
 * Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic
 * longitude with respect to the Sun's geocentric ecliptic longitude.
 * When the Moon and the Sun have the same longitude, that is defined as a new moon.
 * When their longitudes are 180 degrees apart, that is defined as a full moon.
 *
 * This function searches for any value of the lunar phase expressed as an
 * angle in degrees in the range [0, 360).
 *
 * If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter)
 * it is much easier to call the functions [searchMoonQuarter] and [nextMoonQuarter].
 * This function is useful for finding general phase angles outside those four quarters.
 *
 * @param targetLon
 *      The difference in geocentric longitude between the Sun and Moon
 *      that specifies the lunar phase being sought. This can be any value
 *      in the range [0, 360).  Certain values have conventional names:
 *      0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter.
 *
 * @param startTime
 *      The beginning of the time window in which to search for the Moon reaching the specified phase.
 *
 * @param limitDays
 *      The number of days after `startTime` that limits the time window for the search.
 *
 * @return
 * If successful, returns the date and time the moon reaches the phase specified by
 * `targetlon`. This function will return `null` if the phase does not
 * occur within `limitDays` of `startTime`; that is, if the search window is too small.
 */
fun searchMoonPhase(targetLon: Double, startTime: AstroTime, limitDays: Double): AstroTime? {
    // To avoid discontinuities in the moonOffset function causing problems,
    // we need to approximate when that function will next return 0.
    // We probe it with the start time and take advantage of the fact
    // that every lunar phase repeats roughly every 29.5 days.
    // There is a surprising uncertainty in the quarter timing,
    // due to the eccentricity of the moon's orbit.
    // I have seen more than 0.9 days away from the simple prediction.
    // To be safe, we take the predicted time of the event and search
    // +/-1.5 days around it (a 3-day wide window).
    class Context(val targetLon : Double) : SearchContext {
        override fun eval(time: AstroTime) = longitudeOffset(moonPhase(time) - targetLon)
    }
    val moonOffset = Context(targetLon)
    var ya = moonOffset.eval(startTime)
    if (ya > 0.0) ya -= 360.0  // force searching forward in time, not backward
    val uncertainty = 1.5
    val estDt = -(MEAN_SYNODIC_MONTH * ya) / 360.0
    val dt1 = estDt - uncertainty
    if (dt1 > limitDays)
        return null    // not possible for moon phase to occur within specified window (too short)
    val dt2 = min(limitDays, estDt + uncertainty)
    val t1 = startTime.addDays(dt1)
    val t2 = startTime.addDays(dt2)
    return search(moonOffset, t1, t2, 1.0)
}

/**
 * Finds the first lunar quarter after the specified date and time.
 * A lunar quarter is one of the following four lunar phase events:
 * new moon, first quarter, full moon, third quarter.
 * This function finds the lunar quarter that happens soonest
 * after the specified date and time.
 *
 * To continue iterating through consecutive lunar quarters, call this function once,
 * followed by calls to #NextMoonQuarter as many times as desired.
 *
 * @param startTime The date and time at which to start the search.
 * @return A [MoonQuarterInfo] object reporting the next quarter phase and the time it will occur.
 */
fun searchMoonQuarter(startTime: AstroTime): MoonQuarterInfo {
    val currentPhaseAngle = moonPhase(startTime)
    val quarter: Int = (1 + floor(currentPhaseAngle / 90.0).toInt()) % 4
    val quarterTime = searchMoonPhase(90.0 * quarter, startTime, 10.0) ?:
        throw InternalError("Unable to find moon quarter $quarter for startTime=$startTime")
    return MoonQuarterInfo(quarter, quarterTime)
}

/**
 * Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.
 *
 * The *hour angle* of a celestial body indicates its position in the sky with respect
 * to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
 * The hour angle is 0 when the body reaches its highest angle above the horizon in a given day.
 * The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
 * to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
 * the number of hours that have passed since the most recent time that the body has culminated,
 * or reached its highest point.
 *
 * This function searches for the next time a celestial body reaches the given hour angle
 * after the date and time specified by `startTime`.
 * To find when a body culminates, pass 0 for `hourAngle`.
 * To find when a body reaches its lowest point in the sky, pass 12 for `hourAngle`.
 *
 * Note that, especially close to the Earth's poles, a body as seen on a given day
 * may always be above the horizon or always below the horizon, so the caller cannot
 * assume that a culminating object is visible nor that an object is below the horizon
 * at its minimum altitude.
 *
 * On success, the function reports the date and time, along with the horizontal coordinates
 * of the body at that time, as seen by the given observer.
 *
 * @param body
 *      The celestial body, which can the Sun, the Moon, or any planet other than the Earth.
 * @param observer
 *      A location on or near the surface of the Earth where the observer is located.
 * @param hourAngle
 *      An hour angle value in the range [0, 24) indicating the number of sidereal hours after the
 *      body's most recent culmination.
 * @param startTime
 *      The date and time at which to start the search.
 * @return The time when the body reaches the hour angle, and the horizontal coordinates of the body at that time.
 */
fun searchHourAngle(
    body: Body,
    observer: Observer,
    hourAngle: Double,
    startTime: AstroTime
): HourAngleInfo {
    if (body == Body.Earth)
        throw EarthNotAllowedException()

    if (hourAngle < 0.0 || hourAngle >= 24.0)
        throw IllegalArgumentException("hourAngle=$hourAngle is out of the allowed range [0, 24).")

    var time = startTime
    var iter = 0
    while (true) {
        ++iter

        // Calculate Greenwich Apparent Sidereal Time (GAST) at the given time.
        val gast = siderealTime(time)

        // Obtain equatorial coordinates of date for the body.
        val ofdate = equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected)

        // Calculate the adjustment needed in sidereal time
        // to bring the hour angle to the desired value.
        var deltaSiderealHours = ((hourAngle + ofdate.ra - observer.longitude/15.0) - gast) % 24.0
        if (iter == 1) {
            // On the first iteration, always search forward in time.
            if (deltaSiderealHours < 0.0)
                deltaSiderealHours += 24.0
        } else {
            // On subsequent iterations, we make the smallest possible adjustment,
            // either forward or backward in time.
            if (deltaSiderealHours < -12.0)
                deltaSiderealHours += 24.0
            else if (deltaSiderealHours > +12.0)
                deltaSiderealHours -= 24.0
        }

        // If the error is tolerable (less than 0.1 seconds), the search has succeeded.
        if (deltaSiderealHours.absoluteValue * 3600.0 < 0.1) {
            val hor = horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.Normal)
            return HourAngleInfo(time, hor)
        }

        // We need to loop another time to get more accuracy.
        // Update the terrestrial time (in solar days) by adjusting sidereal time.
        time = time.addDays((deltaSiderealHours / 24.0) * SOLAR_DAYS_PER_SIDEREAL_DAY)
    }
}

private fun internalSearchAltitude(
    body: Body,
    observer: Observer,
    direction: Direction,
    startTime: AstroTime,
    limitDays: Double,
    context: SearchContext
): AstroTime? {
    if (body == Body.Earth)
        throw EarthNotAllowedException()

    // Find the pair of hour angles that bound the desired event.
    // When a body's hour angle is 0, it means it is at its highest point
    // in the observer's sky, called culmination.
    // If the body's hour angle is 12, it means it is at its lowest point in the sky.
    // (Note that it is possible for a body to be above OR below the horizon in either case.)
    // If the caller wants a rising event, we want the pair haBefore=12, haAfter=0.
    // If the caller wants a setting event, the desired pair is haBefore=0, haAfter=12.
    val haBefore: Double = when (direction) {
        Direction.Rise -> 12.0      // minimum altitude (bottom) happens before the body rises
        Direction.Set  ->  0.0      // culmination happens before the body sets
    }
    val haAfter: Double = 12.0 - haBefore

    // See if the body is currently above/below the horizon.
    // If we are looking for next rise time and the body is below the horizon,
    // we use the current time as the lower time bound and the next culmination
    // as the upper bound.
    // If the body is above the horizon, we search for the next bottom and use it
    // as the lower bound and the next culmination after that bottom as the upper bound.
    // The same logic applies for finding set times, only we swap the hour angles.

    var altBefore = context.eval(startTime)
    var timeBefore: AstroTime
    if (altBefore > 0.0) {
        // We are past the sought event, so we have to wait for the next "before" event (culm/bottom).
        timeBefore = searchHourAngle(body, observer, haBefore, startTime).time
        altBefore = context.eval(timeBefore)
    } else {
        // We are before or at the sought event, so we find the next "after" event,
        // and use the current time as the "before" event.
        timeBefore = startTime
    }

    var timeAfter = searchHourAngle(body, observer, haAfter, timeBefore).time
    var altAfter = context.eval(timeAfter)

    while (true) {
        if (altBefore <= 0.0 && altAfter > 0.0) {
            // The body crosses the horizon during the time interval.
            // Search between evtBefore and evtAfter for the desired event.
            val time = search(context, timeBefore, timeAfter, 1.0)
            if (time != null)
                return time
        }

        // If we didn't find the desired event, find the next hour angle bracket and try again.
        val evtBefore = searchHourAngle(body, observer, haBefore, timeAfter)
        val evtAfter = searchHourAngle(body, observer, haAfter, timeBefore)

        if (evtBefore.time.ut >= startTime.ut + limitDays)
            return null

        timeBefore = evtBefore.time
        timeAfter = evtAfter.time
        altBefore = context.eval(timeBefore)
        altAfter = context.eval(timeAfter)
    }
}

private class SearchContextPeakAltitude(
    private val body: Body,
    private val direction: Direction,
    private val observer: Observer
): SearchContext {
    private val bodyRadiusAu: Double

    init {
        bodyRadiusAu = when(body) {
            Body.Sun -> SUN_RADIUS_AU
            Body.Moon -> MOON_EQUATORIAL_RADIUS_AU
            else -> 0.0
        }
    }

    override fun eval(time: AstroTime): Double {
        // Return the angular altitude above or below the horizon
        // of the highest part (the peak) of the given object.
        // This is defined as the apparent altitude of the center of the body plus
        // the body's angular radius.
        // The 'direction' parameter controls whether the angle is measured
        // positive above the horizon or positive below the horizon,
        // depending on whether the caller wants rise times or set times, respectively.

        val ofdate: Equatorial = equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected)
        val hor: Topocentric = horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None)
        return direction.sign * (hor.altitude + (bodyRadiusAu / ofdate.dist).radiansToDegrees() + REFRACTION_NEAR_HORIZON)
    }
}

/**
 * Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.
 *
 * This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
 * Rise time is when the body first starts to be visible above the horizon.
 * For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
 * Set time is the moment when the body appears to vanish below the horizon.
 *
 * This function corrects for typical atmospheric refraction, which causes celestial
 * bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
 * It also adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).
 *
 * Note that rise or set may not occur in every 24 hour period.
 * For example, near the Earth's poles, there are long periods of time where
 * the Sun stays below the horizon, never rising.
 * Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day.
 * This is because the Moon sets nearly an hour later each day due to orbiting the Earth a
 * significant amount during each rotation of the Earth.
 * Therefore callers must not assume that the function will always succeed.
 *
 * @param body
 *      The Sun, Moon, or any planet other than the Earth.
 *
 * @param observer
 *      The location where observation takes place.
 *
 * @param direction
 *      Either [Direction.Rise] to find a rise time or [Direction.Set] to find a set time.
 *
 * @param startTime
 *      The date and time at which to start the search.
 *
 * @param limitDays
 *      Limits how many days to search for a rise or set time.
 *      To limit a rise or set time to the same day, you can use a value of 1 day.
 *      In cases where you want to find the next rise or set time no matter how far
 *      in the future (for example, for an observer near the south pole), you can
 *      pass in a larger value like 365.
 *
 * @return
 * On success, returns the date and time of the rise or set time as requested.
 * If the function returns `null`, it means the rise or set event does not occur
 * within `limitDays` days of `startTime`. This is a normal condition,
 * not an error.
 */
fun searchRiseSet(
    body: Body,
    observer: Observer,
    direction: Direction,
    startTime: AstroTime,
    limitDays: Double
): AstroTime? {
    val context = SearchContextPeakAltitude(body, direction, observer)
    return internalSearchAltitude(body, observer, direction, startTime, limitDays, context)
}

private class SearchContextAltitudeError(
    private val body: Body,
    private val direction: Direction,
    private val observer: Observer,
    private val altitude: Double
): SearchContext {
    override fun eval(time: AstroTime): Double {
        val ofdate = equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected)
        val hor = horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None)
        return direction.sign * (hor.altitude - altitude)
    }
}

/**
 * Finds the next time a body reaches a given altitude.
 *
 * Finds when the given body ascends or descends through a given
 * altitude angle, as seen by an observer at the specified location on the Earth.
 * By using the appropriate combination of `direction` and `altitude` parameters,
 * this function can be used to find when civil, nautical, or astronomical twilight
 * begins (dawn) or ends (dusk).
 *
 * Civil dawn begins before sunrise when the Sun ascends through 6 degrees below
 * the horizon. To find civil dawn, pass `Direction.Rise` for `direction` and -6 for `altitude`.
 *
 * Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon.
 * To find civil dusk, pass `Direction.Set` for `direction` and -6 for `altitude`.
 *
 * Nautical twilight is similar to civil twilight, only the `altitude` value should be -12 degrees.
 *
 * Astronomical twilight uses -18 degrees as the `altitude` value.
 *
 * @param body The Sun, Moon, or any planet other than the Earth.
 * @param observer The location where observation takes place.
 * @param direction Either `Direction.Rise` to find an ascending altitude event or `Direction.Set` to find a descending altitude event.
 * @param startTime The date and time at which to start the search.
 * @param limitDays The fractional number of days after `dateStart` that limits when the altitude event is to be found. Must be a positive number.
 * @param altitude The desired altitude angle of the body's center above (positive) or below (negative) the observer's local horizon, expressed in degrees. Must be in the range [-90, +90].
 * @return The date and time of the altitude event, or `null` if no such event occurs within the specified time window.
 */
fun searchAltitude(
    body: Body,
    observer: Observer,
    direction: Direction,
    startTime: AstroTime,
    limitDays: Double,
    altitude: Double
): AstroTime? {
    val context = SearchContextAltitudeError(body, direction, observer, altitude)
    return internalSearchAltitude(body, observer, direction, startTime, limitDays, context)
}

/**
 * Continues searching for lunar quarters from a previous search.
 *
 * After calling [searchMoonQuarter], this function can be called
 * one or more times to continue finding consecutive lunar quarters.
 * This function finds the next consecutive moon quarter event after
 * the one passed in as the parameter `mq`.
 *
 * @param The previous moon quarter found by a call to [searchMoonQuarter] or `nextMoonQuarter`.
 * @return The moon quarter that occurs next in time after the one passed in `mq`.
 */
fun nextMoonQuarter(mq: MoonQuarterInfo): MoonQuarterInfo {
    // Skip 6 days past the previous found moon quarter to find the next one.
    // This is less than the minimum possible increment.
    // So far I have seen the interval well contained by the range (6.5, 8.3) days.
    val time = mq.time.addDays(6.0)
    val nextMoonQuarter = searchMoonQuarter(time)
    // Verify that we found the expected moon quarter.
    val expected = (1 + mq.quarter) % 4
    if (nextMoonQuarter.quarter != expected)
        throw InternalError("Expected to find next quarter $expected, but found ${nextMoonQuarter.quarter}")
    return nextMoonQuarter
}

/**
 * Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using equator at J2000 epoch.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 */
fun rotationEqjEcl(): RotationMatrix {
    // ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
    val c = 0.9174821430670688    // cos(ob)
    val s = 0.3977769691083922    // sin(ob)
    return RotationMatrix(
        1.0, 0.0, 0.0,
        0.0,  +c,  -s,
        0.0,  +s,  +c
    )
}

/**
 * Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: EQJ = equatorial system, using equator at J2000 epoch.
 */
fun rotationEclEqj(): RotationMatrix {
    // ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
    val c = 0.9174821430670688    // cos(ob)
    val s = 0.3977769691083922    // sin(ob)
    return RotationMatrix(
        1.0, 0.0, 0.0,
        0.0,  +c,  +s,
        0.0,  -s,  +c
    )
}

/**
 * Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using equator at J2000 epoch.
 * Target: EQD = equatorial system, using equator of the specified date/time.
 *
 * @param time
 *      The date and time at which the Earth's equator defines the target orientation.
 */
fun rotationEqjEqd(time: AstroTime): RotationMatrix =
    precessionRot(time, PrecessDirection.From2000) combine
    nutationRot(time, PrecessDirection.From2000)

/**
 * Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of the specified date/time.
 * Target: EQJ = equatorial system, using equator at J2000 epoch.
 *
 * @param time
 *      The date and time at which the Earth's equator defines the source orientation.
 */
fun rotationEqdEqj(time: AstroTime): RotationMatrix =
    nutationRot(time, PrecessDirection.Into2000) combine
    precessionRot(time, PrecessDirection.Into2000)

/**
 * Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of the specified date/time.
 * Target: HOR = horizontal system.
 *
 * @param time
 *      The date and time at which the Earth's equator applies.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return
 * A rotation matrix that converts EQD to HOR at `time` and for `observer`.
 * The components of the horizontal vector are:
 * x = north, y = west, z = zenith (straight up from the observer).
 * These components are chosen so that the "right-hand rule" works for the vector
 * and so that north represents the direction where azimuth = 0.
 */
fun rotationEqdHor(time: AstroTime, observer: Observer): RotationMatrix {
    // See the `horizon` function for more explanation of how this works.

    val sinlat = dsin(observer.latitude)
    val coslat = dcos(observer.latitude)
    val sinlon = dsin(observer.longitude)
    val coslon = dcos(observer.longitude)

    val uze = AstroVector(coslat * coslon, coslat * sinlon, sinlat, time)
    val une = AstroVector(-sinlat * coslon, -sinlat * sinlon, coslat, time)
    val uwe = AstroVector(sinlon, -coslon, 0.0, time)

    // Multiply sidereal hours by -15 to convert to degrees and flip eastward
    // rotation of the Earth to westward apparent movement of objects with time.

    val angle = -15.0 * siderealTime(time)
    val uz = spin(angle, uze)
    val un = spin(angle, une)
    val uw = spin(angle, uwe)

    return RotationMatrix(
        un.x, uw.x, uz.x,
        un.y, uw.y, uz.y,
        un.z, uw.z, uz.z
    )
}

/**
 * Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system (x=North, y=West, z=Zenith).
 * Target: EQD = equatorial system, using equator of the specified date/time.
 *
 * @param time
 *      The date and time at which the Earth's equator applies.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return A rotation matrix that converts HOR to EQD at `time` and for `observer`.
 */
fun rotationHorEqd(time: AstroTime, observer: Observer): RotationMatrix =
    rotationEqdHor(time, observer).inverse()

/**
 * Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system (x=North, y=West, z=Zenith).
 * Target: EQJ = equatorial system, using equator at the J2000 epoch.
 *
 * @param time
 *      The date and time of the observation.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return A rotation matrix that converts HOR to EQJ at `time` and for `observer`.
 */
fun rotationHorEqj(time: AstroTime, observer: Observer): RotationMatrix =
    rotationHorEqd(time, observer) combine
    rotationEqdEqj(time)

/**
 * Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using the equator at the J2000 epoch.
 * Target: HOR = horizontal system.
 *
 * @param time
 *      The date and time of the observation.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return
 * A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
 * The components of the horizontal vector are:
 * x = north, y = west, z = zenith (straight up from the observer).
 * These components are chosen so that the "right-hand rule" works for the vector
 * and so that north represents the direction where azimuth = 0.
 */
fun rotationEqjHor(time: AstroTime, observer: Observer): RotationMatrix =
    rotationHorEqj(time, observer).inverse()

/**
 * Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of date.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @param time
 *      The date and time of the source equator.
 *
 * @return A rotation matrix that converts EQD to ECL.
 */
fun rotationEqdEcl(time: AstroTime): RotationMatrix =
    rotationEqdEqj(time) combine
    rotationEqjEcl()

/**
 * Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: EQD = equatorial system, using equator of date.
 *
 * @param time
 *      The date and time of the desired equator.
 *
 * @return A rotation matrix that converts ECL to EQD.
 */
fun rotationEclEqd(time: AstroTime): RotationMatrix =
    rotationEqdEcl(time).inverse()

/**
 * Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: HOR = horizontal system.
 *
 * @param time
 *      The date and time of the desired horizontal orientation.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return
 * A rotation matrix that converts ECL to HOR at `time` and for `observer`.
 * The components of the horizontal vector are:
 * x = north, y = west, z = zenith (straight up from the observer).
 * These components are chosen so that the "right-hand rule" works for the vector
 * and so that north represents the direction where azimuth = 0.
 */
fun rotationEclHor(time: AstroTime, observer: Observer): RotationMatrix =
    rotationEclEqd(time) combine
    rotationEqdHor(time, observer)

/**
 * Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @param time
 *      The date and time of the horizontal observation.
 *
 * @param observer
 *      The location of the horizontal observer.
 *
 * @return A rotation matrix that converts HOR to ECL.
 */
fun rotationHorEcl(time: AstroTime, observer: Observer): RotationMatrix =
    rotationEclHor(time, observer).inverse()

/**
 * Calculates a rotation matrix from galactic (GAL) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: GAL = galactic system (IAU 1958 definition).
 * Target: EQJ = equatorial system, using the equator at the J2000 epoch.
 *
 * @return A rotation matrix that converts GAL to EQJ.
 */
fun rotationEqjGal() =
    // This rotation matrix was calculated by the following script:
    // demo/python/galeqj_matrix.py
    RotationMatrix(
        -0.0548624779711344, +0.4941095946388765, -0.8676668813529025,
        -0.8734572784246782, -0.4447938112296831, -0.1980677870294097,
        -0.4838000529948520, +0.7470034631630423, +0.4559861124470794
    )

/**
 * Calculates a rotation matrix from galactic (GAL) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: GAL = galactic system (IAU 1958 definition).
 * Target: EQJ = equatorial system, using the equator at the J2000 epoch.
 *
 * @return A rotation matrix that converts GAL to EQJ.
 */
fun rotationGalEqj() =
    // This rotation matrix was calculated by the following script:
    // demo/python/galeqj_matrix.py
    RotationMatrix(
        -0.0548624779711344, -0.8734572784246782, -0.4838000529948520,
        +0.4941095946388765, -0.4447938112296831, +0.7470034631630423,
        -0.8676668813529025, -0.1980677870294097, +0.4559861124470794
    )

/**
 * Determines the constellation that contains the given point in the sky.
 *
 * Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the
 * constellation that contains that point.
 *
 * @param ra
 *      The right ascension (RA) of a point in the sky, using the J2000 equatorial system (EQJ).
 *
 * @param dec
 *      The declination (DEC) of a point in the sky, using the J2000 equatorial system (EQJ).
 *
 * @return
 * A structure that contains the 3-letter abbreviation and full name
 * of the constellation that contains the given (ra,dec), along with
 * the converted B1875 (ra,dec) for that point.
 */
fun constellation(ra: Double, dec: Double): ConstellationInfo {
    if (dec < -90.0 || dec > +90.0)
        throw IllegalArgumentException("Invalid declination angle $dec. Must be -90..+90.")

    // Convert right ascension to degrees, and normalize to the range [0, 360).
    val raDeg = (ra * 15.0).withMinDegreeValue(0.0)

    // Convert coordinates from J2000 to B1875.
    var sph2000 = Spherical(dec, raDeg, 1.0)
    val vec2000 = sph2000.toVector(epoch2000)
    val vec1875 = constelRot.rotate(vec2000)
    val equ1875 = vec1875.toEquatorial()

    // Convert DEC from degrees and RA from hours,
    // into compact angle units used in the constelBounds table.
    val xDec = 24.0 * equ1875.dec
    val xRa = (24.0 * 15.0) * equ1875.ra

    // Search for the constellation using the B1875 coordinates.
    for (b in constelBounds)
        if ((b.decLo <= xDec) && (b.raHi > xRa) && (b.raLo <= xRa))
            return ConstellationInfo(constelNames[b.index].symbol, constelNames[b.index].name, equ1875.ra, equ1875.dec)

    // This should never happen! If it does, please report to: https://github.com/cosinekitty/astronomy/issues
    throw InternalError("Unable to find constellation for coordinates RA=$ra, DEC=$dec")
}

//=======================================================================================
// Generated code goes to the bottom of the source file,
// so that most line numbers are consistent between template code and target code.

//---------------------------------------------------------------------------------------
// Table of coefficients for calculating nutation using the IAU2000b model.

private class IauRow(
    val nals0: Int,
    val nals1: Int,
    val nals2: Int,
    val nals3: Int,
    val nals4: Int,
    val cls0: Double,
    val cls1: Double,
    val cls2: Double,
    val cls3: Double,
    val cls4: Double,
    val cls5: Double
)

private val iauRow: Array<IauRow> = arrayOf($ASTRO_IAU_DATA())

//---------------------------------------------------------------------------------------
// VSOP87 coefficients, for calculating major planet state vectors.

$ASTRO_KOTLIN_VSOP(Mercury)
$ASTRO_KOTLIN_VSOP(Venus)
$ASTRO_KOTLIN_VSOP(Earth)
$ASTRO_KOTLIN_VSOP(Mars)
$ASTRO_KOTLIN_VSOP(Jupiter)
$ASTRO_KOTLIN_VSOP(Saturn)
$ASTRO_KOTLIN_VSOP(Uranus)
$ASTRO_KOTLIN_VSOP(Neptune)

//---------------------------------------------------------------------------------------
// Geocentric Moon

private fun addSolarTerms(context: MoonContext) {$ASTRO_ADDSOL()}

//---------------------------------------------------------------------------------------
// Pluto state table

$ASTRO_PLUTO_TABLE()
private val plutoCache = HashMap<Int, List<BodyGravCalc>>()

//---------------------------------------------------------------------------------------
// Models for Jupiter's four largest moons.

$ASTRO_JUPITER_MOONS()

//---------------------------------------------------------------------------------------
// Constellation lookup table.

$ASTRO_CONSTEL()

//---------------------------------------------------------------------------------------

// Create a rotation matrix for converting J2000 to B1875.
// Need to calculate the B1875 epoch. Based on this:
// https://en.wikipedia.org/wiki/Epoch_(astronomy)#Besselian_years
// B = 1900 + (JD - 2415020.31352) / 365.242198781
// I'm interested in using TT instead of JD, giving:
// B = 1900 + ((TT+2451545) - 2415020.31352) / 365.242198781
// B = 1900 + (TT + 36524.68648) / 365.242198781
// TT = 365.242198781*(B - 1900) - 36524.68648 = -45655.741449525
// But the AstroTime constructor wants UT, not TT.
// Near that date, I get a historical correction of ut-tt = 3.2 seconds.
// That gives UT = -45655.74141261017 for the B1875 epoch,
// or 1874-12-31T18:12:21.950Z.
private val constelRot: RotationMatrix = rotationEqjEqd(AstroTime(-45655.74141261017))
