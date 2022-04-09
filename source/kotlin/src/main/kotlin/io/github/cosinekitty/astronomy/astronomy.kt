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

private val iauRow: Array<IauRow> = arrayOf(
    IauRow( 0,  0,  0,  0,  1, -172064161.0,    -174666.0,      33386.0,   92052331.0,       9086.0,      15377.0),
    IauRow( 0,  0,  2, -2,  2,  -13170906.0,      -1675.0,     -13696.0,    5730336.0,      -3015.0,      -4587.0),
    IauRow( 0,  0,  2,  0,  2,   -2276413.0,       -234.0,       2796.0,     978459.0,       -485.0,       1374.0),
    IauRow( 0,  0,  0,  0,  2,    2074554.0,        207.0,       -698.0,    -897492.0,        470.0,       -291.0),
    IauRow( 0,  1,  0,  0,  0,    1475877.0,      -3633.0,      11817.0,      73871.0,       -184.0,      -1924.0),
    IauRow( 0,  1,  2, -2,  2,    -516821.0,       1226.0,       -524.0,     224386.0,       -677.0,       -174.0),
    IauRow( 1,  0,  0,  0,  0,     711159.0,         73.0,       -872.0,      -6750.0,          0.0,        358.0),
    IauRow( 0,  0,  2,  0,  1,    -387298.0,       -367.0,        380.0,     200728.0,         18.0,        318.0),
    IauRow( 1,  0,  2,  0,  2,    -301461.0,        -36.0,        816.0,     129025.0,        -63.0,        367.0),
    IauRow( 0, -1,  2, -2,  2,     215829.0,       -494.0,        111.0,     -95929.0,        299.0,        132.0),
    IauRow( 0,  0,  2, -2,  1,     128227.0,        137.0,        181.0,     -68982.0,         -9.0,         39.0),
    IauRow(-1,  0,  2,  0,  2,     123457.0,         11.0,         19.0,     -53311.0,         32.0,         -4.0),
    IauRow(-1,  0,  0,  2,  0,     156994.0,         10.0,       -168.0,      -1235.0,          0.0,         82.0),
    IauRow( 1,  0,  0,  0,  1,      63110.0,         63.0,         27.0,     -33228.0,          0.0,         -9.0),
    IauRow(-1,  0,  0,  0,  1,     -57976.0,        -63.0,       -189.0,      31429.0,          0.0,        -75.0),
    IauRow(-1,  0,  2,  2,  2,     -59641.0,        -11.0,        149.0,      25543.0,        -11.0,         66.0),
    IauRow( 1,  0,  2,  0,  1,     -51613.0,        -42.0,        129.0,      26366.0,          0.0,         78.0),
    IauRow(-2,  0,  2,  0,  1,      45893.0,         50.0,         31.0,     -24236.0,        -10.0,         20.0),
    IauRow( 0,  0,  0,  2,  0,      63384.0,         11.0,       -150.0,      -1220.0,          0.0,         29.0),
    IauRow( 0,  0,  2,  2,  2,     -38571.0,         -1.0,        158.0,      16452.0,        -11.0,         68.0),
    IauRow( 0, -2,  2, -2,  2,      32481.0,          0.0,          0.0,     -13870.0,          0.0,          0.0),
    IauRow(-2,  0,  0,  2,  0,     -47722.0,          0.0,        -18.0,        477.0,          0.0,        -25.0),
    IauRow( 2,  0,  2,  0,  2,     -31046.0,         -1.0,        131.0,      13238.0,        -11.0,         59.0),
    IauRow( 1,  0,  2, -2,  2,      28593.0,          0.0,         -1.0,     -12338.0,         10.0,         -3.0),
    IauRow(-1,  0,  2,  0,  1,      20441.0,         21.0,         10.0,     -10758.0,          0.0,         -3.0),
    IauRow( 2,  0,  0,  0,  0,      29243.0,          0.0,        -74.0,       -609.0,          0.0,         13.0),
    IauRow( 0,  0,  2,  0,  0,      25887.0,          0.0,        -66.0,       -550.0,          0.0,         11.0),
    IauRow( 0,  1,  0,  0,  1,     -14053.0,        -25.0,         79.0,       8551.0,         -2.0,        -45.0),
    IauRow(-1,  0,  0,  2,  1,      15164.0,         10.0,         11.0,      -8001.0,          0.0,         -1.0),
    IauRow( 0,  2,  2, -2,  2,     -15794.0,         72.0,        -16.0,       6850.0,        -42.0,         -5.0),
    IauRow( 0,  0, -2,  2,  0,      21783.0,          0.0,         13.0,       -167.0,          0.0,         13.0),
    IauRow( 1,  0,  0, -2,  1,     -12873.0,        -10.0,        -37.0,       6953.0,          0.0,        -14.0),
    IauRow( 0, -1,  0,  0,  1,     -12654.0,         11.0,         63.0,       6415.0,          0.0,         26.0),
    IauRow(-1,  0,  2,  2,  1,     -10204.0,          0.0,         25.0,       5222.0,          0.0,         15.0),
    IauRow( 0,  2,  0,  0,  0,      16707.0,        -85.0,        -10.0,        168.0,         -1.0,         10.0),
    IauRow( 1,  0,  2,  2,  2,      -7691.0,          0.0,         44.0,       3268.0,          0.0,         19.0),
    IauRow(-2,  0,  2,  0,  0,     -11024.0,          0.0,        -14.0,        104.0,          0.0,          2.0),
    IauRow( 0,  1,  2,  0,  2,       7566.0,        -21.0,        -11.0,      -3250.0,          0.0,         -5.0),
    IauRow( 0,  0,  2,  2,  1,      -6637.0,        -11.0,         25.0,       3353.0,          0.0,         14.0),
    IauRow( 0, -1,  2,  0,  2,      -7141.0,         21.0,          8.0,       3070.0,          0.0,          4.0),
    IauRow( 0,  0,  0,  2,  1,      -6302.0,        -11.0,          2.0,       3272.0,          0.0,          4.0),
    IauRow( 1,  0,  2, -2,  1,       5800.0,         10.0,          2.0,      -3045.0,          0.0,         -1.0),
    IauRow( 2,  0,  2, -2,  2,       6443.0,          0.0,         -7.0,      -2768.0,          0.0,         -4.0),
    IauRow(-2,  0,  0,  2,  1,      -5774.0,        -11.0,        -15.0,       3041.0,          0.0,         -5.0),
    IauRow( 2,  0,  2,  0,  1,      -5350.0,          0.0,         21.0,       2695.0,          0.0,         12.0),
    IauRow( 0, -1,  2, -2,  1,      -4752.0,        -11.0,         -3.0,       2719.0,          0.0,         -3.0),
    IauRow( 0,  0,  0, -2,  1,      -4940.0,        -11.0,        -21.0,       2720.0,          0.0,         -9.0),
    IauRow(-1, -1,  0,  2,  0,       7350.0,          0.0,         -8.0,        -51.0,          0.0,          4.0),
    IauRow( 2,  0,  0, -2,  1,       4065.0,          0.0,          6.0,      -2206.0,          0.0,          1.0),
    IauRow( 1,  0,  0,  2,  0,       6579.0,          0.0,        -24.0,       -199.0,          0.0,          2.0),
    IauRow( 0,  1,  2, -2,  1,       3579.0,          0.0,          5.0,      -1900.0,          0.0,          1.0),
    IauRow( 1, -1,  0,  0,  0,       4725.0,          0.0,         -6.0,        -41.0,          0.0,          3.0),
    IauRow(-2,  0,  2,  0,  2,      -3075.0,          0.0,         -2.0,       1313.0,          0.0,         -1.0),
    IauRow( 3,  0,  2,  0,  2,      -2904.0,          0.0,         15.0,       1233.0,          0.0,          7.0),
    IauRow( 0, -1,  0,  2,  0,       4348.0,          0.0,        -10.0,        -81.0,          0.0,          2.0),
    IauRow( 1, -1,  2,  0,  2,      -2878.0,          0.0,          8.0,       1232.0,          0.0,          4.0),
    IauRow( 0,  0,  0,  1,  0,      -4230.0,          0.0,          5.0,        -20.0,          0.0,         -2.0),
    IauRow(-1, -1,  2,  2,  2,      -2819.0,          0.0,          7.0,       1207.0,          0.0,          3.0),
    IauRow(-1,  0,  2,  0,  0,      -4056.0,          0.0,          5.0,         40.0,          0.0,         -2.0),
    IauRow( 0, -1,  2,  2,  2,      -2647.0,          0.0,         11.0,       1129.0,          0.0,          5.0),
    IauRow(-2,  0,  0,  0,  1,      -2294.0,          0.0,        -10.0,       1266.0,          0.0,         -4.0),
    IauRow( 1,  1,  2,  0,  2,       2481.0,          0.0,         -7.0,      -1062.0,          0.0,         -3.0),
    IauRow( 2,  0,  0,  0,  1,       2179.0,          0.0,         -2.0,      -1129.0,          0.0,         -2.0),
    IauRow(-1,  1,  0,  1,  0,       3276.0,          0.0,          1.0,         -9.0,          0.0,          0.0),
    IauRow( 1,  1,  0,  0,  0,      -3389.0,          0.0,          5.0,         35.0,          0.0,         -2.0),
    IauRow( 1,  0,  2,  0,  0,       3339.0,          0.0,        -13.0,       -107.0,          0.0,          1.0),
    IauRow(-1,  0,  2, -2,  1,      -1987.0,          0.0,         -6.0,       1073.0,          0.0,         -2.0),
    IauRow( 1,  0,  0,  0,  2,      -1981.0,          0.0,          0.0,        854.0,          0.0,          0.0),
    IauRow(-1,  0,  0,  1,  0,       4026.0,          0.0,       -353.0,       -553.0,          0.0,       -139.0),
    IauRow( 0,  0,  2,  1,  2,       1660.0,          0.0,         -5.0,       -710.0,          0.0,         -2.0),
    IauRow(-1,  0,  2,  4,  2,      -1521.0,          0.0,          9.0,        647.0,          0.0,          4.0),
    IauRow(-1,  1,  0,  1,  1,       1314.0,          0.0,          0.0,       -700.0,          0.0,          0.0),
    IauRow( 0, -2,  2, -2,  1,      -1283.0,          0.0,          0.0,        672.0,          0.0,          0.0),
    IauRow( 1,  0,  2,  2,  1,      -1331.0,          0.0,          8.0,        663.0,          0.0,          4.0),
    IauRow(-2,  0,  2,  2,  2,       1383.0,          0.0,         -2.0,       -594.0,          0.0,         -2.0),
    IauRow(-1,  0,  0,  0,  2,       1405.0,          0.0,          4.0,       -610.0,          0.0,          2.0),
    IauRow( 1,  1,  2, -2,  2,       1290.0,          0.0,          0.0,       -556.0,          0.0,          0.0)
)

//---------------------------------------------------------------------------------------
// VSOP87 coefficients, for calculating major planet state vectors.

private val vsopLonMercury0 = VsopSeries(arrayOf(
    VsopTerm(4.40250710144, 0.00000000000, 0.00000000000),
    VsopTerm(0.40989414977, 1.48302034195, 26087.90314157420),
    VsopTerm(0.05046294200, 4.47785489551, 52175.80628314840),
    VsopTerm(0.00855346844, 1.16520322459, 78263.70942472259),
    VsopTerm(0.00165590362, 4.11969163423, 104351.61256629678),
    VsopTerm(0.00034561897, 0.77930768443, 130439.51570787099),
    VsopTerm(0.00007583476, 3.71348404924, 156527.41884944518)
))

private val vsopLonMercury1 = VsopSeries(arrayOf(
    VsopTerm(26087.90313685529, 0.00000000000, 0.00000000000),
    VsopTerm(0.01131199811, 6.21874197797, 26087.90314157420),
    VsopTerm(0.00292242298, 3.04449355541, 52175.80628314840),
    VsopTerm(0.00075775081, 6.08568821653, 78263.70942472259),
    VsopTerm(0.00019676525, 2.80965111777, 104351.61256629678)
))

private val vsopLonMercury = VsopFormula(arrayOf(
    vsopLonMercury0,
    vsopLonMercury1
))

private val vsopLatMercury0 = VsopSeries(arrayOf(
    VsopTerm(0.11737528961, 1.98357498767, 26087.90314157420),
    VsopTerm(0.02388076996, 5.03738959686, 52175.80628314840),
    VsopTerm(0.01222839532, 3.14159265359, 0.00000000000),
    VsopTerm(0.00543251810, 1.79644363964, 78263.70942472259),
    VsopTerm(0.00129778770, 4.83232503958, 104351.61256629678),
    VsopTerm(0.00031866927, 1.58088495658, 130439.51570787099),
    VsopTerm(0.00007963301, 4.60972126127, 156527.41884944518)
))

private val vsopLatMercury1 = VsopSeries(arrayOf(
    VsopTerm(0.00274646065, 3.95008450011, 26087.90314157420),
    VsopTerm(0.00099737713, 3.14159265359, 0.00000000000)
))

private val vsopLatMercury = VsopFormula(arrayOf(
    vsopLatMercury0,
    vsopLatMercury1
))

private val vsopRadMercury0 = VsopSeries(arrayOf(
    VsopTerm(0.39528271651, 0.00000000000, 0.00000000000),
    VsopTerm(0.07834131818, 6.19233722598, 26087.90314157420),
    VsopTerm(0.00795525558, 2.95989690104, 52175.80628314840),
    VsopTerm(0.00121281764, 6.01064153797, 78263.70942472259),
    VsopTerm(0.00021921969, 2.77820093972, 104351.61256629678),
    VsopTerm(0.00004354065, 5.82894543774, 130439.51570787099)
))

private val vsopRadMercury1 = VsopSeries(arrayOf(
    VsopTerm(0.00217347740, 4.65617158665, 26087.90314157420),
    VsopTerm(0.00044141826, 1.42385544001, 52175.80628314840)
))

private val vsopRadMercury = VsopFormula(arrayOf(
    vsopRadMercury0,
    vsopRadMercury1
))


private val vsopLonVenus0 = VsopSeries(arrayOf(
    VsopTerm(3.17614666774, 0.00000000000, 0.00000000000),
    VsopTerm(0.01353968419, 5.59313319619, 10213.28554621100),
    VsopTerm(0.00089891645, 5.30650047764, 20426.57109242200),
    VsopTerm(0.00005477194, 4.41630661466, 7860.41939243920),
    VsopTerm(0.00003455741, 2.69964447820, 11790.62908865880),
    VsopTerm(0.00002372061, 2.99377542079, 3930.20969621960),
    VsopTerm(0.00001317168, 5.18668228402, 26.29831979980),
    VsopTerm(0.00001664146, 4.25018630147, 1577.34354244780),
    VsopTerm(0.00001438387, 4.15745084182, 9683.59458111640),
    VsopTerm(0.00001200521, 6.15357116043, 30639.85663863300)
))

private val vsopLonVenus1 = VsopSeries(arrayOf(
    VsopTerm(10213.28554621638, 0.00000000000, 0.00000000000),
    VsopTerm(0.00095617813, 2.46406511110, 10213.28554621100),
    VsopTerm(0.00007787201, 0.62478482220, 20426.57109242200)
))

private val vsopLonVenus = VsopFormula(arrayOf(
    vsopLonVenus0,
    vsopLonVenus1
))

private val vsopLatVenus0 = VsopSeries(arrayOf(
    VsopTerm(0.05923638472, 0.26702775812, 10213.28554621100),
    VsopTerm(0.00040107978, 1.14737178112, 20426.57109242200),
    VsopTerm(0.00032814918, 3.14159265359, 0.00000000000)
))

private val vsopLatVenus1 = VsopSeries(arrayOf(
    VsopTerm(0.00287821243, 1.88964962838, 10213.28554621100)
))

private val vsopLatVenus = VsopFormula(arrayOf(
    vsopLatVenus0,
    vsopLatVenus1
))

private val vsopRadVenus0 = VsopSeries(arrayOf(
    VsopTerm(0.72334820891, 0.00000000000, 0.00000000000),
    VsopTerm(0.00489824182, 4.02151831717, 10213.28554621100),
    VsopTerm(0.00001658058, 4.90206728031, 20426.57109242200),
    VsopTerm(0.00001378043, 1.12846591367, 11790.62908865880),
    VsopTerm(0.00001632096, 2.84548795207, 7860.41939243920),
    VsopTerm(0.00000498395, 2.58682193892, 9683.59458111640),
    VsopTerm(0.00000221985, 2.01346696541, 19367.18916223280),
    VsopTerm(0.00000237454, 2.55136053886, 15720.83878487840)
))

private val vsopRadVenus1 = VsopSeries(arrayOf(
    VsopTerm(0.00034551041, 0.89198706276, 10213.28554621100)
))

private val vsopRadVenus = VsopFormula(arrayOf(
    vsopRadVenus0,
    vsopRadVenus1
))


private val vsopLonEarth0 = VsopSeries(arrayOf(
    VsopTerm(1.75347045673, 0.00000000000, 0.00000000000),
    VsopTerm(0.03341656453, 4.66925680415, 6283.07584999140),
    VsopTerm(0.00034894275, 4.62610242189, 12566.15169998280),
    VsopTerm(0.00003417572, 2.82886579754, 3.52311834900),
    VsopTerm(0.00003497056, 2.74411783405, 5753.38488489680),
    VsopTerm(0.00003135899, 3.62767041756, 77713.77146812050),
    VsopTerm(0.00002676218, 4.41808345438, 7860.41939243920),
    VsopTerm(0.00002342691, 6.13516214446, 3930.20969621960),
    VsopTerm(0.00001273165, 2.03709657878, 529.69096509460),
    VsopTerm(0.00001324294, 0.74246341673, 11506.76976979360),
    VsopTerm(0.00000901854, 2.04505446477, 26.29831979980),
    VsopTerm(0.00001199167, 1.10962946234, 1577.34354244780),
    VsopTerm(0.00000857223, 3.50849152283, 398.14900340820),
    VsopTerm(0.00000779786, 1.17882681962, 5223.69391980220),
    VsopTerm(0.00000990250, 5.23268072088, 5884.92684658320),
    VsopTerm(0.00000753141, 2.53339052847, 5507.55323866740),
    VsopTerm(0.00000505267, 4.58292599973, 18849.22754997420),
    VsopTerm(0.00000492392, 4.20505711826, 775.52261132400),
    VsopTerm(0.00000356672, 2.91954114478, 0.06731030280),
    VsopTerm(0.00000284125, 1.89869240932, 796.29800681640),
    VsopTerm(0.00000242879, 0.34481445893, 5486.77784317500),
    VsopTerm(0.00000317087, 5.84901948512, 11790.62908865880),
    VsopTerm(0.00000271112, 0.31486255375, 10977.07880469900),
    VsopTerm(0.00000206217, 4.80646631478, 2544.31441988340),
    VsopTerm(0.00000205478, 1.86953770281, 5573.14280143310),
    VsopTerm(0.00000202318, 2.45767790232, 6069.77675455340),
    VsopTerm(0.00000126225, 1.08295459501, 20.77539549240),
    VsopTerm(0.00000155516, 0.83306084617, 213.29909543800)
))

private val vsopLonEarth1 = VsopSeries(arrayOf(
    VsopTerm(6283.07584999140, 0.00000000000, 0.00000000000),
    VsopTerm(0.00206058863, 2.67823455808, 6283.07584999140),
    VsopTerm(0.00004303419, 2.63512233481, 12566.15169998280)
))

private val vsopLonEarth2 = VsopSeries(arrayOf(
    VsopTerm(0.00008721859, 1.07253635559, 6283.07584999140)
))

private val vsopLonEarth = VsopFormula(arrayOf(
    vsopLonEarth0,
    vsopLonEarth1,
    vsopLonEarth2
))

private val vsopLatEarth0 = VsopSeries(arrayOf(
))

private val vsopLatEarth1 = VsopSeries(arrayOf(
    VsopTerm(0.00227777722, 3.41376620530, 6283.07584999140),
    VsopTerm(0.00003805678, 3.37063423795, 12566.15169998280)
))

private val vsopLatEarth = VsopFormula(arrayOf(
    vsopLatEarth0,
    vsopLatEarth1
))

private val vsopRadEarth0 = VsopSeries(arrayOf(
    VsopTerm(1.00013988784, 0.00000000000, 0.00000000000),
    VsopTerm(0.01670699632, 3.09846350258, 6283.07584999140),
    VsopTerm(0.00013956024, 3.05524609456, 12566.15169998280),
    VsopTerm(0.00003083720, 5.19846674381, 77713.77146812050),
    VsopTerm(0.00001628463, 1.17387558054, 5753.38488489680),
    VsopTerm(0.00001575572, 2.84685214877, 7860.41939243920),
    VsopTerm(0.00000924799, 5.45292236722, 11506.76976979360),
    VsopTerm(0.00000542439, 4.56409151453, 3930.20969621960),
    VsopTerm(0.00000472110, 3.66100022149, 5884.92684658320),
    VsopTerm(0.00000085831, 1.27079125277, 161000.68573767410),
    VsopTerm(0.00000057056, 2.01374292245, 83996.84731811189),
    VsopTerm(0.00000055736, 5.24159799170, 71430.69561812909),
    VsopTerm(0.00000174844, 3.01193636733, 18849.22754997420),
    VsopTerm(0.00000243181, 4.27349530790, 11790.62908865880)
))

private val vsopRadEarth1 = VsopSeries(arrayOf(
    VsopTerm(0.00103018607, 1.10748968172, 6283.07584999140),
    VsopTerm(0.00001721238, 1.06442300386, 12566.15169998280)
))

private val vsopRadEarth2 = VsopSeries(arrayOf(
    VsopTerm(0.00004359385, 5.78455133808, 6283.07584999140)
))

private val vsopRadEarth = VsopFormula(arrayOf(
    vsopRadEarth0,
    vsopRadEarth1,
    vsopRadEarth2
))


private val vsopLonMars0 = VsopSeries(arrayOf(
    VsopTerm(6.20347711581, 0.00000000000, 0.00000000000),
    VsopTerm(0.18656368093, 5.05037100270, 3340.61242669980),
    VsopTerm(0.01108216816, 5.40099836344, 6681.22485339960),
    VsopTerm(0.00091798406, 5.75478744667, 10021.83728009940),
    VsopTerm(0.00027744987, 5.97049513147, 3.52311834900),
    VsopTerm(0.00010610235, 2.93958560338, 2281.23049651060),
    VsopTerm(0.00012315897, 0.84956094002, 2810.92146160520),
    VsopTerm(0.00008926784, 4.15697846427, 0.01725365220),
    VsopTerm(0.00008715691, 6.11005153139, 13362.44970679920),
    VsopTerm(0.00006797556, 0.36462229657, 398.14900340820),
    VsopTerm(0.00007774872, 3.33968761376, 5621.84292321040),
    VsopTerm(0.00003575078, 1.66186505710, 2544.31441988340),
    VsopTerm(0.00004161108, 0.22814971327, 2942.46342329160),
    VsopTerm(0.00003075252, 0.85696614132, 191.44826611160),
    VsopTerm(0.00002628117, 0.64806124465, 3337.08930835080),
    VsopTerm(0.00002937546, 6.07893711402, 0.06731030280),
    VsopTerm(0.00002389414, 5.03896442664, 796.29800681640),
    VsopTerm(0.00002579844, 0.02996736156, 3344.13554504880),
    VsopTerm(0.00001528141, 1.14979301996, 6151.53388830500),
    VsopTerm(0.00001798806, 0.65634057445, 529.69096509460),
    VsopTerm(0.00001264357, 3.62275122593, 5092.15195811580),
    VsopTerm(0.00001286228, 3.06796065034, 2146.16541647520),
    VsopTerm(0.00001546404, 2.91579701718, 1751.53953141600),
    VsopTerm(0.00001024902, 3.69334099279, 8962.45534991020),
    VsopTerm(0.00000891566, 0.18293837498, 16703.06213349900),
    VsopTerm(0.00000858759, 2.40093811940, 2914.01423582380),
    VsopTerm(0.00000832715, 2.46418619474, 3340.59517304760),
    VsopTerm(0.00000832720, 4.49495782139, 3340.62968035200),
    VsopTerm(0.00000712902, 3.66335473479, 1059.38193018920),
    VsopTerm(0.00000748723, 3.82248614017, 155.42039943420),
    VsopTerm(0.00000723861, 0.67497311481, 3738.76143010800),
    VsopTerm(0.00000635548, 2.92182225127, 8432.76438481560),
    VsopTerm(0.00000655162, 0.48864064125, 3127.31333126180),
    VsopTerm(0.00000550474, 3.81001042328, 0.98032106820),
    VsopTerm(0.00000552750, 4.47479317037, 1748.01641306700),
    VsopTerm(0.00000425966, 0.55364317304, 6283.07584999140),
    VsopTerm(0.00000415131, 0.49662285038, 213.29909543800),
    VsopTerm(0.00000472167, 3.62547124025, 1194.44701022460),
    VsopTerm(0.00000306551, 0.38052848348, 6684.74797174860),
    VsopTerm(0.00000312141, 0.99853944405, 6677.70173505060),
    VsopTerm(0.00000293198, 4.22131299634, 20.77539549240),
    VsopTerm(0.00000302375, 4.48618007156, 3532.06069281140),
    VsopTerm(0.00000274027, 0.54222167059, 3340.54511639700),
    VsopTerm(0.00000281079, 5.88163521788, 1349.86740965880),
    VsopTerm(0.00000231183, 1.28242156993, 3870.30339179440),
    VsopTerm(0.00000283602, 5.76885434940, 3149.16416058820),
    VsopTerm(0.00000236117, 5.75503217933, 3333.49887969900),
    VsopTerm(0.00000274033, 0.13372524985, 3340.67973700260),
    VsopTerm(0.00000299395, 2.78323740866, 6254.62666252360)
))

private val vsopLonMars1 = VsopSeries(arrayOf(
    VsopTerm(3340.61242700512, 0.00000000000, 0.00000000000),
    VsopTerm(0.01457554523, 3.60433733236, 3340.61242669980),
    VsopTerm(0.00168414711, 3.92318567804, 6681.22485339960),
    VsopTerm(0.00020622975, 4.26108844583, 10021.83728009940),
    VsopTerm(0.00003452392, 4.73210393190, 3.52311834900),
    VsopTerm(0.00002586332, 4.60670058555, 13362.44970679920),
    VsopTerm(0.00000841535, 4.45864030426, 2281.23049651060)
))

private val vsopLonMars2 = VsopSeries(arrayOf(
    VsopTerm(0.00058152577, 2.04961712429, 3340.61242669980),
    VsopTerm(0.00013459579, 2.45738706163, 6681.22485339960)
))

private val vsopLonMars = VsopFormula(arrayOf(
    vsopLonMars0,
    vsopLonMars1,
    vsopLonMars2
))

private val vsopLatMars0 = VsopSeries(arrayOf(
    VsopTerm(0.03197134986, 3.76832042431, 3340.61242669980),
    VsopTerm(0.00298033234, 4.10616996305, 6681.22485339960),
    VsopTerm(0.00289104742, 0.00000000000, 0.00000000000),
    VsopTerm(0.00031365539, 4.44651053090, 10021.83728009940),
    VsopTerm(0.00003484100, 4.78812549260, 13362.44970679920)
))

private val vsopLatMars1 = VsopSeries(arrayOf(
    VsopTerm(0.00217310991, 6.04472194776, 3340.61242669980),
    VsopTerm(0.00020976948, 3.14159265359, 0.00000000000),
    VsopTerm(0.00012834709, 1.60810667915, 6681.22485339960)
))

private val vsopLatMars = VsopFormula(arrayOf(
    vsopLatMars0,
    vsopLatMars1
))

private val vsopRadMars0 = VsopSeries(arrayOf(
    VsopTerm(1.53033488271, 0.00000000000, 0.00000000000),
    VsopTerm(0.14184953160, 3.47971283528, 3340.61242669980),
    VsopTerm(0.00660776362, 3.81783443019, 6681.22485339960),
    VsopTerm(0.00046179117, 4.15595316782, 10021.83728009940),
    VsopTerm(0.00008109733, 5.55958416318, 2810.92146160520),
    VsopTerm(0.00007485318, 1.77239078402, 5621.84292321040),
    VsopTerm(0.00005523191, 1.36436303770, 2281.23049651060),
    VsopTerm(0.00003825160, 4.49407183687, 13362.44970679920),
    VsopTerm(0.00002306537, 0.09081579001, 2544.31441988340),
    VsopTerm(0.00001999396, 5.36059617709, 3337.08930835080),
    VsopTerm(0.00002484394, 4.92545639920, 2942.46342329160),
    VsopTerm(0.00001960195, 4.74249437639, 3344.13554504880),
    VsopTerm(0.00001167119, 2.11260868341, 5092.15195811580),
    VsopTerm(0.00001102816, 5.00908403998, 398.14900340820),
    VsopTerm(0.00000899066, 4.40791133207, 529.69096509460),
    VsopTerm(0.00000992252, 5.83861961952, 6151.53388830500),
    VsopTerm(0.00000807354, 2.10217065501, 1059.38193018920),
    VsopTerm(0.00000797915, 3.44839203899, 796.29800681640),
    VsopTerm(0.00000740975, 1.49906336885, 2146.16541647520)
))

private val vsopRadMars1 = VsopSeries(arrayOf(
    VsopTerm(0.01107433345, 2.03250524857, 3340.61242669980),
    VsopTerm(0.00103175887, 2.37071847807, 6681.22485339960),
    VsopTerm(0.00012877200, 0.00000000000, 0.00000000000),
    VsopTerm(0.00010815880, 2.70888095665, 10021.83728009940)
))

private val vsopRadMars2 = VsopSeries(arrayOf(
    VsopTerm(0.00044242249, 0.47930604954, 3340.61242669980),
    VsopTerm(0.00008138042, 0.86998389204, 6681.22485339960)
))

private val vsopRadMars = VsopFormula(arrayOf(
    vsopRadMars0,
    vsopRadMars1,
    vsopRadMars2
))


private val vsopLonJupiter0 = VsopSeries(arrayOf(
    VsopTerm(0.59954691494, 0.00000000000, 0.00000000000),
    VsopTerm(0.09695898719, 5.06191793158, 529.69096509460),
    VsopTerm(0.00573610142, 1.44406205629, 7.11354700080),
    VsopTerm(0.00306389205, 5.41734730184, 1059.38193018920),
    VsopTerm(0.00097178296, 4.14264726552, 632.78373931320),
    VsopTerm(0.00072903078, 3.64042916389, 522.57741809380),
    VsopTerm(0.00064263975, 3.41145165351, 103.09277421860),
    VsopTerm(0.00039806064, 2.29376740788, 419.48464387520),
    VsopTerm(0.00038857767, 1.27231755835, 316.39186965660),
    VsopTerm(0.00027964629, 1.78454591820, 536.80451209540),
    VsopTerm(0.00013589730, 5.77481040790, 1589.07289528380),
    VsopTerm(0.00008246349, 3.58227925840, 206.18554843720),
    VsopTerm(0.00008768704, 3.63000308199, 949.17560896980),
    VsopTerm(0.00007368042, 5.08101194270, 735.87651353180),
    VsopTerm(0.00006263150, 0.02497628807, 213.29909543800),
    VsopTerm(0.00006114062, 4.51319998626, 1162.47470440780),
    VsopTerm(0.00004905396, 1.32084470588, 110.20632121940),
    VsopTerm(0.00005305285, 1.30671216791, 14.22709400160),
    VsopTerm(0.00005305441, 4.18625634012, 1052.26838318840),
    VsopTerm(0.00004647248, 4.69958103684, 3.93215326310),
    VsopTerm(0.00003045023, 4.31676431084, 426.59819087600),
    VsopTerm(0.00002609999, 1.56667394063, 846.08283475120),
    VsopTerm(0.00002028191, 1.06376530715, 3.18139373770),
    VsopTerm(0.00001764763, 2.14148655117, 1066.49547719000),
    VsopTerm(0.00001722972, 3.88036268267, 1265.56747862640),
    VsopTerm(0.00001920945, 0.97168196472, 639.89728631400),
    VsopTerm(0.00001633223, 3.58201833555, 515.46387109300),
    VsopTerm(0.00001431999, 4.29685556046, 625.67019231240),
    VsopTerm(0.00000973272, 4.09764549134, 95.97922721780)
))

private val vsopLonJupiter1 = VsopSeries(arrayOf(
    VsopTerm(529.69096508814, 0.00000000000, 0.00000000000),
    VsopTerm(0.00489503243, 4.22082939470, 529.69096509460),
    VsopTerm(0.00228917222, 6.02646855621, 7.11354700080),
    VsopTerm(0.00030099479, 4.54540782858, 1059.38193018920),
    VsopTerm(0.00020720920, 5.45943156902, 522.57741809380),
    VsopTerm(0.00012103653, 0.16994816098, 536.80451209540),
    VsopTerm(0.00006067987, 4.42422292017, 103.09277421860),
    VsopTerm(0.00005433968, 3.98480737746, 419.48464387520),
    VsopTerm(0.00004237744, 5.89008707199, 14.22709400160)
))

private val vsopLonJupiter2 = VsopSeries(arrayOf(
    VsopTerm(0.00047233601, 4.32148536482, 7.11354700080),
    VsopTerm(0.00030649436, 2.92977788700, 529.69096509460),
    VsopTerm(0.00014837605, 3.14159265359, 0.00000000000)
))

private val vsopLonJupiter = VsopFormula(arrayOf(
    vsopLonJupiter0,
    vsopLonJupiter1,
    vsopLonJupiter2
))

private val vsopLatJupiter0 = VsopSeries(arrayOf(
    VsopTerm(0.02268615702, 3.55852606721, 529.69096509460),
    VsopTerm(0.00109971634, 3.90809347197, 1059.38193018920),
    VsopTerm(0.00110090358, 0.00000000000, 0.00000000000),
    VsopTerm(0.00008101428, 3.60509572885, 522.57741809380),
    VsopTerm(0.00006043996, 4.25883108339, 1589.07289528380),
    VsopTerm(0.00006437782, 0.30627119215, 536.80451209540)
))

private val vsopLatJupiter1 = VsopSeries(arrayOf(
    VsopTerm(0.00078203446, 1.52377859742, 529.69096509460)
))

private val vsopLatJupiter = VsopFormula(arrayOf(
    vsopLatJupiter0,
    vsopLatJupiter1
))

private val vsopRadJupiter0 = VsopSeries(arrayOf(
    VsopTerm(5.20887429326, 0.00000000000, 0.00000000000),
    VsopTerm(0.25209327119, 3.49108639871, 529.69096509460),
    VsopTerm(0.00610599976, 3.84115365948, 1059.38193018920),
    VsopTerm(0.00282029458, 2.57419881293, 632.78373931320),
    VsopTerm(0.00187647346, 2.07590383214, 522.57741809380),
    VsopTerm(0.00086792905, 0.71001145545, 419.48464387520),
    VsopTerm(0.00072062974, 0.21465724607, 536.80451209540),
    VsopTerm(0.00065517248, 5.97995884790, 316.39186965660),
    VsopTerm(0.00029134542, 1.67759379655, 103.09277421860),
    VsopTerm(0.00030135335, 2.16132003734, 949.17560896980),
    VsopTerm(0.00023453271, 3.54023522184, 735.87651353180),
    VsopTerm(0.00022283743, 4.19362594399, 1589.07289528380),
    VsopTerm(0.00023947298, 0.27458037480, 7.11354700080),
    VsopTerm(0.00013032614, 2.96042965363, 1162.47470440780),
    VsopTerm(0.00009703360, 1.90669633585, 206.18554843720),
    VsopTerm(0.00012749023, 2.71550286592, 1052.26838318840),
    VsopTerm(0.00007057931, 2.18184839926, 1265.56747862640),
    VsopTerm(0.00006137703, 6.26418240033, 846.08283475120),
    VsopTerm(0.00002616976, 2.00994012876, 1581.95934828300)
))

private val vsopRadJupiter1 = VsopSeries(arrayOf(
    VsopTerm(0.01271801520, 2.64937512894, 529.69096509460),
    VsopTerm(0.00061661816, 3.00076460387, 1059.38193018920),
    VsopTerm(0.00053443713, 3.89717383175, 522.57741809380),
    VsopTerm(0.00031185171, 4.88276958012, 536.80451209540),
    VsopTerm(0.00041390269, 0.00000000000, 0.00000000000)
))

private val vsopRadJupiter = VsopFormula(arrayOf(
    vsopRadJupiter0,
    vsopRadJupiter1
))


private val vsopLonSaturn0 = VsopSeries(arrayOf(
    VsopTerm(0.87401354025, 0.00000000000, 0.00000000000),
    VsopTerm(0.11107659762, 3.96205090159, 213.29909543800),
    VsopTerm(0.01414150957, 4.58581516874, 7.11354700080),
    VsopTerm(0.00398379389, 0.52112032699, 206.18554843720),
    VsopTerm(0.00350769243, 3.30329907896, 426.59819087600),
    VsopTerm(0.00206816305, 0.24658372002, 103.09277421860),
    VsopTerm(0.00079271300, 3.84007056878, 220.41264243880),
    VsopTerm(0.00023990355, 4.66976924553, 110.20632121940),
    VsopTerm(0.00016573588, 0.43719228296, 419.48464387520),
    VsopTerm(0.00014906995, 5.76903183869, 316.39186965660),
    VsopTerm(0.00015820290, 0.93809155235, 632.78373931320),
    VsopTerm(0.00014609559, 1.56518472000, 3.93215326310),
    VsopTerm(0.00013160301, 4.44891291899, 14.22709400160),
    VsopTerm(0.00015053543, 2.71669915667, 639.89728631400),
    VsopTerm(0.00013005299, 5.98119023644, 11.04570026390),
    VsopTerm(0.00010725067, 3.12939523827, 202.25339517410),
    VsopTerm(0.00005863206, 0.23656938524, 529.69096509460),
    VsopTerm(0.00005227757, 4.20783365759, 3.18139373770),
    VsopTerm(0.00006126317, 1.76328667907, 277.03499374140),
    VsopTerm(0.00005019687, 3.17787728405, 433.71173787680),
    VsopTerm(0.00004592550, 0.61977744975, 199.07200143640),
    VsopTerm(0.00004005867, 2.24479718502, 63.73589830340),
    VsopTerm(0.00002953796, 0.98280366998, 95.97922721780),
    VsopTerm(0.00003873670, 3.22283226966, 138.51749687070),
    VsopTerm(0.00002461186, 2.03163875071, 735.87651353180),
    VsopTerm(0.00003269484, 0.77492638211, 949.17560896980),
    VsopTerm(0.00001758145, 3.26580109940, 522.57741809380),
    VsopTerm(0.00001640172, 5.50504453050, 846.08283475120),
    VsopTerm(0.00001391327, 4.02333150505, 323.50541665740),
    VsopTerm(0.00001580648, 4.37265307169, 309.27832265580),
    VsopTerm(0.00001123498, 2.83726798446, 415.55249061210),
    VsopTerm(0.00001017275, 3.71700135395, 227.52618943960),
    VsopTerm(0.00000848642, 3.19150170830, 209.36694217490)
))

private val vsopLonSaturn1 = VsopSeries(arrayOf(
    VsopTerm(213.29909521690, 0.00000000000, 0.00000000000),
    VsopTerm(0.01297370862, 1.82834923978, 213.29909543800),
    VsopTerm(0.00564345393, 2.88499717272, 7.11354700080),
    VsopTerm(0.00093734369, 1.06311793502, 426.59819087600),
    VsopTerm(0.00107674962, 2.27769131009, 206.18554843720),
    VsopTerm(0.00040244455, 2.04108104671, 220.41264243880),
    VsopTerm(0.00019941774, 1.27954390470, 103.09277421860),
    VsopTerm(0.00010511678, 2.74880342130, 14.22709400160),
    VsopTerm(0.00006416106, 0.38238295041, 639.89728631400),
    VsopTerm(0.00004848994, 2.43037610229, 419.48464387520),
    VsopTerm(0.00004056892, 2.92133209468, 110.20632121940),
    VsopTerm(0.00003768635, 3.64965330780, 3.93215326310)
))

private val vsopLonSaturn2 = VsopSeries(arrayOf(
    VsopTerm(0.00116441330, 1.17988132879, 7.11354700080),
    VsopTerm(0.00091841837, 0.07325195840, 213.29909543800),
    VsopTerm(0.00036661728, 0.00000000000, 0.00000000000),
    VsopTerm(0.00015274496, 4.06493179167, 206.18554843720)
))

private val vsopLonSaturn = VsopFormula(arrayOf(
    vsopLonSaturn0,
    vsopLonSaturn1,
    vsopLonSaturn2
))

private val vsopLatSaturn0 = VsopSeries(arrayOf(
    VsopTerm(0.04330678039, 3.60284428399, 213.29909543800),
    VsopTerm(0.00240348302, 2.85238489373, 426.59819087600),
    VsopTerm(0.00084745939, 0.00000000000, 0.00000000000),
    VsopTerm(0.00030863357, 3.48441504555, 220.41264243880),
    VsopTerm(0.00034116062, 0.57297307557, 206.18554843720),
    VsopTerm(0.00014734070, 2.11846596715, 639.89728631400),
    VsopTerm(0.00009916667, 5.79003188904, 419.48464387520),
    VsopTerm(0.00006993564, 4.73604689720, 7.11354700080),
    VsopTerm(0.00004807588, 5.43305312061, 316.39186965660)
))

private val vsopLatSaturn1 = VsopSeries(arrayOf(
    VsopTerm(0.00198927992, 4.93901017903, 213.29909543800),
    VsopTerm(0.00036947916, 3.14159265359, 0.00000000000),
    VsopTerm(0.00017966989, 0.51979431110, 426.59819087600)
))

private val vsopLatSaturn = VsopFormula(arrayOf(
    vsopLatSaturn0,
    vsopLatSaturn1
))

private val vsopRadSaturn0 = VsopSeries(arrayOf(
    VsopTerm(9.55758135486, 0.00000000000, 0.00000000000),
    VsopTerm(0.52921382865, 2.39226219573, 213.29909543800),
    VsopTerm(0.01873679867, 5.23549604660, 206.18554843720),
    VsopTerm(0.01464663929, 1.64763042902, 426.59819087600),
    VsopTerm(0.00821891141, 5.93520042303, 316.39186965660),
    VsopTerm(0.00547506923, 5.01532618980, 103.09277421860),
    VsopTerm(0.00371684650, 2.27114821115, 220.41264243880),
    VsopTerm(0.00361778765, 3.13904301847, 7.11354700080),
    VsopTerm(0.00140617506, 5.70406606781, 632.78373931320),
    VsopTerm(0.00108974848, 3.29313390175, 110.20632121940),
    VsopTerm(0.00069006962, 5.94099540992, 419.48464387520),
    VsopTerm(0.00061053367, 0.94037691801, 639.89728631400),
    VsopTerm(0.00048913294, 1.55733638681, 202.25339517410),
    VsopTerm(0.00034143772, 0.19519102597, 277.03499374140),
    VsopTerm(0.00032401773, 5.47084567016, 949.17560896980),
    VsopTerm(0.00020936596, 0.46349251129, 735.87651353180),
    VsopTerm(0.00009796004, 5.20477537945, 1265.56747862640),
    VsopTerm(0.00011993338, 5.98050967385, 846.08283475120),
    VsopTerm(0.00020839300, 1.52102476129, 433.71173787680),
    VsopTerm(0.00015298404, 3.05943814940, 529.69096509460),
    VsopTerm(0.00006465823, 0.17732249942, 1052.26838318840),
    VsopTerm(0.00011380257, 1.73105427040, 522.57741809380),
    VsopTerm(0.00003419618, 4.94550542171, 1581.95934828300)
))

private val vsopRadSaturn1 = VsopSeries(arrayOf(
    VsopTerm(0.06182981340, 0.25843511480, 213.29909543800),
    VsopTerm(0.00506577242, 0.71114625261, 206.18554843720),
    VsopTerm(0.00341394029, 5.79635741658, 426.59819087600),
    VsopTerm(0.00188491195, 0.47215589652, 220.41264243880),
    VsopTerm(0.00186261486, 3.14159265359, 0.00000000000),
    VsopTerm(0.00143891146, 1.40744822888, 7.11354700080)
))

private val vsopRadSaturn2 = VsopSeries(arrayOf(
    VsopTerm(0.00436902572, 4.78671677509, 213.29909543800)
))

private val vsopRadSaturn = VsopFormula(arrayOf(
    vsopRadSaturn0,
    vsopRadSaturn1,
    vsopRadSaturn2
))


private val vsopLonUranus0 = VsopSeries(arrayOf(
    VsopTerm(5.48129294297, 0.00000000000, 0.00000000000),
    VsopTerm(0.09260408234, 0.89106421507, 74.78159856730),
    VsopTerm(0.01504247898, 3.62719260920, 1.48447270830),
    VsopTerm(0.00365981674, 1.89962179044, 73.29712585900),
    VsopTerm(0.00272328168, 3.35823706307, 149.56319713460),
    VsopTerm(0.00070328461, 5.39254450063, 63.73589830340),
    VsopTerm(0.00068892678, 6.09292483287, 76.26607127560),
    VsopTerm(0.00061998615, 2.26952066061, 2.96894541660),
    VsopTerm(0.00061950719, 2.85098872691, 11.04570026390),
    VsopTerm(0.00026468770, 3.14152083966, 71.81265315070),
    VsopTerm(0.00025710476, 6.11379840493, 454.90936652730),
    VsopTerm(0.00021078850, 4.36059339067, 148.07872442630),
    VsopTerm(0.00017818647, 1.74436930289, 36.64856292950),
    VsopTerm(0.00014613507, 4.73732166022, 3.93215326310),
    VsopTerm(0.00011162509, 5.82681796350, 224.34479570190),
    VsopTerm(0.00010997910, 0.48865004018, 138.51749687070),
    VsopTerm(0.00009527478, 2.95516862826, 35.16409022120),
    VsopTerm(0.00007545601, 5.23626582400, 109.94568878850),
    VsopTerm(0.00004220241, 3.23328220918, 70.84944530420),
    VsopTerm(0.00004051900, 2.27755017300, 151.04766984290),
    VsopTerm(0.00003354596, 1.06549007380, 4.45341812490),
    VsopTerm(0.00002926718, 4.62903718891, 9.56122755560),
    VsopTerm(0.00003490340, 5.48306144511, 146.59425171800),
    VsopTerm(0.00003144069, 4.75199570434, 77.75054398390),
    VsopTerm(0.00002922333, 5.35235361027, 85.82729883120),
    VsopTerm(0.00002272788, 4.36600400036, 70.32818044240),
    VsopTerm(0.00002051219, 1.51773566586, 0.11187458460),
    VsopTerm(0.00002148602, 0.60745949945, 38.13303563780),
    VsopTerm(0.00001991643, 4.92437588682, 277.03499374140),
    VsopTerm(0.00001376226, 2.04283539351, 65.22037101170),
    VsopTerm(0.00001666902, 3.62744066769, 380.12776796000),
    VsopTerm(0.00001284107, 3.11347961505, 202.25339517410),
    VsopTerm(0.00001150429, 0.93343589092, 3.18139373770),
    VsopTerm(0.00001533221, 2.58594681212, 52.69019803950),
    VsopTerm(0.00001281604, 0.54271272721, 222.86032299360),
    VsopTerm(0.00001372139, 4.19641530878, 111.43016149680),
    VsopTerm(0.00001221029, 0.19900650030, 108.46121608020),
    VsopTerm(0.00000946181, 1.19253165736, 127.47179660680),
    VsopTerm(0.00001150989, 4.17898916639, 33.67961751290)
))

private val vsopLonUranus1 = VsopSeries(arrayOf(
    VsopTerm(74.78159860910, 0.00000000000, 0.00000000000),
    VsopTerm(0.00154332863, 5.24158770553, 74.78159856730),
    VsopTerm(0.00024456474, 1.71260334156, 1.48447270830),
    VsopTerm(0.00009258442, 0.42829732350, 11.04570026390),
    VsopTerm(0.00008265977, 1.50218091379, 63.73589830340),
    VsopTerm(0.00009150160, 1.41213765216, 149.56319713460)
))

private val vsopLonUranus = VsopFormula(arrayOf(
    vsopLonUranus0,
    vsopLonUranus1
))

private val vsopLatUranus0 = VsopSeries(arrayOf(
    VsopTerm(0.01346277648, 2.61877810547, 74.78159856730),
    VsopTerm(0.00062341400, 5.08111189648, 149.56319713460),
    VsopTerm(0.00061601196, 3.14159265359, 0.00000000000),
    VsopTerm(0.00009963722, 1.61603805646, 76.26607127560),
    VsopTerm(0.00009926160, 0.57630380333, 73.29712585900)
))

private val vsopLatUranus1 = VsopSeries(arrayOf(
    VsopTerm(0.00034101978, 0.01321929936, 74.78159856730)
))

private val vsopLatUranus = VsopFormula(arrayOf(
    vsopLatUranus0,
    vsopLatUranus1
))

private val vsopRadUranus0 = VsopSeries(arrayOf(
    VsopTerm(19.21264847206, 0.00000000000, 0.00000000000),
    VsopTerm(0.88784984413, 5.60377527014, 74.78159856730),
    VsopTerm(0.03440836062, 0.32836099706, 73.29712585900),
    VsopTerm(0.02055653860, 1.78295159330, 149.56319713460),
    VsopTerm(0.00649322410, 4.52247285911, 76.26607127560),
    VsopTerm(0.00602247865, 3.86003823674, 63.73589830340),
    VsopTerm(0.00496404167, 1.40139935333, 454.90936652730),
    VsopTerm(0.00338525369, 1.58002770318, 138.51749687070),
    VsopTerm(0.00243509114, 1.57086606044, 71.81265315070),
    VsopTerm(0.00190522303, 1.99809394714, 1.48447270830),
    VsopTerm(0.00161858838, 2.79137786799, 148.07872442630),
    VsopTerm(0.00143706183, 1.38368544947, 11.04570026390),
    VsopTerm(0.00093192405, 0.17437220467, 36.64856292950),
    VsopTerm(0.00071424548, 4.24509236074, 224.34479570190),
    VsopTerm(0.00089806014, 3.66105364565, 109.94568878850),
    VsopTerm(0.00039009723, 1.66971401684, 70.84944530420),
    VsopTerm(0.00046677296, 1.39976401694, 35.16409022120),
    VsopTerm(0.00039025624, 3.36234773834, 277.03499374140),
    VsopTerm(0.00036755274, 3.88649278513, 146.59425171800),
    VsopTerm(0.00030348723, 0.70100838798, 151.04766984290),
    VsopTerm(0.00029156413, 3.18056336700, 77.75054398390),
    VsopTerm(0.00022637073, 0.72518687029, 529.69096509460),
    VsopTerm(0.00011959076, 1.75043392140, 984.60033162190),
    VsopTerm(0.00025620756, 5.25656086672, 380.12776796000)
))

private val vsopRadUranus1 = VsopSeries(arrayOf(
    VsopTerm(0.01479896629, 3.67205697578, 74.78159856730)
))

private val vsopRadUranus = VsopFormula(arrayOf(
    vsopRadUranus0,
    vsopRadUranus1
))


private val vsopLonNeptune0 = VsopSeries(arrayOf(
    VsopTerm(5.31188633046, 0.00000000000, 0.00000000000),
    VsopTerm(0.01798475530, 2.90101273890, 38.13303563780),
    VsopTerm(0.01019727652, 0.48580922867, 1.48447270830),
    VsopTerm(0.00124531845, 4.83008090676, 36.64856292950),
    VsopTerm(0.00042064466, 5.41054993053, 2.96894541660),
    VsopTerm(0.00037714584, 6.09221808686, 35.16409022120),
    VsopTerm(0.00033784738, 1.24488874087, 76.26607127560),
    VsopTerm(0.00016482741, 0.00007727998, 491.55792945680),
    VsopTerm(0.00009198584, 4.93747051954, 39.61750834610),
    VsopTerm(0.00008994250, 0.27462171806, 175.16605980020)
))

private val vsopLonNeptune1 = VsopSeries(arrayOf(
    VsopTerm(38.13303563957, 0.00000000000, 0.00000000000),
    VsopTerm(0.00016604172, 4.86323329249, 1.48447270830),
    VsopTerm(0.00015744045, 2.27887427527, 38.13303563780)
))

private val vsopLonNeptune = VsopFormula(arrayOf(
    vsopLonNeptune0,
    vsopLonNeptune1
))

private val vsopLatNeptune0 = VsopSeries(arrayOf(
    VsopTerm(0.03088622933, 1.44104372644, 38.13303563780),
    VsopTerm(0.00027780087, 5.91271884599, 76.26607127560),
    VsopTerm(0.00027623609, 0.00000000000, 0.00000000000),
    VsopTerm(0.00015355489, 2.52123799551, 36.64856292950),
    VsopTerm(0.00015448133, 3.50877079215, 39.61750834610)
))

private val vsopLatNeptune = VsopFormula(arrayOf(
    vsopLatNeptune0
))

private val vsopRadNeptune0 = VsopSeries(arrayOf(
    VsopTerm(30.07013205828, 0.00000000000, 0.00000000000),
    VsopTerm(0.27062259632, 1.32999459377, 38.13303563780),
    VsopTerm(0.01691764014, 3.25186135653, 36.64856292950),
    VsopTerm(0.00807830553, 5.18592878704, 1.48447270830),
    VsopTerm(0.00537760510, 4.52113935896, 35.16409022120),
    VsopTerm(0.00495725141, 1.57105641650, 491.55792945680),
    VsopTerm(0.00274571975, 1.84552258866, 175.16605980020),
    VsopTerm(0.00012012320, 1.92059384991, 1021.24889455140),
    VsopTerm(0.00121801746, 5.79754470298, 76.26607127560),
    VsopTerm(0.00100896068, 0.37702724930, 73.29712585900),
    VsopTerm(0.00135134092, 3.37220609835, 39.61750834610),
    VsopTerm(0.00007571796, 1.07149207335, 388.46515523820)
))

private val vsopRadNeptune = VsopFormula(arrayOf(
    vsopRadNeptune0
))



//---------------------------------------------------------------------------------------
// Geocentric Moon

private fun addSolarTerms(context: MoonContext) {
    context.addSol(    13.9020,    14.0600,    -0.0010,     0.2607, 0, 0, 0, 4)
    context.addSol(     0.4030,    -4.0100,     0.3940,     0.0023, 0, 0, 0, 3)
    context.addSol(  2369.9120,  2373.3600,     0.6010,    28.2333, 0, 0, 0, 2)
    context.addSol(  -125.1540,  -112.7900,    -0.7250,    -0.9781, 0, 0, 0, 1)
    context.addSol(     1.9790,     6.9800,    -0.4450,     0.0433, 1, 0, 0, 4)
    context.addSol(   191.9530,   192.7200,     0.0290,     3.0861, 1, 0, 0, 2)
    context.addSol(    -8.4660,   -13.5100,     0.4550,    -0.1093, 1, 0, 0, 1)
    context.addSol( 22639.5000, 22609.0700,     0.0790,   186.5398, 1, 0, 0, 0)
    context.addSol(    18.6090,     3.5900,    -0.0940,     0.0118, 1, 0, 0,-1)
    context.addSol( -4586.4650, -4578.1300,    -0.0770,    34.3117, 1, 0, 0,-2)
    context.addSol(     3.2150,     5.4400,     0.1920,    -0.0386, 1, 0, 0,-3)
    context.addSol(   -38.4280,   -38.6400,     0.0010,     0.6008, 1, 0, 0,-4)
    context.addSol(    -0.3930,    -1.4300,    -0.0920,     0.0086, 1, 0, 0,-6)
    context.addSol(    -0.2890,    -1.5900,     0.1230,    -0.0053, 0, 1, 0, 4)
    context.addSol(   -24.4200,   -25.1000,     0.0400,    -0.3000, 0, 1, 0, 2)
    context.addSol(    18.0230,    17.9300,     0.0070,     0.1494, 0, 1, 0, 1)
    context.addSol(  -668.1460,  -126.9800,    -1.3020,    -0.3997, 0, 1, 0, 0)
    context.addSol(     0.5600,     0.3200,    -0.0010,    -0.0037, 0, 1, 0,-1)
    context.addSol(  -165.1450,  -165.0600,     0.0540,     1.9178, 0, 1, 0,-2)
    context.addSol(    -1.8770,    -6.4600,    -0.4160,     0.0339, 0, 1, 0,-4)
    context.addSol(     0.2130,     1.0200,    -0.0740,     0.0054, 2, 0, 0, 4)
    context.addSol(    14.3870,    14.7800,    -0.0170,     0.2833, 2, 0, 0, 2)
    context.addSol(    -0.5860,    -1.2000,     0.0540,    -0.0100, 2, 0, 0, 1)
    context.addSol(   769.0160,   767.9600,     0.1070,    10.1657, 2, 0, 0, 0)
    context.addSol(     1.7500,     2.0100,    -0.0180,     0.0155, 2, 0, 0,-1)
    context.addSol(  -211.6560,  -152.5300,     5.6790,    -0.3039, 2, 0, 0,-2)
    context.addSol(     1.2250,     0.9100,    -0.0300,    -0.0088, 2, 0, 0,-3)
    context.addSol(   -30.7730,   -34.0700,    -0.3080,     0.3722, 2, 0, 0,-4)
    context.addSol(    -0.5700,    -1.4000,    -0.0740,     0.0109, 2, 0, 0,-6)
    context.addSol(    -2.9210,   -11.7500,     0.7870,    -0.0484, 1, 1, 0, 2)
    context.addSol(     1.2670,     1.5200,    -0.0220,     0.0164, 1, 1, 0, 1)
    context.addSol(  -109.6730,  -115.1800,     0.4610,    -0.9490, 1, 1, 0, 0)
    context.addSol(  -205.9620,  -182.3600,     2.0560,     1.4437, 1, 1, 0,-2)
    context.addSol(     0.2330,     0.3600,     0.0120,    -0.0025, 1, 1, 0,-3)
    context.addSol(    -4.3910,    -9.6600,    -0.4710,     0.0673, 1, 1, 0,-4)
    context.addSol(     0.2830,     1.5300,    -0.1110,     0.0060, 1,-1, 0, 4)
    context.addSol(    14.5770,    31.7000,    -1.5400,     0.2302, 1,-1, 0, 2)
    context.addSol(   147.6870,   138.7600,     0.6790,     1.1528, 1,-1, 0, 0)
    context.addSol(    -1.0890,     0.5500,     0.0210,     0.0000, 1,-1, 0,-1)
    context.addSol(    28.4750,    23.5900,    -0.4430,    -0.2257, 1,-1, 0,-2)
    context.addSol(    -0.2760,    -0.3800,    -0.0060,    -0.0036, 1,-1, 0,-3)
    context.addSol(     0.6360,     2.2700,     0.1460,    -0.0102, 1,-1, 0,-4)
    context.addSol(    -0.1890,    -1.6800,     0.1310,    -0.0028, 0, 2, 0, 2)
    context.addSol(    -7.4860,    -0.6600,    -0.0370,    -0.0086, 0, 2, 0, 0)
    context.addSol(    -8.0960,   -16.3500,    -0.7400,     0.0918, 0, 2, 0,-2)
    context.addSol(    -5.7410,    -0.0400,     0.0000,    -0.0009, 0, 0, 2, 2)
    context.addSol(     0.2550,     0.0000,     0.0000,     0.0000, 0, 0, 2, 1)
    context.addSol(  -411.6080,    -0.2000,     0.0000,    -0.0124, 0, 0, 2, 0)
    context.addSol(     0.5840,     0.8400,     0.0000,     0.0071, 0, 0, 2,-1)
    context.addSol(   -55.1730,   -52.1400,     0.0000,    -0.1052, 0, 0, 2,-2)
    context.addSol(     0.2540,     0.2500,     0.0000,    -0.0017, 0, 0, 2,-3)
    context.addSol(     0.0250,    -1.6700,     0.0000,     0.0031, 0, 0, 2,-4)
    context.addSol(     1.0600,     2.9600,    -0.1660,     0.0243, 3, 0, 0, 2)
    context.addSol(    36.1240,    50.6400,    -1.3000,     0.6215, 3, 0, 0, 0)
    context.addSol(   -13.1930,   -16.4000,     0.2580,    -0.1187, 3, 0, 0,-2)
    context.addSol(    -1.1870,    -0.7400,     0.0420,     0.0074, 3, 0, 0,-4)
    context.addSol(    -0.2930,    -0.3100,    -0.0020,     0.0046, 3, 0, 0,-6)
    context.addSol(    -0.2900,    -1.4500,     0.1160,    -0.0051, 2, 1, 0, 2)
    context.addSol(    -7.6490,   -10.5600,     0.2590,    -0.1038, 2, 1, 0, 0)
    context.addSol(    -8.6270,    -7.5900,     0.0780,    -0.0192, 2, 1, 0,-2)
    context.addSol(    -2.7400,    -2.5400,     0.0220,     0.0324, 2, 1, 0,-4)
    context.addSol(     1.1810,     3.3200,    -0.2120,     0.0213, 2,-1, 0, 2)
    context.addSol(     9.7030,    11.6700,    -0.1510,     0.1268, 2,-1, 0, 0)
    context.addSol(    -0.3520,    -0.3700,     0.0010,    -0.0028, 2,-1, 0,-1)
    context.addSol(    -2.4940,    -1.1700,    -0.0030,    -0.0017, 2,-1, 0,-2)
    context.addSol(     0.3600,     0.2000,    -0.0120,    -0.0043, 2,-1, 0,-4)
    context.addSol(    -1.1670,    -1.2500,     0.0080,    -0.0106, 1, 2, 0, 0)
    context.addSol(    -7.4120,    -6.1200,     0.1170,     0.0484, 1, 2, 0,-2)
    context.addSol(    -0.3110,    -0.6500,    -0.0320,     0.0044, 1, 2, 0,-4)
    context.addSol(     0.7570,     1.8200,    -0.1050,     0.0112, 1,-2, 0, 2)
    context.addSol(     2.5800,     2.3200,     0.0270,     0.0196, 1,-2, 0, 0)
    context.addSol(     2.5330,     2.4000,    -0.0140,    -0.0212, 1,-2, 0,-2)
    context.addSol(    -0.3440,    -0.5700,    -0.0250,     0.0036, 0, 3, 0,-2)
    context.addSol(    -0.9920,    -0.0200,     0.0000,     0.0000, 1, 0, 2, 2)
    context.addSol(   -45.0990,    -0.0200,     0.0000,    -0.0010, 1, 0, 2, 0)
    context.addSol(    -0.1790,    -9.5200,     0.0000,    -0.0833, 1, 0, 2,-2)
    context.addSol(    -0.3010,    -0.3300,     0.0000,     0.0014, 1, 0, 2,-4)
    context.addSol(    -6.3820,    -3.3700,     0.0000,    -0.0481, 1, 0,-2, 2)
    context.addSol(    39.5280,    85.1300,     0.0000,    -0.7136, 1, 0,-2, 0)
    context.addSol(     9.3660,     0.7100,     0.0000,    -0.0112, 1, 0,-2,-2)
    context.addSol(     0.2020,     0.0200,     0.0000,     0.0000, 1, 0,-2,-4)
    context.addSol(     0.4150,     0.1000,     0.0000,     0.0013, 0, 1, 2, 0)
    context.addSol(    -2.1520,    -2.2600,     0.0000,    -0.0066, 0, 1, 2,-2)
    context.addSol(    -1.4400,    -1.3000,     0.0000,     0.0014, 0, 1,-2, 2)
    context.addSol(     0.3840,    -0.0400,     0.0000,     0.0000, 0, 1,-2,-2)
    context.addSol(     1.9380,     3.6000,    -0.1450,     0.0401, 4, 0, 0, 0)
    context.addSol(    -0.9520,    -1.5800,     0.0520,    -0.0130, 4, 0, 0,-2)
    context.addSol(    -0.5510,    -0.9400,     0.0320,    -0.0097, 3, 1, 0, 0)
    context.addSol(    -0.4820,    -0.5700,     0.0050,    -0.0045, 3, 1, 0,-2)
    context.addSol(     0.6810,     0.9600,    -0.0260,     0.0115, 3,-1, 0, 0)
    context.addSol(    -0.2970,    -0.2700,     0.0020,    -0.0009, 2, 2, 0,-2)
    context.addSol(     0.2540,     0.2100,    -0.0030,     0.0000, 2,-2, 0,-2)
    context.addSol(    -0.2500,    -0.2200,     0.0040,     0.0014, 1, 3, 0,-2)
    context.addSol(    -3.9960,     0.0000,     0.0000,     0.0004, 2, 0, 2, 0)
    context.addSol(     0.5570,    -0.7500,     0.0000,    -0.0090, 2, 0, 2,-2)
    context.addSol(    -0.4590,    -0.3800,     0.0000,    -0.0053, 2, 0,-2, 2)
    context.addSol(    -1.2980,     0.7400,     0.0000,     0.0004, 2, 0,-2, 0)
    context.addSol(     0.5380,     1.1400,     0.0000,    -0.0141, 2, 0,-2,-2)
    context.addSol(     0.2630,     0.0200,     0.0000,     0.0000, 1, 1, 2, 0)
    context.addSol(     0.4260,     0.0700,     0.0000,    -0.0006, 1, 1,-2,-2)
    context.addSol(    -0.3040,     0.0300,     0.0000,     0.0003, 1,-1, 2, 0)
    context.addSol(    -0.3720,    -0.1900,     0.0000,    -0.0027, 1,-1,-2, 2)
    context.addSol(     0.4180,     0.0000,     0.0000,     0.0000, 0, 0, 4, 0)
    context.addSol(    -0.3300,    -0.0400,     0.0000,     0.0000, 3, 0, 2, 0)
}

//---------------------------------------------------------------------------------------
// Pluto state table

private const val PLUTO_NUM_STATES = 51
private const val PLUTO_TIME_STEP  = 29200
private const val PLUTO_DT         = 146
private const val PLUTO_NSTEPS     = 201

private val plutoStateTable: Array<BodyState> = arrayOf(
    BodyState( -730000.0, TerseVector(-26.118207232108, -14.376168177825,   3.384402515299), TerseVector( 1.6339372163656e-03, -2.7861699588508e-03, -1.3585880229445e-03))
,   BodyState( -700800.0, TerseVector( 41.974905202127,  -0.448502952929, -12.770351505989), TerseVector( 7.3458569351457e-04,  2.2785014891658e-03,  4.8619778602049e-04))
,   BodyState( -671600.0, TerseVector( 14.706930780744,  44.269110540027,   9.353698474772), TerseVector(-2.1000147999800e-03,  2.2295915939915e-04,  7.0143443551414e-04))
,   BodyState( -642400.0, TerseVector(-29.441003929957,  -6.430161530570,   6.858481011305), TerseVector( 8.4495803960544e-04, -3.0783914758711e-03, -1.2106305981192e-03))
,   BodyState( -613200.0, TerseVector( 39.444396946234,  -6.557989760571, -13.913760296463), TerseVector( 1.1480029005873e-03,  2.2400006880665e-03,  3.5168075922288e-04))
,   BodyState( -584000.0, TerseVector( 20.230380950700,  43.266966657189,   7.382966091923), TerseVector(-1.9754081700585e-03,  5.3457141292226e-04,  7.5929169129793e-04))
,   BodyState( -554800.0, TerseVector(-30.658325364620,   2.093818874552,   9.880531138071), TerseVector( 6.1010603013347e-05, -3.1326500935382e-03, -9.9346125151067e-04))
,   BodyState( -525600.0, TerseVector( 35.737703251673, -12.587706024764, -14.677847247563), TerseVector( 1.5802939375649e-03,  2.1347678412429e-03,  1.9074436384343e-04))
,   BodyState( -496400.0, TerseVector( 25.466295188546,  41.367478338417,   5.216476873382), TerseVector(-1.8054401046468e-03,  8.3283083599510e-04,  8.0260156912107e-04))
,   BodyState( -467200.0, TerseVector(-29.847174904071,  10.636426313081,  12.297904180106), TerseVector(-6.3257063052907e-04, -2.9969577578221e-03, -7.4476074151596e-04))
,   BodyState( -438000.0, TerseVector( 30.774692107687, -18.236637015304, -14.945535879896), TerseVector( 2.0113162005465e-03,  1.9353827024189e-03, -2.0937793168297e-06))
,   BodyState( -408800.0, TerseVector( 30.243153324028,  38.656267888503,   2.938501750218), TerseVector(-1.6052508674468e-03,  1.1183495337525e-03,  8.3333973416824e-04))
,   BodyState( -379600.0, TerseVector(-27.288984772533,  18.643162147874,  14.023633623329), TerseVector(-1.1856388898191e-03, -2.7170609282181e-03, -4.9015526126399e-04))
,   BodyState( -350400.0, TerseVector( 24.519605196774, -23.245756064727, -14.626862367368), TerseVector( 2.4322321483154e-03,  1.6062008146048e-03, -2.3369181613312e-04))
,   BodyState( -321200.0, TerseVector( 34.505274805875,  35.125338586954,   0.557361475637), TerseVector(-1.3824391637782e-03,  1.3833397561817e-03,  8.4823598806262e-04))
,   BodyState( -292000.0, TerseVector(-23.275363915119,  25.818514298769,  15.055381588598), TerseVector(-1.6062295460975e-03, -2.3395961498533e-03, -2.4377362639479e-04))
,   BodyState( -262800.0, TerseVector( 17.050384798092, -27.180376290126, -13.608963321694), TerseVector( 2.8175521080578e-03,  1.1358749093955e-03, -4.9548725258825e-04))
,   BodyState( -233600.0, TerseVector( 38.093671910285,  30.880588383337,  -1.843688067413), TerseVector(-1.1317697153459e-03,  1.6128814698472e-03,  8.4177586176055e-04))
,   BodyState( -204400.0, TerseVector(-18.197852930878,  31.932869934309,  15.438294826279), TerseVector(-1.9117272501813e-03, -1.9146495909842e-03, -1.9657304369835e-05))
,   BodyState( -175200.0, TerseVector(  8.528924039997, -29.618422200048, -11.805400994258), TerseVector( 3.1034370787005e-03,  5.1393633292430e-04, -7.7293066202546e-04))
,   BodyState( -146000.0, TerseVector( 40.946857258640,  25.904973592021,  -4.256336240499), TerseVector(-8.3652705194051e-04,  1.8129497136404e-03,  8.1564228273060e-04))
,   BodyState( -116800.0, TerseVector(-12.326958895325,  36.881883446292,  15.217158258711), TerseVector(-2.1166103705038e-03, -1.4814420035990e-03,  1.7401209844705e-04))
,   BodyState(  -87600.0, TerseVector( -0.633258375909, -30.018759794709,  -9.171932874950), TerseVector( 3.2016994581737e-03, -2.5279858672148e-04, -1.0411088271861e-03))
,   BodyState(  -58400.0, TerseVector( 42.936048423883,  20.344685584452,  -6.588027007912), TerseVector(-5.0525450073192e-04,  1.9910074335507e-03,  7.7440196540269e-04))
,   BodyState(  -29200.0, TerseVector( -5.975910552974,  40.611809958460,  14.470131723673), TerseVector(-2.2184202156107e-03, -1.0562361130164e-03,  3.3652250216211e-04))
,   BodyState(       0.0, TerseVector( -9.875369580774, -27.978926224737,  -5.753711824704), TerseVector( 3.0287533248818e-03, -1.1276087003636e-03, -1.2651326732361e-03))
,   BodyState(   29200.0, TerseVector( 43.958831986165,  14.214147973292,  -8.808306227163), TerseVector(-1.4717608981871e-04,  2.1404187242141e-03,  7.1486567806614e-04))
,   BodyState(   58400.0, TerseVector(  0.678136763520,  43.094461639362,  13.243238780721), TerseVector(-2.2358226110718e-03, -6.3233636090933e-04,  4.7664798895648e-04))
,   BodyState(   87600.0, TerseVector(-18.282602096834, -23.305039586660,  -1.766620508028), TerseVector( 2.5567245263557e-03, -1.9902940754171e-03, -1.3943491701082e-03))
,   BodyState(  116800.0, TerseVector( 43.873338744526,   7.700705617215, -10.814273666425), TerseVector( 2.3174803055677e-04,  2.2402163127924e-03,  6.2988756452032e-04))
,   BodyState(  146000.0, TerseVector(  7.392949027906,  44.382678951534,  11.629500214854), TerseVector(-2.1932815453830e-03, -2.1751799585364e-04,  5.9556516201114e-04))
,   BodyState(  175200.0, TerseVector(-24.981690229261, -16.204012851426,   2.466457544298), TerseVector( 1.8193989149580e-03, -2.6765419531201e-03, -1.3848283502247e-03))
,   BodyState(  204400.0, TerseVector( 42.530187039511,   0.845935508021, -12.554907527683), TerseVector( 6.5059779150669e-04,  2.2725657282262e-03,  5.1133743202822e-04))
,   BodyState(  233600.0, TerseVector( 13.999526486822,  44.462363044894,   9.669418486465), TerseVector(-2.1079296569252e-03,  1.7533423831993e-04,  6.9128485798076e-04))
,   BodyState(  262800.0, TerseVector(-29.184024803031,  -7.371243995762,   6.493275957928), TerseVector( 9.3581363109681e-04, -3.0610357109184e-03, -1.2364201089345e-03))
,   BodyState(  292000.0, TerseVector( 39.831980671753,  -6.078405766765, -13.909815358656), TerseVector( 1.1117769689167e-03,  2.2362097830152e-03,  3.6230548231153e-04))
,   BodyState(  321200.0, TerseVector( 20.294955108476,  43.417190420251,   7.450091985932), TerseVector(-1.9742157451535e-03,  5.3102050468554e-04,  7.5938408813008e-04))
,   BodyState(  350400.0, TerseVector(-30.669992302160,   2.318743558955,   9.973480913858), TerseVector( 4.5605107450676e-05, -3.1308219926928e-03, -9.9066533301924e-04))
,   BodyState(  379600.0, TerseVector( 35.626122155983, -12.897647509224, -14.777586508444), TerseVector( 1.6015684949743e-03,  2.1171931182284e-03,  1.8002516202204e-04))
,   BodyState(  408800.0, TerseVector( 26.133186148561,  41.232139187599,   5.006401326220), TerseVector(-1.7857704419579e-03,  8.6046232702817e-04,  8.0614690298954e-04))
,   BodyState(  438000.0, TerseVector(-29.576740229230,  11.863535943587,  12.631323039872), TerseVector(-7.2292830060955e-04, -2.9587820140709e-03, -7.0824296450300e-04))
,   BodyState(  467200.0, TerseVector( 29.910805787391, -19.159019294000, -15.013363865194), TerseVector( 2.0871080437997e-03,  1.8848372554514e-03, -3.8528655083926e-05))
,   BodyState(  496400.0, TerseVector( 31.375957451819,  38.050372720763,   2.433138343754), TerseVector(-1.5546055556611e-03,  1.1699815465629e-03,  8.3565439266001e-04))
,   BodyState(  525600.0, TerseVector(-26.360071336928,  20.662505904952,  14.414696258958), TerseVector(-1.3142373118349e-03, -2.6236647854842e-03, -4.2542017598193e-04))
,   BodyState(  554800.0, TerseVector( 22.599441488648, -24.508879898306, -14.484045731468), TerseVector( 2.5454108304806e-03,  1.4917058755191e-03, -3.0243665086079e-04))
,   BodyState(  584000.0, TerseVector( 35.877864013014,  33.894226366071,  -0.224524636277), TerseVector(-1.2941245730845e-03,  1.4560427668319e-03,  8.4762160640137e-04))
,   BodyState(  613200.0, TerseVector(-21.538149762417,  28.204068269761,  15.321973799534), TerseVector(-1.7312117409010e-03, -2.1939631314577e-03, -1.6316913275180e-04))
,   BodyState(  642400.0, TerseVector( 13.971521374415, -28.339941764789, -13.083792871886), TerseVector( 2.9334630526035e-03,  9.1860931752944e-04, -5.9939422488627e-04))
,   BodyState(  671600.0, TerseVector( 39.526942044143,  28.939897360110,  -2.872799527539), TerseVector(-1.0068481658095e-03,  1.7021132888090e-03,  8.3578230511981e-04))
,   BodyState(  700800.0, TerseVector(-15.576200701394,  34.399412961275,  15.466033737854), TerseVector(-2.0098814612884e-03, -1.7191109825989e-03,  7.0414782780416e-05))
,   BodyState(  730000.0, TerseVector(  4.243252837090, -30.118201690825, -10.707441231349), TerseVector( 3.1725847067411e-03,  1.6098461202270e-04, -9.0672150593868e-04))
)
private val plutoCache = HashMap<Int, List<BodyGravCalc>>()

//---------------------------------------------------------------------------------------
// Models for Jupiter's four largest moons.

private val rotationJupEqj = RotationMatrix(
     9.99432765338654e-01, -3.36771074697641e-02,  0.00000000000000e+00,
     3.03959428906285e-02,  9.02057912352809e-01,  4.30543388542295e-01,
    -1.44994559663353e-02, -4.30299169409101e-01,  9.02569881273754e-01
)

private val jupiterMoonModel: Array<JupiterMoon> = arrayOf(
    // [0] Io
    JupiterMoon(
         2.8248942843381399e-07,
         1.4462132960212239e+00,
         3.5515522861824000e+00,
        arrayOf(  // a
            VsopTerm( 0.0028210960212903,  0.0000000000000000e+00,  0.0000000000000000e+00)
        ),
        arrayOf(  // l
            VsopTerm(-0.0001925258348666,  4.9369589722644998e+00,  1.3584836583050000e-02),
            VsopTerm(-0.0000970803596076,  4.3188796477322002e+00,  1.3034138432430000e-02),
            VsopTerm(-0.0000898817416500,  1.9080016428616999e+00,  3.0506486715799999e-03),
            VsopTerm(-0.0000553101050262,  1.4936156681568999e+00,  1.2938928911549999e-02)
        ),
        arrayOf(  // z
            VsopTerm( 0.0041510849668155,  4.0899396355450000e+00, -1.2906864146660001e-02),
            VsopTerm( 0.0006260521444113,  1.4461888986270000e+00,  3.5515522949801999e+00),
            VsopTerm( 0.0000352747346169,  2.1256287034577999e+00,  1.2727416566999999e-04)
        ),
        arrayOf(  // zeta
            VsopTerm( 0.0003142172466014,  2.7964219722923001e+00, -2.3150960980000000e-03),
            VsopTerm( 0.0000904169207946,  1.0477061879627001e+00, -5.6920638196000003e-04)
        )
    ),

    // [1] Europa
    JupiterMoon(
         2.8248327439289299e-07,
        -3.7352634374713622e-01,
         1.7693227111234699e+00,
        arrayOf(  // a
            VsopTerm( 0.0044871037804314,  0.0000000000000000e+00,  0.0000000000000000e+00),
            VsopTerm( 0.0000004324367498,  1.8196456062910000e+00,  1.7822295777568000e+00)
        ),
        arrayOf(  // l
            VsopTerm( 0.0008576433172936,  4.3188693178264002e+00,  1.3034138308049999e-02),
            VsopTerm( 0.0004549582875086,  1.4936531751079001e+00,  1.2938928819619999e-02),
            VsopTerm( 0.0003248939825174,  1.8196494533458001e+00,  1.7822295777568000e+00),
            VsopTerm(-0.0003074250079334,  4.9377037005910998e+00,  1.3584832867240000e-02),
            VsopTerm( 0.0001982386144784,  1.9079869054759999e+00,  3.0510121286900001e-03),
            VsopTerm( 0.0001834063551804,  2.1402853388529000e+00,  1.4500978933800000e-03),
            VsopTerm(-0.0001434383188452,  5.6222140366630002e+00,  8.9111478887838003e-01),
            VsopTerm(-0.0000771939140944,  4.3002724372349999e+00,  2.6733443704265998e+00)
        ),
        arrayOf(  // z
            VsopTerm(-0.0093589104136341,  4.0899396509038999e+00, -1.2906864146660001e-02),
            VsopTerm( 0.0002988994545555,  5.9097265185595003e+00,  1.7693227079461999e+00),
            VsopTerm( 0.0002139036390350,  2.1256289300016000e+00,  1.2727418406999999e-04),
            VsopTerm( 0.0001980963564781,  2.7435168292649998e+00,  6.7797343008999997e-04),
            VsopTerm( 0.0001210388158965,  5.5839943711203004e+00,  3.2056614899999997e-05),
            VsopTerm( 0.0000837042048393,  1.6094538368039000e+00, -9.0402165808846002e-01),
            VsopTerm( 0.0000823525166369,  1.4461887708689001e+00,  3.5515522949801999e+00)
        ),
        arrayOf(  // zeta
            VsopTerm( 0.0040404917832303,  1.0477063169425000e+00, -5.6920640539999997e-04),
            VsopTerm( 0.0002200421034564,  3.3368857864364001e+00, -1.2491307306999999e-04),
            VsopTerm( 0.0001662544744719,  2.4134862374710999e+00,  0.0000000000000000e+00),
            VsopTerm( 0.0000590282470983,  5.9719930968366004e+00, -3.0561602250000000e-05)
        )
    ),

    // [2] Ganymede
    JupiterMoon(
         2.8249818418472298e-07,
         2.8740893911433479e-01,
         8.7820792358932798e-01,
        arrayOf(  // a
            VsopTerm( 0.0071566594572575,  0.0000000000000000e+00,  0.0000000000000000e+00),
            VsopTerm( 0.0000013930299110,  1.1586745884981000e+00,  2.6733443704265998e+00)
        ),
        arrayOf(  // l
            VsopTerm( 0.0002310797886226,  2.1402987195941998e+00,  1.4500978438400001e-03),
            VsopTerm(-0.0001828635964118,  4.3188672736968003e+00,  1.3034138282630000e-02),
            VsopTerm( 0.0001512378778204,  4.9373102372298003e+00,  1.3584834812520000e-02),
            VsopTerm(-0.0001163720969778,  4.3002659861490002e+00,  2.6733443704265998e+00),
            VsopTerm(-0.0000955478069846,  1.4936612842567001e+00,  1.2938928798570001e-02),
            VsopTerm( 0.0000815246854464,  5.6222137132535002e+00,  8.9111478887838003e-01),
            VsopTerm(-0.0000801219679602,  1.2995922951532000e+00,  1.0034433456728999e+00),
            VsopTerm(-0.0000607017260182,  6.4978769669238001e-01,  5.0172167043264004e-01)
        ),
        arrayOf(  // z
            VsopTerm( 0.0014289811307319,  2.1256295942738999e+00,  1.2727413029000001e-04),
            VsopTerm( 0.0007710931226760,  5.5836330003496002e+00,  3.2064341100000001e-05),
            VsopTerm( 0.0005925911780766,  4.0899396636447998e+00, -1.2906864146660001e-02),
            VsopTerm( 0.0002045597496146,  5.2713683670371996e+00, -1.2523544076106000e-01),
            VsopTerm( 0.0001785118648258,  2.8743156721063001e-01,  8.7820792442520001e-01),
            VsopTerm( 0.0001131999784893,  1.4462127277818000e+00,  3.5515522949801999e+00),
            VsopTerm(-0.0000658778169210,  2.2702423990985001e+00, -1.7951364394536999e+00),
            VsopTerm( 0.0000497058888328,  5.9096792204858000e+00,  1.7693227129285001e+00)
        ),
        arrayOf(  // zeta
            VsopTerm( 0.0015932721570848,  3.3368862796665000e+00, -1.2491307058000000e-04),
            VsopTerm( 0.0008533093128905,  2.4133881688166001e+00,  0.0000000000000000e+00),
            VsopTerm( 0.0003513347911037,  5.9720789850126996e+00, -3.0561017709999999e-05),
            VsopTerm(-0.0001441929255483,  1.0477061764435001e+00, -5.6920632124000004e-04)
        )
    ),

    // [3] Callisto
    JupiterMoon(
         2.8249214488990899e-07,
        -3.6203412913757038e-01,
         3.7648623343382798e-01,
        arrayOf(  // a
            VsopTerm( 0.0125879701715314,  0.0000000000000000e+00,  0.0000000000000000e+00),
            VsopTerm( 0.0000035952049470,  6.4965776007116005e-01,  5.0172168165034003e-01),
            VsopTerm( 0.0000027580210652,  1.8084235781510001e+00,  3.1750660413359002e+00)
        ),
        arrayOf(  // l
            VsopTerm( 0.0005586040123824,  2.1404207189814999e+00,  1.4500979323100001e-03),
            VsopTerm(-0.0003805813868176,  2.7358844897852999e+00,  2.9729650620000000e-05),
            VsopTerm( 0.0002205152863262,  6.4979652596399995e-01,  5.0172167243580001e-01),
            VsopTerm( 0.0001877895151158,  1.8084787604004999e+00,  3.1750660413359002e+00),
            VsopTerm( 0.0000766916975242,  6.2720114319754998e+00,  1.3928364636651001e+00),
            VsopTerm( 0.0000747056855106,  1.2995916202344000e+00,  1.0034433456728999e+00)
        ),
        arrayOf(  // z
            VsopTerm( 0.0073755808467977,  5.5836071576083999e+00,  3.2065099140000001e-05),
            VsopTerm( 0.0002065924169942,  5.9209831565786004e+00,  3.7648624194703001e-01),
            VsopTerm( 0.0001589869764021,  2.8744006242622999e-01,  8.7820792442520001e-01),
            VsopTerm(-0.0001561131605348,  2.1257397865089001e+00,  1.2727441285000001e-04),
            VsopTerm( 0.0001486043380971,  1.4462134301023000e+00,  3.5515522949801999e+00),
            VsopTerm( 0.0000635073108731,  5.9096803285953996e+00,  1.7693227129285001e+00),
            VsopTerm( 0.0000599351698525,  4.1125517584797997e+00, -2.7985797954588998e+00),
            VsopTerm( 0.0000540660842731,  5.5390350845569003e+00,  2.8683408228299999e-03),
            VsopTerm(-0.0000489596900866,  4.6218149483337996e+00, -6.2695712529518999e-01)
        ),
        arrayOf(  // zeta
            VsopTerm( 0.0038422977898495,  2.4133922085556998e+00,  0.0000000000000000e+00),
            VsopTerm( 0.0022453891791894,  5.9721736773277003e+00, -3.0561255249999997e-05),
            VsopTerm(-0.0002604479450559,  3.3368746306408998e+00, -1.2491309972000001e-04),
            VsopTerm( 0.0000332112143230,  5.5604137742336999e+00,  2.9003768850700000e-03)
        )
    )
)

//---------------------------------------------------------------------------------------
// Constellation lookup table.

internal val constelNames: Array<ConstellationText> = arrayOf(
    ConstellationText("And", "Andromeda"           )    //  0
,   ConstellationText("Ant", "Antila"              )    //  1
,   ConstellationText("Aps", "Apus"                )    //  2
,   ConstellationText("Aql", "Aquila"              )    //  3
,   ConstellationText("Aqr", "Aquarius"            )    //  4
,   ConstellationText("Ara", "Ara"                 )    //  5
,   ConstellationText("Ari", "Aries"               )    //  6
,   ConstellationText("Aur", "Auriga"              )    //  7
,   ConstellationText("Boo", "Bootes"              )    //  8
,   ConstellationText("Cae", "Caelum"              )    //  9
,   ConstellationText("Cam", "Camelopardis"        )    // 10
,   ConstellationText("Cap", "Capricornus"         )    // 11
,   ConstellationText("Car", "Carina"              )    // 12
,   ConstellationText("Cas", "Cassiopeia"          )    // 13
,   ConstellationText("Cen", "Centaurus"           )    // 14
,   ConstellationText("Cep", "Cepheus"             )    // 15
,   ConstellationText("Cet", "Cetus"               )    // 16
,   ConstellationText("Cha", "Chamaeleon"          )    // 17
,   ConstellationText("Cir", "Circinus"            )    // 18
,   ConstellationText("CMa", "Canis Major"         )    // 19
,   ConstellationText("CMi", "Canis Minor"         )    // 20
,   ConstellationText("Cnc", "Cancer"              )    // 21
,   ConstellationText("Col", "Columba"             )    // 22
,   ConstellationText("Com", "Coma Berenices"      )    // 23
,   ConstellationText("CrA", "Corona Australis"    )    // 24
,   ConstellationText("CrB", "Corona Borealis"     )    // 25
,   ConstellationText("Crt", "Crater"              )    // 26
,   ConstellationText("Cru", "Crux"                )    // 27
,   ConstellationText("Crv", "Corvus"              )    // 28
,   ConstellationText("CVn", "Canes Venatici"      )    // 29
,   ConstellationText("Cyg", "Cygnus"              )    // 30
,   ConstellationText("Del", "Delphinus"           )    // 31
,   ConstellationText("Dor", "Dorado"              )    // 32
,   ConstellationText("Dra", "Draco"               )    // 33
,   ConstellationText("Equ", "Equuleus"            )    // 34
,   ConstellationText("Eri", "Eridanus"            )    // 35
,   ConstellationText("For", "Fornax"              )    // 36
,   ConstellationText("Gem", "Gemini"              )    // 37
,   ConstellationText("Gru", "Grus"                )    // 38
,   ConstellationText("Her", "Hercules"            )    // 39
,   ConstellationText("Hor", "Horologium"          )    // 40
,   ConstellationText("Hya", "Hydra"               )    // 41
,   ConstellationText("Hyi", "Hydrus"              )    // 42
,   ConstellationText("Ind", "Indus"               )    // 43
,   ConstellationText("Lac", "Lacerta"             )    // 44
,   ConstellationText("Leo", "Leo"                 )    // 45
,   ConstellationText("Lep", "Lepus"               )    // 46
,   ConstellationText("Lib", "Libra"               )    // 47
,   ConstellationText("LMi", "Leo Minor"           )    // 48
,   ConstellationText("Lup", "Lupus"               )    // 49
,   ConstellationText("Lyn", "Lynx"                )    // 50
,   ConstellationText("Lyr", "Lyra"                )    // 51
,   ConstellationText("Men", "Mensa"               )    // 52
,   ConstellationText("Mic", "Microscopium"        )    // 53
,   ConstellationText("Mon", "Monoceros"           )    // 54
,   ConstellationText("Mus", "Musca"               )    // 55
,   ConstellationText("Nor", "Norma"               )    // 56
,   ConstellationText("Oct", "Octans"              )    // 57
,   ConstellationText("Oph", "Ophiuchus"           )    // 58
,   ConstellationText("Ori", "Orion"               )    // 59
,   ConstellationText("Pav", "Pavo"                )    // 60
,   ConstellationText("Peg", "Pegasus"             )    // 61
,   ConstellationText("Per", "Perseus"             )    // 62
,   ConstellationText("Phe", "Phoenix"             )    // 63
,   ConstellationText("Pic", "Pictor"              )    // 64
,   ConstellationText("PsA", "Pisces Austrinus"    )    // 65
,   ConstellationText("Psc", "Pisces"              )    // 66
,   ConstellationText("Pup", "Puppis"              )    // 67
,   ConstellationText("Pyx", "Pyxis"               )    // 68
,   ConstellationText("Ret", "Reticulum"           )    // 69
,   ConstellationText("Scl", "Sculptor"            )    // 70
,   ConstellationText("Sco", "Scorpius"            )    // 71
,   ConstellationText("Sct", "Scutum"              )    // 72
,   ConstellationText("Ser", "Serpens"             )    // 73
,   ConstellationText("Sex", "Sextans"             )    // 74
,   ConstellationText("Sge", "Sagitta"             )    // 75
,   ConstellationText("Sgr", "Sagittarius"         )    // 76
,   ConstellationText("Tau", "Taurus"              )    // 77
,   ConstellationText("Tel", "Telescopium"         )    // 78
,   ConstellationText("TrA", "Triangulum Australe" )    // 79
,   ConstellationText("Tri", "Triangulum"          )    // 80
,   ConstellationText("Tuc", "Tucana"              )    // 81
,   ConstellationText("UMa", "Ursa Major"          )    // 82
,   ConstellationText("UMi", "Ursa Minor"          )    // 83
,   ConstellationText("Vel", "Vela"                )    // 84
,   ConstellationText("Vir", "Virgo"               )    // 85
,   ConstellationText("Vol", "Volans"              )    // 86
,   ConstellationText("Vul", "Vulpecula"           )    // 87
)

internal val constelBounds: Array<ConstellationBoundary> = arrayOf(
    ConstellationBoundary(83,    0.0, 8640.0,  2112.0)    // UMi
,   ConstellationBoundary(83, 2880.0, 5220.0,  2076.0)    // UMi
,   ConstellationBoundary(83, 7560.0, 8280.0,  2068.0)    // UMi
,   ConstellationBoundary(83, 6480.0, 7560.0,  2064.0)    // UMi
,   ConstellationBoundary(15,    0.0, 2880.0,  2040.0)    // Cep
,   ConstellationBoundary(10, 3300.0, 3840.0,  1968.0)    // Cam
,   ConstellationBoundary(15,    0.0, 1800.0,  1920.0)    // Cep
,   ConstellationBoundary(10, 3840.0, 5220.0,  1920.0)    // Cam
,   ConstellationBoundary(83, 6300.0, 6480.0,  1920.0)    // UMi
,   ConstellationBoundary(33, 7260.0, 7560.0,  1920.0)    // Dra
,   ConstellationBoundary(15,    0.0, 1263.0,  1848.0)    // Cep
,   ConstellationBoundary(10, 4140.0, 4890.0,  1848.0)    // Cam
,   ConstellationBoundary(83, 5952.0, 6300.0,  1800.0)    // UMi
,   ConstellationBoundary(15, 7260.0, 7440.0,  1800.0)    // Cep
,   ConstellationBoundary(10, 2868.0, 3300.0,  1764.0)    // Cam
,   ConstellationBoundary(33, 3300.0, 4080.0,  1764.0)    // Dra
,   ConstellationBoundary(83, 4680.0, 5952.0,  1680.0)    // UMi
,   ConstellationBoundary(13, 1116.0, 1230.0,  1632.0)    // Cas
,   ConstellationBoundary(33, 7350.0, 7440.0,  1608.0)    // Dra
,   ConstellationBoundary(33, 4080.0, 4320.0,  1596.0)    // Dra
,   ConstellationBoundary(15,    0.0,  120.0,  1584.0)    // Cep
,   ConstellationBoundary(83, 5040.0, 5640.0,  1584.0)    // UMi
,   ConstellationBoundary(15, 8490.0, 8640.0,  1584.0)    // Cep
,   ConstellationBoundary(33, 4320.0, 4860.0,  1536.0)    // Dra
,   ConstellationBoundary(33, 4860.0, 5190.0,  1512.0)    // Dra
,   ConstellationBoundary(15, 8340.0, 8490.0,  1512.0)    // Cep
,   ConstellationBoundary(10, 2196.0, 2520.0,  1488.0)    // Cam
,   ConstellationBoundary(33, 7200.0, 7350.0,  1476.0)    // Dra
,   ConstellationBoundary(15, 7393.2, 7416.0,  1462.0)    // Cep
,   ConstellationBoundary(10, 2520.0, 2868.0,  1440.0)    // Cam
,   ConstellationBoundary(82, 2868.0, 3030.0,  1440.0)    // UMa
,   ConstellationBoundary(33, 7116.0, 7200.0,  1428.0)    // Dra
,   ConstellationBoundary(15, 7200.0, 7393.2,  1428.0)    // Cep
,   ConstellationBoundary(15, 8232.0, 8340.0,  1418.0)    // Cep
,   ConstellationBoundary(13,    0.0,  876.0,  1404.0)    // Cas
,   ConstellationBoundary(33, 6990.0, 7116.0,  1392.0)    // Dra
,   ConstellationBoundary(13,  612.0,  687.0,  1380.0)    // Cas
,   ConstellationBoundary(13,  876.0, 1116.0,  1368.0)    // Cas
,   ConstellationBoundary(10, 1116.0, 1140.0,  1368.0)    // Cam
,   ConstellationBoundary(15, 8034.0, 8232.0,  1350.0)    // Cep
,   ConstellationBoundary(10, 1800.0, 2196.0,  1344.0)    // Cam
,   ConstellationBoundary(82, 5052.0, 5190.0,  1332.0)    // UMa
,   ConstellationBoundary(33, 5190.0, 6990.0,  1332.0)    // Dra
,   ConstellationBoundary(10, 1140.0, 1200.0,  1320.0)    // Cam
,   ConstellationBoundary(15, 7968.0, 8034.0,  1320.0)    // Cep
,   ConstellationBoundary(15, 7416.0, 7908.0,  1316.0)    // Cep
,   ConstellationBoundary(13,    0.0,  612.0,  1296.0)    // Cas
,   ConstellationBoundary(50, 2196.0, 2340.0,  1296.0)    // Lyn
,   ConstellationBoundary(82, 4350.0, 4860.0,  1272.0)    // UMa
,   ConstellationBoundary(33, 5490.0, 5670.0,  1272.0)    // Dra
,   ConstellationBoundary(15, 7908.0, 7968.0,  1266.0)    // Cep
,   ConstellationBoundary(10, 1200.0, 1800.0,  1260.0)    // Cam
,   ConstellationBoundary(13, 8232.0, 8400.0,  1260.0)    // Cas
,   ConstellationBoundary(33, 5670.0, 6120.0,  1236.0)    // Dra
,   ConstellationBoundary(62,  735.0,  906.0,  1212.0)    // Per
,   ConstellationBoundary(33, 6120.0, 6564.0,  1212.0)    // Dra
,   ConstellationBoundary(13,    0.0,  492.0,  1200.0)    // Cas
,   ConstellationBoundary(62,  492.0,  600.0,  1200.0)    // Per
,   ConstellationBoundary(50, 2340.0, 2448.0,  1200.0)    // Lyn
,   ConstellationBoundary(13, 8400.0, 8640.0,  1200.0)    // Cas
,   ConstellationBoundary(82, 4860.0, 5052.0,  1164.0)    // UMa
,   ConstellationBoundary(13,    0.0,  402.0,  1152.0)    // Cas
,   ConstellationBoundary(13, 8490.0, 8640.0,  1152.0)    // Cas
,   ConstellationBoundary(39, 6543.0, 6564.0,  1140.0)    // Her
,   ConstellationBoundary(33, 6564.0, 6870.0,  1140.0)    // Dra
,   ConstellationBoundary(30, 6870.0, 6900.0,  1140.0)    // Cyg
,   ConstellationBoundary(62,  600.0,  735.0,  1128.0)    // Per
,   ConstellationBoundary(82, 3030.0, 3300.0,  1128.0)    // UMa
,   ConstellationBoundary(13,   60.0,  312.0,  1104.0)    // Cas
,   ConstellationBoundary(82, 4320.0, 4350.0,  1080.0)    // UMa
,   ConstellationBoundary(50, 2448.0, 2652.0,  1068.0)    // Lyn
,   ConstellationBoundary(30, 7887.0, 7908.0,  1056.0)    // Cyg
,   ConstellationBoundary(30, 7875.0, 7887.0,  1050.0)    // Cyg
,   ConstellationBoundary(30, 6900.0, 6984.0,  1044.0)    // Cyg
,   ConstellationBoundary(82, 3300.0, 3660.0,  1008.0)    // UMa
,   ConstellationBoundary(82, 3660.0, 3882.0,   960.0)    // UMa
,   ConstellationBoundary( 8, 5556.0, 5670.0,   960.0)    // Boo
,   ConstellationBoundary(39, 5670.0, 5880.0,   960.0)    // Her
,   ConstellationBoundary(50, 3330.0, 3450.0,   954.0)    // Lyn
,   ConstellationBoundary( 0,    0.0,  906.0,   882.0)    // And
,   ConstellationBoundary(62,  906.0,  924.0,   882.0)    // Per
,   ConstellationBoundary(51, 6969.0, 6984.0,   876.0)    // Lyr
,   ConstellationBoundary(62, 1620.0, 1689.0,   864.0)    // Per
,   ConstellationBoundary(30, 7824.0, 7875.0,   864.0)    // Cyg
,   ConstellationBoundary(44, 7875.0, 7920.0,   864.0)    // Lac
,   ConstellationBoundary( 7, 2352.0, 2652.0,   852.0)    // Aur
,   ConstellationBoundary(50, 2652.0, 2790.0,   852.0)    // Lyn
,   ConstellationBoundary( 0,    0.0,  720.0,   840.0)    // And
,   ConstellationBoundary(44, 7920.0, 8214.0,   840.0)    // Lac
,   ConstellationBoundary(44, 8214.0, 8232.0,   828.0)    // Lac
,   ConstellationBoundary( 0, 8232.0, 8460.0,   828.0)    // And
,   ConstellationBoundary(62,  924.0,  978.0,   816.0)    // Per
,   ConstellationBoundary(82, 3882.0, 3960.0,   816.0)    // UMa
,   ConstellationBoundary(29, 4320.0, 4440.0,   816.0)    // CVn
,   ConstellationBoundary(50, 2790.0, 3330.0,   804.0)    // Lyn
,   ConstellationBoundary(48, 3330.0, 3558.0,   804.0)    // LMi
,   ConstellationBoundary( 0,  258.0,  507.0,   792.0)    // And
,   ConstellationBoundary( 8, 5466.0, 5556.0,   792.0)    // Boo
,   ConstellationBoundary( 0, 8460.0, 8550.0,   770.0)    // And
,   ConstellationBoundary(29, 4440.0, 4770.0,   768.0)    // CVn
,   ConstellationBoundary( 0, 8550.0, 8640.0,   752.0)    // And
,   ConstellationBoundary(29, 5025.0, 5052.0,   738.0)    // CVn
,   ConstellationBoundary(80,  870.0,  978.0,   736.0)    // Tri
,   ConstellationBoundary(62,  978.0, 1620.0,   736.0)    // Per
,   ConstellationBoundary( 7, 1620.0, 1710.0,   720.0)    // Aur
,   ConstellationBoundary(51, 6543.0, 6969.0,   720.0)    // Lyr
,   ConstellationBoundary(82, 3960.0, 4320.0,   696.0)    // UMa
,   ConstellationBoundary(30, 7080.0, 7530.0,   696.0)    // Cyg
,   ConstellationBoundary( 7, 1710.0, 2118.0,   684.0)    // Aur
,   ConstellationBoundary(48, 3558.0, 3780.0,   684.0)    // LMi
,   ConstellationBoundary(29, 4770.0, 5025.0,   684.0)    // CVn
,   ConstellationBoundary( 0,    0.0,   24.0,   672.0)    // And
,   ConstellationBoundary(80,  507.0,  600.0,   672.0)    // Tri
,   ConstellationBoundary( 7, 2118.0, 2352.0,   672.0)    // Aur
,   ConstellationBoundary(37, 2838.0, 2880.0,   672.0)    // Gem
,   ConstellationBoundary(30, 7530.0, 7824.0,   672.0)    // Cyg
,   ConstellationBoundary(30, 6933.0, 7080.0,   660.0)    // Cyg
,   ConstellationBoundary(80,  690.0,  870.0,   654.0)    // Tri
,   ConstellationBoundary(25, 5820.0, 5880.0,   648.0)    // CrB
,   ConstellationBoundary( 8, 5430.0, 5466.0,   624.0)    // Boo
,   ConstellationBoundary(25, 5466.0, 5820.0,   624.0)    // CrB
,   ConstellationBoundary(51, 6612.0, 6792.0,   624.0)    // Lyr
,   ConstellationBoundary(48, 3870.0, 3960.0,   612.0)    // LMi
,   ConstellationBoundary(51, 6792.0, 6933.0,   612.0)    // Lyr
,   ConstellationBoundary(80,  600.0,  690.0,   600.0)    // Tri
,   ConstellationBoundary(66,  258.0,  306.0,   570.0)    // Psc
,   ConstellationBoundary(48, 3780.0, 3870.0,   564.0)    // LMi
,   ConstellationBoundary(87, 7650.0, 7710.0,   564.0)    // Vul
,   ConstellationBoundary(77, 2052.0, 2118.0,   548.0)    // Tau
,   ConstellationBoundary( 0,   24.0,   51.0,   528.0)    // And
,   ConstellationBoundary(73, 5730.0, 5772.0,   528.0)    // Ser
,   ConstellationBoundary(37, 2118.0, 2238.0,   516.0)    // Gem
,   ConstellationBoundary(87, 7140.0, 7290.0,   510.0)    // Vul
,   ConstellationBoundary(87, 6792.0, 6930.0,   506.0)    // Vul
,   ConstellationBoundary( 0,   51.0,  306.0,   504.0)    // And
,   ConstellationBoundary(87, 7290.0, 7404.0,   492.0)    // Vul
,   ConstellationBoundary(37, 2811.0, 2838.0,   480.0)    // Gem
,   ConstellationBoundary(87, 7404.0, 7650.0,   468.0)    // Vul
,   ConstellationBoundary(87, 6930.0, 7140.0,   460.0)    // Vul
,   ConstellationBoundary( 6, 1182.0, 1212.0,   456.0)    // Ari
,   ConstellationBoundary(75, 6792.0, 6840.0,   444.0)    // Sge
,   ConstellationBoundary(59, 2052.0, 2076.0,   432.0)    // Ori
,   ConstellationBoundary(37, 2238.0, 2271.0,   420.0)    // Gem
,   ConstellationBoundary(75, 6840.0, 7140.0,   388.0)    // Sge
,   ConstellationBoundary(77, 1788.0, 1920.0,   384.0)    // Tau
,   ConstellationBoundary(39, 5730.0, 5790.0,   384.0)    // Her
,   ConstellationBoundary(75, 7140.0, 7290.0,   378.0)    // Sge
,   ConstellationBoundary(77, 1662.0, 1788.0,   372.0)    // Tau
,   ConstellationBoundary(77, 1920.0, 2016.0,   372.0)    // Tau
,   ConstellationBoundary(23, 4620.0, 4860.0,   360.0)    // Com
,   ConstellationBoundary(39, 6210.0, 6570.0,   344.0)    // Her
,   ConstellationBoundary(23, 4272.0, 4620.0,   336.0)    // Com
,   ConstellationBoundary(37, 2700.0, 2811.0,   324.0)    // Gem
,   ConstellationBoundary(39, 6030.0, 6210.0,   308.0)    // Her
,   ConstellationBoundary(61,    0.0,   51.0,   300.0)    // Peg
,   ConstellationBoundary(77, 2016.0, 2076.0,   300.0)    // Tau
,   ConstellationBoundary(37, 2520.0, 2700.0,   300.0)    // Gem
,   ConstellationBoundary(61, 7602.0, 7680.0,   300.0)    // Peg
,   ConstellationBoundary(37, 2271.0, 2496.0,   288.0)    // Gem
,   ConstellationBoundary(39, 6570.0, 6792.0,   288.0)    // Her
,   ConstellationBoundary(31, 7515.0, 7578.0,   284.0)    // Del
,   ConstellationBoundary(61, 7578.0, 7602.0,   284.0)    // Peg
,   ConstellationBoundary(45, 4146.0, 4272.0,   264.0)    // Leo
,   ConstellationBoundary(59, 2247.0, 2271.0,   240.0)    // Ori
,   ConstellationBoundary(37, 2496.0, 2520.0,   240.0)    // Gem
,   ConstellationBoundary(21, 2811.0, 2853.0,   240.0)    // Cnc
,   ConstellationBoundary(61, 8580.0, 8640.0,   240.0)    // Peg
,   ConstellationBoundary( 6,  600.0, 1182.0,   238.0)    // Ari
,   ConstellationBoundary(31, 7251.0, 7308.0,   204.0)    // Del
,   ConstellationBoundary( 8, 4860.0, 5430.0,   192.0)    // Boo
,   ConstellationBoundary(61, 8190.0, 8580.0,   180.0)    // Peg
,   ConstellationBoundary(21, 2853.0, 3330.0,   168.0)    // Cnc
,   ConstellationBoundary(45, 3330.0, 3870.0,   168.0)    // Leo
,   ConstellationBoundary(58, 6570.0, 6718.4,   150.0)    // Oph
,   ConstellationBoundary( 3, 6718.4, 6792.0,   150.0)    // Aql
,   ConstellationBoundary(31, 7500.0, 7515.0,   144.0)    // Del
,   ConstellationBoundary(20, 2520.0, 2526.0,   132.0)    // CMi
,   ConstellationBoundary(73, 6570.0, 6633.0,   108.0)    // Ser
,   ConstellationBoundary(39, 5790.0, 6030.0,    96.0)    // Her
,   ConstellationBoundary(58, 6570.0, 6633.0,    72.0)    // Oph
,   ConstellationBoundary(61, 7728.0, 7800.0,    66.0)    // Peg
,   ConstellationBoundary(66,    0.0,  720.0,    48.0)    // Psc
,   ConstellationBoundary(73, 6690.0, 6792.0,    48.0)    // Ser
,   ConstellationBoundary(31, 7308.0, 7500.0,    48.0)    // Del
,   ConstellationBoundary(34, 7500.0, 7680.0,    48.0)    // Equ
,   ConstellationBoundary(61, 7680.0, 7728.0,    48.0)    // Peg
,   ConstellationBoundary(61, 7920.0, 8190.0,    48.0)    // Peg
,   ConstellationBoundary(61, 7800.0, 7920.0,    42.0)    // Peg
,   ConstellationBoundary(20, 2526.0, 2592.0,    36.0)    // CMi
,   ConstellationBoundary(77, 1290.0, 1662.0,     0.0)    // Tau
,   ConstellationBoundary(59, 1662.0, 1680.0,     0.0)    // Ori
,   ConstellationBoundary(20, 2592.0, 2910.0,     0.0)    // CMi
,   ConstellationBoundary(85, 5280.0, 5430.0,     0.0)    // Vir
,   ConstellationBoundary(58, 6420.0, 6570.0,     0.0)    // Oph
,   ConstellationBoundary(16,  954.0, 1182.0,   -42.0)    // Cet
,   ConstellationBoundary(77, 1182.0, 1290.0,   -42.0)    // Tau
,   ConstellationBoundary(73, 5430.0, 5856.0,   -78.0)    // Ser
,   ConstellationBoundary(59, 1680.0, 1830.0,   -96.0)    // Ori
,   ConstellationBoundary(59, 2100.0, 2247.0,   -96.0)    // Ori
,   ConstellationBoundary(73, 6420.0, 6468.0,   -96.0)    // Ser
,   ConstellationBoundary(73, 6570.0, 6690.0,   -96.0)    // Ser
,   ConstellationBoundary( 3, 6690.0, 6792.0,   -96.0)    // Aql
,   ConstellationBoundary(66, 8190.0, 8580.0,   -96.0)    // Psc
,   ConstellationBoundary(45, 3870.0, 4146.0,  -144.0)    // Leo
,   ConstellationBoundary(85, 4146.0, 4260.0,  -144.0)    // Vir
,   ConstellationBoundary(66,    0.0,  120.0,  -168.0)    // Psc
,   ConstellationBoundary(66, 8580.0, 8640.0,  -168.0)    // Psc
,   ConstellationBoundary(85, 5130.0, 5280.0,  -192.0)    // Vir
,   ConstellationBoundary(58, 5730.0, 5856.0,  -192.0)    // Oph
,   ConstellationBoundary( 3, 7200.0, 7392.0,  -216.0)    // Aql
,   ConstellationBoundary( 4, 7680.0, 7872.0,  -216.0)    // Aqr
,   ConstellationBoundary(58, 6180.0, 6468.0,  -240.0)    // Oph
,   ConstellationBoundary(54, 2100.0, 2910.0,  -264.0)    // Mon
,   ConstellationBoundary(35, 1770.0, 1830.0,  -264.0)    // Eri
,   ConstellationBoundary(59, 1830.0, 2100.0,  -264.0)    // Ori
,   ConstellationBoundary(41, 2910.0, 3012.0,  -264.0)    // Hya
,   ConstellationBoundary(74, 3450.0, 3870.0,  -264.0)    // Sex
,   ConstellationBoundary(85, 4260.0, 4620.0,  -264.0)    // Vir
,   ConstellationBoundary(58, 6330.0, 6360.0,  -280.0)    // Oph
,   ConstellationBoundary( 3, 6792.0, 7200.0,  -288.8)    // Aql
,   ConstellationBoundary(35, 1740.0, 1770.0,  -348.0)    // Eri
,   ConstellationBoundary( 4, 7392.0, 7680.0,  -360.0)    // Aqr
,   ConstellationBoundary(73, 6180.0, 6570.0,  -384.0)    // Ser
,   ConstellationBoundary(72, 6570.0, 6792.0,  -384.0)    // Sct
,   ConstellationBoundary(41, 3012.0, 3090.0,  -408.0)    // Hya
,   ConstellationBoundary(58, 5856.0, 5895.0,  -438.0)    // Oph
,   ConstellationBoundary(41, 3090.0, 3270.0,  -456.0)    // Hya
,   ConstellationBoundary(26, 3870.0, 3900.0,  -456.0)    // Crt
,   ConstellationBoundary(71, 5856.0, 5895.0,  -462.0)    // Sco
,   ConstellationBoundary(47, 5640.0, 5730.0,  -480.0)    // Lib
,   ConstellationBoundary(28, 4530.0, 4620.0,  -528.0)    // Crv
,   ConstellationBoundary(85, 4620.0, 5130.0,  -528.0)    // Vir
,   ConstellationBoundary(41, 3270.0, 3510.0,  -576.0)    // Hya
,   ConstellationBoundary(16,  600.0,  954.0,  -585.2)    // Cet
,   ConstellationBoundary(35,  954.0, 1350.0,  -585.2)    // Eri
,   ConstellationBoundary(26, 3900.0, 4260.0,  -588.0)    // Crt
,   ConstellationBoundary(28, 4260.0, 4530.0,  -588.0)    // Crv
,   ConstellationBoundary(47, 5130.0, 5370.0,  -588.0)    // Lib
,   ConstellationBoundary(58, 5856.0, 6030.0,  -590.0)    // Oph
,   ConstellationBoundary(16,    0.0,  600.0,  -612.0)    // Cet
,   ConstellationBoundary(11, 7680.0, 7872.0,  -612.0)    // Cap
,   ConstellationBoundary( 4, 7872.0, 8580.0,  -612.0)    // Aqr
,   ConstellationBoundary(16, 8580.0, 8640.0,  -612.0)    // Cet
,   ConstellationBoundary(41, 3510.0, 3690.0,  -636.0)    // Hya
,   ConstellationBoundary(35, 1692.0, 1740.0,  -654.0)    // Eri
,   ConstellationBoundary(46, 1740.0, 2202.0,  -654.0)    // Lep
,   ConstellationBoundary(11, 7200.0, 7680.0,  -672.0)    // Cap
,   ConstellationBoundary(41, 3690.0, 3810.0,  -700.0)    // Hya
,   ConstellationBoundary(41, 4530.0, 5370.0,  -708.0)    // Hya
,   ConstellationBoundary(47, 5370.0, 5640.0,  -708.0)    // Lib
,   ConstellationBoundary(71, 5640.0, 5760.0,  -708.0)    // Sco
,   ConstellationBoundary(35, 1650.0, 1692.0,  -720.0)    // Eri
,   ConstellationBoundary(58, 6030.0, 6336.0,  -720.0)    // Oph
,   ConstellationBoundary(76, 6336.0, 6420.0,  -720.0)    // Sgr
,   ConstellationBoundary(41, 3810.0, 3900.0,  -748.0)    // Hya
,   ConstellationBoundary(19, 2202.0, 2652.0,  -792.0)    // CMa
,   ConstellationBoundary(41, 4410.0, 4530.0,  -792.0)    // Hya
,   ConstellationBoundary(41, 3900.0, 4410.0,  -840.0)    // Hya
,   ConstellationBoundary(36, 1260.0, 1350.0,  -864.0)    // For
,   ConstellationBoundary(68, 3012.0, 3372.0,  -882.0)    // Pyx
,   ConstellationBoundary(35, 1536.0, 1650.0,  -888.0)    // Eri
,   ConstellationBoundary(76, 6420.0, 6900.0,  -888.0)    // Sgr
,   ConstellationBoundary(65, 7680.0, 8280.0,  -888.0)    // PsA
,   ConstellationBoundary(70, 8280.0, 8400.0,  -888.0)    // Scl
,   ConstellationBoundary(36, 1080.0, 1260.0,  -950.0)    // For
,   ConstellationBoundary( 1, 3372.0, 3960.0,  -954.0)    // Ant
,   ConstellationBoundary(70,    0.0,  600.0,  -960.0)    // Scl
,   ConstellationBoundary(36,  600.0, 1080.0,  -960.0)    // For
,   ConstellationBoundary(35, 1392.0, 1536.0,  -960.0)    // Eri
,   ConstellationBoundary(70, 8400.0, 8640.0,  -960.0)    // Scl
,   ConstellationBoundary(14, 5100.0, 5370.0, -1008.0)    // Cen
,   ConstellationBoundary(49, 5640.0, 5760.0, -1008.0)    // Lup
,   ConstellationBoundary(71, 5760.0, 5911.5, -1008.0)    // Sco
,   ConstellationBoundary( 9, 1740.0, 1800.0, -1032.0)    // Cae
,   ConstellationBoundary(22, 1800.0, 2370.0, -1032.0)    // Col
,   ConstellationBoundary(67, 2880.0, 3012.0, -1032.0)    // Pup
,   ConstellationBoundary(35, 1230.0, 1392.0, -1056.0)    // Eri
,   ConstellationBoundary(71, 5911.5, 6420.0, -1092.0)    // Sco
,   ConstellationBoundary(24, 6420.0, 6900.0, -1092.0)    // CrA
,   ConstellationBoundary(76, 6900.0, 7320.0, -1092.0)    // Sgr
,   ConstellationBoundary(53, 7320.0, 7680.0, -1092.0)    // Mic
,   ConstellationBoundary(35, 1080.0, 1230.0, -1104.0)    // Eri
,   ConstellationBoundary( 9, 1620.0, 1740.0, -1116.0)    // Cae
,   ConstellationBoundary(49, 5520.0, 5640.0, -1152.0)    // Lup
,   ConstellationBoundary(63,    0.0,  840.0, -1156.0)    // Phe
,   ConstellationBoundary(35,  960.0, 1080.0, -1176.0)    // Eri
,   ConstellationBoundary(40, 1470.0, 1536.0, -1176.0)    // Hor
,   ConstellationBoundary( 9, 1536.0, 1620.0, -1176.0)    // Cae
,   ConstellationBoundary(38, 7680.0, 7920.0, -1200.0)    // Gru
,   ConstellationBoundary(67, 2160.0, 2880.0, -1218.0)    // Pup
,   ConstellationBoundary(84, 2880.0, 2940.0, -1218.0)    // Vel
,   ConstellationBoundary(35,  870.0,  960.0, -1224.0)    // Eri
,   ConstellationBoundary(40, 1380.0, 1470.0, -1224.0)    // Hor
,   ConstellationBoundary(63,    0.0,  660.0, -1236.0)    // Phe
,   ConstellationBoundary(12, 2160.0, 2220.0, -1260.0)    // Car
,   ConstellationBoundary(84, 2940.0, 3042.0, -1272.0)    // Vel
,   ConstellationBoundary(40, 1260.0, 1380.0, -1276.0)    // Hor
,   ConstellationBoundary(32, 1380.0, 1440.0, -1276.0)    // Dor
,   ConstellationBoundary(63,    0.0,  570.0, -1284.0)    // Phe
,   ConstellationBoundary(35,  780.0,  870.0, -1296.0)    // Eri
,   ConstellationBoundary(64, 1620.0, 1800.0, -1296.0)    // Pic
,   ConstellationBoundary(49, 5418.0, 5520.0, -1296.0)    // Lup
,   ConstellationBoundary(84, 3042.0, 3180.0, -1308.0)    // Vel
,   ConstellationBoundary(12, 2220.0, 2340.0, -1320.0)    // Car
,   ConstellationBoundary(14, 4260.0, 4620.0, -1320.0)    // Cen
,   ConstellationBoundary(49, 5100.0, 5418.0, -1320.0)    // Lup
,   ConstellationBoundary(56, 5418.0, 5520.0, -1320.0)    // Nor
,   ConstellationBoundary(32, 1440.0, 1560.0, -1356.0)    // Dor
,   ConstellationBoundary(84, 3180.0, 3960.0, -1356.0)    // Vel
,   ConstellationBoundary(14, 3960.0, 4050.0, -1356.0)    // Cen
,   ConstellationBoundary( 5, 6300.0, 6480.0, -1368.0)    // Ara
,   ConstellationBoundary(78, 6480.0, 7320.0, -1368.0)    // Tel
,   ConstellationBoundary(38, 7920.0, 8400.0, -1368.0)    // Gru
,   ConstellationBoundary(40, 1152.0, 1260.0, -1380.0)    // Hor
,   ConstellationBoundary(64, 1800.0, 1980.0, -1380.0)    // Pic
,   ConstellationBoundary(12, 2340.0, 2460.0, -1392.0)    // Car
,   ConstellationBoundary(63,    0.0,  480.0, -1404.0)    // Phe
,   ConstellationBoundary(35,  480.0,  780.0, -1404.0)    // Eri
,   ConstellationBoundary(63, 8400.0, 8640.0, -1404.0)    // Phe
,   ConstellationBoundary(32, 1560.0, 1650.0, -1416.0)    // Dor
,   ConstellationBoundary(56, 5520.0, 5911.5, -1440.0)    // Nor
,   ConstellationBoundary(43, 7320.0, 7680.0, -1440.0)    // Ind
,   ConstellationBoundary(64, 1980.0, 2160.0, -1464.0)    // Pic
,   ConstellationBoundary(18, 5460.0, 5520.0, -1464.0)    // Cir
,   ConstellationBoundary( 5, 5911.5, 5970.0, -1464.0)    // Ara
,   ConstellationBoundary(18, 5370.0, 5460.0, -1526.0)    // Cir
,   ConstellationBoundary( 5, 5970.0, 6030.0, -1526.0)    // Ara
,   ConstellationBoundary(64, 2160.0, 2460.0, -1536.0)    // Pic
,   ConstellationBoundary(12, 2460.0, 3252.0, -1536.0)    // Car
,   ConstellationBoundary(14, 4050.0, 4260.0, -1536.0)    // Cen
,   ConstellationBoundary(27, 4260.0, 4620.0, -1536.0)    // Cru
,   ConstellationBoundary(14, 4620.0, 5232.0, -1536.0)    // Cen
,   ConstellationBoundary(18, 4860.0, 4920.0, -1560.0)    // Cir
,   ConstellationBoundary( 5, 6030.0, 6060.0, -1560.0)    // Ara
,   ConstellationBoundary(40,  780.0, 1152.0, -1620.0)    // Hor
,   ConstellationBoundary(69, 1152.0, 1650.0, -1620.0)    // Ret
,   ConstellationBoundary(18, 5310.0, 5370.0, -1620.0)    // Cir
,   ConstellationBoundary( 5, 6060.0, 6300.0, -1620.0)    // Ara
,   ConstellationBoundary(60, 6300.0, 6480.0, -1620.0)    // Pav
,   ConstellationBoundary(81, 7920.0, 8400.0, -1620.0)    // Tuc
,   ConstellationBoundary(32, 1650.0, 2370.0, -1680.0)    // Dor
,   ConstellationBoundary(18, 4920.0, 5310.0, -1680.0)    // Cir
,   ConstellationBoundary(79, 5310.0, 6120.0, -1680.0)    // TrA
,   ConstellationBoundary(81,    0.0,  480.0, -1800.0)    // Tuc
,   ConstellationBoundary(42, 1260.0, 1650.0, -1800.0)    // Hyi
,   ConstellationBoundary(86, 2370.0, 3252.0, -1800.0)    // Vol
,   ConstellationBoundary(12, 3252.0, 4050.0, -1800.0)    // Car
,   ConstellationBoundary(55, 4050.0, 4920.0, -1800.0)    // Mus
,   ConstellationBoundary(60, 6480.0, 7680.0, -1800.0)    // Pav
,   ConstellationBoundary(43, 7680.0, 8400.0, -1800.0)    // Ind
,   ConstellationBoundary(81, 8400.0, 8640.0, -1800.0)    // Tuc
,   ConstellationBoundary(81,  270.0,  480.0, -1824.0)    // Tuc
,   ConstellationBoundary(42,    0.0, 1260.0, -1980.0)    // Hyi
,   ConstellationBoundary(17, 2760.0, 4920.0, -1980.0)    // Cha
,   ConstellationBoundary( 2, 4920.0, 6480.0, -1980.0)    // Aps
,   ConstellationBoundary(52, 1260.0, 2760.0, -2040.0)    // Men
,   ConstellationBoundary(57,    0.0, 8640.0, -2160.0)    // Oct
)

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
