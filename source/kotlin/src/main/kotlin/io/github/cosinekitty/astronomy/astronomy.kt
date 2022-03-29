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
import kotlin.math.acos
import kotlin.math.atan2
import kotlin.math.cos
import kotlin.math.floor
import kotlin.math.hypot
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
private const val ASEC180 = 180.0 * 60.0 * 60.0;         // arcseconds per 180 degrees (or pi radians)
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
enum class Body(internal val massProduct: Double?) {
    /**
     * A placeholder value representing an invalid or unknown celestial body.
     */
    Invalid(massProduct = null),

    /**
     * The planet Mercury.
     */
    Mercury(massProduct = MERCURY_GM),

    /**
     * The planet Venus.
     */
    Venus(massProduct = VENUS_GM),

    /**
     * The planet Earth.
     * Some functions that accept a `Body` parameter will fail if passed this value
     * because they assume that an observation is being made from the Earth,
     * and therefore the Earth is not a target of observation.
     */
    Earth(massProduct = EARTH_GM),

    /**
     * The planet Mars.
     */
    Mars(massProduct = MARS_GM),

    /**
     * The planet Jupiter.
     */
    Jupiter(massProduct = JUPITER_GM),

    /**
     * The planet Saturn.
     */
    Saturn(massProduct = SATURN_GM),

    /**
     * The planet Uranus.
     */
    Uranus(massProduct = URANUS_GM),

    /**
     * The planet Neptune.
     */
    Neptune(massProduct = NEPTUNE_GM),

    /**
     * The planet Pluto.
     */
    Pluto(massProduct = PLUTO_GM),

    /**
     * The Sun.
     */
    Sun(massProduct = SUN_GM),

    /**
     * The Earth's natural satellite, the Moon.
     */
    Moon(massProduct = MOON_GM),

    /**
     * The Earth/Moon Barycenter.
     */
    EMB(massProduct = EARTH_GM + MOON_GM),

    /**
     * The Solar System Barycenter.
     */
    SSB(massProduct = null),
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

    val quadrature get() = (x * x) + (y * y) + (z * z)
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
    fun length() = sqrt((x * x) + (y * y) + (z * z))

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

    infix fun dot(other: AstroVector): Double {
        verifyIdenticalTimes(other.t)
        return x*other.x + y*other.y + z*other.z
    }

    fun angleWith(other: AstroVector): Double {
        val r = length() * other.length()
        val d = (this dot other) / r
        return (
            if (d <= -1.0)
                180.0
            else if (d >= +1.0)
                0.0
            else
                acos(d).radiansToDegrees()
        )
    }

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
     * *IMPORTANT:* This function differs from #AstroVector.toSpherical in two ways:
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
            sphere.lat + Astronomy.refractionAngle(refraction, sphere.lat),
            toggleAzimuthDirection(sphere.lon),
            sphere.dist
        )
    }
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
        rotate(state.position),
        rotate(state.velocity),
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
        val j = (axis + 2) % 3;
        val k = axis

        val piv = arrayOf(DoubleArray(3), DoubleArray(3), DoubleArray(3))

        piv[i][i] = c*rot[i][i] - s*rot[i][j];
        piv[i][j] = s*rot[i][i] + c*rot[i][j];
        piv[i][k] = rot[i][k];
        piv[j][i] = c*rot[j][i] - s*rot[j][j];
        piv[j][j] = s*rot[j][i] + c*rot[j][j];
        piv[j][k] = rot[j][k];
        piv[k][i] = c*rot[k][i] - s*rot[k][j];
        piv[k][j] = s*rot[k][i] + c*rot[k][j];
        piv[k][k] = rot[k][k];

        return RotationMatrix(piv)
    }

    companion object {
        /**
         * The identity rotation matrix.
         *
         * A matrix that has no effect on orientation.
         * This matrix can be the starting point for other operations,
         * such as calling a series of #RotationMatrix.combine or #RotationMatrix.pivot.
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
     * includes the time, as required by the type #AstroVector.
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
     *      #AstroVector requires a valid time value when passed to certain other functions.
     *
     * @param refraction
     *      The refraction option used to model atmospheric lensing. See #Astronomy.refractionAngle.
     *      This specifies how refraction is to be removed from the altitude stored in `this.lat`.
     *
     * @returns
     *      A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).
     */
    fun toVectorFromHorizon(time: AstroTime, refraction: Refraction): AstroVector =
        Spherical(
            lat + Astronomy.inverseRefractionAngle(refraction, lat),
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
 * Call #Astronomy.seasons to calculate this data structure for a given year.
 */
class SeasonsInfo(
    /**
     * The date and time of the March equinox for the specified year.
     */
    val marEquinox: AstroTime,

    /**
     * The date and time of the June soltice for the specified year.
     */
    val junSolstice: AstroTime,

    /**
     * The date and time of the September equinox for the specified year.
     */
    val sepEquinox: AstroTime,

    /**
     * The date and time of the December solstice for the specified year.
     */
    val decSolstice: AstroTime
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
 * Lunar libration angles, returned by #Astronomy.libration.
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
 * Returned by the function #Astronomy.searchHourAngle to report information about
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
 * See #Astronomy.elongation for more detailed information about the members of this class.
 * See also #Astronomy.searchMaxElongation for how to search for maximum elongation events.
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
 * This data structure is returned by #Astronomy.searchLunarApsis and #Astronomy.nextLunarApsis
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
 * Returned by #Astronomy.searchLunarEclipse or #Astronomy.nextLunarEclipse
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
 * Returned by #Astronomy.searchGlobalSolarEclipse or #Astronomy.nextGlobalSolarEclipse
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
 * Returned by #Astronomy.searchLocalSolarEclipse or #Astronomy.nextLocalSolarEclipse
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
 * Returned by #Astronomy.searchTransit or #Astronomy.nextTransit to report
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
 * Returned by the functions #Astronomy.illumination and #Astronomy.searchPeakMagnitude
 * to report the visual magnitude and illuminated fraction of a celestial body at a given date and time.
 */
class IllumInfo(
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
 * This structure is returned by #Astronomy.rotationAxis to report
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
 * This object is returned by #Astronomy.searchMoonNode and #Astronomy.nextMoonNode
 * to report information about the center of the Moon passing through the ecliptic plane.
 */
class NodeEventInfo(
    val time: AstroTime,
    val kind: NodeEventKind
)


interface SearchContext {
    /**
     * Evaluates a scalar function at a given time.
     *
     * @param time
     *      The time at which to evaluate the function.
     *
     * @returns
     *      The floating point value of the scalar function at the given time.
     */
    fun Eval(time: AstroTime): Double
}


/**
 * Reports the constellation that a given celestial point lies within.
 *
 * The #Astronomy.constellation function returns this object
 * to report which constellation corresponds with a given point in the sky.
 * Constellations are defined with respect to the B1875 equatorial system
 * per IAU standard. Although `Astronomy.constellation` requires J2000 equatorial
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


//---------------------------------------------------------------------------------------
// VSOP87: semi-analytic model of major planet positions

private class VsopTerm(
    val amplitude: Double,
    val phase: Double,
    val frequency: Double
)

private class VsopSeries(
    val term: Array<VsopTerm>
)

private class VsopFormula(
    val series: Array<VsopSeries>
)

private class VsopModel(
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

private class BodyState(
    val tt: Double,
    val r: TerseVector,
    val v: TerseVector
)

private fun calcVsopPosVel(model: VsopModel, tt: Double): BodyState {
    val t = tt / DAYS_PER_MILLENNIUM

    // Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates.
    val lon = vsopFormulaCalc(model.lon, t, true)
    val lat = vsopFormulaCalc(model.lat, t, false)
    val rad = vsopFormulaCalc(model.rad, t, false)

    val eclipPos = vsopSphereToRect(lon, lat, rad)

    // Calculate derivatives of spherical coordinates with respect to time.
    val dlon = vsopDerivCalc(model.lon, t);     // dlon = d(lon) / dt
    val dlat = vsopDerivCalc(model.lat, t);     // dlat = d(lat) / dt
    val drad = vsopDerivCalc(model.rad, t);     // drad = d(rad) / dt

    // Use spherical coords and spherical derivatives to calculate
    // the velocity vector in rectangular coordinates.

    val coslon = cos(lon);
    val sinlon = sin(lon);
    val coslat = cos(lat);
    val sinlat = sin(lat);

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
    // for easy of porting and consistent code generation.
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
                CO[-J,I] =  CO[J,I];
                SI[-J,I] = -SI[J,I];
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
        val lat_seconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*sin(S)-6.24*sin(3*S) + N
        return Spherical(
            lat_seconds / 3600.0,
            360.0 * frac((L0+DLAM/ARC) / PI2),
            (ARC * EARTH_EQUATORIAL_RADIUS_AU) / (0.999953253 * SINPI)
        )
    }
}

//---------------------------------------------------------------------------------------

/**
 * The main container of astronomy calculation functions.
 */
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
     * @returns
     *      The mass product of the given body in au^3/day^2.
     */
    fun massProduct(body: Body) = body.massProduct ?: throw InvalidBodyException(body)

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
     * @returns
     *      An EQD state vector that holds the geocentric position and velocity
     *      of the observer at the given time.
     */
    private fun terra(observer: Observer, time: AstroTime): StateVector {
        val st = siderealTime(time)
        val phi = observer.latitude.degreesToRadians()
        val sinphi = sin(phi)
        val cosphi = cos(phi)
        val c = 1.0 / hypot(cosphi,  EARTH_FLATTENING * sinphi)
        val s = c * EARTH_FLATTENING_SQUARED
        val ht_km = observer.height / 1000.0
        val ach = (EARTH_EQUATORIAL_RADIUS_KM * c) + ht_km
        val ash = (EARTH_EQUATORIAL_RADIUS_KM * s) + ht_km
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
     * @returns
     *      The location on or near the Earth's surface corresponding to
     *      the given position vector and time.
     */
    private fun inverseTerra(ovec: AstroVector): Observer {
        var lon_deg: Double
        var lat_deg: Double
        var height_km: Double

        // Convert from AU to kilometers.
        val x = ovec.x * KM_PER_AU
        val y = ovec.y * KM_PER_AU
        val z = ovec.z * KM_PER_AU
        val p = hypot(x, y)
        if (p < 1.0e-6) {
            // Special case: within 1 millimeter of a pole!
            // Use arbitrary longitude, and latitude determined by polarity of z.
            lon_deg = 0.0
            lat_deg = if (z > 0.0) +90.0 else -90.0
            // Elevation is calculated directly from z.
            height_km = z.absoluteValue - EARTH_POLAR_RADIUS_KM
        } else {
            // Calculate exact longitude in the half-open range (-180, +180].
            val stlocl = atan2(y, x)
            lon_deg = longitudeOffset(stlocl.radiansToDegrees() - (15 * siderealTime(ovec.t)))
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
            lat_deg = lat.radiansToDegrees()
            // Solve for exact height in kilometers.
            // There are two formulas I can use. Use whichever has the less risky denominator.
            val adjust = EARTH_EQUATORIAL_RADIUS_KM / denom
            height_km =
                if (s.absoluteValue > c.absoluteValue)
                    z/s - F*adjust
                else
                    p/c - adjust
        }

        return Observer(lat_deg, lon_deg, 1000.0 * height_km)
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
     * @returns
     *      Angular coordinates expressed in the same equatorial system as `vector`.
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
     * See #AxisInfo for more detailed information.
     *
     * @param body
     *      One of the following values:
     *      `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`,
     *      `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
     *
     * @param time
     *      The time at which to calculate the body's rotation axis.
     *
     * @returns
     *      North pole orientation and body spin angle.
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
    * To calculate an equatorial J2000 vector instead, use #Astronomy.geoMoon.
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
     * @returns
     *      The Moon's position vector in J2000 equatorial coordinates (EQJ).
     */
    fun geoMoon(time: AstroTime): AstroVector {
        val eclSphere = eclipticGeoMoon(time)
        val eclVec = eclSphere.toVector(time)
        val equVec = eclipticToEquatorial(eclVec)
        return precession(equVec, PrecessDirection.Into2000)
    }

    private fun helioEarth(time: AstroTime) = calcVsop(vsopTable[2], time)

    /**
     * Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.
     *
     * This function calculates the position of the given celestial body as a vector,
     * using the center of the Sun as the origin.  The result is expressed as a Cartesian
     * vector in the J2000 equatorial system: the coordinates are based on the mean equator
     * of the Earth at noon UTC on 1 January 2000.
     *
     * The position is not corrected for light travel time or aberration.
     * This is different from the behavior of #Astronomy.geoVector.
     *
     * If given an invalid value for `body`, this function will throw an #InvalidBodyException.
     *
     * @param body
     *      A body for which to calculate a heliocentric position:
     *      the Sun, Moon, EMB, SSB, or any of the planets.
     *
     * @param time
     *      The date and time for which to calculate the position.
     *
     * @returns
     *      The heliocentric position vector of the center of the given body.
     */
    fun helioVector(body: Body, time: AstroTime): AstroVector =
        when (body) {
            Body.Sun     -> AstroVector(0.0, 0.0, 0.0, time)
            Body.Mercury -> calcVsop(vsopTable[0], time)
            Body.Venus   -> calcVsop(vsopTable[1], time)
            Body.Earth   -> calcVsop(vsopTable[2], time)
            Body.Mars    -> calcVsop(vsopTable[3], time)
            Body.Jupiter -> calcVsop(vsopTable[4], time)
            Body.Saturn  -> calcVsop(vsopTable[5], time)
            Body.Uranus  -> calcVsop(vsopTable[6], time)
            Body.Neptune -> calcVsop(vsopTable[7], time)
            Body.Moon    -> helioEarth(time) + geoMoon(time)
            Body.EMB     -> helioEarth(time) + (geoMoon(time) / (1.0 + EARTH_MOON_MASS_RATIO))
            // FIXFIXFIX: add Pluto
            // FIXFIXFIX: add Solar System Barycenter
            else -> throw InvalidBodyException(body)
        }

    /**
     * Calculates the distance between a body and the Sun at a given time.
     *
     * Given a date and time, this function calculates the distance between
     * the center of `body` and the center of the Sun, expressed in AU.
     * For the planets Mercury through Neptune, this function is significantly
     * more efficient than calling #Astronomy.helioVector followed by taking the length
     * of the resulting vector.
     *
     * @param body
     *      A body for which to calculate a heliocentric distance:
     *      the Sun, Moon, EMB, SSB, or any of the planets.
     *
     * @param time
     *      The date and time for which to calculate the distance.
     *
     * @returns
     *      The heliocentric distance in AU.
     */
    fun helioDistance(body: Body, time: AstroTime): Double =
        when (body) {
            Body.Sun     -> 0.0
            Body.Mercury -> vsopDistance(vsopTable[0], time)
            Body.Venus   -> vsopDistance(vsopTable[1], time)
            Body.Earth   -> vsopDistance(vsopTable[2], time)
            Body.Mars    -> vsopDistance(vsopTable[3], time)
            Body.Jupiter -> vsopDistance(vsopTable[4], time)
            Body.Saturn  -> vsopDistance(vsopTable[5], time)
            Body.Uranus  -> vsopDistance(vsopTable[6], time)
            Body.Neptune -> vsopDistance(vsopTable[7], time)
            else -> helioVector(body, time).length()
        }
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



private val vsopTable: Array<VsopModel> = arrayOf (
    VsopModel(vsopLonMercury,  vsopLatMercury,   vsopRadMercury),   // 0
    VsopModel(vsopLonVenus,    vsopLatVenus,     vsopRadVenus  ),   // 1
    VsopModel(vsopLonEarth,    vsopLatEarth,     vsopRadEarth  ),   // 2
    VsopModel(vsopLonMars,     vsopLatMars,      vsopRadMars   ),   // 3
    VsopModel(vsopLonJupiter,  vsopLatJupiter,   vsopRadJupiter),   // 4
    VsopModel(vsopLonSaturn,   vsopLatSaturn,    vsopRadSaturn ),   // 5
    VsopModel(vsopLonUranus,   vsopLatUranus,    vsopRadUranus ),   // 6
    VsopModel(vsopLonNeptune,  vsopLatNeptune,   vsopRadNeptune)    // 7
)

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
