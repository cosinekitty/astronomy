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
import kotlin.math.atan2
import kotlin.math.cos
import kotlin.math.roundToLong
import kotlin.math.sin
import kotlin.math.sqrt
import kotlin.math.tan
import kotlin.math.PI

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

private const val ASEC360 = 1296000.0
private const val ASEC2RAD = 4.848136811095359935899141e-6
private const val PI2 = 2.0 * PI
private const val ARC = 3600.0 * 180.0 / PI       // arcseconds per radian
private const val SUN_RADIUS_KM  = 695700.0
private const val SUN_RADIUS_AU  = SUN_RADIUS_KM / KM_PER_AU
private const val EARTH_FLATTENING = 0.996647180302104
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

    internal fun julianCenturies() = tt / 36525.0

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

    infix fun dot(other: AstroVector): Double {
        verifyIdenticalTimes(other.t)
        return x*other.x + y*other.y + z*other.z
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
            angle = (1.02 / tan((hd+10.3/(hd+5.11)).degreesToRadians())) / 60.0

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
    fun massProduct(body: Body) =
        when (body) {
            Body.Sun     -> SUN_GM
            Body.Mercury -> MERCURY_GM
            Body.Venus   -> VENUS_GM
            Body.Earth   -> EARTH_GM
            Body.Moon    -> MOON_GM
            Body.EMB     -> EARTH_GM + MOON_GM
            Body.Mars    -> MARS_GM
            Body.Jupiter -> JUPITER_GM
            Body.Saturn  -> SATURN_GM
            Body.Uranus  -> URANUS_GM
            Body.Neptune -> NEPTUNE_GM
            Body.Pluto   -> PLUTO_GM
            else -> throw InvalidBodyException(body)
        }

    private enum class PrecessDirection {
        From2000,
        Into2000,
    }

    private fun precessionRot(time: AstroTime, dir: PrecessDirection): RotationMatrix {
        val t = time.julianCenturies()
        val eps0 = 84381.406 * ASEC2RAD

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

        val sa = sin(eps0)
        val ca = cos(eps0)
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

    private fun precession(pos: AstroVector, time: AstroTime, dir: PrecessDirection) =
        precessionRot(time, dir).rotate(pos)

    private fun precessionPosVel(state: StateVector, time: AstroTime, dir: PrecessDirection) =
        precessionRot(time, dir).rotate(state)

    private class EarthTilt(
        val tt: Double,
        val dpsi: Double,
        val deps: Double,
        val ee: Double,
        val mobl: Double,
        val tobl: Double
    )

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
        val ee = time.psi * cos(mobl.degreesToRadians()) / 15.0
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

    //==================================================================================================
    // Generated code goes to the bottom of the source file,
    // so that most line numbers are consistent between template code and target code.

    // Table of coefficients for calculating nutation using the IAU2000b model.
    private val iauRow: Array<IauRow> = arrayOf($ASTRO_IAU_DATA())
}
