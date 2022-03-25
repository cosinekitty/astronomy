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
 * Convert an angle expressed in radians to an angle expressed in degrees.
 */
fun Double.radiansToDegrees() = this * RAD2DEG

private val TimeZoneUtc = TimeZone.getTimeZone("UTC")
private const val DAYS_PER_TROPICAL_YEAR = 365.24217

private fun Double.withMinDegreeValue(min: Double): Double {
    var deg = this
    while (deg < min)
        deg += 360.0
    while (deg >= min + 360.0)
        deg -= 360.0
    return deg
}

private fun toggleAzimuthDirection(az: Double) = (360.0 - az).withMinDegreeValue(0.0)

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
    mar_equinox: AstroTime,

    /**
     * The date and time of the June soltice for the specified year.
     */
    jun_solstice: AstroTime,

    /**
     * The date and time of the September equinox for the specified year.
     */
    sep_equinox: AstroTime,

    /**
     * The date and time of the December solstice for the specified year.
     */
    dec_solstice: AstroTime
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
    val dist_km: Double,

    /**
     * The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth.
     */
    val diam_deg: Double
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
    val dist_au: Double,

    /**
     * The distance between the centers of the bodies in kilometers.
     */
    val dist_km: Double
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
 * Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
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
    val sd_penum: Double,

    /**
     * The semi-duration of the partial phase in minutes, or 0.0 if none.
     */
    val sd_partial: Double,

    /**
     * The semi-duration of the total phase in minutes, or 0.0 if none.
     */
    val sd_total: Double
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
}
