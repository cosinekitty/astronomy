/**
    @preserve

    Astronomy library for JavaScript (browser and Node.js).
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019-2021 Don Cross <cosinekitty@gmail.com>

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
export declare type FlexibleDateTime = Date | number | AstroTime;
/**
 * @brief Calculates the angle in degrees between two vectors.
 *
 * The angle is measured in the plane that contains both vectors.
 *
 * @param {Vector} a
 *      The first of a pair of vectors between which to measure an angle.
 *
 * @param {Vector} b
 *      The second of a pair of vectors between which to measure an angle.
 *
 * @returns {number}
 *      The angle between the two vectors expressed in degrees.
 *      The value is in the range [0, 180].
 */
export declare function AngleBetween(a: Vector, b: Vector): number;
/**
 * @constant {string[]} Bodies
 *      An array of strings, each a name of a supported astronomical body.
 *      Not all bodies are valid for all functions, but any string not in this
 *      list is not supported at all.
 */
export declare const Bodies: string[];
export declare function DeltaT_EspenakMeeus(ut: number): number;
export declare type DeltaTimeFunction = (ut: number) => number;
export declare function DeltaT_JplHorizons(ut: number): number;
export declare function SetDeltaTFunction(func: DeltaTimeFunction): void;
/**
 * @brief The date and time of an astronomical observation.
 *
 * Objects of type `AstroTime` are used throughout the internals
 * of the Astronomy library, and are included in certain return objects.
 * Use the constructor or the {@link MakeTime} function to create an `AstroTime` object.
 *
 * @property {Date} date
 *      The JavaScript Date object for the given date and time.
 *      This Date corresponds to the numeric day value stored in the `ut` property.
 *
 * @property {number} ut
 *      Universal Time (UT1/UTC) in fractional days since the J2000 epoch.
 *      Universal Time represents time measured with respect to the Earth's rotation,
 *      tracking mean solar days.
 *      The Astronomy library approximates UT1 and UTC as being the same thing.
 *      This gives sufficient accuracy for the precision requirements of this project.
 *
 * @property {number} tt
 *      Terrestrial Time in fractional days since the J2000 epoch.
 *      TT represents a continuously flowing ephemeris timescale independent of
 *      any variations of the Earth's rotation, and is adjusted from UT
 *      using historical and predictive models of those variations.
 */
export declare class AstroTime {
    date: Date;
    ut: number;
    tt: number;
    /**
     * @param {FlexibleDateTime} date
     *      A JavaScript Date object, a numeric UTC value expressed in J2000 days, or another AstroTime object.
     */
    constructor(date: FlexibleDateTime);
    /**
     * Formats an `AstroTime` object as an [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601)
     * date/time string in UTC, to millisecond resolution.
     * Example: `2018-08-17T17:22:04.050Z`
     * @returns {string}
     */
    toString(): string;
    /**
     * Returns a new `AstroTime` object adjusted by the floating point number of days.
     * Does NOT modify the original `AstroTime` object.
     *
     * @param {number} days
     *      The floating point number of days by which to adjust the given date and time.
     *      Positive values adjust the date toward the future, and
     *      negative values adjust the date toward the past.
     *
     * @returns {AstroTime}
     */
    AddDays(days: number): AstroTime;
}
/**
 * @brief A `Date`, `number`, or `AstroTime` value that specifies the date and time of an astronomical event.
 *
 * `FlexibleDateTime` is a placeholder type that represents three different types
 * that may be passed to many Astronomy Engine functions: a JavaScript `Date` object,
 * a number representing the real-valued number of UT days since the J2000 epoch,
 * or an {@link AstroTime} object.
 *
 * This flexibility is for convenience of outside callers.
 * Internally, Astronomy Engine always converts a `FlexibleTime` parameter
 * to an `AstroTime` object by calling {@link MakeTime}.
 *
 * @typedef {Date | number | AstroTime} FlexibleDateTime
 */
/**
 * @brief Converts multiple date/time formats to `AstroTime` format.
 *
 * Given a Date object or a number days since noon (12:00) on January 1, 2000 (UTC),
 * this function creates an {@link AstroTime} object.
 *
 * Given an {@link AstroTime} object, returns the same object unmodified.
 * Use of this function is not required for any of the other exposed functions in this library,
 * because they all guarantee converting date/time parameters to `AstroTime`
 * as needed. However, it may be convenient for callers who need to understand
 * the difference between UTC and TT (Terrestrial Time). In some use cases,
 * converting once to `AstroTime` format and passing the result into multiple
 * function calls may be more efficient than passing in native JavaScript Date objects.
 *
 * @param {FlexibleDateTime} date
 *      A Date object, a number of UTC days since the J2000 epoch (noon on January 1, 2000),
 *      or an AstroTime object. See remarks above.
 *
 * @returns {AstroTime}
 */
export declare function MakeTime(date: FlexibleDateTime): AstroTime;
export declare let CalcMoonCount: number;
/**
 * @brief A 3D Cartesian vector with a time attached to it.
 *
 * Holds the Cartesian coordinates of a vector in 3D space,
 * along with the time at which the vector is valid.
 *
 * @property {number} x        The x-coordinate expressed in astronomical units (AU).
 * @property {number} y        The y-coordinate expressed in astronomical units (AU).
 * @property {number} z        The z-coordinate expressed in astronomical units (AU).
 * @property {AstroTime} t     The time at which the vector is valid.
 */
export declare class Vector {
    x: number;
    y: number;
    z: number;
    t: AstroTime;
    constructor(x: number, y: number, z: number, t: AstroTime);
    /**
     * Returns the length of the vector in astronomical units (AU).
     * @returns {number}
     */
    Length(): number;
}
/**
 * @brief Holds spherical coordinates: latitude, longitude, distance.
 *
 * Spherical coordinates represent the location of
 * a point using two angles and a distance.
 *
 * @property {number} lat       The latitude angle: -90..+90 degrees.
 * @property {number} lon       The longitude angle: 0..360 degrees.
 * @property {number} dist      Distance in AU.
 */
export declare class Spherical {
    lat: number;
    lon: number;
    dist: number;
    constructor(lat: number, lon: number, dist: number);
}
/**
 * @brief Holds right ascension, declination, and distance of a celestial object.
 *
 * @property {number} ra
 *      Right ascension in sidereal hours: [0, 24).
 *
 * @property {number} dec
 *      Declination in degrees: [-90, +90].
 *
 * @property {number} dist
 *      Distance to the celestial object expressed in
 *      <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU).
 */
export declare class EquatorialCoordinates {
    ra: number;
    dec: number;
    dist: number;
    constructor(ra: number, dec: number, dist: number);
}
/**
 * @brief Contains a rotation matrix that can be used to transform one coordinate system to another.
 *
 * @property {number[][]} rot
 *      A normalized 3x3 rotation matrix. For example, the identity matrix is represented
 *      as `[[1, 0, 0], [0, 1, 0], [0, 0, 1]]`.
 */
export declare class RotationMatrix {
    rot: number[][];
    constructor(rot: number[][]);
}
/**
 * @brief Creates a rotation matrix that can be used to transform one coordinate system to another.
 *
 * This function verifies that the `rot` parameter is of the correct format:
 * a number[3][3] array. It throws an exception if `rot` is not of that shape.
 * Otherwise it creates a new {@link RotationMatrix} object based on `rot`.
 *
 * @param {number[][]} rot
 *      An array [3][3] of numbers. Defines a rotation matrix used to premultiply
 *      a 3D vector to reorient it into another coordinate system.
 *
 * @returns {RotationMatrix}
 */
export declare function MakeRotation(rot: number[][]): RotationMatrix;
/**
 * @brief Represents the location of an object seen by an observer on the Earth.
 *
 * Holds azimuth (compass direction) and altitude (angle above/below the horizon)
 * of a celestial object as seen by an observer at a particular location on the Earth's surface.
 * Also holds right ascension and declination of the same object.
 * All of these coordinates are optionally adjusted for atmospheric refraction;
 * therefore the right ascension and declination values may not exactly match
 * those found inside a corresponding {@link EquatorialCoordinates} object.
 *
 * @property {number} azimuth
 *      A horizontal compass direction angle in degrees measured starting at north
 *      and increasing positively toward the east.
 *      The value is in the range [0, 360).
 *      North = 0, east = 90, south = 180, west = 270.
 *
 * @property {number} altitude
 *      A vertical angle in degrees above (positive) or below (negative) the horizon.
 *      The value is in the range [-90, +90].
 *      The altitude angle is optionally adjusted upward due to atmospheric refraction.
 *
 * @property {number} ra
 *      The right ascension of the celestial body in sidereal hours.
 *      The value is in the reange [0, 24).
 *      If `altitude` was adjusted for atmospheric reaction, `ra`
 *      is likewise adjusted.
 *
 * @property {number} dec
 *      The declination of of the celestial body in degrees.
 *      The value in the range [-90, +90].
 *      If `altitude` was adjusted for atmospheric reaction, `dec`
 *      is likewise adjusted.
 */
export declare class HorizontalCoordinates {
    azimuth: number;
    altitude: number;
    ra: number;
    dec: number;
    constructor(azimuth: number, altitude: number, ra: number, dec: number);
}
/**
 * @brief Ecliptic coordinates of a celestial body.
 *
 * The origin and date of the coordinate system may vary depending on the caller's usage.
 * In general, ecliptic coordinates are measured with respect to the mean plane of the Earth's
 * orbit around the Sun.
 * Includes Cartesian coordinates `(ex, ey, ez)` measured in
 * <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU)
 * and spherical coordinates `(elon, elat)` measured in degrees.
 *
 * @property {number} ex
 *      The Cartesian x-coordinate of the body in astronomical units (AU).
 *      The x-axis is within the ecliptic plane and is oriented in the direction of the
 *      <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>.
 *
 * @property {number} ey
 *      The Cartesian y-coordinate of the body in astronomical units (AU).
 *      The y-axis is within the ecliptic plane and is oriented 90 degrees
 *      counterclockwise from the equinox, as seen from above the Sun's north pole.
 *
 * @property {number} ez
 *      The Cartesian z-coordinate of the body in astronomical units (AU).
 *      The z-axis is oriented perpendicular to the ecliptic plane,
 *      along the direction of the Sun's north pole.
 *
 * @property {number} elat
 *      The ecliptic latitude of the body in degrees.
 *      This is the angle north or south of the ecliptic plane.
 *      The value is in the range [-90, +90].
 *      Positive values are north and negative values are south.
 *
 * @property {number} elon
 *      The ecliptic longitude of the body in degrees.
 *      This is the angle measured counterclockwise around the ecliptic plane,
 *      as seen from above the Sun's north pole.
 *      This is the same direction that the Earth orbits around the Sun.
 *      The angle is measured starting at 0 from the equinox and increases
 *      up to 360 degrees.
 */
export declare class EclipticCoordinates {
    ex: number;
    ey: number;
    ez: number;
    elat: number;
    elon: number;
    constructor(ex: number, ey: number, ez: number, elat: number, elon: number);
}
/**
 * @brief Converts equatorial coordinates to horizontal coordinates.
 *
 * Given a date and time, a geographic location of an observer on the Earth, and
 * equatorial coordinates (right ascension and declination) of a celestial body,
 * returns horizontal coordinates (azimuth and altitude angles) for that body
 * as seen by that observer. Allows optional correction for atmospheric refraction.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to find horizontal coordinates.
 *
 * @param {Observer} observer
 *      The location of the observer for which to find horizontal coordinates.
 *
 * @param {number} ra
 *      Right ascension in sidereal hours of the celestial object,
 *      referred to the mean equinox of date for the J2000 epoch.
 *
 * @param {number} dec
 *      Declination in degrees of the celestial object,
 *      referred to the mean equator of date for the J2000 epoch.
 *      Positive values are north of the celestial equator and negative values are south.
 *
 * @param {string} refraction
 *      If omitted or has a false-like value (false, null, undefined, etc.)
 *      the calculations are performed without any correction for atmospheric
 *      refraction. If the value is the string `"normal"`,
 *      uses the recommended refraction correction based on Meeus "Astronomical Algorithms"
 *      with a linear taper more than 1 degree below the horizon. The linear
 *      taper causes the refraction to linearly approach 0 as the altitude of the
 *      body approaches the nadir (-90 degrees).
 *      If the value is the string `"jplhor"`, uses a JPL Horizons
 *      compatible formula. This is the same algorithm as `"normal"`,
 *      only without linear tapering; this can result in physically impossible
 *      altitudes of less than -90 degrees, which may cause problems for some applications.
 *      (The `"jplhor"` option was created for unit testing against data
 *      generated by JPL Horizons, and is otherwise not recommended for use.)
 *
 * @returns {HorizontalCoordinates}
 */
export declare function Horizon(date: FlexibleDateTime, observer: Observer, ra: number, dec: number, refraction?: string): HorizontalCoordinates;
/**
 * @brief Represents the geographic location of an observer on the surface of the Earth.
 *
 * @property {number} latitude
 *      The observer's geographic latitude in degrees north of the Earth's equator.
 *      The value is negative for observers south of the equator.
 *      Must be in the range -90 to +90.
 *
 * @property {number} longitude
 *      The observer's geographic longitude in degrees east of the prime meridian
 *      passing through Greenwich, England.
 *      The value is negative for observers west of the prime meridian.
 *      The value should be kept in the range -180 to +180 to minimize floating point errors.
 *
 * @property {number} height
 *      The observer's elevation above mean sea level, expressed in meters.
 */
export declare class Observer {
    latitude: number;
    longitude: number;
    height: number;
    constructor(latitude: number, longitude: number, height: number);
}
/**
 * @brief Returns apparent geocentric true ecliptic coordinates of date for the Sun.
 *
 * This function is used for calculating the times of equinoxes and solstices.
 *
 * <i>Geocentric</i> means coordinates as the Sun would appear to a hypothetical observer
 * at the center of the Earth.
 * <i>Ecliptic coordinates of date</i> are measured along the plane of the Earth's mean
 * orbit around the Sun, using the
 * <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a>
 * of the Earth as adjusted for precession and nutation of the Earth's
 * axis of rotation on the given date.
 *
 * @param {FlexibleDateTime} date
 *      The date and time at which to calculate the Sun's apparent location as seen from
 *      the center of the Earth.
 *
 * @returns {EclipticCoordinates}
 */
export declare function SunPosition(date: FlexibleDateTime): EclipticCoordinates;
/**
 * @brief Calculates equatorial coordinates of a Solar System body at a given time.
 *
 * Returns topocentric equatorial coordinates (right ascension and declination)
 * in one of two different systems: J2000 or true-equator-of-date.
 * Allows optional correction for aberration.
 * Always corrects for light travel time (represents the object as seen by the observer
 * with light traveling to the Earth at finite speed, not where the object is right now).
 * <i>Topocentric</i> refers to a position as seen by an observer on the surface of the Earth.
 * This function corrects for
 * <a href="https://en.wikipedia.org/wiki/Parallax">parallax</a>
 * of the object between a geocentric observer and a topocentric observer.
 * This is most significant for the Moon, because it is so close to the Earth.
 * However, it can have a small effect on the apparent positions of other bodies.
 *
 * @param {string} body
 *      The name of the body for which to find equatorial coordinates.
 *      Not allowed to be `"Earth"`.
 *
 * @param {FlexibleDateTime} date
 *      Specifies the date and time at which the body is to be observed.
 *
 * @param {Observer} observer
 *      The location on the Earth of the observer.
 *
 * @param {bool} ofdate
 *      Pass `true` to return equatorial coordinates of date,
 *      i.e. corrected for precession and nutation at the given date.
 *      This is needed to get correct horizontal coordinates when you call
 *      {@link Horizon}.
 *      Pass `false` to return equatorial coordinates in the J2000 system.
 *
 * @param {bool} aberration
 *      Pass `true` to correct for
 *      <a href="https://en.wikipedia.org/wiki/Aberration_of_light">aberration</a>,
 *      or `false` to leave uncorrected.
 *
 * @returns {EquatorialCoordinates}
 *      The topocentric coordinates of the body as adjusted for the given observer.
 */
export declare function Equator(body: string, date: FlexibleDateTime, observer: Observer, ofdate: boolean, aberration: boolean): EquatorialCoordinates;
/**
 * @brief Converts equatorial Cartesian coordinates to ecliptic Cartesian and angular coordinates.
 *
 * Given J2000 equatorial Cartesian coordinates,
 * returns J2000 ecliptic latitude, longitude, and cartesian coordinates.
 * You can call {@link GeoVector} and use its (x, y, z) return values
 * to pass into this function.
 *
 * @param {number} gx
 *      The x-coordinate of a 3D vector in the J2000 equatorial coordinate system.
 *
 * @param {number} gy
 *      The y-coordinate of a 3D vector in the J2000 equatorial coordinate system.
 *
 * @param {number} gz
 *      The z-coordinate of a 3D vector in the J2000 equatorial coordinate system.
 *
 * @returns {EclipticCoordinates}
 */
export declare function Ecliptic(gx: number, gy: number, gz: number): EclipticCoordinates;
/**
 * @brief Calculates the geocentric Cartesian coordinates for the Moon in the J2000 equatorial system.
 *
 * Based on the Nautical Almanac Office's <i>Improved Lunar Ephemeris</i> of 1954,
 * which in turn derives from E. W. Brown's lunar theories.
 * Adapted from Turbo Pascal code from the book
 * <a href="https://www.springer.com/us/book/9783540672210">Astronomy on the Personal Computer</a>
 * by Montenbruck and Pfleger.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the Moon's geocentric position.
 *
 * @returns {Vector}
 */
export declare function GeoMoon(date: FlexibleDateTime): Vector;
/**
 * @brief Calculates a vector from the center of the Sun to the given body at the given time.
 *
 * Calculates heliocentric (i.e., with respect to the center of the Sun)
 * Cartesian coordinates in the J2000 equatorial system of a celestial
 * body at a specified time. The position is not corrected for light travel time or aberration.
 *
 * @param {string} body
 *      One of the strings
 *      `"Sun"`, `"Moon"`, `"Mercury"`, `"Venus"`,
 *      `"Earth"`, `"Mars"`, `"Jupiter"`, `"Saturn"`,
 *      `"Uranus"`, `"Neptune"`, `"Pluto"`,
 *      `"SSB"`, or `"EMB"`.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which the body's position is to be calculated.
 *
 * @returns {Vector}
 */
export declare function HelioVector(body: string, date: FlexibleDateTime): Vector;
/**
 * @brief Calculates the distance between a body and the Sun at a given time.
 *
 * Given a date and time, this function calculates the distance between
 * the center of `body` and the center of the Sun.
 * For the planets Mercury through Neptune, this function is significantly
 * more efficient than calling {@link HelioVector} followed by taking the length
 * of the resulting vector.
 *
 * @param {string} body
 *      A body for which to calculate a heliocentric distance:
 *      the Sun, Moon, or any of the planets.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the heliocentric distance.
 *
 * @returns {number}
 *      The heliocentric distance in AU.
 */
export declare function HelioDistance(body: string, date: FlexibleDateTime): number;
/**
 * @brief Calculates a vector from the center of the Earth to the given body at the given time.
 *
 * Calculates geocentric (i.e., with respect to the center of the Earth)
 * Cartesian coordinates in the J2000 equatorial system of a celestial
 * body at a specified time. The position is always corrected for light travel time:
 * this means the position of the body is "back-dated" based on how long it
 * takes light to travel from the body to an observer on the Earth.
 * Also, the position can optionally be corrected for aberration, an effect
 * causing the apparent direction of the body to be shifted based on
 * transverse movement of the Earth with respect to the rays of light
 * coming from that body.
 *
 * @param {string} body
 *      One of the strings
 *      `"Sun"`, `"Moon"`, `"Mercury"`, `"Venus"`,
 *      `"Earth"`, `"Mars"`, `"Jupiter"`, `"Saturn"`,
 *      `"Uranus"`, `"Neptune"`, or `"Pluto"`.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which the body's position is to be calculated.
 *
 * @param {bool} aberration
 *      Pass `true` to correct for
 *      <a href="https://en.wikipedia.org/wiki/Aberration_of_light">aberration</a>,
 *      or `false` to leave uncorrected.
 *
 * @returns {Vector}
 */
export declare function GeoVector(body: string, date: FlexibleDateTime, aberration: boolean): Vector;
export interface SearchOptions {
    dt_tolerance_seconds?: number;
    init_f1?: number;
    init_f2?: number;
    iter_limit?: number;
}
/**
 * @brief Options for the {@link Search} function.
 *
 * @typedef {object} SearchOptions
 *
 * @property {number | undefined} dt_tolerance_seconds
 *      The number of seconds for a time window smaller than which the search
 *      is considered successful.  Using too large a tolerance can result in
 *      an inaccurate time estimate.  Using too small a tolerance can cause
 *      excessive computation, or can even cause the search to fail because of
 *      limited floating-point resolution.  Defaults to 1 second.
 *
 * @property {number | undefined} init_f1
 *      As an optimization, if the caller of {@link Search}
 *      has already calculated the value of the function being searched (the parameter `func`)
 *      at the time coordinate `t1`, it can pass in that value as `init_f1`.
 *      For very expensive calculations, this can measurably improve performance.
 *
 * @property {number | undefined} init_f2
 *      The same as `init_f1`, except this is the optional initial value of `func(t2)`
 *      instead of `func(t1)`.
 *
 * @property {number | undefined} iter_limit
 */
/**
 * @brief Finds the time when a function ascends through zero.
 *
 * Search for next time <i>t</i> (such that <i>t</i> is between `t1` and `t2`)
 * that `func(t)` crosses from a negative value to a non-negative value.
 * The given function must have "smooth" behavior over the entire inclusive range [`t1`, `t2`],
 * meaning that it behaves like a continuous differentiable function.
 * It is not required that `t1` &lt; `t2`; `t1` &gt; `t2`
 * allows searching backward in time.
 * Note: `t1` and `t2` must be chosen such that there is no possibility
 * of more than one zero-crossing (ascending or descending), or it is possible
 * that the "wrong" event will be found (i.e. not the first event after t1)
 * or even that the function will return `null`, indicating that no event was found.
 *
 * @param {function(AstroTime): number} func
 *      The function to find an ascending zero crossing for.
 *      The function must accept a single parameter of type {@link AstroTime}
 *      and return a numeric value.
 *
 * @param {AstroTime} t1
 *      The lower time bound of a search window.
 *
 * @param {AstroTime} t2
 *      The upper time bound of a search window.
 *
 * @param {SearchOptions | undefined} options
 *      Options that can tune the behavior of the search.
 *      Most callers can omit this argument.
 *
 * @returns {AstroTime | null}
 *      If the search is successful, returns the date and time of the solution.
 *      If the search fails, returns `null`.
 */
export declare function Search(f: (t: AstroTime) => number, t1: AstroTime, t2: AstroTime, options?: SearchOptions): AstroTime | null;
/**
 * @brief Searches for when the Sun reaches a given ecliptic longitude.
 *
 * Searches for the moment in time when the center of the Sun reaches a given apparent
 * ecliptic longitude, as seen from the center of the Earth, within a given range of dates.
 * This function can be used to determine equinoxes and solstices.
 * However, it is usually more convenient and efficient to call {@link Seasons}
 * to calculate equinoxes and solstices for a given calendar year.
 * `SearchSunLongitude` is more general in that it allows searching for arbitrary longitude values.
 *
 * @param {number} targetLon
 *      The desired ecliptic longitude of date in degrees.
 *      This may be any value in the range [0, 360), although certain
 *      values have conventional meanings:
 *
 *      When `targetLon` is 0, finds the March equinox,
 *      which is the moment spring begins in the northern hemisphere
 *      and the beginning of autumn in the southern hemisphere.
 *
 *      When `targetLon` is 180, finds the September equinox,
 *      which is the moment autumn begins in the northern hemisphere and
 *      spring begins in the southern hemisphere.
 *
 *      When `targetLon` is 90, finds the northern solstice, which is the
 *      moment summer begins in the northern hemisphere and winter
 *      begins in the southern hemisphere.
 *
 *      When `targetLon` is 270, finds the southern solstice, which is the
 *      moment winter begins in the northern hemisphere and summer
 *      begins in the southern hemisphere.
 *
 * @param {FlexibleDateTime} dateStart
 *      A date and time known to be earlier than the desired longitude event.
 *
 * @param {number} limitDays
 *      A floating point number of days, which when added to `dateStart`,
 *      yields a date and time known to be after the desired longitude event.
 *
 * @returns {AstroTime | null}
 *      The date and time when the Sun reaches the apparent ecliptic longitude `targetLon`
 *      within the range of times specified by `dateStart` and `limitDays`.
 *      If the Sun does not reach the target longitude within the specified time range, or the
 *      time range is excessively wide, the return value is `null`.
 *      To avoid a `null` return value, the caller must pick a time window around
 *      the event that is within a few days but not so small that the event might fall outside the window.
 */
export declare function SearchSunLongitude(targetLon: number, dateStart: FlexibleDateTime, limitDays: number): AstroTime | null;
/**
 * @brief Calculates the longitude separation between the Sun and the given body.
 *
 * Calculates the ecliptic longitude difference
 * between the given body and the Sun as seen from
 * the Earth at a given moment in time.
 * The returned value ranges [0, 360) degrees.
 * By definition, the Earth and the Sun are both in the plane of the ecliptic.
 * Ignores the height of the `body` above or below the ecliptic plane;
 * the resulting angle is measured around the ecliptic plane for the "shadow"
 * of the body onto that plane.
 *
 * Use {@link AngleFromSun} instead, if you wish to calculate the full angle
 * between the Sun and a body, instead of just their longitude difference.
 *
 * @param {string} body
 *      The name of a supported celestial body other than the Earth.
 *
 * @param {FlexibleDateTime} date
 *      The time at which the relative longitude is to be found.
 *
 * @returns {number}
 *      An angle in degrees in the range [0, 360).
 *      Values less than 180 indicate that the body is to the east
 *      of the Sun as seen from the Earth; that is, the body sets after
 *      the Sun does and is visible in the evening sky.
 *      Values greater than 180 indicate that the body is to the west of
 *      the Sun and is visible in the morning sky.
 */
export declare function LongitudeFromSun(body: string, date: FlexibleDateTime): number;
/**
 * @brief Calculates the angular separation between the Sun and the given body.
 *
 * Returns the full angle seen from
 * the Earth, between the given body and the Sun.
 * Unlike {@link LongitudeFromSun}, this function does not
 * project the body's "shadow" onto the ecliptic;
 * the angle is measured in 3D space around the plane that
 * contains the centers of the Earth, the Sun, and `body`.
 *
 * @param {string} body
 *      The name of a supported celestial body other than the Earth.
 *
 * @param {FlexibleDateTime} date
 *      The time at which the angle from the Sun is to be found.
 *
 * @returns {number}
 *      An angle in degrees in the range [0, 180].
 */
export declare function AngleFromSun(body: string, date: FlexibleDateTime): number;
/**
 * @brief Calculates heliocentric ecliptic longitude based on the J2000 equinox.
 *
 * @param {string} body
 *      The name of a celestial body other than the Sun.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the ecliptic longitude.
 *
 * @returns {number}
 *      The ecliptic longitude angle of the body in degrees measured counterclockwise around the mean
 *      plane of the Earth's orbit, as seen from above the Sun's north pole.
 *      Ecliptic longitude starts at 0 at the J2000
 *      <a href="https://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)">equinox</a> and
 *      increases in the same direction the Earth orbits the Sun.
 *      The returned value is always in the range [0, 360).
 */
export declare function EclipticLongitude(body: string, date: FlexibleDateTime): number;
/**
 * @brief Information about the apparent brightness and sunlit phase of a celestial object.
 *
 * @property {AstroTime} time
 *      The date and time pertaining to the other calculated values in this object.
 *
 * @property {number} mag
 *      The <a href="https://en.wikipedia.org/wiki/Apparent_magnitude">apparent visual magnitude</a> of the celestial body.
 *
 * @property {number} phase_angle
 *      The angle in degrees as seen from the center of the celestial body between the Sun and the Earth.
 *      The value is always in the range 0 to 180.
 *      The phase angle provides a measure of what fraction of the body's face appears
 *      illuminated by the Sun as seen from the Earth.
 *      When the observed body is the Sun, the `phase` property is set to 0,
 *      although this has no physical meaning because the Sun emits, rather than reflects, light.
 *      When the phase is near 0 degrees, the body appears "full".
 *      When it is 90 degrees, the body appears "half full".
 *      And when it is 180 degrees, the body appears "new" and is very difficult to see
 *      because it is both dim and lost in the Sun's glare as seen from the Earth.
 *
 * @property {number} phase_fraction
 *      The fraction of the body's face that is illuminated by the Sun, as seen from the Earth.
 *      Calculated from `phase_angle` for convenience.
 *      This value ranges from 0 to 1.
 *
 * @property {number} helio_dist
 *      The distance between the center of the Sun and the center of the body in
 *      <a href="https://en.wikipedia.org/wiki/Astronomical_unit">astronomical units</a> (AU).
 *
 * @property {number} geo_dist
 *      The distance between the center of the Earth and the center of the body in AU.
 *
 * @property {Vector} gc
 *      Geocentric coordinates: the 3D vector from the center of the Earth to the center of the body.
 *      The components are in expressed in AU and are oriented with respect to the J2000 equatorial plane.
 *
 * @property {Vector} hc
 *      Heliocentric coordinates: The 3D vector from the center of the Sun to the center of the body.
 *      Like `gc`, `hc` is expressed in AU and oriented with respect
 *      to the J2000 equatorial plane.
 *
 * @property {number | undefined} ring_tilt
 *      For Saturn, this is the angular tilt of the planet's rings in degrees away
 *      from the line of sight from the Earth. When the value is near 0, the rings
 *      appear edge-on from the Earth and are therefore difficult to see.
 *      When `ring_tilt` approaches its maximum value (about 27 degrees),
 *      the rings appear widest and brightest from the Earth.
 *      Unlike the <a href="https://ssd.jpl.nasa.gov/horizons.cgi">JPL Horizons</a> online tool,
 *      this library includes the effect of the ring tilt angle in the calculated value
 *      for Saturn's visual magnitude.
 *      For all bodies other than Saturn, the value of `ring_tilt` is `undefined`.
 */
export declare class IlluminationInfo {
    time: AstroTime;
    mag: number;
    phase_angle: number;
    helio_dist: number;
    geo_dist: number;
    gc: Vector;
    hc: Vector;
    ring_tilt?: number | undefined;
    phase_fraction: number;
    constructor(time: AstroTime, mag: number, phase_angle: number, helio_dist: number, geo_dist: number, gc: Vector, hc: Vector, ring_tilt?: number | undefined);
}
/**
 * @brief Calculates visual magnitude and related information about a body.
 *
 * Calculates the phase angle, visual magnitude,
 * and other values relating to the body's illumination
 * at the given date and time, as seen from the Earth.
 *
 * @param {string} body
 *      The name of the celestial body being observed.
 *      Not allowed to be `"Earth"`.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the illumination data for the given body.
 *
 * @returns {IlluminationInfo}
 */
export declare function Illumination(body: string, date: FlexibleDateTime): IlluminationInfo;
/**
 * @brief Searches for when the Earth and a given body reach a relative ecliptic longitude separation.
 *
 * Searches for the date and time the relative ecliptic longitudes of
 * the specified body and the Earth, as seen from the Sun, reach a certain
 * difference. This function is useful for finding conjunctions and oppositions
 * of the planets. For the opposition of a superior planet (Mars, Jupiter, ..., Pluto),
 * or the inferior conjunction of an inferior planet (Mercury, Venus),
 * call with `targetRelLon` = 0. The 0 value indicates that both
 * planets are on the same ecliptic longitude line, ignoring the other planet's
 * distance above or below the plane of the Earth's orbit.
 * For superior conjunctions, call with `targetRelLon` = 180.
 * This means the Earth and the other planet are on opposite sides of the Sun.
 *
 * @param {string} body
 *      The name of a planet other than the Earth.
 *
 * @param {number} targetRelLon
 *      The desired angular difference in degrees between the ecliptic longitudes
 *      of `body` and the Earth. Must be in the range (-180, +180].
 *
 * @param {FlexibleDateTime} startDate
 *      The date and time after which to find the next occurrence of the
 *      body and the Earth reaching the desired relative longitude.
 *
 * @returns {AstroTime}
 *      The time when the Earth and the body next reach the specified relative longitudes.
 */
export declare function SearchRelativeLongitude(body: string, targetRelLon: number, startDate: FlexibleDateTime): AstroTime;
/**
 * @brief Determines the moon's phase expressed as an ecliptic longitude.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for which to calculate the moon's phase.
 *
 * @returns {number}
 *      A value in the range [0, 360) indicating the difference
 *      in ecliptic longitude between the center of the Sun and the
 *      center of the Moon, as seen from the center of the Earth.
 *      Certain longitude values have conventional meanings:
 *
 * * 0 = new moon
 * * 90 = first quarter
 * * 180 = full moon
 * * 270 = third quarter
 */
export declare function MoonPhase(date: FlexibleDateTime): number;
/**
 * @brief Searches for the date and time that the Moon reaches a specified phase.
 *
 * Lunar phases are defined in terms of geocentric ecliptic longitudes
 * with respect to the Sun.  When the Moon and the Sun have the same ecliptic
 * longitude, that is defined as a new moon. When the two ecliptic longitudes
 * are 180 degrees apart, that is defined as a full moon.
 * To enumerate quarter lunar phases, it is simpler to call
 * {@link SearchMoonQuarter} once, followed by repeatedly calling
 * {@link NextMoonQuarter}. `SearchMoonPhase` is only
 * necessary for finding other lunar phases than the usual quarter phases.
 *
 * @param {number} targetLon
 *      The difference in geocentric ecliptic longitude between the Sun and Moon
 *      that specifies the lunar phase being sought. This can be any value
 *      in the range [0, 360). Here are some helpful examples:
 *      0 = new moon,
 *      90 = first quarter,
 *      180 = full moon,
 *      270 = third quarter.
 *
 * @param {FlexibleDateTime} dateStart
 *      The beginning of the window of time in which to search.
 *
 * @param {number} limitDays
 *      The floating point number of days after `dateStart`
 *      that limits the window of time in which to search.
 *
 * @returns {AstroTime | null}
 *      If the specified lunar phase occurs after `dateStart`
 *      and before `limitDays` days after `dateStart`,
 *      this function returns the date and time of the first such occurrence.
 *      Otherwise, it returns `null`.
 */
export declare function SearchMoonPhase(targetLon: number, dateStart: FlexibleDateTime, limitDays: number): AstroTime | null;
/**
 * @brief A quarter lunar phase, along with when it occurs.
 *
 * @property {number} quarter
 *      An integer as follows:
 *      0 = new moon,
 *      1 = first quarter,
 *      2 = full moon,
 *      3 = third quarter.
 *
 * @property {AstroTime} time
 *      The date and time of the quarter lunar phase.
 */
export declare class MoonQuarter {
    quarter: number;
    time: AstroTime;
    constructor(quarter: number, time: AstroTime);
}
/**
 * @brief Finds the first quarter lunar phase after the specified date and time.
 *
 * The quarter lunar phases are: new moon, first quarter, full moon, and third quarter.
 * To enumerate quarter lunar phases, call `SearchMoonQuarter` once,
 * then pass its return value to {@link NextMoonQuarter} to find the next
 * `MoonQuarter`. Keep calling `NextMoonQuarter` in a loop,
 * passing the previous return value as the argument to the next call.
 *
 * @param {FlexibleDateTime} dateStart
 *      The date and time after which to find the first quarter lunar phase.
 *
 * @returns {MoonQuarter}
 */
export declare function SearchMoonQuarter(dateStart: FlexibleDateTime): MoonQuarter;
/**
 * @brief Finds the next quarter lunar phase in a series.
 *
 * Given a {@link MoonQuarter} object, finds the next consecutive
 * quarter lunar phase. See remarks in {@link SearchMoonQuarter}
 * for explanation of usage.
 *
 * @param {MoonQuarter} mq
 *      The return value of a prior call to {@link MoonQuarter} or `NextMoonQuarter`.
 */
export declare function NextMoonQuarter(mq: MoonQuarter): MoonQuarter;
/**
 * @brief Finds the next rise or set time for a body.
 *
 * Finds a rise or set time for the given body as
 * seen by an observer at the specified location on the Earth.
 * Rise time is defined as the moment when the top of the body
 * is observed to first appear above the horizon in the east.
 * Set time is defined as the moment the top of the body
 * is observed to sink below the horizon in the west.
 * The times are adjusted for typical atmospheric refraction conditions.
 *
 * @param {string} body
 *      The name of the body to find the rise or set time for.
 *
 * @param {Observer} observer
 *      Specifies the geographic coordinates and elevation above sea level of the observer.
 *
 * @param {number} direction
 *      Either +1 to find rise time or -1 to find set time.
 *      Any other value will cause an exception to be thrown.
 *
 * @param {FlexibleDateTime} dateStart
 *      The date and time after which the specified rise or set time is to be found.
 *
 * @param {number} limitDays
 *      The fractional number of days after `dateStart` that limits
 *      when the rise or set time is to be found.
 *
 * @returns {AstroTime | null}
 *      The date and time of the rise or set event, or null if no such event
 *      occurs within the specified time window.
 */
export declare function SearchRiseSet(body: string, observer: Observer, direction: number, dateStart: FlexibleDateTime, limitDays: number): AstroTime | null;
/**
 * @brief Horizontal position of a body upon reaching an hour angle.
 *
 * Returns information about an occurrence of a celestial body
 * reaching a given hour angle as seen by an observer at a given
 * location on the surface of the Earth.
 *
 * @property {AstroTime} time
 *      The date and time of the celestial body reaching the hour angle.
 *
 * @property {HorizontalCoordinates} hor
 *      Topocentric horizontal coordinates for the body
 *      at the time indicated by the `time` property.
 */
export declare class HourAngleEvent {
    time: AstroTime;
    hor: HorizontalCoordinates;
    constructor(time: AstroTime, hor: HorizontalCoordinates);
}
/**
 * @brief Finds when a body will reach a given hour angle.
 *
 * Finds the next time the given body is seen to reach the specified
 * <a href="https://en.wikipedia.org/wiki/Hour_angle">hour angle</a>
 * by the given observer.
 * Providing `hourAngle` = 0 finds the next maximum altitude event (culmination).
 * Providing `hourAngle` = 12 finds the next minimum altitude event.
 * Note that, especially close to the Earth's poles, a body as seen on a given day
 * may always be above the horizon or always below the horizon, so the caller cannot
 * assume that a culminating object is visible nor that an object is below the horizon
 * at its minimum altitude.
 *
 * @param {string} body
 *      The name of a celestial body other than the Earth.
 *
 * @param {Observer} observer
 *      Specifies the geographic coordinates and elevation above sea level of the observer.
 *
 * @param {number} hourAngle
 *      The hour angle expressed in
 *      <a href="https://en.wikipedia.org/wiki/Sidereal_time">sidereal</a>
 *      hours for which the caller seeks to find the body attain.
 *      The value must be in the range [0, 24).
 *      The hour angle represents the number of sidereal hours that have
 *      elapsed since the most recent time the body crossed the observer's local
 *      <a href="https://en.wikipedia.org/wiki/Meridian_(astronomy)">meridian</a>.
 *      This specifying `hourAngle` = 0 finds the moment in time
 *      the body reaches the highest angular altitude in a given sidereal day.
 *
 * @param {FlexibleDateTime} dateStart
 *      The date and time after which the desired hour angle crossing event
 *      is to be found.
 *
 * @returns {HourAngleEvent}
 */
export declare function SearchHourAngle(body: string, observer: Observer, hourAngle: number, dateStart: FlexibleDateTime): HourAngleEvent;
/**
 * @brief When the seasons change for a given calendar year.
 *
 * Represents the dates and times of the two solstices
 * and the two equinoxes in a given calendar year.
 * These four events define the changing of the seasons on the Earth.
 *
 * @property {AstroTime} mar_equinox
 *      The date and time of the March equinox in the given calendar year.
 *      This is the moment in March that the plane of the Earth's equator passes
 *      through the center of the Sun; thus the Sun's declination
 *      changes from a negative number to a positive number.
 *      The March equinox defines
 *      the beginning of spring in the northern hemisphere and
 *      the beginning of autumn in the southern hemisphere.
 *
 * @property {AstroTime} jun_solstice
 *      The date and time of the June solstice in the given calendar year.
 *      This is the moment in June that the Sun reaches its most positive
 *      declination value.
 *      At this moment the Earth's north pole is most tilted most toward the Sun.
 *      The June solstice defines
 *      the beginning of summer in the northern hemisphere and
 *      the beginning of winter in the southern hemisphere.
 *
 * @property {AstroTime} sep_equinox
 *      The date and time of the September equinox in the given calendar year.
 *      This is the moment in September that the plane of the Earth's equator passes
 *      through the center of the Sun; thus the Sun's declination
 *      changes from a positive number to a negative number.
 *      The September equinox defines
 *      the beginning of autumn in the northern hemisphere and
 *      the beginning of spring in the southern hemisphere.
 *
 * @property {AstroTime} dec_solstice
 *      The date and time of the December solstice in the given calendar year.
 *      This is the moment in December that the Sun reaches its most negative
 *      declination value.
 *      At this moment the Earth's south pole is tilted most toward the Sun.
 *      The December solstice defines
 *      the beginning of winter in the northern hemisphere and
 *      the beginning of summer in the southern hemisphere.
 */
export declare class SeasonInfo {
    mar_equinox: AstroTime;
    jun_solstice: AstroTime;
    sep_equinox: AstroTime;
    dec_solstice: AstroTime;
    constructor(mar_equinox: AstroTime, jun_solstice: AstroTime, sep_equinox: AstroTime, dec_solstice: AstroTime);
}
/**
 * @brief Finds the equinoxes and solstices for a given calendar year.
 *
 * @param {number | AstroTime} year
 *      The integer value or `AstroTime` object that specifies
 *      the UTC calendar year for which to find equinoxes and solstices.
 *
 * @returns {SeasonInfo}
 */
export declare function Seasons(year: (number | AstroTime)): SeasonInfo;
/**
 * @brief The viewing conditions of a body relative to the Sun.
 *
 * Represents the angular separation of a body from the Sun as seen from the Earth
 * and the relative ecliptic longitudes between that body and the Earth as seen from the Sun.
 *
 * @property {AstroTime} time
 *      The date and time of the observation.
 *
 * @property {string}  visibility
 *      Either `"morning"` or `"evening"`,
 *      indicating when the body is most easily seen.
 *
 * @property {number}  elongation
 *      The angle in degrees, as seen from the center of the Earth,
 *      of the apparent separation between the body and the Sun.
 *      This angle is measured in 3D space and is not projected onto the ecliptic plane.
 *      When `elongation` is less than a few degrees, the body is very
 *      difficult to see from the Earth because it is lost in the Sun's glare.
 *      The elongation is always in the range [0, 180].
 *
 * @property {number}  ecliptic_separation
 *      The absolute value of the difference between the body's ecliptic longitude
 *      and the Sun's ecliptic longitude, both as seen from the center of the Earth.
 *      This angle measures around the plane of the Earth's orbit (the ecliptic),
 *      and ignores how far above or below that plane the body is.
 *      The ecliptic separation is measured in degrees and is always in the range [0, 180].
 *
 * @see {@link Elongation}
 */
export declare class ElongationEvent {
    time: AstroTime;
    visibility: string;
    elongation: number;
    ecliptic_separation: number;
    constructor(time: AstroTime, visibility: string, elongation: number, ecliptic_separation: number);
}
/**
 * @brief Calculates the viewing conditions of a body relative to the Sun.
 *
 * Calculates angular separation of a body from the Sun as seen from the Earth
 * and the relative ecliptic longitudes between that body and the Earth as seen from the Sun.
 * See the return type {@link ElongationEvent} for details.
 *
 * This function is helpful for determining how easy
 * it is to view a planet away from the Sun's glare on a given date.
 * It also determines whether the object is visible in the morning or evening;
 * this is more important the smaller the elongation is.
 * It is also used to determine how far a planet is from opposition, conjunction, or quadrature.
 *
 * @param {string} body
 *      The name of the observed body. Not allowed to be `"Earth"`.
 *
 * @returns {ElongationEvent}
 */
export declare function Elongation(body: string, date: FlexibleDateTime): ElongationEvent;
/**
 * @brief Finds the next time Mercury or Venus reaches maximum elongation.
 *
 * Searches for the next maximum elongation event for Mercury or Venus
 * that occurs after the given start date. Calling with other values
 * of `body` will result in an exception.
 * Maximum elongation occurs when the body has the greatest
 * angular separation from the Sun, as seen from the Earth.
 * Returns an `ElongationEvent` object containing the date and time of the next
 * maximum elongation, the elongation in degrees, and whether
 * the body is visible in the morning or evening.
 *
 * @param {string} body     Either `"Mercury"` or `"Venus"`.
 * @param {FlexibleDateTime} startDate  The date and time after which to search for the next maximum elongation event.
 *
 * @returns {ElongationEvent}
 */
export declare function SearchMaxElongation(body: string, startDate: FlexibleDateTime): ElongationEvent;
/**
 * @brief Searches for the date and time Venus will next appear brightest as seen from the Earth.
 *
 * @param {string} body
 *      Currently only `"Venus"` is supported.
 *      Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from Earth,
 *      so peak magnitude events have little practical value for that planet.
 *      The Moon reaches peak magnitude very close to full moon, which can be found using
 *      {@link SearchMoonQuarter} or {@link SearchMoonPhase}.
 *      The other planets reach peak magnitude very close to opposition,
 *      which can be found using {@link SearchRelativeLongitude}.
 *
 * @param {FlexibleDateTime} startDate
 *      The date and time after which to find the next peak magnitude event.
 *
 * @returns {IlluminationInfo}
 */
export declare function SearchPeakMagnitude(body: string, startDate: FlexibleDateTime): IlluminationInfo;
/**
 * @brief A closest or farthest point in a body's orbit around its primary.
 *
 * For a planet orbiting the Sun, apsis is a perihelion or aphelion, respectively.
 * For the Moon orbiting the Earth, apsis is a perigee or apogee, respectively.
 *
 * @property {AstroTime} time
 *      The date and time of the apsis.
 *
 * @property {number} kind
 *      For a closest approach (perigee or perihelion), `kind` is 0.
 *      For a farthest distance event (apogee or aphelion), `kind` is 1.
 *
 * @property {number} dist_au
 *      The distance between the centers of the two bodies in astronomical units (AU).
 *
 * @property {number} dist_km
 *      The distance between the centers of the two bodies in kilometers.
 *
 * @see {@link SearchLunarApsis}
 * @see {@link NextLunarApsis}
 */
export declare class Apsis {
    time: AstroTime;
    kind: number;
    dist_au: number;
    dist_km: number;
    constructor(time: AstroTime, kind: number, dist_au: number);
}
/**
 * @brief Finds the next perigee or apogee of the Moon.
 *
 * Finds the next perigee (closest approach) or apogee (farthest remove) of the Moon
 * that occurs after the specified date and time.
 *
 * @param {FlexibleDateTime} startDate
 *      The date and time after which to find the next perigee or apogee.
 *
 * @returns {Apsis}
 */
export declare function SearchLunarApsis(startDate: FlexibleDateTime): Apsis;
/**
 * @brief Finds the next lunar apsis (perigee or apogee) in a series.
 *
 * Given a lunar apsis returned by an initial call to {@link SearchLunarApsis},
 * or a previous call to `NextLunarApsis`, finds the next lunar apsis.
 * If the given apsis is a perigee, this function finds the next apogee, and vice versa.
 *
 * @param {Apsis} apsis
 *      A lunar perigee or apogee event.
 *
 * @returns {Apsis}
 *      The successor apogee for the given perigee, or the successor perigee for the given apogee.
 */
export declare function NextLunarApsis(apsis: Apsis): Apsis;
/**
 * @brief Finds the next perihelion or aphelion of a planet.
 *
 * Finds the date and time of a planet's perihelion (closest approach to the Sun)
 * or aphelion (farthest distance from the Sun) after a given time.
 *
 * Given a date and time to start the search in `startTime`, this function finds the
 * next date and time that the center of the specified planet reaches the closest or farthest point
 * in its orbit with respect to the center of the Sun, whichever comes first
 * after `startTime`.
 *
 * The closest point is called *perihelion* and the farthest point is called *aphelion*.
 * The word *apsis* refers to either event.
 *
 * To iterate through consecutive alternating perihelion and aphelion events,
 * call `SearchPlanetApsis` once, then use the return value to call
 * {@link NextPlanetApsis}. After that, keep feeding the previous return value
 * from `NextPlanetApsis` into another call of `NextPlanetApsis`
 * as many times as desired.
 *
 * @param {string} body
 *      The planet for which to find the next perihelion/aphelion event.
 *      Not allowed to be `"Sun"` or `"Moon"`.
 *
 * @param {AstroTime} startTime
 *      The date and time at which to start searching for the next perihelion or aphelion.
 *
 * @returns {Apsis}
 *      The next perihelion or aphelion that occurs after `startTime`.
 */
export declare function SearchPlanetApsis(body: string, startTime: AstroTime): Apsis;
/**
 * @brief Finds the next planetary perihelion or aphelion event in a series.
 *
 * This function requires an {@link Apsis} value obtained from a call
 * to {@link SearchPlanetApsis} or `NextPlanetApsis`.
 * Given an aphelion event, this function finds the next perihelion event, and vice versa.
 * See {@link SearchPlanetApsis} for more details.
 *
 * @param {string} body
 *      The planet for which to find the next perihelion/aphelion event.
 *      Not allowed to be `"Sun"` or `"Moon"`.
 *      Must match the body passed into the call that produced the `apsis` parameter.
 *
 * @param {Apsis} apsis
 *      An apsis event obtained from a call to {@link SearchPlanetApsis} or `NextPlanetApsis`.
 *
 * @returns {Apsis}
 *      Same as the return value for {@link SearchPlanetApsis}.
 */
export declare function NextPlanetApsis(body: string, apsis: Apsis): Apsis;
/**
 * @brief Calculates the inverse of a rotation matrix.
 *
 * Given a rotation matrix that performs some coordinate transform,
 * this function returns the matrix that reverses that trasnform.
 *
 * @param {RotationMatrix} rotation
 *      The rotation matrix to be inverted.
 *
 * @returns {RotationMatrix}
 *      The inverse rotation matrix.
 */
export declare function InverseRotation(rotation: RotationMatrix): RotationMatrix;
/**
 * @brief Creates a rotation based on applying one rotation followed by another.
 *
 * Given two rotation matrices, returns a combined rotation matrix that is
 * equivalent to rotating based on the first matrix, followed by the second.
 *
 * @param {RotationMatrix} a
 *      The first rotation to apply.
 *
 * @param {RotationMatrix} b
 *      The second rotation to apply.
 *
 * @returns {RotationMatrix}
 *      The combined rotation matrix.
 */
export declare function CombineRotation(a: RotationMatrix, b: RotationMatrix): RotationMatrix;
/**
 * @brief Converts spherical coordinates to Cartesian coordinates.
 *
 * Given spherical coordinates and a time at which they are valid,
 * returns a vector of Cartesian coordinates. The returned value
 * includes the time, as required by `AstroTime`.
 *
 * @param {Spherical} sphere
 *      Spherical coordinates to be converted.
 *
 * @param {AstroTime} time
 *      The time that should be included in the returned vector.
 *
 * @returns {Vector}
 *      The vector form of the supplied spherical coordinates.
 */
export declare function VectorFromSphere(sphere: Spherical, time: AstroTime): Vector;
/**
 * @brief Given angular equatorial coordinates, calculates the equatorial vector.
 *
 * @param {EquatorialCoordinates} equ
 *      An object that contains angular equatorial coordinates to be converted to a vector.
 *
 * @param {AstroTime} time
 *      The date and time of the observation. This is needed because the returned
 *      vector object requires a valid time value when passed to certain other functions.
 *
 * @returns {Vector}
 *      A vector in the equatorial system.
 */
export declare function VectorFromEquator(equ: EquatorialCoordinates, time: AstroTime): Vector;
/**
 * @brief Given an equatorial vector, calculates equatorial angular coordinates.
 *
 * @param {Vector} vec
 *      A vector in an equatorial coordinate system.
 *
 * @returns {EquatorialCoordinates}
 *      Angular coordinates expressed in the same equatorial system as `vec`.
 */
export declare function EquatorFromVector(vec: Vector): EquatorialCoordinates;
/**
 * @brief Converts Cartesian coordinates to spherical coordinates.
 *
 * Given a Cartesian vector, returns latitude, longitude, and distance.
 *
 * @param {Vector} vector
 *      Cartesian vector to be converted to spherical coordinates.
 *
 * @returns {Spherical}
 *      Spherical coordinates that are equivalent to the given vector.
 */
export declare function SphereFromVector(vector: Vector): Spherical;
/**
 * @brief Converts Cartesian coordinates to horizontal coordinates.
 *
 * Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.
 *
 * *IMPORTANT:* This function differs from {@link SphereFromVector} in two ways:
 * - `SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise
 *   from north (e.g., west = +90), but this function represents a clockwise rotation
 *   (e.g., east = +90). The difference is because `SphereFromVector` is intended
 *   to preserve the vector "right-hand rule", while this function defines azimuth in a more
 *   traditional way as used in navigation and cartography.
 * - This function optionally corrects for atmospheric refraction, while `SphereFromVector` does not.
 *
 * The returned object contains the azimuth in `lon`.
 * It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.
 *
 * The altitude is stored in `lat`.
 *
 * The distance to the observed object is stored in `dist`,
 * and is expressed in astronomical units (AU).
 *
 * @param {Vector} vector
 *      Cartesian vector to be converted to horizontal coordinates.
 *
 * @param {string} refraction
 *      `"normal"`: correct altitude for atmospheric refraction (recommended).
 *      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.
 *      `null`: no atmospheric refraction correction is performed.
 *
 * @returns {Spherical}
 */
export declare function HorizonFromVector(vector: Vector, refraction: string): Spherical;
/**
 * @brief Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.
 *
 * @param {Spherical} sphere
 *      A structure that contains apparent horizontal coordinates:
 *      `lat` holds the refracted azimuth angle,
 *      `lon` holds the azimuth in degrees clockwise from north,
 *      and `dist` holds the distance from the observer to the object in AU.
 *
 * @param {AstroTime} time
 *      The date and time of the observation. This is needed because the returned
 *      vector object requires a valid time value when passed to certain other functions.
 *
 * @param {string} refraction
 *      `"normal"`: correct altitude for atmospheric refraction (recommended).
 *      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.
 *      `null`: no atmospheric refraction correction is performed.
 *
 * @returns {Vector}
 *      A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).
 */
export declare function VectorFromHorizon(sphere: Spherical, time: AstroTime, refraction: string): Vector;
/**
 * @brief Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.
 *
 * Given an altitude angle and a refraction option, calculates
 * the amount of "lift" caused by atmospheric refraction.
 * This is the number of degrees higher in the sky an object appears
 * due to the lensing of the Earth's atmosphere.
 *
 * @param {string} refraction
 *      `"normal"`: correct altitude for atmospheric refraction (recommended).
 *      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.
 *      `null`: no atmospheric refraction correction is performed.
 *
 * @param {number} altitude
 *      An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90.
 *
 * @returns {number}
 *      The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing.
 */
export declare function Refraction(refraction: string, altitude: number): number;
/**
 * @brief Calculates the inverse of an atmospheric refraction angle.
 *
 * Given an observed altitude angle that includes atmospheric refraction,
 * calculate the negative angular correction to obtain the unrefracted
 * altitude. This is useful for cases where observed horizontal
 * coordinates are to be converted to another orientation system,
 * but refraction first must be removed from the observed position.
 *
 * @param {string} refraction
 *      `"normal"`: correct altitude for atmospheric refraction (recommended).
 *      `"jplhor"`: for JPL Horizons compatibility testing only; not recommended for normal use.
 *      `null`: no atmospheric refraction correction is performed.
 *
 * @param {number} bent_altitude
 *      The apparent altitude that includes atmospheric refraction.
 *
 * @returns {number}
 *      The angular adjustment in degrees to be added to the
 *      altitude angle to correct for atmospheric lensing.
 *      This will be less than or equal to zero.
 */
export declare function InverseRefraction(refraction: string, bent_altitude: number): number;
/**
 * @brief Applies a rotation to a vector, yielding a rotated vector.
 *
 * This function transforms a vector in one orientation to a vector
 * in another orientation.
 *
 * @param {RotationMatrix} rotation
 *      A rotation matrix that specifies how the orientation of the vector is to be changed.
 *
 * @param {Vector} vector
 *      The vector whose orientation is to be changed.
 *
 * @returns {Vector}
 *      A vector in the orientation specified by `rotation`.
 */
export declare function RotateVector(rotation: RotationMatrix, vector: Vector): Vector;
/**
 * @brief Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using equator at J2000 epoch.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQJ to ECL.
 */
export declare function Rotation_EQJ_ECL(): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: EQJ = equatorial system, using equator at J2000 epoch.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts ECL to EQJ.
 */
export declare function Rotation_ECL_EQJ(): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using equator at J2000 epoch.
 * Target: EQD = equatorial system, using equator of the specified date/time.
 *
 * @param {AstroTime} time
 *      The date and time at which the Earth's equator defines the target orientation.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQJ to EQD at `time`.
 */
export declare function Rotation_EQJ_EQD(time: AstroTime): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of the specified date/time.
 * Target: EQJ = equatorial system, using equator at J2000 epoch.
 *
 * @param {AstroTime} time
 *      The date and time at which the Earth's equator defines the source orientation.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQD at `time` to EQJ.
 */
export declare function Rotation_EQD_EQJ(time: AstroTime): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of the specified date/time.
 * Target: HOR = horizontal system.
 *
 * Use `HorizonFromVector` to convert the return value
 * to a traditional altitude/azimuth pair.
 *
 * @param {AstroTime} time
 *      The date and time at which the Earth's equator applies.
 *
 * @param {Observer} observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQD to HOR at `time` and for `observer`.
 *      The components of the horizontal vector are:
 *      x = north, y = west, z = zenith (straight up from the observer).
 *      These components are chosen so that the "right-hand rule" works for the vector
 *      and so that north represents the direction where azimuth = 0.
 */
export declare function Rotation_EQD_HOR(time: AstroTime, observer: Observer): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system (x=North, y=West, z=Zenith).
 * Target: EQD = equatorial system, using equator of the specified date/time.
 *
 * @param {AstroTime} time
 *      The date and time at which the Earth's equator applies.
 *
 * @param {Observer} observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts HOR to EQD at `time` and for `observer`.
 */
export declare function Rotation_HOR_EQD(time: AstroTime, observer: Observer): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system (x=North, y=West, z=Zenith).
 * Target: EQJ = equatorial system, using equator at the J2000 epoch.
 *
 * @param {AstroTime} time
 *      The date and time of the observation.
 *
 * @param {Observer} observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts HOR to EQD at `time` and for `observer`.
 */
export declare function Rotation_HOR_EQJ(time: AstroTime, observer: Observer): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQJ = equatorial system, using the equator at the J2000 epoch.
 * Target: HOR = horizontal system.
 *
 * Use {@link HorizonFromVector} to convert the return value
 * to a traditional altitude/azimuth pair.
 *
 * @param time
 *      The date and time of the desired horizontal orientation.
 *
 * @param observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @return
 *      A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
 *      The components of the horizontal vector are:
 *      x = north, y = west, z = zenith (straight up from the observer).
 *      These components are chosen so that the "right-hand rule" works for the vector
 *      and so that north represents the direction where azimuth = 0.
 */
export declare function Rotation_EQJ_HOR(time: AstroTime, observer: Observer): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: EQD = equatorial system, using equator of date.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @param {AstroTime} time
 *      The date and time of the source equator.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts EQD to ECL.
 */
export declare function Rotation_EQD_ECL(time: AstroTime): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: EQD = equatorial system, using equator of date.
 *
 * @param {AstroTime} time
 *      The date and time of the desired equator.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts ECL to EQD.
 */
export declare function Rotation_ECL_EQD(time: AstroTime): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: ECL = ecliptic system, using equator at J2000 epoch.
 * Target: HOR = horizontal system.
 *
 * Use {@link HorizonFromVector} to convert the return value
 * to a traditional altitude/azimuth pair.
 *
 * @param {AstroTime} time
 *      The date and time of the desired horizontal orientation.
 *
 * @param {Observer} observer
 *      A location near the Earth's mean sea level that defines the observer's horizon.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts ECL to HOR at `time` and for `observer`.
 *      The components of the horizontal vector are:
 *      x = north, y = west, z = zenith (straight up from the observer).
 *      These components are chosen so that the "right-hand rule" works for the vector
 *      and so that north represents the direction where azimuth = 0.
 */
export declare function Rotation_ECL_HOR(time: AstroTime, observer: Observer): RotationMatrix;
/**
 * @brief Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).
 *
 * This is one of the family of functions that returns a rotation matrix
 * for converting from one orientation to another.
 * Source: HOR = horizontal system.
 * Target: ECL = ecliptic system, using equator at J2000 epoch.
 *
 * @param {AstroTime} time
 *      The date and time of the horizontal observation.
 *
 * @param {Observer} observer
 *      The location of the horizontal observer.
 *
 * @returns {RotationMatrix}
 *      A rotation matrix that converts HOR to ECL.
 */
export declare function Rotation_HOR_ECL(time: AstroTime, observer: Observer): RotationMatrix;
/**
 * @brief Reports the constellation that a given celestial point lies within.
 *
 * @property {string} symbol
 *      3-character mnemonic symbol for the constellation, e.g. "Ori".
 *
 * @property {string} name
 *      Full name of constellation, e.g. "Orion".
 *
 * @property {number} ra1875
 *      Right ascension expressed in B1875 coordinates.
 *
 * @property {number} dec1875
 *      Declination expressed in B1875 coordinates.
 */
export declare class ConstellationInfo {
    symbol: string;
    name: string;
    ra1875: number;
    dec1875: number;
    constructor(symbol: string, name: string, ra1875: number, dec1875: number);
}
/**
 * @brief Determines the constellation that contains the given point in the sky.
 *
 * Given J2000 equatorial (EQJ) coordinates of a point in the sky,
 * determines the constellation that contains that point.
 *
 * @param {number} ra
 *      The right ascension (RA) of a point in the sky, using the J2000 equatorial system.
 *
 * @param {number} dec
 *      The declination (DEC) of a point in the sky, using the J2000 equatorial system.
 *
 * @returns {ConstellationInfo}
 *      An object that contains the 3-letter abbreviation and full name
 *      of the constellation that contains the given (ra,dec), along with
 *      the converted B1875 (ra,dec) for that point.
 */
export declare function Constellation(ra: number, dec: number): ConstellationInfo;
/**
 * @brief Returns information about a lunar eclipse.
 *
 * Returned by {@link SearchLunarEclipse} or {@link NextLunarEclipse}
 * to report information about a lunar eclipse event.
 * When a lunar eclipse is found, it is classified as penumbral, partial, or total.
 * Penumbral eclipses are difficult to observe, because the moon is only slightly dimmed
 * by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
 * Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
 * Total eclipses occur when the entire Moon passes into the Earth's umbra.
 *
 * The `kind` field thus holds one of the strings `"penumbral"`, `"partial"`,
 * or `"total"`, depending on the kind of lunar eclipse found.
 *
 * Field `peak` holds the date and time of the peak of the eclipse, when it is at its peak.
 *
 * Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
 * of the eclipse, which is half of the amount of time the eclipse spends in each
 * phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
 * By converting from minutes to days, and subtracting/adding with `peak`, the caller
 * may determine the date and time of the beginning/end of each eclipse phase.
 *
 * @property {string} kind
 *      The type of lunar eclipse found.
 *
 * @property {AstroTime} peak
 *      The time of the eclipse at its peak.
 *
 * @property {number} sd_penum
 *      The semi-duration of the penumbral phase in minutes.
 *
 * @property {number} sd_partial
 *      The semi-duration of the penumbral phase in minutes, or 0.0 if none.
 *
 * @property {number} sd_total
 *      The semi-duration of the penumbral phase in minutes, or 0.0 if none.
 *
 */
export declare class LunarEclipseInfo {
    kind: string;
    peak: AstroTime;
    sd_penum: number;
    sd_partial: number;
    sd_total: number;
    constructor(kind: string, peak: AstroTime, sd_penum: number, sd_partial: number, sd_total: number);
}
/**
 * @brief Searches for a lunar eclipse.
 *
 * This function finds the first lunar eclipse that occurs after `startTime`.
 * A lunar eclipse may be penumbral, partial, or total.
 * See {@link LunarEclipseInfo} for more information.
 * To find a series of lunar eclipses, call this function once,
 * then keep calling {@link NextLunarEclipse} as many times as desired,
 * passing in the `center` value returned from the previous call.
 *
 * @param {FlexibleDateTime} date
 *      The date and time for starting the search for a lunar eclipse.
 *
 * @returns {LunarEclipseInfo}
 */
export declare function SearchLunarEclipse(date: FlexibleDateTime): LunarEclipseInfo;
/**
    @brief Reports the time and geographic location of the peak of a solar eclipse.

    Returned by {@link SearchGlobalSolarEclipse} or {@link NextGlobalSolarEclipse}
    to report information about a solar eclipse event.

    Field `peak` holds the date and time of the peak of the eclipse, defined as
    the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.

    The eclipse is classified as partial, annular, or total, depending on the
    maximum amount of the Sun's disc obscured, as seen at the peak location
    on the surface of the Earth.

    The `kind` field thus holds one of the strings `"partial"`, `"annular"`, or `"total"`.
    A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
    An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
    to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
    A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
    observes either a total or annular eclipse.

    If `kind` is `"total"` or `"annular"`, the `latitude` and `longitude`
    fields give the geographic coordinates of the center of the Moon's shadow projected
    onto the daytime side of the Earth at the instant of the eclipse's peak.
    If `kind` has any other value, `latitude` and `longitude` are undefined and should
    not be used.

    @property {string} kind
        One of the following string values: `"partial"`, `"annular"`, `"total"`.

    @property {AstroTime} peak
        The date and time of the peak of the eclipse, defined as the instant
        when the axis of the Moon's shadow cone passes closest to the Earth's center.

    @property {number} distance
        The distance in kilometers between the axis of the Moon's shadow cone
        and the center of the Earth at the time indicated by `peak`.

    @property {number | undefined} latitude
        If `kind` holds `"total"`, the geographic latitude in degrees
        where the center of the Moon's shadow falls on the Earth at the
        time indicated by `peak`; otherwise, `latitude` holds `undefined`.

    @property {number | undefined} longitude
        If `kind` holds `"total"`, the geographic longitude in degrees
        where the center of the Moon's shadow falls on the Earth at the
        time indicated by `peak`; otherwise, `longitude` holds `undefined`.
*/
export declare class GlobalSolarEclipseInfo {
    kind: string;
    peak: AstroTime;
    distance: number;
    latitude?: number | undefined;
    longitude?: number | undefined;
    constructor(kind: string, peak: AstroTime, distance: number, latitude?: number | undefined, longitude?: number | undefined);
}
/**
 * @brief Searches for the next lunar eclipse in a series.
 *
 * After using {@link SearchLunarEclipse} to find the first lunar eclipse
 * in a series, you can call this function to find the next consecutive lunar eclipse.
 * Pass in the `center` value from the {@link LunarEclipseInfo} returned by the
 * previous call to `SearchLunarEclipse` or `NextLunarEclipse`
 * to find the next lunar eclipse.
 *
 * @param {AstroTime} prevEclipseTime
 *      A date and time near a full moon. Lunar eclipse search will start at the next full moon.
 *
 * @returns {LunarEclipseInfo}
 */
export declare function NextLunarEclipse(prevEclipseTime: AstroTime): LunarEclipseInfo;
/**
 * @brief Searches for a solar eclipse visible anywhere on the Earth's surface.
 *
 * This function finds the first solar eclipse that occurs after `startTime`.
 * A solar eclipse may be partial, annular, or total.
 * See {@link GlobalSolarEclipseInfo} for more information.
 * To find a series of solar eclipses, call this function once,
 * then keep calling {@link NextGlobalSolarEclipse} as many times as desired,
 * passing in the `peak` value returned from the previous call.
 *
 * @param {AstroTime} startTime
 *      The date and time for starting the search for a solar eclipse.
 *
 * @returns {GlobalSolarEclipseInfo}
 */
export declare function SearchGlobalSolarEclipse(startTime: AstroTime): GlobalSolarEclipseInfo;
/**
 * @brief Searches for the next global solar eclipse in a series.
 *
 * After using {@link SearchGlobalSolarEclipse} to find the first solar eclipse
 * in a series, you can call this function to find the next consecutive solar eclipse.
 * Pass in the `peak` value from the {@link GlobalSolarEclipseInfo} returned by the
 * previous call to `SearchGlobalSolarEclipse` or `NextGlobalSolarEclipse`
 * to find the next solar eclipse.
 *
 * @param {AstroTime} prevEclipseTime
 *      A date and time near a new moon. Solar eclipse search will start at the next new moon.
 *
 * @returns {GlobalSolarEclipseInfo}
 */
export declare function NextGlobalSolarEclipse(prevEclipseTime: AstroTime): GlobalSolarEclipseInfo;
/**
 * @brief Holds a time and the observed altitude of the Sun at that time.
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
 *
 * @property {AstroTime} time
 *      The date and time of the event.
 *
 * @property {number} altitude
 *      The angular altitude of the center of the Sun above/below the horizon, at `time`,
 *      corrected for atmospheric refraction and expressed in degrees.
 */
export declare class EclipseEvent {
    time: AstroTime;
    altitude: number;
    constructor(time: AstroTime, altitude: number);
}
/**
 * @brief Information about a solar eclipse as seen by an observer at a given time and geographic location.
 *
 * Returned by {@link SearchLocalSolarEclipse} or {@link NextLocalSolarEclipse}
 * to report information about a solar eclipse as seen at a given geographic location.
 *
 * When a solar eclipse is found, it is classified by setting `kind`
 * to `"partial"`, `"annular"`, or `"total"`.
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
 * The fields `partial_begin` and `partial_end` are always set, and indicate when
 * the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
 * `total_begin` and `total_end` indicate when the total/annular phase begins/ends.
 * When an event field is valid, the caller must also check its `altitude` field to
 * see whether the Sun is above the horizon at the time indicated by the `time` field.
 * See {@link EclipseEvent} for more information.
 *
 * @property {string} kind
 *      The type of solar eclipse found: `"partial"`, `"annular"`, or `"total"`.
 *
 * @property {EclipseEvent} partial_begin
 *      The time and Sun altitude at the beginning of the eclipse.
 *
 * @property {EclipseEvent | undefined} total_begin
 *      If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise undefined.
 *
 * @property {EclipseEvent} peak
 *      The time and Sun altitude when the eclipse reaches its peak.
 *
 * @property {EclipseEvent | undefined} total_end
 *      If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise undefined.
 *
 * @property {EclipseEvent} partial_end
 *      The time and Sun altitude at the end of the eclipse.
 */
export declare class LocalSolarEclipseInfo {
    kind: string;
    partial_begin: EclipseEvent;
    total_begin: EclipseEvent | undefined;
    peak: EclipseEvent;
    total_end: EclipseEvent | undefined;
    partial_end: EclipseEvent;
    constructor(kind: string, partial_begin: EclipseEvent, total_begin: EclipseEvent | undefined, peak: EclipseEvent, total_end: EclipseEvent | undefined, partial_end: EclipseEvent);
}
/**
 * @brief Searches for a solar eclipse visible at a specific location on the Earth's surface.
 *
 * This function finds the first solar eclipse that occurs after `startTime`.
 * A solar eclipse may be partial, annular, or total.
 * See {@link LocalSolarEclipseInfo} for more information.
 *
 * To find a series of solar eclipses, call this function once,
 * then keep calling {@link NextLocalSolarEclipse} as many times as desired,
 * passing in the `peak` value returned from the previous call.
 *
 * IMPORTANT: An eclipse reported by this function might be partly or
 * completely invisible to the observer due to the time of day.
 * See {@link LocalSolarEclipseInfo} for more information about this topic.
 *
 * @param {AstroTime} startTime
 *      The date and time for starting the search for a solar eclipse.
 *
 * @param {Observer} observer
 *      The geographic location of the observer.
 *
 * @returns {LocalSolarEclipseInfo}
 */
export declare function SearchLocalSolarEclipse(startTime: AstroTime, observer: Observer): LocalSolarEclipseInfo;
/**
 * @brief Searches for the next local solar eclipse in a series.
 *
 * After using {@link SearchLocalSolarEclipse} to find the first solar eclipse
 * in a series, you can call this function to find the next consecutive solar eclipse.
 * Pass in the `peak` value from the {@link LocalSolarEclipseInfo} returned by the
 * previous call to `SearchLocalSolarEclipse` or `NextLocalSolarEclipse`
 * to find the next solar eclipse.
 * This function finds the first solar eclipse that occurs after `startTime`.
 * A solar eclipse may be partial, annular, or total.
 * See {@link LocalSolarEclipseInfo} for more information.
 *
 * @param {AstroTime} prevEclipseTime
 *      The date and time for starting the search for a solar eclipse.
 *
 * @param {Observer} observer
 *      The geographic location of the observer.
 *
 * @returns {LocalSolarEclipseInfo}
 */
export declare function NextLocalSolarEclipse(prevEclipseTime: AstroTime, observer: Observer): LocalSolarEclipseInfo;
/**
 * @brief Information about a transit of Mercury or Venus, as seen from the Earth.
 *
 * Returned by {@link SearchTransit} or {@link NextTransit} to report
 * information about a transit of Mercury or Venus.
 * A transit is when Mercury or Venus passes between the Sun and Earth so that
 * the other planet is seen in silhouette against the Sun.
 *
 * The calculations are performed from the point of view of a geocentric observer.
 *
 * @property {AstroTime} start
 *      The date and time at the beginning of the transit.
 *      This is the moment the planet first becomes visible against the Sun in its background.
 *
 * @property {AstroTime} peak
 *      When the planet is most aligned with the Sun, as seen from the Earth.
 *
 * @property {AstroTime} finish
 *      The date and time at the end of the transit.
 *      This is the moment the planet is last seen against the Sun in its background.
 *
 * @property {number} separation
 *      The minimum angular separation, in arcminutes, between the centers of the Sun and the planet.
 *      This angle pertains to the time stored in `peak`.
 */
export declare class TransitInfo {
    start: AstroTime;
    peak: AstroTime;
    finish: AstroTime;
    separation: number;
    constructor(start: AstroTime, peak: AstroTime, finish: AstroTime, separation: number);
}
/**
 * @brief Searches for the first transit of Mercury or Venus after a given date.
 *
 * Finds the first transit of Mercury or Venus after a specified date.
 * A transit is when an inferior planet passes between the Sun and the Earth
 * so that the silhouette of the planet is visible against the Sun in the background.
 * To continue the search, pass the `finish` time in the returned structure to
 * {@link NextTransit}.
 *
 * @param {string} body
 *      The planet whose transit is to be found. Must be `"Mercury"` or `"Venus"`.
 *
 * @param {AstroTime} startTime
 *      The date and time for starting the search for a transit.
 *
 * @returns {TransitInfo}
 */
export declare function SearchTransit(body: string, startTime: AstroTime): TransitInfo;
/**
 * @brief Searches for the next transit of Mercury or Venus in a series.
 *
 * After calling {@link SearchTransit} to find a transit of Mercury or Venus,
 * this function finds the next transit after that.
 * Keep calling this function as many times as you want to keep finding more transits.
 *
 * @param {string} body
 *      The planet whose transit is to be found. Must be `"Mercury"` or `"Venus"`.
 *
 * @param {AstroTime} prevTransitTime
 *      A date and time near the previous transit.
 *
 * @returns {TransitInfo}
 */
export declare function NextTransit(body: string, prevTransitTime: AstroTime): TransitInfo;
