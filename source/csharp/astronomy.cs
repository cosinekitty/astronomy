/*
    Astronomy Engine for C# / .NET.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019-2020 Don Cross <cosinekitty@gmail.com>

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

using System;

namespace CosineKitty
{
    /// <summary>
    /// This exception is thrown by certain Astronomy Engine functions
    /// when an invalid attempt is made to use the Earth as the observed
    /// celestial body. Usually this happens for cases where the Earth itself
    /// is the location of the observer.
    /// </summary>
    public class EarthNotAllowedException: ArgumentException
    {
        /// <summary>Creates an exception indicating that the Earth is not allowed as a target body.</summary>
        public EarthNotAllowedException():
            base("The Earth is not allowed as the body parameter.")
            {}
    }

    /// <summary>
    /// This exception is thrown by certain Astronomy Engine functions
    /// when a body is specified that is not appropriate for the given operation.
    /// </summary>
    public class InvalidBodyException: ArgumentException
    {
        /// <summary>Creates an exception indicating that the given body is not valid for this operation.</summary>
        public InvalidBodyException(Body body):
            base(string.Format("Invalid body: {0}", body))
            {}
    }

    /// <summary>
    /// Pluto supports calculations only within the year range Astronomy.MinYear to Astronomy.MaxYear.
    /// </summary>
    public class BadTimeException : ArgumentException
    {
        /// <summary>Creates an exception indicating Astronomy Engine cannot calculate for this time.</summary>
        public BadTimeException(AstroTime time):
            base(string.Format("Pluto supports calculations only within the year range {0}..{1}. Given time {2} is out of bounds.",
                Astronomy.MinYear, Astronomy.MaxYear, time))
            {}
    }

    /// <summary>Defines a function type for calculating Delta T.</summary>
    /// <remarks>
    /// Delta T is the discrepancy between times measured using an atomic clock
    /// and times based on observations of the Earth's rotation, which is gradually
    /// slowing down over time. Delta T = TT - UT, where
    /// TT = Terrestrial Time, based on atomic time, and
    /// UT = Universal Time, civil time based on the Earth's rotation.
    /// Astronomy Engine defaults to using a Delta T function defined by
    /// Espenak and Meeus in their "Five Millennium Canon of Solar Eclipses".
    /// See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    /// </remarks>
    public delegate double DeltaTimeFunc(double ut);

    /// <summary>
    /// The enumeration of celestial bodies supported by Astronomy Engine.
    /// </summary>
    public enum Body
    {
        /// <summary>
        /// A placeholder value representing an invalid or unknown celestial body.
        /// </summary>
        Invalid = -1,

        /// <summary>
        /// The planet Mercury.
        /// </summary>
        Mercury,

        /// <summary>
        /// The planet Venus.
        /// </summary>
        Venus,

        /// <summary>
        /// The planet Earth.
        /// Some functions that accept a `Body` parameter will fail if passed this value
        /// because they assume that an observation is being made from the Earth,
        /// and therefore the Earth is not a target of observation.
        /// </summary>
        Earth,

        /// <summary>
        /// The planet Mars.
        /// </summary>
        Mars,

        /// <summary>
        /// The planet Jupiter.
        /// </summary>
        Jupiter,

        /// <summary>
        /// The planet Saturn.
        /// </summary>
        Saturn,

        /// <summary>
        /// The planet Uranus.
        /// </summary>
        Uranus,

        /// <summary>
        /// The planet Neptune.
        /// </summary>
        Neptune,

        /// <summary>
        /// The planet Pluto.
        /// </summary>
        Pluto,

        /// <summary>
        /// The Sun.
        /// </summary>
        Sun,

        /// <summary>
        /// The Earth's natural satellite, the Moon.
        /// </summary>
        Moon,

        /// <summary>
        /// The Earth/Moon Barycenter.
        /// </summary>
        EMB,

        /// <summary>
        /// The Solar System Barycenter.
        /// </summary>
        SSB,
    }

    /// <summary>
    /// A date and time used for astronomical calculations.
    /// </summary>
    public class AstroTime
    {
        private static readonly DateTime Origin = new DateTime(2000, 1, 1, 12, 0, 0, DateTimeKind.Utc);

        /// <summary>
        /// UT1/UTC number of days since noon on January 1, 2000.
        /// </summary>
        /// <remarks>
        /// The floating point number of days of Universal Time since noon UTC January 1, 2000.
        /// Astronomy Engine approximates UTC and UT1 as being the same thing, although they are
        /// not exactly equivalent; UTC and UT1 can disagree by up to plus or minus 0.9 seconds.
        /// This approximation is sufficient for the accuracy requirements of Astronomy Engine.
        ///
        /// Universal Time Coordinate (UTC) is the international standard for legal and civil
        /// timekeeping and replaces the older Greenwich Mean Time (GMT) standard.
        /// UTC is kept in sync with unpredictable observed changes in the Earth's rotation
        /// by occasionally adding leap seconds as needed.
        ///
        /// UT1 is an idealized time scale based on observed rotation of the Earth, which
        /// gradually slows down in an unpredictable way over time, due to tidal drag by the Moon and Sun,
        /// large scale weather events like hurricanes, and internal seismic and convection effects.
        /// Conceptually, UT1 drifts from atomic time continuously and erratically, whereas UTC
        /// is adjusted by a scheduled whole number of leap seconds as needed.
        ///
        /// The value in `ut` is appropriate for any calculation involving the Earth's rotation,
        /// such as calculating rise/set times, culumination, and anything involving apparent
        /// sidereal time.
        ///
        /// Before the era of atomic timekeeping, days based on the Earth's rotation
        /// were often known as *mean solar days*.
        /// </remarks>
        public readonly double ut;

        /// <summary>
        /// Terrestrial Time days since noon on January 1, 2000.
        /// </summary>
        /// <remarks>
        /// Terrestrial Time is an atomic time scale defined as a number of days since noon on January 1, 2000.
        /// In this system, days are not based on Earth rotations, but instead by
        /// the number of elapsed [SI seconds](https://physics.nist.gov/cuu/Units/second.html)
        /// divided by 86400. Unlike `ut`, `tt` increases uniformly without adjustments
        /// for changes in the Earth's rotation.
        ///
        /// The value in `tt` is used for calculations of movements not involving the Earth's rotation,
        /// such as the orbits of planets around the Sun, or the Moon around the Earth.
        ///
        /// Historically, Terrestrial Time has also been known by the term *Ephemeris Time* (ET).
        /// </remarks>
        public readonly double tt;

        internal double psi;    // For internal use only. Used to optimize Earth tilt calculations.
        internal double eps;    // For internal use only. Used to optimize Earth tilt calculations.

        /// <summary>
        /// Creates an `AstroTime` object from a Universal Time day value.
        /// </summary>
        /// <param name="ut">The number of days after the J2000 epoch.</param>
        public AstroTime(double ut)
        {
            this.ut = ut;
            this.tt = Astronomy.TerrestrialTime(ut);
            this.psi = this.eps = double.NaN;
        }

        /// <summary>
        /// Creates an `AstroTime` object from a .NET `DateTime` object.
        /// </summary>
        /// <param name="d">The date and time to be converted to AstroTime format.</param>
        public AstroTime(DateTime d)
            : this((d.ToUniversalTime() - Origin).TotalDays)
        {
        }

        /// <summary>
        /// Creates an `AstroTime` object from a UTC year, month, day, hour, minute and second.
        /// </summary>
        /// <param name="year">The UTC year value.</param>
        /// <param name="month">The UTC month value 1..12.</param>
        /// <param name="day">The UTC day of the month 1..31.</param>
        /// <param name="hour">The UTC hour value 0..23.</param>
        /// <param name="minute">The UTC minute value 0..59.</param>
        /// <param name="second">The UTC second value 0..59.</param>
        public AstroTime(int year, int month, int day, int hour, int minute, int second)
            : this(new DateTime(year, month, day, hour, minute, second, DateTimeKind.Utc))
        {
        }

        /// <summary>
        /// Converts this object to .NET `DateTime` format.
        /// </summary>
        /// <returns>a UTC `DateTime` object for this `AstroTime` value.</returns>
        public DateTime ToUtcDateTime()
        {
            return Origin.AddDays(ut).ToUniversalTime();
        }

        /// <summary>
        /// Converts this `AstroTime` to ISO 8601 format, expressed in UTC with millisecond resolution.
        /// </summary>
        /// <returns>Example: "2019-08-30T17:45:22.763".</returns>
        public override string ToString()
        {
            return ToUtcDateTime().ToString("yyyy-MM-ddTHH:mm:ss.fffZ");
        }

        /// <summary>
        /// Calculates the sum or difference of an #AstroTime with a specified floating point number of days.
        /// </summary>
        /// <remarks>
        /// Sometimes we need to adjust a given #AstroTime value by a certain amount of time.
        /// This function adds the given real number of days in `days` to the date and time in this object.
        ///
        /// More precisely, the result's Universal Time field `ut` is exactly adjusted by `days` and
        /// the Terrestrial Time field `tt` is adjusted correctly for the resulting UTC date and time,
        /// according to the historical and predictive Delta-T model provided by the
        /// [United States Naval Observatory](http://maia.usno.navy.mil/ser7/).
        /// </remarks>
        /// <param name="days">A floating point number of days by which to adjust `time`. May be negative, 0, or positive.</param>
        /// <returns>A date and time that is conceptually equal to `time + days`.</returns>
        public AstroTime AddDays(double days)
        {
            return new AstroTime(this.ut + days);
        }
    }

    /// <summary>
    /// A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).
    /// </summary>
    public struct AstroVector
    {
        /// <summary>
        /// The Cartesian x-coordinate of the vector in AU.
        /// </summary>
        public double x;

        /// <summary>
        /// The Cartesian y-coordinate of the vector in AU.
        /// </summary>
        public double y;

        /// <summary>
        /// The Cartesian z-coordinate of the vector in AU.
        /// </summary>
        public double z;

        /// <summary>
        /// The date and time at which this vector is valid.
        /// </summary>
        public AstroTime t;

        /// <summary>
        /// Creates an AstroVector.
        /// </summary>
        /// <param name="x">A Cartesian x-coordinate expressed in AU.</param>
        /// <param name="y">A Cartesian y-coordinate expressed in AU.</param>
        /// <param name="z">A Cartesian z-coordinate expressed in AU.</param>
        /// <param name="t">The date and time at which this vector is valid.</param>
        public AstroVector(double x, double y, double z, AstroTime t)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.t = t;
        }

        /// <summary>
        /// Calculates the total distance in AU represented by this vector.
        /// </summary>
        /// <returns>The nonnegative length of the Cartisian vector in AU.</returns>
        public double Length()
        {
            return Math.Sqrt(x*x + y*y + z*z);
        }
    }

    /// <summary>
    /// Contains a rotation matrix that can be used to transform one coordinate system to another.
    /// </summary>
    public struct RotationMatrix
    {
        /// <summary>A normalized 3x3 rotation matrix.</summary>
        public readonly double[,] rot;

        /// <summary>Creates a rotation matrix.</summary>
        /// <param name="rot">A 3x3 array of floating point numbers defining the rotation matrix.</param>
        public RotationMatrix(double[,] rot)
        {
            if (rot == null || rot.GetLength(0) != 3 || rot.GetLength(1) != 3)
                throw new ArgumentException("Rotation matrix must be given a 3x3 array.");

            this.rot = rot;
        }
    }

    /// <summary>
    /// Spherical coordinates: latitude, longitude, distance.
    /// </summary>
    public struct Spherical
    {
        /// <summary>The latitude angle: -90..+90 degrees.</summary>
        public readonly double lat;

        /// <summary>The longitude angle: 0..360 degrees.</summary>
        public readonly double lon;

        /// <summary>Distance in AU.</summary>
        public readonly double dist;

        /// <summary>
        /// Creates a set of spherical coordinates.
        /// </summary>
        /// <param name="lat">The latitude angle: -90..+90 degrees.</param>
        /// <param name="lon">The longitude angle: 0..360 degrees.</param>
        /// <param name="dist">Distance in AU.</param>
        public Spherical(double lat, double lon, double dist)
        {
            this.lat = lat;
            this.lon = lon;
            this.dist = dist;
        }
    }

    /// <summary>
    /// The location of an observer on (or near) the surface of the Earth.
    /// </summary>
    /// <remarks>
    /// This structure is passed to functions that calculate phenomena as observed
    /// from a particular place on the Earth.
    /// </remarks>
    public struct Observer
    {
        /// <summary>
        /// Geographic latitude in degrees north (positive) or south (negative) of the equator.
        /// </summary>
        public readonly double latitude;

        /// <summary>
        /// Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.
        /// </summary>
        public readonly double longitude;

        /// <summary>
        /// The height above (positive) or below (negative) sea level, expressed in meters.
        /// </summary>
        public readonly double height;

        /// <summary>
        /// Creates an Observer object.
        /// </summary>
        /// <param name="latitude">Geographic latitude in degrees north (positive) or south (negative) of the equator.</param>
        /// <param name="longitude">Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.</param>
        /// <param name="height">The height above (positive) or below (negative) sea level, expressed in meters.</param>
        public Observer(double latitude, double longitude, double height)
        {
            this.latitude = latitude;
            this.longitude = longitude;
            this.height = height;
        }
    }

    /// <summary>
    /// Selects the date for which the Earth's equator is to be used for representing equatorial coordinates.
    /// </summary>
    /// <remarks>
    /// The Earth's equator is not always in the same plane due to precession and nutation.
    ///
    /// Sometimes it is useful to have a fixed plane of reference for equatorial coordinates
    /// across different calendar dates.  In these cases, a fixed *epoch*, or reference time,
    /// is helpful. Astronomy Engine provides the J2000 epoch for such cases.  This refers
    /// to the plane of the Earth's orbit as it was on noon UTC on 1 January 2000.
    ///
    /// For some other purposes, it is more helpful to represent coordinates using the Earth's
    /// equator exactly as it is on that date. For example, when calculating rise/set times
    /// or horizontal coordinates, it is most accurate to use the orientation of the Earth's
    /// equator at that same date and time. For these uses, Astronomy Engine allows *of-date*
    /// calculations.
    /// </remarks>
    public enum EquatorEpoch
    {
        /// <summary>
        /// Represent equatorial coordinates in the J2000 epoch.
        /// </summary>
        J2000,

        /// <summary>
        /// Represent equatorial coordinates using the Earth's equator at the given date and time.
        /// </summary>
        OfDate,
    }

    /// <summary>
    /// Aberration calculation options.
    /// </summary>
    /// <remarks>
    /// [Aberration](https://en.wikipedia.org/wiki/Aberration_of_light) is an effect
    /// causing the apparent direction of an observed body to be shifted due to transverse
    /// movement of the Earth with respect to the rays of light coming from that body.
    /// This angular correction can be anywhere from 0 to about 20 arcseconds,
    /// depending on the position of the observed body relative to the instantaneous
    /// velocity vector of the Earth.
    ///
    /// Some Astronomy Engine functions allow optional correction for aberration by
    /// passing in a value of this enumerated type.
    ///
    /// Aberration correction is useful to improve accuracy of coordinates of
    /// apparent locations of bodies seen from the Earth.
    /// However, because aberration affects not only the observed body (such as a planet)
    /// but the surrounding stars, aberration may be unhelpful (for example)
    /// for determining exactly when a planet crosses from one constellation to another.
    /// </remarks>
    public enum Aberration
    {
        /// <summary>
        /// Request correction for aberration.
        /// </summary>
        Corrected,

        /// <summary>
        /// Do not correct for aberration.
        /// </summary>
        None,
    }

    /// <summary>
    /// Selects whether to correct for atmospheric refraction, and if so, how.
    /// </summary>
    public enum Refraction
    {
        /// <summary>
        /// No atmospheric refraction correction (airless).
        /// </summary>
        None,

        /// <summary>
        /// Recommended correction for standard atmospheric refraction.
        /// </summary>
        Normal,

        /// <summary>
        /// Used only for compatibility testing with JPL Horizons online tool.
        /// </summary>
        JplHor,
    }

    /// <summary>
    /// Selects whether to search for a rising event or a setting event for a celestial body.
    /// </summary>
    public enum Direction
    {
        /// <summary>
        /// Indicates a rising event: a celestial body is observed to rise above the horizon by an observer on the Earth.
        /// </summary>
        Rise = +1,

        /// <summary>
        /// Indicates a setting event: a celestial body is observed to sink below the horizon by an observer on the Earth.
        /// </summary>
        Set = -1,
    }

    /// <summary>
    /// Indicates whether a body (especially Mercury or Venus) is best seen in the morning or evening.
    /// </summary>
    public enum Visibility
    {
        /// <summary>
        /// The body is best visible in the morning, before sunrise.
        /// </summary>
        Morning,

        /// <summary>
        /// The body is best visible in the evening, after sunset.
        /// </summary>
        Evening,
    }

    /// <summary>
    /// Equatorial angular coordinates.
    /// </summary>
    /// <remarks>
    /// Coordinates of a celestial body as seen from the Earth
    /// (geocentric or topocentric, depending on context),
    /// oriented with respect to the projection of the Earth's equator onto the sky.
    /// </remarks>
    public struct Equatorial
    {
        /// <summary>
        /// Right ascension in sidereal hours.
        /// </summary>
        public readonly double ra;

        /// <summary>
        /// Declination in degrees.
        /// </summary>
        public readonly double dec;

        /// <summary>
        /// Distance to the celestial body in AU.
        /// </summary>
        public readonly double dist;

        /// <summary>
        /// Creates an equatorial coordinates object.
        /// </summary>
        /// <param name="ra">Right ascension in sidereal hours.</param>
        /// <param name="dec">Declination in degrees.</param>
        /// <param name="dist">Distance to the celestial body in AU.</param>
        public Equatorial(double ra, double dec, double dist)
        {
            this.ra = ra;
            this.dec = dec;
            this.dist = dist;
        }
    }

    /// <summary>
    /// Ecliptic angular and Cartesian coordinates.
    /// </summary>
    /// <remarks>
    /// Coordinates of a celestial body as seen from the center of the Sun (heliocentric),
    /// oriented with respect to the plane of the Earth's orbit around the Sun (the ecliptic).
    /// </remarks>
    public struct Ecliptic
    {
        /// <summary>
        /// Cartesian x-coordinate: in the direction of the equinox along the ecliptic plane.
        /// </summary>
        public readonly double ex;

        /// <summary>
        /// Cartesian y-coordinate: in the ecliptic plane 90 degrees prograde from the equinox.
        /// </summary>
        public readonly double ey;

        /// <summary>
        /// Cartesian z-coordinate: perpendicular to the ecliptic plane. Positive is north.
        /// </summary>
        public readonly double ez;

        /// <summary>
        /// Latitude in degrees north (positive) or south (negative) of the ecliptic plane.
        /// </summary>
        public readonly double elat;

        /// <summary>
        /// Longitude in degrees around the ecliptic plane prograde from the equinox.
        /// </summary>
        public readonly double elon;

        /// <summary>
        /// Creates an object that holds Cartesian and angular ecliptic coordinates.
        /// </summary>
        /// <param name="ex">x-coordinate of the ecliptic position</param>
        /// <param name="ey">y-coordinate of the ecliptic position</param>
        /// <param name="ez">z-coordinate of the ecliptic position</param>
        /// <param name="elat">ecliptic latitude</param>
        /// <param name="elon">ecliptic longitude</param>
        public Ecliptic(double ex, double ey, double ez, double elat, double elon)
        {
            this.ex = ex;
            this.ey = ey;
            this.ez = ez;
            this.elat = elat;
            this.elon = elon;
        }
    }

    /// <summary>
    /// Coordinates of a celestial body as seen by a topocentric observer.
    /// </summary>
    /// <remarks>
    /// Contains horizontal and equatorial coordinates seen by an observer on or near
    /// the surface of the Earth (a topocentric observer).
    /// Optionally corrected for atmospheric refraction.
    /// </remarks>
    public struct Topocentric
    {
        /// <summary>
        /// Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West.
        /// </summary>
        public readonly double azimuth;

        /// <summary>
        /// Angle in degrees above (positive) or below (negative) the observer's horizon.
        /// </summary>
        public readonly double altitude;

        /// <summary>
        /// Right ascension in sidereal hours.
        /// </summary>
        public readonly double ra;

        /// <summary>
        /// Declination in degrees.
        /// </summary>
        public readonly double dec;

        /// <summary>
        /// Creates a topocentric position object.
        /// </summary>
        /// <param name="azimuth">Compass direction around the horizon in degrees. 0=North, 90=East, 180=South, 270=West.</param>
        /// <param name="altitude">Angle in degrees above (positive) or below (negative) the observer's horizon.</param>
        /// <param name="ra">Right ascension in sidereal hours.</param>
        /// <param name="dec">Declination in degrees.</param>
        public Topocentric(double azimuth, double altitude, double ra, double dec)
        {
            this.azimuth = azimuth;
            this.altitude = altitude;
            this.ra = ra;
            this.dec = dec;
        }
    }

    /// <summary>
    /// The dates and times of changes of season for a given calendar year.
    /// Call #Astronomy.Seasons to calculate this data structure for a given year.
    /// </summary>
    public struct SeasonsInfo
    {
        /// <summary>
        /// The date and time of the March equinox for the specified year.
        /// </summary>
        public readonly AstroTime mar_equinox;

        /// <summary>
        /// The date and time of the June soltice for the specified year.
        /// </summary>
        public readonly AstroTime jun_solstice;

        /// <summary>
        /// The date and time of the September equinox for the specified year.
        /// </summary>
        public readonly AstroTime sep_equinox;

        /// <summary>
        /// The date and time of the December solstice for the specified year.
        /// </summary>
        public readonly AstroTime dec_solstice;

        internal SeasonsInfo(AstroTime mar_equinox, AstroTime jun_solstice, AstroTime sep_equinox, AstroTime dec_solstice)
        {
            this.mar_equinox = mar_equinox;
            this.jun_solstice = jun_solstice;
            this.sep_equinox = sep_equinox;
            this.dec_solstice = dec_solstice;
        }
    }

    /// <summary>
    /// A lunar quarter event (new moon, first quarter, full moon, or third quarter) along with its date and time.
    /// </summary>
    public struct MoonQuarterInfo
    {
        /// <summary>
        /// 0=new moon, 1=first quarter, 2=full moon, 3=third quarter.
        /// </summary>
        public readonly int quarter;

        /// <summary>
        /// The date and time of the lunar quarter.
        /// </summary>
        public readonly AstroTime time;

        internal MoonQuarterInfo(int quarter, AstroTime time)
        {
            this.quarter = quarter;
            this.time = time;
        }
    }

    /// <summary>
    /// Information about a celestial body crossing a specific hour angle.
    /// </summary>
    /// <remarks>
    /// Returned by the function #Astronomy.SearchHourAngle to report information about
    /// a celestial body crossing a certain hour angle as seen by a specified topocentric observer.
    /// </remarks>
    public struct HourAngleInfo
    {
        /// <summary>The date and time when the body crosses the specified hour angle.</summary>
        public readonly AstroTime time;

        /// <summary>Apparent coordinates of the body at the time it crosses the specified hour angle.</summary>
        public readonly Topocentric hor;

        /// <summary>
        /// Creates a struct that represents a celestial body crossing a specific hour angle.
        /// </summary>
        /// <param name="time">The date and time when the body crosses the specified hour angle.</param>
        /// <param name="hor">Apparent coordinates of the body at the time it crosses the specified hour angle.</param>
        public HourAngleInfo(AstroTime time, Topocentric hor)
        {
            this.time = time;
            this.hor = hor;
        }
    }

    /// <summary>
    /// Contains information about the visibility of a celestial body at a given date and time.
    /// See #Astronomy.Elongation for more detailed information about the members of this structure.
    /// See also #Astronomy.SearchMaxElongation for how to search for maximum elongation events.
    /// </summary>
    public struct ElongationInfo
    {
        /// <summary>The date and time of the observation.</summary>
        public readonly AstroTime time;

        /// <summary>Whether the body is best seen in the morning or the evening.</summary>
        public readonly Visibility visibility;

        /// <summary>The angle in degrees between the body and the Sun, as seen from the Earth.</summary>
        public readonly double elongation;

        /// <summary>The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth.</summary>
        public readonly double ecliptic_separation;

        /// <summary>
        /// Creates a structure that represents an elongation event.
        /// </summary>
        /// <param name="time">The date and time of the observation.</param>
        /// <param name="visibility">Whether the body is best seen in the morning or the evening.</param>
        /// <param name="elongation">The angle in degrees between the body and the Sun, as seen from the Earth.</param>
        /// <param name="ecliptic_separation">The difference between the ecliptic longitudes of the body and the Sun, as seen from the Earth.</param>
        public ElongationInfo(AstroTime time, Visibility visibility, double elongation, double ecliptic_separation)
        {
            this.time = time;
            this.visibility = visibility;
            this.elongation = elongation;
            this.ecliptic_separation = ecliptic_separation;
        }
    }

    /// <summary>
    /// The type of apsis: pericenter (closest approach) or apocenter (farthest distance).
    /// </summary>
    public enum ApsisKind
    {
        /// <summary>The body is at its closest approach to the object it orbits.</summary>
        Pericenter,

        /// <summary>The body is at its farthest distance from the object it orbits.</summary>
        Apocenter,
    }

    /// <summary>
    /// An apsis event: pericenter (closest approach) or apocenter (farthest distance).
    /// </summary>
    /// <remarks>
    /// For the Moon orbiting the Earth, or a planet orbiting the Sun, an *apsis* is an
    /// event where the orbiting body reaches its closest or farthest point from the primary body.
    /// The closest approach is called *pericenter* and the farthest point is *apocenter*.
    ///
    /// More specific terminology is common for particular orbiting bodies.
    /// The Moon's closest approach to the Earth is called *perigee* and its farthest
    /// point is called *apogee*. The closest approach of a planet to the Sun is called
    /// *perihelion* and the furthest point is called *aphelion*.
    ///
    /// This data structure is returned by #Astronomy.SearchLunarApsis and #Astronomy.NextLunarApsis
    /// to iterate through consecutive alternating perigees and apogees.
    /// </remarks>
    public struct ApsisInfo
    {
        /// <summary>The date and time of the apsis.</summary>
        public readonly AstroTime time;

        /// <summary>Whether this is a pericenter or apocenter event.</summary>
        public readonly ApsisKind kind;

        /// <summary>The distance between the centers of the bodies in astronomical units.</summary>
        public readonly double dist_au;

        /// <summary>The distance between the centers of the bodies in kilometers.</summary>
        public readonly double dist_km;

        internal ApsisInfo(AstroTime time, ApsisKind kind, double dist_au)
        {
            this.time = time;
            this.kind = kind;
            this.dist_au = dist_au;
            this.dist_km = dist_au * Astronomy.KM_PER_AU;
        }
    }

    /// <summary>different kinds of lunar/solar eclipses.</summary>
    public enum EclipseKind
    {
        /// <summary>No eclipse found.</summary>
        None,

        /// <summary>A penumbral lunar eclipse. (Never used for a solar eclipse.)</summary>
        Penumbral,

        /// <summary>A partial lunar/solar eclipse.</summary>
        Partial,

        /// <summary>An annular solar eclipse. (Never used for a lunar eclipse.)</summary>
        Annular,

        /// <summary>A total lunar/solar eclipse.</summary>
        Total,
    }

    /// <summary>
    /// Information about a lunar eclipse.
    /// </summary>
    /// <remarks>
    /// Returned by #Astronomy.SearchLunarEclipse or #Astronomy.NextLunarEclipse
    /// to report information about a lunar eclipse event.
    /// When a lunar eclipse is found, it is classified as penumbral, partial, or total.
    /// Penumbral eclipses are difficult to observe, because the moon is only slightly dimmed
    /// by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
    /// Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
    /// Total eclipses occur when the entire Moon passes into the Earth's umbra.
    ///
    /// The `kind` field thus holds `EclipseKind.Penumbral`, `EclipseKind.Partial`,
    /// or `EclipseKind.Total`, depending on the kind of lunar eclipse found.
    ///
    /// Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
    ///
    /// Fields `sd_penum`, `sd_partial`, and `sd_total` hold the semi-duration of each phase
    /// of the eclipse, which is half of the amount of time the eclipse spends in each
    /// phase (expressed in minutes), or 0 if the eclipse never reaches that phase.
    /// By converting from minutes to days, and subtracting/adding with `peak`, the caller
    /// may determine the date and time of the beginning/end of each eclipse phase.
    /// </remarks>
    public struct LunarEclipseInfo
    {
        /// <summary>The type of lunar eclipse found.</summary>
        public EclipseKind kind;

        /// <summary>The time of the eclipse at its peak.</summary>
        public AstroTime peak;

        /// <summary>The semi-duration of the penumbral phase in minutes.</summary>
        public double sd_penum;

        /// <summary>The semi-duration of the partial phase in minutes, or 0.0 if none.</summary>
        public double sd_partial;

        /// <summary>The semi-duration of the total phase in minutes, or 0.0 if none.</summary>
        public double sd_total;

        internal LunarEclipseInfo(EclipseKind kind, AstroTime peak, double sd_penum, double sd_partial, double sd_total)
        {
            this.kind = kind;
            this.peak = peak;
            this.sd_penum = sd_penum;
            this.sd_partial = sd_partial;
            this.sd_total = sd_total;
        }
    }


    /// <summary>
    /// Reports the time and geographic location of the peak of a solar eclipse.
    /// </summary>
    /// <remarks>
    /// Returned by #Astronomy.SearchGlobalSolarEclipse or #Astronomy.NextGlobalSolarEclipse
    /// to report information about a solar eclipse event.
    ///
    /// Field `peak` holds the date and time of the peak of the eclipse, defined as
    /// the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.
    ///
    /// The eclipse is classified as partial, annular, or total, depending on the
    /// maximum amount of the Sun's disc obscured, as seen at the peak location
    /// on the surface of the Earth.
    ///
    /// The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    /// A total eclipse is when the peak observer sees the Sun completely blocked by the Moon.
    /// An annular eclipse is like a total eclipse, but the Moon is too far from the Earth's surface
    /// to completely block the Sun; instead, the Sun takes on a ring-shaped appearance.
    /// A partial eclipse is when the Moon blocks part of the Sun's disc, but nobody on the Earth
    /// observes either a total or annular eclipse.
    ///
    /// If `kind` is `EclipseKind.Total` or `EclipseKind.Annular`, the `latitude` and `longitude`
    /// fields give the geographic coordinates of the center of the Moon's shadow projected
    /// onto the daytime side of the Earth at the instant of the eclipse's peak.
    /// If `kind` has any other value, `latitude` and `longitude` are undefined and should
    /// not be used.
    /// </remarks>
    public struct GlobalSolarEclipseInfo
    {
        /// <summary>The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.</summary>
        public EclipseKind kind;

        /// <summary>The date and time of the eclipse at its peak.</summary>
        public AstroTime peak;

        /// <summary>The distance between the Sun/Moon shadow axis and the center of the Earth, in kilometers.</summary>
        public double distance;

        /// <summary>The geographic latitude at the center of the peak eclipse shadow.</summary>
        public double latitude;

        /// <summary>The geographic longitude at the center of the peak eclipse shadow.</summary>
        public double longitude;
    }


    /// <summary>
    /// Holds a time and the observed altitude of the Sun at that time.
    /// </summary>
    /// <remarks>
    /// When reporting a solar eclipse observed at a specific location on the Earth
    /// (a "local" solar eclipse), a series of events occur. In addition
    /// to the time of each event, it is important to know the altitude of the Sun,
    /// because each event may be invisible to the observer if the Sun is below
    /// the horizon (i.e. it at night).
    ///
    /// If `altitude` is negative, the event is theoretical only; it would be
    /// visible if the Earth were transparent, but the observer cannot actually see it.
    /// If `altitude` is positive but less than a few degrees, visibility will be impaired by
    /// atmospheric interference (sunrise or sunset conditions).
    /// </remarks>
    public struct EclipseEvent
    {
        /// <summary>The date and time of the event.</summary>
        public AstroTime time;

        /// <summary>
        /// The angular altitude of the center of the Sun above/below the horizon, at `time`,
        /// corrected for atmospheric refraction and expressed in degrees.
        /// </summary>
        public double altitude;
    }


    /// <summary>
    /// Information about a solar eclipse as seen by an observer at a given time and geographic location.
    /// </summary>
    /// <remarks>
    /// Returned by #Astronomy.SearchLocalSolarEclipse or #Astronomy.NextLocalSolarEclipse
    /// to report information about a solar eclipse as seen at a given geographic location.
    ///
    /// When a solar eclipse is found, it is classified as partial, annular, or total.
    /// The `kind` field thus holds `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.
    /// A partial solar eclipse is when the Moon does not line up directly enough with the Sun
    /// to completely block the Sun's light from reaching the observer.
    /// An annular eclipse occurs when the Moon's disc is completely visible against the Sun
    /// but the Moon is too far away to completely block the Sun's light; this leaves the
    /// Sun with a ring-like appearance.
    /// A total eclipse occurs when the Moon is close enough to the Earth and aligned with the
    /// Sun just right to completely block all sunlight from reaching the observer.
    ///
    /// There are 5 "event" fields, each of which contains a time and a solar altitude.
    /// Field `peak` holds the date and time of the center of the eclipse, when it is at its peak.
    /// The fields `partial_begin` and `partial_end` are always set, and indicate when
    /// the eclipse begins/ends. If the eclipse reaches totality or becomes annular,
    /// `total_begin` and `total_end` indicate when the total/annular phase begins/ends.
    /// When an event field is valid, the caller must also check its `altitude` field to
    /// see whether the Sun is above the horizon at the time indicated by the `time` field.
    /// See #EclipseEvent for more information.
    /// </remarks>
    public struct LocalSolarEclipseInfo
    {
        /// <summary>The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.</summary>
        public EclipseKind  kind;

        /// <summary>The time and Sun altitude at the beginning of the eclipse.</summary>
        public EclipseEvent partial_begin;

        /// <summary>If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase begins; otherwise invalid.</summary>
        public EclipseEvent total_begin;

        /// <summary>The time and Sun altitude when the eclipse reaches its peak.</summary>
        public EclipseEvent peak;

        /// <summary>If this is an annular or a total eclipse, the time and Sun altitude when annular/total phase ends; otherwise invalid.</summary>
        public EclipseEvent total_end;

        /// <summary>The time and Sun altitude at the end of the eclipse.</summary>
        public EclipseEvent partial_end;
    }


    /// <summary>
    /// Information about a transit of Mercury or Venus, as seen from the Earth.
    /// </summary>
    /// <remarks>
    /// Returned by #Astronomy.SearchTransit or #Astronomy.NextTransit to report
    /// information about a transit of Mercury or Venus.
    /// A transit is when Mercury or Venus passes between the Sun and Earth so that
    /// the other planet is seen in silhouette against the Sun.
    ///
    /// The `start` field reports the moment in time when the planet first becomes
    /// visible against the Sun in its background.
    /// The `peak` field reports when the planet is most aligned with the Sun,
    /// as seen from the Earth.
    /// The `finish` field reports the last moment when the planet is visible
    /// against the Sun in its background.
    ///
    /// The calculations are performed from the point of view of a geocentric observer.
    /// </remarks>
    public struct TransitInfo
    {
        /// <summary>Date and time at the beginning of the transit.</summary>
        public AstroTime start;

        /// <summary>Date and time of the peak of the transit.</summary>
        public AstroTime peak;

        /// <summary>Date and time at the end of the transit.</summary>
        public AstroTime finish;

        /// <summary>Angular separation in arcminutes between the centers of the Sun and the planet at time `peak`.</summary>
        public double separation;
    }


    internal struct ShadowInfo
    {
        public AstroTime time;
        public double u;    // dot product of (heliocentric earth) and (geocentric moon): defines the shadow plane where the Moon is
        public double r;    // km distance between center of Moon and the line passing through the centers of the Sun and Earth.
        public double k;    // umbra radius in km, at the shadow plane
        public double p;    // penumbra radius in km, at the shadow plane
        public AstroVector target;      // coordinates of target body relative to shadow-casting body at 'time'
        public AstroVector dir;         // heliocentric coordinates of shadow-casting body at 'time'

        public ShadowInfo(AstroTime time, double u, double r, double k, double p, AstroVector target, AstroVector dir)
        {
            this.time = time;
            this.u = u;
            this.r = r;
            this.k = k;
            this.p = p;
            this.target = target;
            this.dir = dir;
        }
    }

    /// <summary>
    /// Information about the brightness and illuminated shape of a celestial body.
    /// </summary>
    /// <remarks>
    /// Returned by the functions #Astronomy.Illumination and #Astronomy.SearchPeakMagnitude
    /// to report the visual magnitude and illuminated fraction of a celestial body at a given date and time.
    /// </remarks>
    public struct IllumInfo
    {
        /// <summary>The date and time of the observation.</summary>
        public readonly AstroTime time;

        /// <summary>The visual magnitude of the body. Smaller values are brighter.</summary>
        public readonly double  mag;

        /// <summary>The angle in degrees between the Sun and the Earth, as seen from the body. Indicates the body's phase as seen from the Earth.</summary>
        public readonly double phase_angle;

        /// <summary>The distance between the Sun and the body at the observation time.</summary>
        public readonly double helio_dist;

        /// <summary>For Saturn, the tilt angle in degrees of its rings as seen from Earth. For all other bodies, 0.</summary>
        public readonly double ring_tilt;

        internal IllumInfo(AstroTime time, double mag, double phase_angle, double helio_dist, double ring_tilt)
        {
            this.time = time;
            this.mag = mag;
            this.phase_angle = phase_angle;
            this.helio_dist = helio_dist;
            this.ring_tilt = ring_tilt;
        }
    }

    /// <summary>
    /// Represents a function whose ascending root is to be found.
    /// See #Astronomy.Search.
    /// </summary>
    public abstract class SearchContext
    {
        /// <summary>
        /// Evaluates the function at a given time
        /// </summary>
        /// <param name="time">The time at which to evaluate the function.</param>
        /// <returns>The floating point value of the function at the specified time.</returns>
        public abstract double Eval(AstroTime time);
    }

    internal class SearchContext_MagnitudeSlope: SearchContext
    {
        private readonly Body body;

        public SearchContext_MagnitudeSlope(Body body)
        {
            this.body = body;
        }

        public override double Eval(AstroTime time)
        {
            /*
                The Search() function finds a transition from negative to positive values.
                The derivative of magnitude y with respect to time t (dy/dt)
                is negative as an object gets brighter, because the magnitude numbers
                get smaller. At peak magnitude dy/dt = 0, then as the object gets dimmer,
                dy/dt > 0.
            */
            const double dt = 0.01;
            AstroTime t1 = time.AddDays(-dt/2);
            AstroTime t2 = time.AddDays(+dt/2);
            IllumInfo y1 = Astronomy.Illumination(body, t1);
            IllumInfo y2 = Astronomy.Illumination(body, t2);
            return (y2.mag - y1.mag) / dt;
        }
    }

    internal class SearchContext_NegElongSlope: SearchContext
    {
        private readonly Body body;

        public SearchContext_NegElongSlope(Body body)
        {
            this.body = body;
        }

        public override double Eval(AstroTime time)
        {
            const double dt = 0.1;
            AstroTime t1 = time.AddDays(-dt/2.0);
            AstroTime t2 = time.AddDays(+dt/2.0);

            double e1 = Astronomy.AngleFromSun(body, t1);
            double e2 = Astronomy.AngleFromSun(body, t2);
            return (e1 - e2)/dt;
        }
    }

    internal class SearchContext_SunOffset: SearchContext
    {
        private readonly double targetLon;

        public SearchContext_SunOffset(double targetLon)
        {
            this.targetLon = targetLon;
        }

        public override double Eval(AstroTime time)
        {
            Ecliptic ecl = Astronomy.SunPosition(time);
            return Astronomy.LongitudeOffset(ecl.elon - targetLon);
        }
    }

    internal class SearchContext_MoonOffset: SearchContext
    {
        private readonly double targetLon;

        public SearchContext_MoonOffset(double targetLon)
        {
            this.targetLon = targetLon;
        }

        public override double Eval(AstroTime time)
        {
            double angle = Astronomy.MoonPhase(time);
            return Astronomy.LongitudeOffset(angle - targetLon);
        }
    }

    internal class SearchContext_PeakAltitude: SearchContext
    {
        private readonly Body body;
        private readonly int direction;
        private readonly Observer observer;
        private readonly double body_radius_au;

        public SearchContext_PeakAltitude(Body body, Direction direction, Observer observer)
        {
            this.body = body;
            this.direction = (int)direction;
            this.observer = observer;

            switch (body)
            {
                case Body.Sun:
                    this.body_radius_au = Astronomy.SUN_RADIUS_AU;
                    break;

                case Body.Moon:
                    this.body_radius_au = Astronomy.MOON_EQUATORIAL_RADIUS_AU;
                    break;

                default:
                    this.body_radius_au = 0.0;
                    break;
            }
        }

        public override double Eval(AstroTime time)
        {
            /*
                Return the angular altitude above or below the horizon
                of the highest part (the peak) of the given object.
                This is defined as the apparent altitude of the center of the body plus
                the body's angular radius.
                The 'direction' parameter controls whether the angle is measured
                positive above the horizon or positive below the horizon,
                depending on whether the caller wants rise times or set times, respectively.
            */

            Equatorial ofdate = Astronomy.Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);

            /* We calculate altitude without refraction, then add fixed refraction near the horizon. */
            /* This gives us the time of rise/set without the extra work. */
            Topocentric hor = Astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None);

            return direction * (hor.altitude + Astronomy.RAD2DEG*(body_radius_au / ofdate.dist) + Astronomy.REFRACTION_NEAR_HORIZON);
        }
    }

    internal class SearchContext_MoonDistanceSlope: SearchContext
    {
        private readonly int direction;

        public SearchContext_MoonDistanceSlope(int direction)
        {
            this.direction = direction;
        }

        public static double MoonDistance(AstroTime time)
        {
            var context = new MoonContext(time.tt / 36525.0);
            MoonResult moon = context.CalcMoon();
            return moon.distance_au;
        }

        public override double Eval(AstroTime time)
        {
            const double dt = 0.001;
            AstroTime t1 = time.AddDays(-dt/2.0);
            AstroTime t2 = time.AddDays(+dt/2.0);
            double dist1 = MoonDistance(t1);
            double dist2 = MoonDistance(t2);
            return direction * (dist2 - dist1)/dt;
        }
    }

    internal class SearchContext_PlanetDistanceSlope: SearchContext
    {
        private readonly double direction;
        private readonly Body body;

        public SearchContext_PlanetDistanceSlope(double direction, Body body)
        {
            this.direction = direction;
            this.body = body;
        }

        public override double Eval(AstroTime time)
        {
            const double dt = 0.001;
            AstroTime t1 = time.AddDays(-dt/2.0);
            AstroTime t2 = time.AddDays(+dt/2.0);
            double r1 = Astronomy.HelioDistance(body, t1);
            double r2 = Astronomy.HelioDistance(body, t2);
            return direction * (r2 - r1) / dt;
        }
    }

    internal class SearchContext_EarthShadow: SearchContext
    {
        private readonly double radius_limit;
        private readonly double direction;

        public SearchContext_EarthShadow(double radius_limit, double direction)
        {
            this.radius_limit = radius_limit;
            this.direction = direction;
        }

        public override double Eval(AstroTime time)
        {
            return direction * (Astronomy.EarthShadow(time).r - radius_limit);
        }
    }

    internal class SearchContext_EarthShadowSlope: SearchContext
    {
        public override double Eval(AstroTime time)
        {
            const double dt = 1.0 / 86400.0;
            AstroTime t1 = time.AddDays(-dt);
            AstroTime t2 = time.AddDays(+dt);
            ShadowInfo shadow1 = Astronomy.EarthShadow(t1);
            ShadowInfo shadow2 = Astronomy.EarthShadow(t2);
            return (shadow2.r - shadow1.r) / dt;
        }
    }

    internal class SearchContext_MoonShadowSlope: SearchContext
    {
        public override double Eval(AstroTime time)
        {
            const double dt = 1.0 / 86400.0;
            AstroTime t1 = time.AddDays(-dt);
            AstroTime t2 = time.AddDays(+dt);
            ShadowInfo shadow1 = Astronomy.MoonShadow(t1);
            ShadowInfo shadow2 = Astronomy.MoonShadow(t2);
            return (shadow2.r - shadow1.r) / dt;
        }
    }

    internal class SearchContext_LocalMoonShadowSlope: SearchContext
    {
        private readonly Observer observer;

        public SearchContext_LocalMoonShadowSlope(Observer observer)
        {
            this.observer = observer;
        }

        public override double Eval(AstroTime time)
        {
            const double dt = 1.0 / 86400.0;
            AstroTime t1 = time.AddDays(-dt);
            AstroTime t2 = time.AddDays(+dt);
            ShadowInfo shadow1 = Astronomy.LocalMoonShadow(t1, observer);
            ShadowInfo shadow2 = Astronomy.LocalMoonShadow(t2, observer);
            return (shadow2.r - shadow1.r) / dt;
        }
    }

    internal class SearchContext_PlanetShadowSlope: SearchContext
    {
        private Body body;
        private double planet_radius_km;

        public SearchContext_PlanetShadowSlope(Body body, double planet_radius_km)
        {
            this.body = body;
            this.planet_radius_km = planet_radius_km;
        }

        public override double Eval(AstroTime time)
        {
            const double dt = 1.0 / 86400.0;
            ShadowInfo shadow1 = Astronomy.PlanetShadow(body, planet_radius_km, time.AddDays(-dt));
            ShadowInfo shadow2 = Astronomy.PlanetShadow(body, planet_radius_km, time.AddDays(+dt));
            return (shadow2.r - shadow1.r) / dt;
        }
    }

    internal class SearchContext_PlanetShadowBoundary: SearchContext
    {
        private Body body;
        private double planet_radius_km;
        private double direction;

        public SearchContext_PlanetShadowBoundary(Body body, double planet_radius_km, double direction)
        {
            this.body = body;
            this.planet_radius_km = planet_radius_km;
            this.direction = direction;
        }

        public override double Eval(AstroTime time)
        {
            ShadowInfo shadow = Astronomy.PlanetShadow(body, planet_radius_km, time);
            return direction * (shadow.r - shadow.p);
        }
    }

    internal class SearchContext_LocalEclipseTransition: SearchContext
    {
        private readonly Func<ShadowInfo,double> func;
        private readonly double direction;
        private readonly Observer observer;

        public SearchContext_LocalEclipseTransition(Func<ShadowInfo,double> func, double direction, Observer observer)
        {
            this.func = func;
            this.direction = direction;
            this.observer = observer;
        }

        public override double Eval(AstroTime time)
        {
            ShadowInfo shadow = Astronomy.LocalMoonShadow(time, observer);
            return direction * func(shadow);
        }
    }


    internal class PascalArray2<ElemType>
    {
        private readonly int xmin;
        private readonly int xmax;
        private readonly int ymin;
        private readonly int ymax;
        private readonly ElemType[,] array;

        public PascalArray2(int xmin, int xmax, int ymin, int ymax)
        {
            this.xmin = xmin;
            this.xmax = xmax;
            this.ymin = ymin;
            this.ymax = ymax;
            this.array = new ElemType[(xmax - xmin) + 1, (ymax - ymin) + 1];
        }

        public ElemType this[int x, int y]
        {
            get { return array[x - xmin, y - ymin]; }
            set { array[x - xmin, y - ymin] = value; }
        }
    }

    internal class MoonContext
    {
        double T;
        double DGAM;
        double DLAM, N, GAM1C, SINPI;
        double L0, L, LS, F, D, S;
        double DL0, DL, DLS, DF, DD, DS;
        PascalArray2<double> CO = new PascalArray2<double>(-6, 6, 1, 4);
        PascalArray2<double> SI = new PascalArray2<double>(-6, 6, 1, 4);

        static double Frac(double x)
        {
            return x - Math.Floor(x);
        }

        static void AddThe(
            double c1, double s1, double c2, double s2,
            out double c, out double s)
        {
            c = c1*c2 - s1*s2;
            s = s1*c2 + c1*s2;
        }

        static double Sine(double phi)
        {
            /* sine, of phi in revolutions, not radians */
            return Math.Sin(2.0 * Math.PI * phi);
        }

        void LongPeriodic()
        {
            double S1 = Sine(0.19833+0.05611*T);
            double S2 = Sine(0.27869+0.04508*T);
            double S3 = Sine(0.16827-0.36903*T);
            double S4 = Sine(0.34734-5.37261*T);
            double S5 = Sine(0.10498-5.37899*T);
            double S6 = Sine(0.42681-0.41855*T);
            double S7 = Sine(0.14943-5.37511*T);

            DL0 = 0.84*S1+0.31*S2+14.27*S3+ 7.26*S4+ 0.28*S5+0.24*S6;
            DL  = 2.94*S1+0.31*S2+14.27*S3+ 9.34*S4+ 1.12*S5+0.83*S6;
            DLS =-6.40*S1                                   -1.89*S6;
            DF  = 0.21*S1+0.31*S2+14.27*S3-88.70*S4-15.30*S5+0.24*S6-1.86*S7;
            DD  = DL0-DLS;
            DGAM  = -3332E-9 * Sine(0.59734-5.37261*T)
                    -539E-9 * Sine(0.35498-5.37899*T)
                    -64E-9 * Sine(0.39943-5.37511*T);
        }

        private readonly int[] I = new int[4];

        void Term(int p, int q, int r, int s, out double x, out double y)
        {
            I[0] = p;
            I[1] = q;
            I[2] = r;
            I[3] = s;
            x = 1.0;
            y = 0.0;

            for (int k=1; k<=4; ++k)
                if (I[k-1] != 0.0)
                    AddThe(x, y, CO[I[k-1], k], SI[I[k-1], k], out x, out y);
        }

        void AddSol(
            double coeffl,
            double coeffs,
            double coeffg,
            double coeffp,
            int p,
            int q,
            int r,
            int s)
        {
            double x, y;
            Term(p, q, r, s, out x, out y);
            DLAM += coeffl*y;
            DS += coeffs*y;
            GAM1C += coeffg*x;
            SINPI += coeffp*x;
        }

        void ADDN(double coeffn, int p, int q, int r, int s, out double x, out double y)
        {
            Term(p, q, r, s, out x, out y);
            N += coeffn * y;
        }

        void SolarN()
        {
            double x, y;

            N = 0.0;
            ADDN(-526.069, 0, 0,1,-2, out x, out y);
            ADDN(  -3.352, 0, 0,1,-4, out x, out y);
            ADDN( +44.297,+1, 0,1,-2, out x, out y);
            ADDN(  -6.000,+1, 0,1,-4, out x, out y);
            ADDN( +20.599,-1, 0,1, 0, out x, out y);
            ADDN( -30.598,-1, 0,1,-2, out x, out y);
            ADDN( -24.649,-2, 0,1, 0, out x, out y);
            ADDN(  -2.000,-2, 0,1,-2, out x, out y);
            ADDN( -22.571, 0,+1,1,-2, out x, out y);
            ADDN( +10.985, 0,-1,1,-2, out x, out y);
        }

        void Planetary()
        {
            DLAM +=
                +0.82*Sine(0.7736  -62.5512*T)+0.31*Sine(0.0466 -125.1025*T)
                +0.35*Sine(0.5785  -25.1042*T)+0.66*Sine(0.4591+1335.8075*T)
                +0.64*Sine(0.3130  -91.5680*T)+1.14*Sine(0.1480+1331.2898*T)
                +0.21*Sine(0.5918+1056.5859*T)+0.44*Sine(0.5784+1322.8595*T)
                +0.24*Sine(0.2275   -5.7374*T)+0.28*Sine(0.2965   +2.6929*T)
                +0.33*Sine(0.3132   +6.3368*T);
        }

        internal MoonContext(double centuries_since_j2000)
        {
            int I, J, MAX;
            double T2, ARG, FAC;
            double c, s;

            T = centuries_since_j2000;
            T2 = T*T;
            DLAM = 0;
            DS = 0;
            GAM1C = 0;
            SINPI = 3422.7000;
            LongPeriodic();
            L0 = Astronomy.PI2*Frac(0.60643382+1336.85522467*T-0.00000313*T2) + DL0/Astronomy.ARC;
            L  = Astronomy.PI2*Frac(0.37489701+1325.55240982*T+0.00002565*T2) + DL /Astronomy.ARC;
            LS = Astronomy.PI2*Frac(0.99312619+  99.99735956*T-0.00000044*T2) + DLS/Astronomy.ARC;
            F  = Astronomy.PI2*Frac(0.25909118+1342.22782980*T-0.00000892*T2) + DF /Astronomy.ARC;
            D  = Astronomy.PI2*Frac(0.82736186+1236.85308708*T-0.00000397*T2) + DD /Astronomy.ARC;
            for (I=1; I<=4; ++I)
            {
                switch(I)
                {
                    case 1:  ARG=L;  MAX=4; FAC=1.000002208;               break;
                    case 2:  ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T; break;
                    case 3:  ARG=F;  MAX=4; FAC=1.000002708+139.978*DGAM;  break;
                    default: ARG=D;  MAX=6; FAC=1.0;                       break;
                }
                CO[0,I] = 1.0;
                CO[1,I] = Math.Cos(ARG)*FAC;
                SI[0,I] = 0.0;
                SI[1,I] = Math.Sin(ARG)*FAC;
                for (J=2; J<=MAX; ++J)
                {
                    AddThe(CO[J-1,I], SI[J-1,I], CO[1,I], SI[1,I], out c, out s);
                    CO[J,I] = c;
                    SI[J,I] = s;
                }

                for (J=1; J<=MAX; ++J)
                {
                    CO[-J,I] =  CO[J,I];
                    SI[-J,I] = -SI[J,I];
                }
            }
        }

        internal MoonResult CalcMoon()
        {
            ++Astronomy.CalcMoonCount;

            AddSol(    13.9020,    14.0600,    -0.0010,     0.2607, 0, 0, 0, 4);
            AddSol(     0.4030,    -4.0100,     0.3940,     0.0023, 0, 0, 0, 3);
            AddSol(  2369.9120,  2373.3600,     0.6010,    28.2333, 0, 0, 0, 2);
            AddSol(  -125.1540,  -112.7900,    -0.7250,    -0.9781, 0, 0, 0, 1);
            AddSol(     1.9790,     6.9800,    -0.4450,     0.0433, 1, 0, 0, 4);
            AddSol(   191.9530,   192.7200,     0.0290,     3.0861, 1, 0, 0, 2);
            AddSol(    -8.4660,   -13.5100,     0.4550,    -0.1093, 1, 0, 0, 1);
            AddSol( 22639.5000, 22609.0700,     0.0790,   186.5398, 1, 0, 0, 0);
            AddSol(    18.6090,     3.5900,    -0.0940,     0.0118, 1, 0, 0,-1);
            AddSol( -4586.4650, -4578.1300,    -0.0770,    34.3117, 1, 0, 0,-2);
            AddSol(     3.2150,     5.4400,     0.1920,    -0.0386, 1, 0, 0,-3);
            AddSol(   -38.4280,   -38.6400,     0.0010,     0.6008, 1, 0, 0,-4);
            AddSol(    -0.3930,    -1.4300,    -0.0920,     0.0086, 1, 0, 0,-6);
            AddSol(    -0.2890,    -1.5900,     0.1230,    -0.0053, 0, 1, 0, 4);
            AddSol(   -24.4200,   -25.1000,     0.0400,    -0.3000, 0, 1, 0, 2);
            AddSol(    18.0230,    17.9300,     0.0070,     0.1494, 0, 1, 0, 1);
            AddSol(  -668.1460,  -126.9800,    -1.3020,    -0.3997, 0, 1, 0, 0);
            AddSol(     0.5600,     0.3200,    -0.0010,    -0.0037, 0, 1, 0,-1);
            AddSol(  -165.1450,  -165.0600,     0.0540,     1.9178, 0, 1, 0,-2);
            AddSol(    -1.8770,    -6.4600,    -0.4160,     0.0339, 0, 1, 0,-4);
            AddSol(     0.2130,     1.0200,    -0.0740,     0.0054, 2, 0, 0, 4);
            AddSol(    14.3870,    14.7800,    -0.0170,     0.2833, 2, 0, 0, 2);
            AddSol(    -0.5860,    -1.2000,     0.0540,    -0.0100, 2, 0, 0, 1);
            AddSol(   769.0160,   767.9600,     0.1070,    10.1657, 2, 0, 0, 0);
            AddSol(     1.7500,     2.0100,    -0.0180,     0.0155, 2, 0, 0,-1);
            AddSol(  -211.6560,  -152.5300,     5.6790,    -0.3039, 2, 0, 0,-2);
            AddSol(     1.2250,     0.9100,    -0.0300,    -0.0088, 2, 0, 0,-3);
            AddSol(   -30.7730,   -34.0700,    -0.3080,     0.3722, 2, 0, 0,-4);
            AddSol(    -0.5700,    -1.4000,    -0.0740,     0.0109, 2, 0, 0,-6);
            AddSol(    -2.9210,   -11.7500,     0.7870,    -0.0484, 1, 1, 0, 2);
            AddSol(     1.2670,     1.5200,    -0.0220,     0.0164, 1, 1, 0, 1);
            AddSol(  -109.6730,  -115.1800,     0.4610,    -0.9490, 1, 1, 0, 0);
            AddSol(  -205.9620,  -182.3600,     2.0560,     1.4437, 1, 1, 0,-2);
            AddSol(     0.2330,     0.3600,     0.0120,    -0.0025, 1, 1, 0,-3);
            AddSol(    -4.3910,    -9.6600,    -0.4710,     0.0673, 1, 1, 0,-4);
            AddSol(     0.2830,     1.5300,    -0.1110,     0.0060, 1,-1, 0, 4);
            AddSol(    14.5770,    31.7000,    -1.5400,     0.2302, 1,-1, 0, 2);
            AddSol(   147.6870,   138.7600,     0.6790,     1.1528, 1,-1, 0, 0);
            AddSol(    -1.0890,     0.5500,     0.0210,     0.0000, 1,-1, 0,-1);
            AddSol(    28.4750,    23.5900,    -0.4430,    -0.2257, 1,-1, 0,-2);
            AddSol(    -0.2760,    -0.3800,    -0.0060,    -0.0036, 1,-1, 0,-3);
            AddSol(     0.6360,     2.2700,     0.1460,    -0.0102, 1,-1, 0,-4);
            AddSol(    -0.1890,    -1.6800,     0.1310,    -0.0028, 0, 2, 0, 2);
            AddSol(    -7.4860,    -0.6600,    -0.0370,    -0.0086, 0, 2, 0, 0);
            AddSol(    -8.0960,   -16.3500,    -0.7400,     0.0918, 0, 2, 0,-2);
            AddSol(    -5.7410,    -0.0400,     0.0000,    -0.0009, 0, 0, 2, 2);
            AddSol(     0.2550,     0.0000,     0.0000,     0.0000, 0, 0, 2, 1);
            AddSol(  -411.6080,    -0.2000,     0.0000,    -0.0124, 0, 0, 2, 0);
            AddSol(     0.5840,     0.8400,     0.0000,     0.0071, 0, 0, 2,-1);
            AddSol(   -55.1730,   -52.1400,     0.0000,    -0.1052, 0, 0, 2,-2);
            AddSol(     0.2540,     0.2500,     0.0000,    -0.0017, 0, 0, 2,-3);
            AddSol(     0.0250,    -1.6700,     0.0000,     0.0031, 0, 0, 2,-4);
            AddSol(     1.0600,     2.9600,    -0.1660,     0.0243, 3, 0, 0, 2);
            AddSol(    36.1240,    50.6400,    -1.3000,     0.6215, 3, 0, 0, 0);
            AddSol(   -13.1930,   -16.4000,     0.2580,    -0.1187, 3, 0, 0,-2);
            AddSol(    -1.1870,    -0.7400,     0.0420,     0.0074, 3, 0, 0,-4);
            AddSol(    -0.2930,    -0.3100,    -0.0020,     0.0046, 3, 0, 0,-6);
            AddSol(    -0.2900,    -1.4500,     0.1160,    -0.0051, 2, 1, 0, 2);
            AddSol(    -7.6490,   -10.5600,     0.2590,    -0.1038, 2, 1, 0, 0);
            AddSol(    -8.6270,    -7.5900,     0.0780,    -0.0192, 2, 1, 0,-2);
            AddSol(    -2.7400,    -2.5400,     0.0220,     0.0324, 2, 1, 0,-4);
            AddSol(     1.1810,     3.3200,    -0.2120,     0.0213, 2,-1, 0, 2);
            AddSol(     9.7030,    11.6700,    -0.1510,     0.1268, 2,-1, 0, 0);
            AddSol(    -0.3520,    -0.3700,     0.0010,    -0.0028, 2,-1, 0,-1);
            AddSol(    -2.4940,    -1.1700,    -0.0030,    -0.0017, 2,-1, 0,-2);
            AddSol(     0.3600,     0.2000,    -0.0120,    -0.0043, 2,-1, 0,-4);
            AddSol(    -1.1670,    -1.2500,     0.0080,    -0.0106, 1, 2, 0, 0);
            AddSol(    -7.4120,    -6.1200,     0.1170,     0.0484, 1, 2, 0,-2);
            AddSol(    -0.3110,    -0.6500,    -0.0320,     0.0044, 1, 2, 0,-4);
            AddSol(     0.7570,     1.8200,    -0.1050,     0.0112, 1,-2, 0, 2);
            AddSol(     2.5800,     2.3200,     0.0270,     0.0196, 1,-2, 0, 0);
            AddSol(     2.5330,     2.4000,    -0.0140,    -0.0212, 1,-2, 0,-2);
            AddSol(    -0.3440,    -0.5700,    -0.0250,     0.0036, 0, 3, 0,-2);
            AddSol(    -0.9920,    -0.0200,     0.0000,     0.0000, 1, 0, 2, 2);
            AddSol(   -45.0990,    -0.0200,     0.0000,    -0.0010, 1, 0, 2, 0);
            AddSol(    -0.1790,    -9.5200,     0.0000,    -0.0833, 1, 0, 2,-2);
            AddSol(    -0.3010,    -0.3300,     0.0000,     0.0014, 1, 0, 2,-4);
            AddSol(    -6.3820,    -3.3700,     0.0000,    -0.0481, 1, 0,-2, 2);
            AddSol(    39.5280,    85.1300,     0.0000,    -0.7136, 1, 0,-2, 0);
            AddSol(     9.3660,     0.7100,     0.0000,    -0.0112, 1, 0,-2,-2);
            AddSol(     0.2020,     0.0200,     0.0000,     0.0000, 1, 0,-2,-4);
            AddSol(     0.4150,     0.1000,     0.0000,     0.0013, 0, 1, 2, 0);
            AddSol(    -2.1520,    -2.2600,     0.0000,    -0.0066, 0, 1, 2,-2);
            AddSol(    -1.4400,    -1.3000,     0.0000,     0.0014, 0, 1,-2, 2);
            AddSol(     0.3840,    -0.0400,     0.0000,     0.0000, 0, 1,-2,-2);
            AddSol(     1.9380,     3.6000,    -0.1450,     0.0401, 4, 0, 0, 0);
            AddSol(    -0.9520,    -1.5800,     0.0520,    -0.0130, 4, 0, 0,-2);
            AddSol(    -0.5510,    -0.9400,     0.0320,    -0.0097, 3, 1, 0, 0);
            AddSol(    -0.4820,    -0.5700,     0.0050,    -0.0045, 3, 1, 0,-2);
            AddSol(     0.6810,     0.9600,    -0.0260,     0.0115, 3,-1, 0, 0);
            AddSol(    -0.2970,    -0.2700,     0.0020,    -0.0009, 2, 2, 0,-2);
            AddSol(     0.2540,     0.2100,    -0.0030,     0.0000, 2,-2, 0,-2);
            AddSol(    -0.2500,    -0.2200,     0.0040,     0.0014, 1, 3, 0,-2);
            AddSol(    -3.9960,     0.0000,     0.0000,     0.0004, 2, 0, 2, 0);
            AddSol(     0.5570,    -0.7500,     0.0000,    -0.0090, 2, 0, 2,-2);
            AddSol(    -0.4590,    -0.3800,     0.0000,    -0.0053, 2, 0,-2, 2);
            AddSol(    -1.2980,     0.7400,     0.0000,     0.0004, 2, 0,-2, 0);
            AddSol(     0.5380,     1.1400,     0.0000,    -0.0141, 2, 0,-2,-2);
            AddSol(     0.2630,     0.0200,     0.0000,     0.0000, 1, 1, 2, 0);
            AddSol(     0.4260,     0.0700,     0.0000,    -0.0006, 1, 1,-2,-2);
            AddSol(    -0.3040,     0.0300,     0.0000,     0.0003, 1,-1, 2, 0);
            AddSol(    -0.3720,    -0.1900,     0.0000,    -0.0027, 1,-1,-2, 2);
            AddSol(     0.4180,     0.0000,     0.0000,     0.0000, 0, 0, 4, 0);
            AddSol(    -0.3300,    -0.0400,     0.0000,     0.0000, 3, 0, 2, 0);

            SolarN();
            Planetary();
            S = F + DS/Astronomy.ARC;

            double lat_seconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*Math.Sin(S)-6.24*Math.Sin(3*S) + N;

            return new MoonResult(
                Astronomy.PI2 * Frac((L0+DLAM/Astronomy.ARC) / Astronomy.PI2),
                lat_seconds * (Astronomy.DEG2RAD / 3600.0),
                (Astronomy.ARC * Astronomy.EARTH_EQUATORIAL_RADIUS_AU) / (0.999953253 * SINPI)
            );
        }
    }

    internal struct MoonResult
    {
        public readonly double geo_eclip_lon;
        public readonly double geo_eclip_lat;
        public readonly double distance_au;

        public MoonResult(double lon, double lat, double dist)
        {
            this.geo_eclip_lon = lon;
            this.geo_eclip_lat = lat;
            this.distance_au = dist;
        }
    }

    /// <summary>
    /// Reports the constellation that a given celestial point lies within.
    /// </summary>
    /// <remarks>
    /// The #Astronomy.Constellation function returns this struct
    /// to report which constellation corresponds with a given point in the sky.
    /// Constellations are defined with respect to the B1875 equatorial system
    /// per IAU standard. Although `Astronomy.Constellation` requires J2000 equatorial
    /// coordinates, the struct contains converted B1875 coordinates for reference.
    /// </remarks>
    public struct ConstellationInfo
    {
        /// <summary>
        /// 3-character mnemonic symbol for the constellation, e.g. "Ori".
        /// </summary>
        public readonly string Symbol;

        /// <summary>
        /// Full name of constellation, e.g. "Orion".
        /// </summary>
        public readonly string Name;

        /// <summary>
        /// Right ascension expressed in B1875 coordinates.
        /// </summary>
        public readonly double Ra1875;

        /// <summary>
        /// Declination expressed in B1875 coordinates.
        /// </summary>
        public readonly double Dec1875;

        internal ConstellationInfo(string symbol, string name, double ra1875, double dec1875)
        {
            this.Symbol = symbol;
            this.Name = name;
            this.Ra1875 = ra1875;
            this.Dec1875 = dec1875;
        }
    }

    /// <summary>
    /// The wrapper class that holds Astronomy Engine functions.
    /// </summary>
    public static class Astronomy
    {
        private const double DAYS_PER_TROPICAL_YEAR = 365.24217;
        internal const double DEG2RAD = 0.017453292519943296;
        internal const double RAD2DEG = 57.295779513082321;
        private const double ASEC360 = 1296000.0;
        private const double ASEC2RAD = 4.848136811095359935899141e-6;
        internal const double PI2 = 2.0 * Math.PI;
        internal const double ARC = 3600.0 * 180.0 / Math.PI;       /* arcseconds per radian */
        private const double C_AUDAY = 173.1446326846693;           /* speed of light in AU/day */
        internal const double KM_PER_AU = 1.4959787069098932e+8;

        internal const double SUN_RADIUS_KM  = 695700.0;
        internal const double SUN_RADIUS_AU  = SUN_RADIUS_KM / KM_PER_AU;

        internal const double EARTH_FLATTENING = 0.996647180302104;
        internal const double EARTH_EQUATORIAL_RADIUS_KM = 6378.1366;
        internal const double EARTH_EQUATORIAL_RADIUS_AU = EARTH_EQUATORIAL_RADIUS_KM / KM_PER_AU;
        internal const double EARTH_MEAN_RADIUS_KM = 6371.0;    /* mean radius of the Earth's geoid, without atmosphere */
        internal const double EARTH_ATMOSPHERE_KM = 88.0;       /* effective atmosphere thickness for lunar eclipses */
        internal const double EARTH_ECLIPSE_RADIUS_KM = EARTH_MEAN_RADIUS_KM + EARTH_ATMOSPHERE_KM;

        internal const double MOON_EQUATORIAL_RADIUS_KM = 1738.1;
        internal const double MOON_MEAN_RADIUS_KM       = 1737.4;
        internal const double MOON_POLAR_RADIUS_KM      = 1736.0;
        internal const double MOON_EQUATORIAL_RADIUS_AU = (MOON_EQUATORIAL_RADIUS_KM / KM_PER_AU);

        private const double ANGVEL = 7.2921150e-5;
        private const double SECONDS_PER_DAY = 24.0 * 3600.0;
        private const double SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592;
        private const double MEAN_SYNODIC_MONTH = 29.530588;     /* average number of days for Moon to return to the same phase */
        private const double EARTH_ORBITAL_PERIOD = 365.256;
        private const double NEPTUNE_ORBITAL_PERIOD = 60189.0;
        internal const double REFRACTION_NEAR_HORIZON = 34.0 / 60.0;   /* degrees of refractive "lift" seen for objects near horizon */
        private const double ASEC180 = 180.0 * 60.0 * 60.0;         /* arcseconds per 180 degrees (or pi radians) */
        private const double AU_PER_PARSEC = (ASEC180 / Math.PI);   /* exact definition of how many AU = one parsec */
        private const double EARTH_MOON_MASS_RATIO = 81.30056;
        private const double SUN_MASS     = 333054.25318;        /* Sun's mass relative to Earth. */
        private const double JUPITER_MASS =    317.84997;        /* Jupiter's mass relative to Earth. */
        private const double SATURN_MASS  =     95.16745;        /* Saturn's mass relative to Earth. */
        private const double URANUS_MASS  =     14.53617;        /* Uranus's mass relative to Earth. */
        private const double NEPTUNE_MASS =     17.14886;        /* Neptune's mass relative to Earth. */

        /// <summary>Counter used for performance testing.</summary>
        public static int CalcMoonCount;

        internal static double LongitudeOffset(double diff)
        {
            double offset = diff;

            while (offset <= -180.0)
                offset += 360.0;

            while (offset > 180.0)
                offset -= 360.0;

            return offset;
        }

        internal static double NormalizeLongitude(double lon)
        {
            while (lon < 0.0)
                lon += 360.0;

            while (lon >= 360.0)
                lon -= 360.0;

            return lon;
        }


        private struct vsop_term_t
        {
            public double amplitude;
            public double phase;
            public double frequency;

            public vsop_term_t(double amplitude, double phase, double frequency)
            {
                this.amplitude = amplitude;
                this.phase = phase;
                this.frequency = frequency;
            }
        }

        private struct vsop_series_t
        {
            public vsop_term_t[] term;

            public vsop_series_t(vsop_term_t[] term)
            {
                this.term = term;
            }
        }

        private struct vsop_formula_t
        {
            public vsop_series_t[] series;

            public vsop_formula_t(vsop_series_t[] series)
            {
                this.series = series;
            }
        }

        private struct vsop_model_t
        {
            public vsop_formula_t lat;
            public vsop_formula_t lon;
            public vsop_formula_t rad;

            public vsop_model_t(vsop_series_t[] lat, vsop_series_t[] lon, vsop_series_t[] rad)
            {
                this.lat = new vsop_formula_t(lat);
                this.lon = new vsop_formula_t(lon);
                this.rad = new vsop_formula_t(rad);
            }
        };

        private static readonly vsop_term_t[] vsop_lat_Mercury_0 = new vsop_term_t[]
        {
            new vsop_term_t(4.40250710144, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.40989414977, 1.48302034195, 26087.90314157420),
            new vsop_term_t(0.05046294200, 4.47785489551, 52175.80628314840),
            new vsop_term_t(0.00855346844, 1.16520322459, 78263.70942472259),
            new vsop_term_t(0.00165590362, 4.11969163423, 104351.61256629678),
            new vsop_term_t(0.00034561897, 0.77930768443, 130439.51570787099),
            new vsop_term_t(0.00007583476, 3.71348404924, 156527.41884944518)
        };

        private static readonly vsop_term_t[] vsop_lat_Mercury_1 = new vsop_term_t[]
        {
            new vsop_term_t(26087.90313685529, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.01131199811, 6.21874197797, 26087.90314157420),
            new vsop_term_t(0.00292242298, 3.04449355541, 52175.80628314840),
            new vsop_term_t(0.00075775081, 6.08568821653, 78263.70942472259),
            new vsop_term_t(0.00019676525, 2.80965111777, 104351.61256629678)
        };

        private static readonly vsop_series_t[] vsop_lat_Mercury = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Mercury_0),
            new vsop_series_t(vsop_lat_Mercury_1)
        };

        private static readonly vsop_term_t[] vsop_lon_Mercury_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.11737528961, 1.98357498767, 26087.90314157420),
            new vsop_term_t(0.02388076996, 5.03738959686, 52175.80628314840),
            new vsop_term_t(0.01222839532, 3.14159265359, 0.00000000000),
            new vsop_term_t(0.00543251810, 1.79644363964, 78263.70942472259),
            new vsop_term_t(0.00129778770, 4.83232503958, 104351.61256629678),
            new vsop_term_t(0.00031866927, 1.58088495658, 130439.51570787099),
            new vsop_term_t(0.00007963301, 4.60972126127, 156527.41884944518)
        };

        private static readonly vsop_term_t[] vsop_lon_Mercury_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00274646065, 3.95008450011, 26087.90314157420),
            new vsop_term_t(0.00099737713, 3.14159265359, 0.00000000000)
        };

        private static readonly vsop_series_t[] vsop_lon_Mercury = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Mercury_0),
            new vsop_series_t(vsop_lon_Mercury_1)
        };

        private static readonly vsop_term_t[] vsop_rad_Mercury_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.39528271651, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.07834131818, 6.19233722598, 26087.90314157420),
            new vsop_term_t(0.00795525558, 2.95989690104, 52175.80628314840),
            new vsop_term_t(0.00121281764, 6.01064153797, 78263.70942472259),
            new vsop_term_t(0.00021921969, 2.77820093972, 104351.61256629678),
            new vsop_term_t(0.00004354065, 5.82894543774, 130439.51570787099)
        };

        private static readonly vsop_term_t[] vsop_rad_Mercury_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00217347740, 4.65617158665, 26087.90314157420),
            new vsop_term_t(0.00044141826, 1.42385544001, 52175.80628314840)
        };

        private static readonly vsop_series_t[] vsop_rad_Mercury = new vsop_series_t[]
        {
            new vsop_series_t(vsop_rad_Mercury_0),
            new vsop_series_t(vsop_rad_Mercury_1)
        };


        private static readonly vsop_term_t[] vsop_lat_Venus_0 = new vsop_term_t[]
        {
            new vsop_term_t(3.17614666774, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.01353968419, 5.59313319619, 10213.28554621100),
            new vsop_term_t(0.00089891645, 5.30650047764, 20426.57109242200),
            new vsop_term_t(0.00005477194, 4.41630661466, 7860.41939243920),
            new vsop_term_t(0.00003455741, 2.69964447820, 11790.62908865880),
            new vsop_term_t(0.00002372061, 2.99377542079, 3930.20969621960),
            new vsop_term_t(0.00001317168, 5.18668228402, 26.29831979980),
            new vsop_term_t(0.00001664146, 4.25018630147, 1577.34354244780),
            new vsop_term_t(0.00001438387, 4.15745084182, 9683.59458111640),
            new vsop_term_t(0.00001200521, 6.15357116043, 30639.85663863300)
        };

        private static readonly vsop_term_t[] vsop_lat_Venus_1 = new vsop_term_t[]
        {
            new vsop_term_t(10213.28554621638, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00095617813, 2.46406511110, 10213.28554621100),
            new vsop_term_t(0.00007787201, 0.62478482220, 20426.57109242200)
        };

        private static readonly vsop_series_t[] vsop_lat_Venus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Venus_0),
            new vsop_series_t(vsop_lat_Venus_1)
        };

        private static readonly vsop_term_t[] vsop_lon_Venus_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.05923638472, 0.26702775812, 10213.28554621100),
            new vsop_term_t(0.00040107978, 1.14737178112, 20426.57109242200),
            new vsop_term_t(0.00032814918, 3.14159265359, 0.00000000000)
        };

        private static readonly vsop_term_t[] vsop_lon_Venus_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00287821243, 1.88964962838, 10213.28554621100)
        };

        private static readonly vsop_series_t[] vsop_lon_Venus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Venus_0),
            new vsop_series_t(vsop_lon_Venus_1)
        };

        private static readonly vsop_term_t[] vsop_rad_Venus_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.72334820891, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00489824182, 4.02151831717, 10213.28554621100),
            new vsop_term_t(0.00001658058, 4.90206728031, 20426.57109242200),
            new vsop_term_t(0.00001378043, 1.12846591367, 11790.62908865880),
            new vsop_term_t(0.00001632096, 2.84548795207, 7860.41939243920),
            new vsop_term_t(0.00000498395, 2.58682193892, 9683.59458111640),
            new vsop_term_t(0.00000221985, 2.01346696541, 19367.18916223280),
            new vsop_term_t(0.00000237454, 2.55136053886, 15720.83878487840)
        };

        private static readonly vsop_term_t[] vsop_rad_Venus_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00034551041, 0.89198706276, 10213.28554621100)
        };

        private static readonly vsop_series_t[] vsop_rad_Venus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_rad_Venus_0),
            new vsop_series_t(vsop_rad_Venus_1)
        };


        private static readonly vsop_term_t[] vsop_lat_Earth_0 = new vsop_term_t[]
        {
            new vsop_term_t(1.75347045673, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.03341656453, 4.66925680415, 6283.07584999140),
            new vsop_term_t(0.00034894275, 4.62610242189, 12566.15169998280),
            new vsop_term_t(0.00003417572, 2.82886579754, 3.52311834900),
            new vsop_term_t(0.00003497056, 2.74411783405, 5753.38488489680),
            new vsop_term_t(0.00003135899, 3.62767041756, 77713.77146812050),
            new vsop_term_t(0.00002676218, 4.41808345438, 7860.41939243920),
            new vsop_term_t(0.00002342691, 6.13516214446, 3930.20969621960),
            new vsop_term_t(0.00001273165, 2.03709657878, 529.69096509460),
            new vsop_term_t(0.00001324294, 0.74246341673, 11506.76976979360),
            new vsop_term_t(0.00000901854, 2.04505446477, 26.29831979980),
            new vsop_term_t(0.00001199167, 1.10962946234, 1577.34354244780),
            new vsop_term_t(0.00000857223, 3.50849152283, 398.14900340820),
            new vsop_term_t(0.00000779786, 1.17882681962, 5223.69391980220),
            new vsop_term_t(0.00000990250, 5.23268072088, 5884.92684658320),
            new vsop_term_t(0.00000753141, 2.53339052847, 5507.55323866740),
            new vsop_term_t(0.00000505267, 4.58292599973, 18849.22754997420),
            new vsop_term_t(0.00000492392, 4.20505711826, 775.52261132400),
            new vsop_term_t(0.00000356672, 2.91954114478, 0.06731030280),
            new vsop_term_t(0.00000284125, 1.89869240932, 796.29800681640),
            new vsop_term_t(0.00000242879, 0.34481445893, 5486.77784317500),
            new vsop_term_t(0.00000317087, 5.84901948512, 11790.62908865880),
            new vsop_term_t(0.00000271112, 0.31486255375, 10977.07880469900),
            new vsop_term_t(0.00000206217, 4.80646631478, 2544.31441988340),
            new vsop_term_t(0.00000205478, 1.86953770281, 5573.14280143310),
            new vsop_term_t(0.00000202318, 2.45767790232, 6069.77675455340),
            new vsop_term_t(0.00000126225, 1.08295459501, 20.77539549240),
            new vsop_term_t(0.00000155516, 0.83306084617, 213.29909543800)
        };

        private static readonly vsop_term_t[] vsop_lat_Earth_1 = new vsop_term_t[]
        {
            new vsop_term_t(6283.07584999140, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00206058863, 2.67823455808, 6283.07584999140),
            new vsop_term_t(0.00004303419, 2.63512233481, 12566.15169998280)
        };

        private static readonly vsop_term_t[] vsop_lat_Earth_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00008721859, 1.07253635559, 6283.07584999140)
        };

        private static readonly vsop_series_t[] vsop_lat_Earth = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Earth_0),
            new vsop_series_t(vsop_lat_Earth_1),
            new vsop_series_t(vsop_lat_Earth_2)
        };

        private static readonly vsop_term_t[] vsop_lon_Earth_0 = new vsop_term_t[]
        {
        };

        private static readonly vsop_term_t[] vsop_lon_Earth_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00227777722, 3.41376620530, 6283.07584999140),
            new vsop_term_t(0.00003805678, 3.37063423795, 12566.15169998280)
        };

        private static readonly vsop_series_t[] vsop_lon_Earth = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Earth_0),
            new vsop_series_t(vsop_lon_Earth_1)
        };

        private static readonly vsop_term_t[] vsop_rad_Earth_0 = new vsop_term_t[]
        {
            new vsop_term_t(1.00013988784, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.01670699632, 3.09846350258, 6283.07584999140),
            new vsop_term_t(0.00013956024, 3.05524609456, 12566.15169998280),
            new vsop_term_t(0.00003083720, 5.19846674381, 77713.77146812050),
            new vsop_term_t(0.00001628463, 1.17387558054, 5753.38488489680),
            new vsop_term_t(0.00001575572, 2.84685214877, 7860.41939243920),
            new vsop_term_t(0.00000924799, 5.45292236722, 11506.76976979360),
            new vsop_term_t(0.00000542439, 4.56409151453, 3930.20969621960),
            new vsop_term_t(0.00000472110, 3.66100022149, 5884.92684658320),
            new vsop_term_t(0.00000085831, 1.27079125277, 161000.68573767410),
            new vsop_term_t(0.00000057056, 2.01374292245, 83996.84731811189),
            new vsop_term_t(0.00000055736, 5.24159799170, 71430.69561812909),
            new vsop_term_t(0.00000174844, 3.01193636733, 18849.22754997420),
            new vsop_term_t(0.00000243181, 4.27349530790, 11790.62908865880)
        };

        private static readonly vsop_term_t[] vsop_rad_Earth_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00103018607, 1.10748968172, 6283.07584999140),
            new vsop_term_t(0.00001721238, 1.06442300386, 12566.15169998280)
        };

        private static readonly vsop_term_t[] vsop_rad_Earth_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00004359385, 5.78455133808, 6283.07584999140)
        };

        private static readonly vsop_series_t[] vsop_rad_Earth = new vsop_series_t[]
        {
            new vsop_series_t(vsop_rad_Earth_0),
            new vsop_series_t(vsop_rad_Earth_1),
            new vsop_series_t(vsop_rad_Earth_2)
        };


        private static readonly vsop_term_t[] vsop_lat_Mars_0 = new vsop_term_t[]
        {
            new vsop_term_t(6.20347711581, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.18656368093, 5.05037100270, 3340.61242669980),
            new vsop_term_t(0.01108216816, 5.40099836344, 6681.22485339960),
            new vsop_term_t(0.00091798406, 5.75478744667, 10021.83728009940),
            new vsop_term_t(0.00027744987, 5.97049513147, 3.52311834900),
            new vsop_term_t(0.00010610235, 2.93958560338, 2281.23049651060),
            new vsop_term_t(0.00012315897, 0.84956094002, 2810.92146160520),
            new vsop_term_t(0.00008926784, 4.15697846427, 0.01725365220),
            new vsop_term_t(0.00008715691, 6.11005153139, 13362.44970679920),
            new vsop_term_t(0.00006797556, 0.36462229657, 398.14900340820),
            new vsop_term_t(0.00007774872, 3.33968761376, 5621.84292321040),
            new vsop_term_t(0.00003575078, 1.66186505710, 2544.31441988340),
            new vsop_term_t(0.00004161108, 0.22814971327, 2942.46342329160),
            new vsop_term_t(0.00003075252, 0.85696614132, 191.44826611160),
            new vsop_term_t(0.00002628117, 0.64806124465, 3337.08930835080),
            new vsop_term_t(0.00002937546, 6.07893711402, 0.06731030280),
            new vsop_term_t(0.00002389414, 5.03896442664, 796.29800681640),
            new vsop_term_t(0.00002579844, 0.02996736156, 3344.13554504880),
            new vsop_term_t(0.00001528141, 1.14979301996, 6151.53388830500),
            new vsop_term_t(0.00001798806, 0.65634057445, 529.69096509460),
            new vsop_term_t(0.00001264357, 3.62275122593, 5092.15195811580),
            new vsop_term_t(0.00001286228, 3.06796065034, 2146.16541647520),
            new vsop_term_t(0.00001546404, 2.91579701718, 1751.53953141600),
            new vsop_term_t(0.00001024902, 3.69334099279, 8962.45534991020),
            new vsop_term_t(0.00000891566, 0.18293837498, 16703.06213349900),
            new vsop_term_t(0.00000858759, 2.40093811940, 2914.01423582380),
            new vsop_term_t(0.00000832715, 2.46418619474, 3340.59517304760),
            new vsop_term_t(0.00000832720, 4.49495782139, 3340.62968035200),
            new vsop_term_t(0.00000712902, 3.66335473479, 1059.38193018920),
            new vsop_term_t(0.00000748723, 3.82248614017, 155.42039943420),
            new vsop_term_t(0.00000723861, 0.67497311481, 3738.76143010800),
            new vsop_term_t(0.00000635548, 2.92182225127, 8432.76438481560),
            new vsop_term_t(0.00000655162, 0.48864064125, 3127.31333126180),
            new vsop_term_t(0.00000550474, 3.81001042328, 0.98032106820),
            new vsop_term_t(0.00000552750, 4.47479317037, 1748.01641306700),
            new vsop_term_t(0.00000425966, 0.55364317304, 6283.07584999140),
            new vsop_term_t(0.00000415131, 0.49662285038, 213.29909543800),
            new vsop_term_t(0.00000472167, 3.62547124025, 1194.44701022460),
            new vsop_term_t(0.00000306551, 0.38052848348, 6684.74797174860),
            new vsop_term_t(0.00000312141, 0.99853944405, 6677.70173505060),
            new vsop_term_t(0.00000293198, 4.22131299634, 20.77539549240),
            new vsop_term_t(0.00000302375, 4.48618007156, 3532.06069281140),
            new vsop_term_t(0.00000274027, 0.54222167059, 3340.54511639700),
            new vsop_term_t(0.00000281079, 5.88163521788, 1349.86740965880),
            new vsop_term_t(0.00000231183, 1.28242156993, 3870.30339179440),
            new vsop_term_t(0.00000283602, 5.76885434940, 3149.16416058820),
            new vsop_term_t(0.00000236117, 5.75503217933, 3333.49887969900),
            new vsop_term_t(0.00000274033, 0.13372524985, 3340.67973700260),
            new vsop_term_t(0.00000299395, 2.78323740866, 6254.62666252360)
        };

        private static readonly vsop_term_t[] vsop_lat_Mars_1 = new vsop_term_t[]
        {
            new vsop_term_t(3340.61242700512, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.01457554523, 3.60433733236, 3340.61242669980),
            new vsop_term_t(0.00168414711, 3.92318567804, 6681.22485339960),
            new vsop_term_t(0.00020622975, 4.26108844583, 10021.83728009940),
            new vsop_term_t(0.00003452392, 4.73210393190, 3.52311834900),
            new vsop_term_t(0.00002586332, 4.60670058555, 13362.44970679920),
            new vsop_term_t(0.00000841535, 4.45864030426, 2281.23049651060)
        };

        private static readonly vsop_term_t[] vsop_lat_Mars_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00058152577, 2.04961712429, 3340.61242669980),
            new vsop_term_t(0.00013459579, 2.45738706163, 6681.22485339960)
        };

        private static readonly vsop_series_t[] vsop_lat_Mars = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Mars_0),
            new vsop_series_t(vsop_lat_Mars_1),
            new vsop_series_t(vsop_lat_Mars_2)
        };

        private static readonly vsop_term_t[] vsop_lon_Mars_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.03197134986, 3.76832042431, 3340.61242669980),
            new vsop_term_t(0.00298033234, 4.10616996305, 6681.22485339960),
            new vsop_term_t(0.00289104742, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00031365539, 4.44651053090, 10021.83728009940),
            new vsop_term_t(0.00003484100, 4.78812549260, 13362.44970679920)
        };

        private static readonly vsop_term_t[] vsop_lon_Mars_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00217310991, 6.04472194776, 3340.61242669980),
            new vsop_term_t(0.00020976948, 3.14159265359, 0.00000000000),
            new vsop_term_t(0.00012834709, 1.60810667915, 6681.22485339960)
        };

        private static readonly vsop_series_t[] vsop_lon_Mars = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Mars_0),
            new vsop_series_t(vsop_lon_Mars_1)
        };

        private static readonly vsop_term_t[] vsop_rad_Mars_0 = new vsop_term_t[]
        {
            new vsop_term_t(1.53033488271, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.14184953160, 3.47971283528, 3340.61242669980),
            new vsop_term_t(0.00660776362, 3.81783443019, 6681.22485339960),
            new vsop_term_t(0.00046179117, 4.15595316782, 10021.83728009940),
            new vsop_term_t(0.00008109733, 5.55958416318, 2810.92146160520),
            new vsop_term_t(0.00007485318, 1.77239078402, 5621.84292321040),
            new vsop_term_t(0.00005523191, 1.36436303770, 2281.23049651060),
            new vsop_term_t(0.00003825160, 4.49407183687, 13362.44970679920),
            new vsop_term_t(0.00002306537, 0.09081579001, 2544.31441988340),
            new vsop_term_t(0.00001999396, 5.36059617709, 3337.08930835080),
            new vsop_term_t(0.00002484394, 4.92545639920, 2942.46342329160),
            new vsop_term_t(0.00001960195, 4.74249437639, 3344.13554504880),
            new vsop_term_t(0.00001167119, 2.11260868341, 5092.15195811580),
            new vsop_term_t(0.00001102816, 5.00908403998, 398.14900340820),
            new vsop_term_t(0.00000899066, 4.40791133207, 529.69096509460),
            new vsop_term_t(0.00000992252, 5.83861961952, 6151.53388830500),
            new vsop_term_t(0.00000807354, 2.10217065501, 1059.38193018920),
            new vsop_term_t(0.00000797915, 3.44839203899, 796.29800681640),
            new vsop_term_t(0.00000740975, 1.49906336885, 2146.16541647520)
        };

        private static readonly vsop_term_t[] vsop_rad_Mars_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.01107433345, 2.03250524857, 3340.61242669980),
            new vsop_term_t(0.00103175887, 2.37071847807, 6681.22485339960),
            new vsop_term_t(0.00012877200, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00010815880, 2.70888095665, 10021.83728009940)
        };

        private static readonly vsop_term_t[] vsop_rad_Mars_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00044242249, 0.47930604954, 3340.61242669980),
            new vsop_term_t(0.00008138042, 0.86998389204, 6681.22485339960)
        };

        private static readonly vsop_series_t[] vsop_rad_Mars = new vsop_series_t[]
        {
            new vsop_series_t(vsop_rad_Mars_0),
            new vsop_series_t(vsop_rad_Mars_1),
            new vsop_series_t(vsop_rad_Mars_2)
        };


        private static readonly vsop_term_t[] vsop_lat_Jupiter_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.59954691494, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.09695898719, 5.06191793158, 529.69096509460),
            new vsop_term_t(0.00573610142, 1.44406205629, 7.11354700080),
            new vsop_term_t(0.00306389205, 5.41734730184, 1059.38193018920),
            new vsop_term_t(0.00097178296, 4.14264726552, 632.78373931320),
            new vsop_term_t(0.00072903078, 3.64042916389, 522.57741809380),
            new vsop_term_t(0.00064263975, 3.41145165351, 103.09277421860),
            new vsop_term_t(0.00039806064, 2.29376740788, 419.48464387520),
            new vsop_term_t(0.00038857767, 1.27231755835, 316.39186965660),
            new vsop_term_t(0.00027964629, 1.78454591820, 536.80451209540),
            new vsop_term_t(0.00013589730, 5.77481040790, 1589.07289528380),
            new vsop_term_t(0.00008246349, 3.58227925840, 206.18554843720),
            new vsop_term_t(0.00008768704, 3.63000308199, 949.17560896980),
            new vsop_term_t(0.00007368042, 5.08101194270, 735.87651353180),
            new vsop_term_t(0.00006263150, 0.02497628807, 213.29909543800),
            new vsop_term_t(0.00006114062, 4.51319998626, 1162.47470440780),
            new vsop_term_t(0.00004905396, 1.32084470588, 110.20632121940),
            new vsop_term_t(0.00005305285, 1.30671216791, 14.22709400160),
            new vsop_term_t(0.00005305441, 4.18625634012, 1052.26838318840),
            new vsop_term_t(0.00004647248, 4.69958103684, 3.93215326310),
            new vsop_term_t(0.00003045023, 4.31676431084, 426.59819087600),
            new vsop_term_t(0.00002609999, 1.56667394063, 846.08283475120),
            new vsop_term_t(0.00002028191, 1.06376530715, 3.18139373770),
            new vsop_term_t(0.00001764763, 2.14148655117, 1066.49547719000),
            new vsop_term_t(0.00001722972, 3.88036268267, 1265.56747862640),
            new vsop_term_t(0.00001920945, 0.97168196472, 639.89728631400),
            new vsop_term_t(0.00001633223, 3.58201833555, 515.46387109300),
            new vsop_term_t(0.00001431999, 4.29685556046, 625.67019231240),
            new vsop_term_t(0.00000973272, 4.09764549134, 95.97922721780)
        };

        private static readonly vsop_term_t[] vsop_lat_Jupiter_1 = new vsop_term_t[]
        {
            new vsop_term_t(529.69096508814, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00489503243, 4.22082939470, 529.69096509460),
            new vsop_term_t(0.00228917222, 6.02646855621, 7.11354700080),
            new vsop_term_t(0.00030099479, 4.54540782858, 1059.38193018920),
            new vsop_term_t(0.00020720920, 5.45943156902, 522.57741809380),
            new vsop_term_t(0.00012103653, 0.16994816098, 536.80451209540),
            new vsop_term_t(0.00006067987, 4.42422292017, 103.09277421860),
            new vsop_term_t(0.00005433968, 3.98480737746, 419.48464387520),
            new vsop_term_t(0.00004237744, 5.89008707199, 14.22709400160)
        };

        private static readonly vsop_term_t[] vsop_lat_Jupiter_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00047233601, 4.32148536482, 7.11354700080),
            new vsop_term_t(0.00030649436, 2.92977788700, 529.69096509460),
            new vsop_term_t(0.00014837605, 3.14159265359, 0.00000000000)
        };

        private static readonly vsop_series_t[] vsop_lat_Jupiter = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Jupiter_0),
            new vsop_series_t(vsop_lat_Jupiter_1),
            new vsop_series_t(vsop_lat_Jupiter_2)
        };

        private static readonly vsop_term_t[] vsop_lon_Jupiter_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.02268615702, 3.55852606721, 529.69096509460),
            new vsop_term_t(0.00109971634, 3.90809347197, 1059.38193018920),
            new vsop_term_t(0.00110090358, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00008101428, 3.60509572885, 522.57741809380),
            new vsop_term_t(0.00006043996, 4.25883108339, 1589.07289528380),
            new vsop_term_t(0.00006437782, 0.30627119215, 536.80451209540)
        };

        private static readonly vsop_term_t[] vsop_lon_Jupiter_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00078203446, 1.52377859742, 529.69096509460)
        };

        private static readonly vsop_series_t[] vsop_lon_Jupiter = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Jupiter_0),
            new vsop_series_t(vsop_lon_Jupiter_1)
        };

        private static readonly vsop_term_t[] vsop_rad_Jupiter_0 = new vsop_term_t[]
        {
            new vsop_term_t(5.20887429326, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.25209327119, 3.49108639871, 529.69096509460),
            new vsop_term_t(0.00610599976, 3.84115365948, 1059.38193018920),
            new vsop_term_t(0.00282029458, 2.57419881293, 632.78373931320),
            new vsop_term_t(0.00187647346, 2.07590383214, 522.57741809380),
            new vsop_term_t(0.00086792905, 0.71001145545, 419.48464387520),
            new vsop_term_t(0.00072062974, 0.21465724607, 536.80451209540),
            new vsop_term_t(0.00065517248, 5.97995884790, 316.39186965660),
            new vsop_term_t(0.00029134542, 1.67759379655, 103.09277421860),
            new vsop_term_t(0.00030135335, 2.16132003734, 949.17560896980),
            new vsop_term_t(0.00023453271, 3.54023522184, 735.87651353180),
            new vsop_term_t(0.00022283743, 4.19362594399, 1589.07289528380),
            new vsop_term_t(0.00023947298, 0.27458037480, 7.11354700080),
            new vsop_term_t(0.00013032614, 2.96042965363, 1162.47470440780),
            new vsop_term_t(0.00009703360, 1.90669633585, 206.18554843720),
            new vsop_term_t(0.00012749023, 2.71550286592, 1052.26838318840),
            new vsop_term_t(0.00007057931, 2.18184839926, 1265.56747862640),
            new vsop_term_t(0.00006137703, 6.26418240033, 846.08283475120),
            new vsop_term_t(0.00002616976, 2.00994012876, 1581.95934828300)
        };

        private static readonly vsop_term_t[] vsop_rad_Jupiter_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.01271801520, 2.64937512894, 529.69096509460),
            new vsop_term_t(0.00061661816, 3.00076460387, 1059.38193018920),
            new vsop_term_t(0.00053443713, 3.89717383175, 522.57741809380),
            new vsop_term_t(0.00031185171, 4.88276958012, 536.80451209540),
            new vsop_term_t(0.00041390269, 0.00000000000, 0.00000000000)
        };

        private static readonly vsop_series_t[] vsop_rad_Jupiter = new vsop_series_t[]
        {
            new vsop_series_t(vsop_rad_Jupiter_0),
            new vsop_series_t(vsop_rad_Jupiter_1)
        };


        private static readonly vsop_term_t[] vsop_lat_Saturn_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.87401354025, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.11107659762, 3.96205090159, 213.29909543800),
            new vsop_term_t(0.01414150957, 4.58581516874, 7.11354700080),
            new vsop_term_t(0.00398379389, 0.52112032699, 206.18554843720),
            new vsop_term_t(0.00350769243, 3.30329907896, 426.59819087600),
            new vsop_term_t(0.00206816305, 0.24658372002, 103.09277421860),
            new vsop_term_t(0.00079271300, 3.84007056878, 220.41264243880),
            new vsop_term_t(0.00023990355, 4.66976924553, 110.20632121940),
            new vsop_term_t(0.00016573588, 0.43719228296, 419.48464387520),
            new vsop_term_t(0.00014906995, 5.76903183869, 316.39186965660),
            new vsop_term_t(0.00015820290, 0.93809155235, 632.78373931320),
            new vsop_term_t(0.00014609559, 1.56518472000, 3.93215326310),
            new vsop_term_t(0.00013160301, 4.44891291899, 14.22709400160),
            new vsop_term_t(0.00015053543, 2.71669915667, 639.89728631400),
            new vsop_term_t(0.00013005299, 5.98119023644, 11.04570026390),
            new vsop_term_t(0.00010725067, 3.12939523827, 202.25339517410),
            new vsop_term_t(0.00005863206, 0.23656938524, 529.69096509460),
            new vsop_term_t(0.00005227757, 4.20783365759, 3.18139373770),
            new vsop_term_t(0.00006126317, 1.76328667907, 277.03499374140),
            new vsop_term_t(0.00005019687, 3.17787728405, 433.71173787680),
            new vsop_term_t(0.00004592550, 0.61977744975, 199.07200143640),
            new vsop_term_t(0.00004005867, 2.24479718502, 63.73589830340),
            new vsop_term_t(0.00002953796, 0.98280366998, 95.97922721780),
            new vsop_term_t(0.00003873670, 3.22283226966, 138.51749687070),
            new vsop_term_t(0.00002461186, 2.03163875071, 735.87651353180),
            new vsop_term_t(0.00003269484, 0.77492638211, 949.17560896980),
            new vsop_term_t(0.00001758145, 3.26580109940, 522.57741809380),
            new vsop_term_t(0.00001640172, 5.50504453050, 846.08283475120),
            new vsop_term_t(0.00001391327, 4.02333150505, 323.50541665740),
            new vsop_term_t(0.00001580648, 4.37265307169, 309.27832265580),
            new vsop_term_t(0.00001123498, 2.83726798446, 415.55249061210),
            new vsop_term_t(0.00001017275, 3.71700135395, 227.52618943960),
            new vsop_term_t(0.00000848642, 3.19150170830, 209.36694217490)
        };

        private static readonly vsop_term_t[] vsop_lat_Saturn_1 = new vsop_term_t[]
        {
            new vsop_term_t(213.29909521690, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.01297370862, 1.82834923978, 213.29909543800),
            new vsop_term_t(0.00564345393, 2.88499717272, 7.11354700080),
            new vsop_term_t(0.00093734369, 1.06311793502, 426.59819087600),
            new vsop_term_t(0.00107674962, 2.27769131009, 206.18554843720),
            new vsop_term_t(0.00040244455, 2.04108104671, 220.41264243880),
            new vsop_term_t(0.00019941774, 1.27954390470, 103.09277421860),
            new vsop_term_t(0.00010511678, 2.74880342130, 14.22709400160),
            new vsop_term_t(0.00006416106, 0.38238295041, 639.89728631400),
            new vsop_term_t(0.00004848994, 2.43037610229, 419.48464387520),
            new vsop_term_t(0.00004056892, 2.92133209468, 110.20632121940),
            new vsop_term_t(0.00003768635, 3.64965330780, 3.93215326310)
        };

        private static readonly vsop_term_t[] vsop_lat_Saturn_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00116441330, 1.17988132879, 7.11354700080),
            new vsop_term_t(0.00091841837, 0.07325195840, 213.29909543800),
            new vsop_term_t(0.00036661728, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00015274496, 4.06493179167, 206.18554843720)
        };

        private static readonly vsop_series_t[] vsop_lat_Saturn = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Saturn_0),
            new vsop_series_t(vsop_lat_Saturn_1),
            new vsop_series_t(vsop_lat_Saturn_2)
        };

        private static readonly vsop_term_t[] vsop_lon_Saturn_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.04330678039, 3.60284428399, 213.29909543800),
            new vsop_term_t(0.00240348302, 2.85238489373, 426.59819087600),
            new vsop_term_t(0.00084745939, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00030863357, 3.48441504555, 220.41264243880),
            new vsop_term_t(0.00034116062, 0.57297307557, 206.18554843720),
            new vsop_term_t(0.00014734070, 2.11846596715, 639.89728631400),
            new vsop_term_t(0.00009916667, 5.79003188904, 419.48464387520),
            new vsop_term_t(0.00006993564, 4.73604689720, 7.11354700080),
            new vsop_term_t(0.00004807588, 5.43305312061, 316.39186965660)
        };

        private static readonly vsop_term_t[] vsop_lon_Saturn_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00198927992, 4.93901017903, 213.29909543800),
            new vsop_term_t(0.00036947916, 3.14159265359, 0.00000000000),
            new vsop_term_t(0.00017966989, 0.51979431110, 426.59819087600)
        };

        private static readonly vsop_series_t[] vsop_lon_Saturn = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Saturn_0),
            new vsop_series_t(vsop_lon_Saturn_1)
        };

        private static readonly vsop_term_t[] vsop_rad_Saturn_0 = new vsop_term_t[]
        {
            new vsop_term_t(9.55758135486, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.52921382865, 2.39226219573, 213.29909543800),
            new vsop_term_t(0.01873679867, 5.23549604660, 206.18554843720),
            new vsop_term_t(0.01464663929, 1.64763042902, 426.59819087600),
            new vsop_term_t(0.00821891141, 5.93520042303, 316.39186965660),
            new vsop_term_t(0.00547506923, 5.01532618980, 103.09277421860),
            new vsop_term_t(0.00371684650, 2.27114821115, 220.41264243880),
            new vsop_term_t(0.00361778765, 3.13904301847, 7.11354700080),
            new vsop_term_t(0.00140617506, 5.70406606781, 632.78373931320),
            new vsop_term_t(0.00108974848, 3.29313390175, 110.20632121940),
            new vsop_term_t(0.00069006962, 5.94099540992, 419.48464387520),
            new vsop_term_t(0.00061053367, 0.94037691801, 639.89728631400),
            new vsop_term_t(0.00048913294, 1.55733638681, 202.25339517410),
            new vsop_term_t(0.00034143772, 0.19519102597, 277.03499374140),
            new vsop_term_t(0.00032401773, 5.47084567016, 949.17560896980),
            new vsop_term_t(0.00020936596, 0.46349251129, 735.87651353180),
            new vsop_term_t(0.00009796004, 5.20477537945, 1265.56747862640),
            new vsop_term_t(0.00011993338, 5.98050967385, 846.08283475120),
            new vsop_term_t(0.00020839300, 1.52102476129, 433.71173787680),
            new vsop_term_t(0.00015298404, 3.05943814940, 529.69096509460),
            new vsop_term_t(0.00006465823, 0.17732249942, 1052.26838318840),
            new vsop_term_t(0.00011380257, 1.73105427040, 522.57741809380),
            new vsop_term_t(0.00003419618, 4.94550542171, 1581.95934828300)
        };

        private static readonly vsop_term_t[] vsop_rad_Saturn_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.06182981340, 0.25843511480, 213.29909543800),
            new vsop_term_t(0.00506577242, 0.71114625261, 206.18554843720),
            new vsop_term_t(0.00341394029, 5.79635741658, 426.59819087600),
            new vsop_term_t(0.00188491195, 0.47215589652, 220.41264243880),
            new vsop_term_t(0.00186261486, 3.14159265359, 0.00000000000),
            new vsop_term_t(0.00143891146, 1.40744822888, 7.11354700080)
        };

        private static readonly vsop_term_t[] vsop_rad_Saturn_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00436902572, 4.78671677509, 213.29909543800)
        };

        private static readonly vsop_series_t[] vsop_rad_Saturn = new vsop_series_t[]
        {
            new vsop_series_t(vsop_rad_Saturn_0),
            new vsop_series_t(vsop_rad_Saturn_1),
            new vsop_series_t(vsop_rad_Saturn_2)
        };


        private static readonly vsop_term_t[] vsop_lat_Uranus_0 = new vsop_term_t[]
        {
            new vsop_term_t(5.48129294297, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.09260408234, 0.89106421507, 74.78159856730),
            new vsop_term_t(0.01504247898, 3.62719260920, 1.48447270830),
            new vsop_term_t(0.00365981674, 1.89962179044, 73.29712585900),
            new vsop_term_t(0.00272328168, 3.35823706307, 149.56319713460),
            new vsop_term_t(0.00070328461, 5.39254450063, 63.73589830340),
            new vsop_term_t(0.00068892678, 6.09292483287, 76.26607127560),
            new vsop_term_t(0.00061998615, 2.26952066061, 2.96894541660),
            new vsop_term_t(0.00061950719, 2.85098872691, 11.04570026390),
            new vsop_term_t(0.00026468770, 3.14152083966, 71.81265315070),
            new vsop_term_t(0.00025710476, 6.11379840493, 454.90936652730),
            new vsop_term_t(0.00021078850, 4.36059339067, 148.07872442630),
            new vsop_term_t(0.00017818647, 1.74436930289, 36.64856292950),
            new vsop_term_t(0.00014613507, 4.73732166022, 3.93215326310),
            new vsop_term_t(0.00011162509, 5.82681796350, 224.34479570190),
            new vsop_term_t(0.00010997910, 0.48865004018, 138.51749687070),
            new vsop_term_t(0.00009527478, 2.95516862826, 35.16409022120),
            new vsop_term_t(0.00007545601, 5.23626582400, 109.94568878850),
            new vsop_term_t(0.00004220241, 3.23328220918, 70.84944530420),
            new vsop_term_t(0.00004051900, 2.27755017300, 151.04766984290),
            new vsop_term_t(0.00003354596, 1.06549007380, 4.45341812490),
            new vsop_term_t(0.00002926718, 4.62903718891, 9.56122755560),
            new vsop_term_t(0.00003490340, 5.48306144511, 146.59425171800),
            new vsop_term_t(0.00003144069, 4.75199570434, 77.75054398390),
            new vsop_term_t(0.00002922333, 5.35235361027, 85.82729883120),
            new vsop_term_t(0.00002272788, 4.36600400036, 70.32818044240),
            new vsop_term_t(0.00002051219, 1.51773566586, 0.11187458460),
            new vsop_term_t(0.00002148602, 0.60745949945, 38.13303563780),
            new vsop_term_t(0.00001991643, 4.92437588682, 277.03499374140),
            new vsop_term_t(0.00001376226, 2.04283539351, 65.22037101170),
            new vsop_term_t(0.00001666902, 3.62744066769, 380.12776796000),
            new vsop_term_t(0.00001284107, 3.11347961505, 202.25339517410),
            new vsop_term_t(0.00001150429, 0.93343589092, 3.18139373770),
            new vsop_term_t(0.00001533221, 2.58594681212, 52.69019803950),
            new vsop_term_t(0.00001281604, 0.54271272721, 222.86032299360),
            new vsop_term_t(0.00001372139, 4.19641530878, 111.43016149680),
            new vsop_term_t(0.00001221029, 0.19900650030, 108.46121608020),
            new vsop_term_t(0.00000946181, 1.19253165736, 127.47179660680),
            new vsop_term_t(0.00001150989, 4.17898916639, 33.67961751290)
        };

        private static readonly vsop_term_t[] vsop_lat_Uranus_1 = new vsop_term_t[]
        {
            new vsop_term_t(74.78159860910, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00154332863, 5.24158770553, 74.78159856730),
            new vsop_term_t(0.00024456474, 1.71260334156, 1.48447270830),
            new vsop_term_t(0.00009258442, 0.42829732350, 11.04570026390),
            new vsop_term_t(0.00008265977, 1.50218091379, 63.73589830340),
            new vsop_term_t(0.00009150160, 1.41213765216, 149.56319713460)
        };

        private static readonly vsop_series_t[] vsop_lat_Uranus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Uranus_0),
            new vsop_series_t(vsop_lat_Uranus_1)
        };

        private static readonly vsop_term_t[] vsop_lon_Uranus_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.01346277648, 2.61877810547, 74.78159856730),
            new vsop_term_t(0.00062341400, 5.08111189648, 149.56319713460),
            new vsop_term_t(0.00061601196, 3.14159265359, 0.00000000000),
            new vsop_term_t(0.00009963722, 1.61603805646, 76.26607127560),
            new vsop_term_t(0.00009926160, 0.57630380333, 73.29712585900)
        };

        private static readonly vsop_term_t[] vsop_lon_Uranus_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00034101978, 0.01321929936, 74.78159856730)
        };

        private static readonly vsop_series_t[] vsop_lon_Uranus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Uranus_0),
            new vsop_series_t(vsop_lon_Uranus_1)
        };

        private static readonly vsop_term_t[] vsop_rad_Uranus_0 = new vsop_term_t[]
        {
            new vsop_term_t(19.21264847206, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.88784984413, 5.60377527014, 74.78159856730),
            new vsop_term_t(0.03440836062, 0.32836099706, 73.29712585900),
            new vsop_term_t(0.02055653860, 1.78295159330, 149.56319713460),
            new vsop_term_t(0.00649322410, 4.52247285911, 76.26607127560),
            new vsop_term_t(0.00602247865, 3.86003823674, 63.73589830340),
            new vsop_term_t(0.00496404167, 1.40139935333, 454.90936652730),
            new vsop_term_t(0.00338525369, 1.58002770318, 138.51749687070),
            new vsop_term_t(0.00243509114, 1.57086606044, 71.81265315070),
            new vsop_term_t(0.00190522303, 1.99809394714, 1.48447270830),
            new vsop_term_t(0.00161858838, 2.79137786799, 148.07872442630),
            new vsop_term_t(0.00143706183, 1.38368544947, 11.04570026390),
            new vsop_term_t(0.00093192405, 0.17437220467, 36.64856292950),
            new vsop_term_t(0.00071424548, 4.24509236074, 224.34479570190),
            new vsop_term_t(0.00089806014, 3.66105364565, 109.94568878850),
            new vsop_term_t(0.00039009723, 1.66971401684, 70.84944530420),
            new vsop_term_t(0.00046677296, 1.39976401694, 35.16409022120),
            new vsop_term_t(0.00039025624, 3.36234773834, 277.03499374140),
            new vsop_term_t(0.00036755274, 3.88649278513, 146.59425171800),
            new vsop_term_t(0.00030348723, 0.70100838798, 151.04766984290),
            new vsop_term_t(0.00029156413, 3.18056336700, 77.75054398390),
            new vsop_term_t(0.00022637073, 0.72518687029, 529.69096509460),
            new vsop_term_t(0.00011959076, 1.75043392140, 984.60033162190),
            new vsop_term_t(0.00025620756, 5.25656086672, 380.12776796000)
        };

        private static readonly vsop_term_t[] vsop_rad_Uranus_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.01479896629, 3.67205697578, 74.78159856730)
        };

        private static readonly vsop_series_t[] vsop_rad_Uranus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_rad_Uranus_0),
            new vsop_series_t(vsop_rad_Uranus_1)
        };


        private static readonly vsop_term_t[] vsop_lat_Neptune_0 = new vsop_term_t[]
        {
            new vsop_term_t(5.31188633046, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.01798475530, 2.90101273890, 38.13303563780),
            new vsop_term_t(0.01019727652, 0.48580922867, 1.48447270830),
            new vsop_term_t(0.00124531845, 4.83008090676, 36.64856292950),
            new vsop_term_t(0.00042064466, 5.41054993053, 2.96894541660),
            new vsop_term_t(0.00037714584, 6.09221808686, 35.16409022120),
            new vsop_term_t(0.00033784738, 1.24488874087, 76.26607127560),
            new vsop_term_t(0.00016482741, 0.00007727998, 491.55792945680),
            new vsop_term_t(0.00009198584, 4.93747051954, 39.61750834610),
            new vsop_term_t(0.00008994250, 0.27462171806, 175.16605980020)
        };

        private static readonly vsop_term_t[] vsop_lat_Neptune_1 = new vsop_term_t[]
        {
            new vsop_term_t(38.13303563957, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00016604172, 4.86323329249, 1.48447270830),
            new vsop_term_t(0.00015744045, 2.27887427527, 38.13303563780)
        };

        private static readonly vsop_series_t[] vsop_lat_Neptune = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Neptune_0),
            new vsop_series_t(vsop_lat_Neptune_1)
        };

        private static readonly vsop_term_t[] vsop_lon_Neptune_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.03088622933, 1.44104372644, 38.13303563780),
            new vsop_term_t(0.00027780087, 5.91271884599, 76.26607127560),
            new vsop_term_t(0.00027623609, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00015355489, 2.52123799551, 36.64856292950),
            new vsop_term_t(0.00015448133, 3.50877079215, 39.61750834610)
        };

        private static readonly vsop_series_t[] vsop_lon_Neptune = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Neptune_0)
        };

        private static readonly vsop_term_t[] vsop_rad_Neptune_0 = new vsop_term_t[]
        {
            new vsop_term_t(30.07013205828, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.27062259632, 1.32999459377, 38.13303563780),
            new vsop_term_t(0.01691764014, 3.25186135653, 36.64856292950),
            new vsop_term_t(0.00807830553, 5.18592878704, 1.48447270830),
            new vsop_term_t(0.00537760510, 4.52113935896, 35.16409022120),
            new vsop_term_t(0.00495725141, 1.57105641650, 491.55792945680),
            new vsop_term_t(0.00274571975, 1.84552258866, 175.16605980020),
            new vsop_term_t(0.00012012320, 1.92059384991, 1021.24889455140),
            new vsop_term_t(0.00121801746, 5.79754470298, 76.26607127560),
            new vsop_term_t(0.00100896068, 0.37702724930, 73.29712585900),
            new vsop_term_t(0.00135134092, 3.37220609835, 39.61750834610),
            new vsop_term_t(0.00007571796, 1.07149207335, 388.46515523820)
        };

        private static readonly vsop_series_t[] vsop_rad_Neptune = new vsop_series_t[]
        {
            new vsop_series_t(vsop_rad_Neptune_0)
        };



        private static readonly vsop_model_t[] vsop = new vsop_model_t[]
        {
            new vsop_model_t(vsop_lat_Mercury,  vsop_lon_Mercury,   vsop_rad_Mercury),
            new vsop_model_t(vsop_lat_Venus,    vsop_lon_Venus,     vsop_rad_Venus  ),
            new vsop_model_t(vsop_lat_Earth,    vsop_lon_Earth,     vsop_rad_Earth  ),
            new vsop_model_t(vsop_lat_Mars,     vsop_lon_Mars,      vsop_rad_Mars   ),
            new vsop_model_t(vsop_lat_Jupiter,  vsop_lon_Jupiter,   vsop_rad_Jupiter),
            new vsop_model_t(vsop_lat_Saturn,   vsop_lon_Saturn,    vsop_rad_Saturn ),
            new vsop_model_t(vsop_lat_Uranus,   vsop_lon_Uranus,    vsop_rad_Uranus ),
            new vsop_model_t(vsop_lat_Neptune,  vsop_lon_Neptune,   vsop_rad_Neptune)
        };

        /// <summary>
        /// The minimum year value supported by Astronomy Engine for calculating Pluto.
        /// All other bodies can be calculated for any year.
        /// </summary>
        public const int MinYear = 1700;

        /// <summary>
        /// The maximum year value supported by Astronomy Engine for calculating Pluto.
        /// All other bodies can be calculated for any year.
        /// </summary>
        public const int MaxYear = 2200;

        /// <summary>The default Delta T function used by Astronomy Engine.</summary>
        /// <remarks>
        /// Espenak and Meeus use a series of piecewise polynomials to
        /// approximate DeltaT of the Earth in their "Five Millennium Canon of Solar Eclipses".
        /// See: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
        /// This is the default Delta T function used by Astronomy Engine.
        /// </remarks>
        /// <param name="ut">The floating point number of days since noon UTC on January 1, 2000.</param>
        /// <returns>The estimated difference TT-UT on the given date, expressed in seconds.</returns>

        public static double DeltaT_EspenakMeeus(double ut)
        {
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
            double u, u2, u3, u4, u5, u6, u7;
            double y = 2000 + ((ut - 14) / DAYS_PER_TROPICAL_YEAR);
            if (y < -500)
            {
                u = (y - 1820)/100;
                return -20 + (32 * u*u);
            }
            if (y < 500)
            {
                u = y/100;
                u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3;
                return 10583.6 - 1014.41*u + 33.78311*u2 - 5.952053*u3 - 0.1798452*u4 + 0.022174192*u5 + 0.0090316521*u6;
            }
            if (y < 1600)
            {
                u = (y - 1000) / 100;
                u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3;
                return 1574.2 - 556.01*u + 71.23472*u2 + 0.319781*u3 - 0.8503463*u4 - 0.005050998*u5 + 0.0083572073*u6;
            }
            if (y < 1700)
            {
                u = y - 1600;
                u2 = u*u; u3 = u*u2;
                return 120 - 0.9808*u - 0.01532*u2 + u3/7129.0;
            }
            if (y < 1800)
            {
                u = y - 1700;
                u2 = u*u; u3 = u*u2; u4 = u2*u2;
                return 8.83 + 0.1603*u - 0.0059285*u2 + 0.00013336*u3 - u4/1174000;
            }
            if (y < 1860)
            {
                u = y - 1800;
                u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3; u6 = u3*u3; u7 = u3*u4;
                return 13.72 - 0.332447*u + 0.0068612*u2 + 0.0041116*u3 - 0.00037436*u4 + 0.0000121272*u5 - 0.0000001699*u6 + 0.000000000875*u7;
            }
            if (y < 1900)
            {
                u = y - 1860;
                u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3;
                return 7.62 + 0.5737*u - 0.251754*u2 + 0.01680668*u3 - 0.0004473624*u4 + u5/233174;
            }
            if (y < 1920)
            {
                u = y - 1900;
                u2 = u*u; u3 = u*u2; u4 = u2*u2;
                return -2.79 + 1.494119*u - 0.0598939*u2 + 0.0061966*u3 - 0.000197*u4;
            }
            if (y < 1941)
            {
                u = y - 1920;
                u2 = u*u; u3 = u*u2;
                return 21.20 + 0.84493*u - 0.076100*u2 + 0.0020936*u3;
            }
            if (y < 1961)
            {
                u = y - 1950;
                u2 = u*u; u3 = u*u2;
                return 29.07 + 0.407*u - u2/233 + u3/2547;
            }
            if (y < 1986)
            {
                u = y - 1975;
                u2 = u*u; u3 = u*u2;
                return 45.45 + 1.067*u - u2/260 - u3/718;
            }
            if (y < 2005)
            {
                u = y - 2000;
                u2 = u*u; u3 = u*u2; u4 = u2*u2; u5 = u2*u3;
                return 63.86 + 0.3345*u - 0.060374*u2 + 0.0017275*u3 + 0.000651814*u4 + 0.00002373599*u5;
            }
            if (y < 2050)
            {
                u = y - 2000;
                return 62.92 + 0.32217*u + 0.005589*u*u;
            }
            if (y < 2150)
            {
                u = (y - 1820) / 100;
                return -20.0 + 32.0*u*u - 0.5628*(2150 - y);
            }

            /* all years after 2150 */
            u = (y - 1820) / 100;
            return -20 + (32 * u*u);
        }

        private static DeltaTimeFunc DeltaT = DeltaT_EspenakMeeus;

        internal static double TerrestrialTime(double ut)
        {
            return ut + DeltaT(ut)/86400.0;
        }

        private static double VsopFormulaCalc(vsop_formula_t formula, double t)
        {
            double coord = 0.0;
            double tpower = 1.0;
            for (int s=0; s < formula.series.Length; ++s)
            {
                double sum = 0.0;
                vsop_series_t series = formula.series[s];
                for (int i=0; i < series.term.Length; ++i)
                {
                    vsop_term_t term = series.term[i];
                    sum += term.amplitude * Math.Cos(term.phase + (t * term.frequency));
                }
                coord += tpower * sum;
                tpower *= t;
            }
            return coord;
        }

        private static AstroVector CalcVsop(vsop_model_t model, AstroTime time)
        {
            double t = time.tt / 365250;    /* millennia since 2000 */

            /* Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates. */
            double sphere0 = VsopFormulaCalc(model.lat, t);
            double sphere1 = VsopFormulaCalc(model.lon, t);
            double sphere2 = VsopFormulaCalc(model.rad, t);

            /* Convert ecliptic spherical coordinates to ecliptic Cartesian coordinates. */
            double r_coslat = sphere2 * Math.Cos(sphere1);
            double eclip0 = r_coslat * Math.Cos(sphere0);
            double eclip1 = r_coslat * Math.Sin(sphere0);
            double eclip2 = sphere2 * Math.Sin(sphere1);

            /* Convert ecliptic Cartesian coordinates to equatorial Cartesian coordinates. */
            double x = eclip0 + 0.000000440360*eclip1 - 0.000000190919*eclip2;
            double y = -0.000000479966*eclip0 + 0.917482137087*eclip1 - 0.397776982902*eclip2;
            double z = 0.397776982902*eclip1 + 0.917482137087*eclip2;

            return new AstroVector(x, y, z, time);
        }

        private static double VsopHelioDistance(vsop_model_t model, AstroTime time)
        {
            /*
                The caller only wants to know the distance between the planet and the Sun.
                So we only need to calculate the radial component of the spherical coordinates.
            */

            double t = time.tt / 365250;    /* millennia since 2000 */
            return VsopFormulaCalc(model.rad, t);
        }

        // TOP2013 model for Pluto

        private struct astro_top_term_t
        {
            public double k;
            public double c;
            public double s;
            public astro_top_term_t(double k, double c, double s)
            {
                this.k = k;
                this.c = c;
                this.s = s;
            }
        }

        private struct top_elliptical_t
        {
            public double a;           /* AU */
            public double lambda;      /* rad */
            public double k;           /* 1 */
            public double h;           /* 1 */
            public double q;           /* 1 */
            public double p;           /* 1 */
        }

        private static readonly astro_top_term_t[][][] top_model_8 = new astro_top_term_t[][][]
        {
            new astro_top_term_t[][]  // f=0
            {
                new astro_top_term_t[]  // f=0, s= 0
                {
                    new astro_top_term_t(        0,  3.9544617144029999e+01,  0.0000000000000000e+00), // f=0, s= 0, t=   0
                    new astro_top_term_t(     1402, -1.8891373533434089e-01, -8.5258197635470073e-02), // f=0, s= 0, t=   1
                    new astro_top_term_t(     1331, -4.1495877833812339e-02, -3.3387415274886263e-02), // f=0, s= 0, t=   2
                    new astro_top_term_t(      522, -4.8502474919249819e-02, -7.3455412272547278e-03), // f=0, s= 0, t=   3
                    new astro_top_term_t(       71,  2.8102132918948221e-02, -5.3468437660152152e-03), // f=0, s= 0, t=   4
                    new astro_top_term_t(     1261, -9.1600251304608978e-03, -1.2309204554431390e-02), // f=0, s= 0, t=   5
                    new astro_top_term_t(      452, -1.2202161344831991e-02, -5.2519856329982890e-03), // f=0, s= 0, t=   6
                    new astro_top_term_t(     2875, -9.7845229475467185e-03, -7.2968560644718955e-04), // f=0, s= 0, t=   7
                    new astro_top_term_t(       35,  4.8494518209585983e-03, -6.8918970374425084e-03), // f=0, s= 0, t=   8
                    new astro_top_term_t(      141,  5.8271375488234932e-03, -3.1946778653436391e-03), // f=0, s= 0, t=   9
                    new astro_top_term_t(      137,  1.5300059509576150e-03, -6.0327729791525954e-03), // f=0, s= 0, t=  10
                    new astro_top_term_t(        4, -1.9494897408412360e-03,  4.5717130739994704e-03), // f=0, s= 0, t=  11
                    new astro_top_term_t(     1190, -1.7524220672664240e-03, -4.2980683502911454e-03), // f=0, s= 0, t=  12
                    new astro_top_term_t(      381, -3.1062775803702681e-03, -2.4728667551542258e-03), // f=0, s= 0, t=  13
                    new astro_top_term_t(        8, -1.7188663433411050e-03,  2.9756270077158122e-03), // f=0, s= 0, t=  14
                    new astro_top_term_t(     1543, -7.7472653128184826e-04,  2.6514626782777680e-03), // f=0, s= 0, t=  15
                    new astro_top_term_t(     1115, -1.5405722125111840e-03, -2.0778390548994150e-03), // f=0, s= 0, t=  16
                    new astro_top_term_t(     2804, -2.0048397209869230e-03, -4.1957951179189120e-04), // f=0, s= 0, t=  17
                    new astro_top_term_t(       67,  9.6850762192148931e-04, -1.5811913714969829e-03), // f=0, s= 0, t=  18
                    new astro_top_term_t(      212,  1.5466715821083480e-03, -9.6836654994834209e-04), // f=0, s= 0, t=  19
                    new astro_top_term_t(     1119, -1.8820367463891121e-04, -1.4293834479379090e-03), // f=0, s= 0, t=  20
                    new astro_top_term_t(    17405, -1.0738845199599739e-03,  9.5985010997943349e-04), // f=0, s= 0, t=  21
                    new astro_top_term_t(    28337,  7.5821211083786067e-04,  1.1416213940389449e-03), // f=0, s= 0, t=  22
                    new astro_top_term_t(      310, -6.9650370158153983e-04, -1.0024667762364200e-03), // f=0, s= 0, t=  23
                    new astro_top_term_t(     1044, -3.2801771454650589e-04, -6.4947140155397116e-04), // f=0, s= 0, t=  24
                    new astro_top_term_t(       63,  4.0424767075158291e-04,  5.7315886355325109e-04), // f=0, s= 0, t=  25
                    new astro_top_term_t(     1614, -1.3238448498587051e-04,  6.7949492229369074e-04), // f=0, s= 0, t=  26
                    new astro_top_term_t(       12, -6.5048480404874143e-04, -1.0430021074697129e-04), // f=0, s= 0, t=  27
                    new astro_top_term_t(      283,  3.3929768009878552e-04, -5.1573678840263412e-04), // f=0, s= 0, t=  28
                    new astro_top_term_t(      133,  4.0138197254613279e-04,  4.5284627712770291e-04), // f=0, s= 0, t=  29
                    new astro_top_term_t(     1421,  5.0823717785117468e-04,  3.1010622577270548e-04), // f=0, s= 0, t=  30
                    new astro_top_term_t(     1383, -5.6813906164891585e-04, -1.7225147327178090e-04), // f=0, s= 0, t=  31
                    new astro_top_term_t(     4348, -5.1094469101100064e-04,  1.4416513132178369e-04), // f=0, s= 0, t=  32
                    new astro_top_term_t(     2733, -4.8192307867155672e-04, -1.8631709444892481e-04), // f=0, s= 0, t=  33
                    new astro_top_term_t(      664,  6.3695948832436563e-05,  5.0429756322537705e-04), // f=0, s= 0, t=  34
                    new astro_top_term_t(      177, -3.0532744995821309e-04,  4.0226535349675440e-04), // f=0, s= 0, t=  35
                    new astro_top_term_t(      204,  3.4677834246970749e-04,  3.2588064363496143e-04), // f=0, s= 0, t=  36
                    new astro_top_term_t(     1048,  5.6219036569383140e-05, -4.5145044715373130e-04), // f=0, s= 0, t=  37
                    new astro_top_term_t(     1406,  3.6764565360298558e-04, -2.4793326161876619e-04), // f=0, s= 0, t=  38
                    new astro_top_term_t(     1398, -6.1399335668996843e-05, -4.0368084856225791e-04), // f=0, s= 0, t=  39
                    new astro_top_term_t(      880,  3.8206123841964649e-04,  1.0727867737651879e-04), // f=0, s= 0, t=  40
                    new astro_top_term_t(      239, -9.2482984334176681e-05, -3.6357760642322443e-04), // f=0, s= 0, t=  41
                    new astro_top_term_t(      541, -3.5492935880988093e-04, -7.2696664262252120e-05), // f=0, s= 0, t=  42
                    new astro_top_term_t(    17334, -3.1608420086033548e-04,  1.6282567766535830e-04), // f=0, s= 0, t=  43
                    new astro_top_term_t(      503,  3.4204624270103872e-04,  6.9800261361389600e-06), // f=0, s= 0, t=  44
                    new astro_top_term_t(    28266,  1.1026801272880990e-04,  3.1928911807394130e-04), // f=0, s= 0, t=  45
                    new astro_top_term_t(      974, -1.5403436032839399e-04, -2.6908883645786191e-04), // f=0, s= 0, t=  46
                    new astro_top_term_t(      345, -2.1738414747494940e-04,  1.5563797862269319e-04), // f=0, s= 0, t=  47
                    new astro_top_term_t(     1924,  2.1398279518720060e-04,  1.4058024849360011e-04), // f=0, s= 0, t=  48
                    new astro_top_term_t(      275,  1.7682640345886169e-04,  1.7533174872185701e-04), // f=0, s= 0, t=  49
                    new astro_top_term_t(      106,  2.0541555697206479e-04, -1.2965626091678131e-04), // f=0, s= 0, t=  50
                    new astro_top_term_t(     1335,  2.1493444414672229e-04, -1.0499061765897920e-04), // f=0, s= 0, t=  51
                    new astro_top_term_t(      271,  2.1059070752496800e-04, -9.6647526271973576e-05)  // f=0, s= 0, t=  52
                },
                new astro_top_term_t[]  // f=0, s= 1
                {
                    new astro_top_term_t(     1402, -2.4541007710854022e-02,  5.3822529675651168e-02), // f=0, s= 1, t=   0
                    new astro_top_term_t(        0,  3.7890000000000000e-02,  0.0000000000000000e+00), // f=0, s= 1, t=   1
                    new astro_top_term_t(     1331, -1.5981733929364941e-02,  1.9623571765255570e-02), // f=0, s= 1, t=   2
                    new astro_top_term_t(      522, -2.1667846285401051e-03,  1.3851523320136261e-02), // f=0, s= 1, t=   3
                    new astro_top_term_t(       71,  8.6723020043621556e-04,  4.7241967094509424e-03), // f=0, s= 1, t=   4
                    new astro_top_term_t(     1261, -3.8008804453790821e-03,  2.7392780267271222e-03), // f=0, s= 1, t=   5
                    new astro_top_term_t(     2875, -5.9015441259769620e-04,  3.3190792202711962e-03), // f=0, s= 1, t=   6
                    new astro_top_term_t(       35, -2.0646844712970550e-03, -2.6221746776036530e-03), // f=0, s= 1, t=   7
                    new astro_top_term_t(        4, -2.7341177070724500e-03,  8.4711655196688536e-04), // f=0, s= 1, t=   8
                    new astro_top_term_t(     1190, -2.1427361371102339e-03,  8.3829844178717706e-04), // f=0, s= 1, t=   9
                    new astro_top_term_t(      452, -6.5191457442425422e-04,  1.3962477909825410e-03), // f=0, s= 1, t=  10
                    new astro_top_term_t(        8, -1.4174353266390050e-03, -4.8800609526062797e-04), // f=0, s= 1, t=  11
                    new astro_top_term_t(      381, -7.7874017402358557e-04,  9.3105764052437356e-04), // f=0, s= 1, t=  12
                    new astro_top_term_t(      137, -1.1495119655719700e-03, -2.9568842990989599e-04), // f=0, s= 1, t=  13
                    new astro_top_term_t(     2804, -2.9320339800219448e-04,  1.0245691529492771e-03), // f=0, s= 1, t=  14
                    new astro_top_term_t(     1119, -9.8344008047702861e-04,  1.2122440012927610e-04), // f=0, s= 1, t=  15
                    new astro_top_term_t(     1543,  7.2662531437559094e-04,  2.2380632879329121e-04), // f=0, s= 1, t=  16
                    new astro_top_term_t(     1115, -4.9059020449034973e-04,  5.4125269506251448e-04), // f=0, s= 1, t=  17
                    new astro_top_term_t(      310, -5.0858314998190228e-04,  3.3071892812221339e-04)  // f=0, s= 1, t=  18
                },
                new astro_top_term_t[]  // f=0, s= 2
                {
                    new astro_top_term_t(     1402,  6.1564734183986881e-03,  6.9727544398706350e-03), // f=0, s= 2, t=   0
                    new astro_top_term_t(     1331,  3.4473718686991932e-03,  5.3547458162748890e-03), // f=0, s= 2, t=   1
                    new astro_top_term_t(        0, -6.0199059744408881e-03,  0.0000000000000000e+00), // f=0, s= 2, t=   2
                    new astro_top_term_t(      522,  1.8506017571807240e-03,  1.2255527089933200e-03), // f=0, s= 2, t=   3
                    new astro_top_term_t(        4, -9.0791196955079877e-04, -1.0597464135497440e-03), // f=0, s= 2, t=   4
                    new astro_top_term_t(     1261, -2.5360705302150178e-04,  1.1006253876392330e-03)  // f=0, s= 2, t=   5
                }
            },
            new astro_top_term_t[][]  // f=1
            {
                new astro_top_term_t[]  // f=1, s= 0
                {
                    new astro_top_term_t(        0,  4.1654711248260003e+00,  0.0000000000000000e+00), // f=1, s= 0, t=   0
                    new astro_top_term_t(     1402,  2.0610516934972331e-03, -4.5672163937433728e-03), // f=1, s= 0, t=   1
                    new astro_top_term_t(        4,  3.4901370659179681e-03,  2.2214931280208532e-03), // f=1, s= 0, t=   2
                    new astro_top_term_t(     1473, -2.5274121559055551e-04,  1.4438005411011950e-03), // f=1, s= 0, t=   3
                    new astro_top_term_t(        8,  1.1811406701290150e-03,  5.6294240646172070e-04), // f=1, s= 0, t=   4
                    new astro_top_term_t(      522,  1.6107717821433110e-04, -1.0623010008063969e-03), // f=1, s= 0, t=   5
                    new astro_top_term_t(     1331,  2.7651217171242138e-04, -3.4287860661004818e-04), // f=1, s= 0, t=   6
                    new astro_top_term_t(      593,  3.3064342491574279e-05,  3.2281574175741223e-04), // f=1, s= 0, t=   7
                    new astro_top_term_t(     2875,  1.7391859631880900e-05, -2.4371046197993790e-04), // f=1, s= 0, t=   8
                    new astro_top_term_t(       71,  4.3540142673103847e-05, -2.3324350638665871e-04), // f=1, s= 0, t=   9
                    new astro_top_term_t(       12, -1.4733740390542089e-05,  1.6458293257830351e-04), // f=1, s= 0, t=  10
                    new astro_top_term_t(       35, -7.7366158840780485e-05, -1.2373453663453069e-04), // f=1, s= 0, t=  11
                    new astro_top_term_t(      137,  1.0874981420582371e-04,  2.9642106239502481e-05), // f=1, s= 0, t=  12
                    new astro_top_term_t(      452,  3.2837879995135691e-05, -8.0083302669255936e-05), // f=1, s= 0, t=  13
                    new astro_top_term_t(      106, -7.1058821662524569e-05, -2.9131217624735831e-05), // f=1, s= 0, t=  14
                    new astro_top_term_t(     2945,  1.2419927037655510e-05,  6.9699905790296052e-05), // f=1, s= 0, t=  15
                    new astro_top_term_t(     1115,  4.9015184839992363e-05, -3.7626776449287022e-05), // f=1, s= 0, t=  16
                    new astro_top_term_t(     1261,  3.7725355301707661e-05, -2.7730601457362580e-05), // f=1, s= 0, t=  17
                    new astro_top_term_t(       63,  3.4620953450636643e-05, -2.5727054077819990e-05), // f=1, s= 0, t=  18
                    new astro_top_term_t(    17405, -2.4737834776937812e-05, -2.7668964560844859e-05), // f=1, s= 0, t=  19
                    new astro_top_term_t(      141, -1.9152627579774910e-05,  3.0608840803230902e-05), // f=1, s= 0, t=  20
                    new astro_top_term_t(     1543,  3.4403440209432860e-05,  8.2302490865683481e-06), // f=1, s= 0, t=  21
                    new astro_top_term_t(    28337, -2.9450824049753070e-05,  1.9509947544013829e-05), // f=1, s= 0, t=  22
                    new astro_top_term_t(      208, -3.3402595486149828e-05,  6.9721024372198099e-07), // f=1, s= 0, t=  23
                    new astro_top_term_t(       16, -2.5708226560219951e-05,  5.4690438392817274e-06), // f=1, s= 0, t=  24
                    new astro_top_term_t(     2804,  8.2434428446966487e-06, -2.2056820101697141e-05), // f=1, s= 0, t=  25
                    new astro_top_term_t(      133,  1.7535496693018641e-05, -1.5317433153649361e-05), // f=1, s= 0, t=  26
                    new astro_top_term_t(      177,  9.4209450622295435e-06,  1.8609791865429379e-05), // f=1, s= 0, t=  27
                    new astro_top_term_t(      204,  1.3418024861975440e-05, -1.5519809537438971e-05), // f=1, s= 0, t=  28
                    new astro_top_term_t(       67,  1.9542040889474939e-05, -2.0438794816549699e-07), // f=1, s= 0, t=  29
                    new astro_top_term_t(     1186, -1.1497747150912189e-05,  1.4072246679624251e-05), // f=1, s= 0, t=  30
                    new astro_top_term_t(     1421, -7.4843280005381670e-06,  1.2303197921178691e-05), // f=1, s= 0, t=  31
                    new astro_top_term_t(     1383,  4.1628765439336240e-06, -1.3719681259667719e-05), // f=1, s= 0, t=  32
                    new astro_top_term_t(      275,  7.1803653660818622e-06, -1.2018565522581090e-05), // f=1, s= 0, t=  33
                    new astro_top_term_t(      212,  1.3642921909731060e-05,  2.1147480912091179e-06), // f=1, s= 0, t=  34
                    new astro_top_term_t(       59, -2.8996742569855378e-06, -1.3093850574478501e-05), // f=1, s= 0, t=  35
                    new astro_top_term_t(     4348, -3.6912767653086940e-06, -1.2768209379424250e-05), // f=1, s= 0, t=  36
                    new astro_top_term_t(    17476,  8.8159150487736961e-06,  5.9364512059773267e-06), // f=1, s= 0, t=  37
                    new astro_top_term_t(       27,  6.2980935945243974e-06, -8.4472797839010716e-06), // f=1, s= 0, t=  38
                    new astro_top_term_t(     1406,  5.7734139509693619e-06,  8.7248712048903149e-06), // f=1, s= 0, t=  39
                    new astro_top_term_t(     1398,  1.0115989678361140e-05, -1.3220473394572090e-06), // f=1, s= 0, t=  40
                    new astro_top_term_t(    28407,  6.7642022362205018e-06, -7.5210649765782616e-06), // f=1, s= 0, t=  41
                    new astro_top_term_t(      664,  7.9736063591689162e-06, -2.2024063975673750e-06), // f=1, s= 0, t=  42
                    new astro_top_term_t(      381,  3.1761072779613261e-06, -7.3448935507901516e-06), // f=1, s= 0, t=  43
                    new astro_top_term_t(      541,  2.5898481389558310e-06, -7.0323742443685127e-06), // f=1, s= 0, t=  44
                    new astro_top_term_t(      503, -1.5979599372675240e-07,  7.4586091247535634e-06), // f=1, s= 0, t=  45
                    new astro_top_term_t(      271, -3.1728511373578081e-06, -6.6603677156048786e-06), // f=1, s= 0, t=  46
                    new astro_top_term_t(      129, -1.8397756711585990e-06, -7.0451043684758007e-06), // f=1, s= 0, t=  47
                    new astro_top_term_t(      247,  3.1533300480278400e-07, -7.1604266714610849e-06), // f=1, s= 0, t=  48
                    new astro_top_term_t(      200, -2.2707299296152400e-06, -6.6783936815674601e-06), // f=1, s= 0, t=  49
                    new astro_top_term_t(       31, -4.2146281840702766e-06, -5.6388188974303981e-06), // f=1, s= 0, t=  50
                    new astro_top_term_t(      341, -3.9640903174182786e-06, -5.5968952746387096e-06), // f=1, s= 0, t=  51
                    new astro_top_term_t(       20, -2.1086487296606532e-06, -4.2993743513743323e-06), // f=1, s= 0, t=  52
                    new astro_top_term_t(      974, -1.9337296464942458e-06,  4.1280092704873532e-06), // f=1, s= 0, t=  53
                    new astro_top_term_t(     1492,  1.2386884071711840e-06, -4.0279260453190260e-06), // f=1, s= 0, t=  54
                    new astro_top_term_t(     1454, -1.8674475925452479e-07,  4.1994903608977473e-06), // f=1, s= 0, t=  55
                    new astro_top_term_t(       55, -4.1054351042334211e-06, -5.3592659852261218e-07), // f=1, s= 0, t=  56
                    new astro_top_term_t(     1044,  3.5262751998119361e-06,  1.9367410649084749e-06), // f=1, s= 0, t=  57
                    new astro_top_term_t(       39, -7.0347175355875052e-07,  3.9228440423731852e-06), // f=1, s= 0, t=  58
                    new astro_top_term_t(      345, -2.7714162211172452e-06, -2.8012028299863999e-06), // f=1, s= 0, t=  59
                    new astro_top_term_t(      283, -3.8855399506750791e-06,  4.2353837211089759e-07), // f=1, s= 0, t=  60
                    new astro_top_term_t(     4418,  1.9380717555480920e-06,  3.3362906758656841e-06), // f=1, s= 0, t=  61
                    new astro_top_term_t(     1708,  3.8363066095252423e-06,  4.0071146882428590e-07), // f=1, s= 0, t=  62
                    new astro_top_term_t(      903, -2.4535489999476891e-06,  2.7713257447896069e-06), // f=1, s= 0, t=  63
                    new astro_top_term_t(      412, -2.4194398659193562e-06, -2.6818837171571071e-06), // f=1, s= 0, t=  64
                    new astro_top_term_t(    17334, -1.5878732753673419e-06, -3.0847459710695832e-06), // f=1, s= 0, t=  65
                    new astro_top_term_t(     1410,  2.1528516528967772e-06,  2.5955495586096718e-06), // f=1, s= 0, t=  66
                    new astro_top_term_t(      318, -1.6007757166478231e-06,  2.9567402302706461e-06), // f=1, s= 0, t=  67
                    new astro_top_term_t(    28266, -3.1282886569147322e-06,  1.0733981547714400e-06), // f=1, s= 0, t=  68
                    new astro_top_term_t(     1394,  3.2570224915190650e-06, -5.9986669416035846e-08), // f=1, s= 0, t=  69
                    new astro_top_term_t(    72490,  7.8958524750142168e-07,  3.1029527324934821e-06), // f=1, s= 0, t=  70
                    new astro_top_term_t(     9220,  2.8325155820290858e-06, -1.4397983817771050e-06), // f=1, s= 0, t=  71
                    new astro_top_term_t(     1614,  3.0748937595146630e-06,  5.0326493238056022e-07), // f=1, s= 0, t=  72
                    new astro_top_term_t(      408, -3.0722440820616450e-06,  4.9345659005564487e-07), // f=1, s= 0, t=  73
                    new astro_top_term_t(      337, -2.9632105019568382e-06,  9.3548023528643972e-08), // f=1, s= 0, t=  74
                    new astro_top_term_t(       79,  2.7845485358443480e-06,  7.4025133333557920e-07), // f=1, s= 0, t=  75
                    new astro_top_term_t(      354,  1.9360680354843408e-06,  2.0552415375649160e-06), // f=1, s= 0, t=  76
                    new astro_top_term_t(      612,  5.3703132984668158e-07,  2.6669817089723940e-06), // f=1, s= 0, t=  77
                    new astro_top_term_t(      279, -2.5413275396432281e-07,  2.6502544396045040e-06), // f=1, s= 0, t=  78
                    new astro_top_term_t(       75, -2.8830829576295098e-07,  2.6106707165907741e-06), // f=1, s= 0, t=  79
                    new astro_top_term_t(      479, -2.4319941706034390e-06,  7.6942886776602959e-07), // f=1, s= 0, t=  80
                    new astro_top_term_t(      416,  2.0096007406999570e-06,  1.5520759224768389e-06), // f=1, s= 0, t=  81
                    new astro_top_term_t(     1096, -2.3257110689311771e-06,  9.7831153778294669e-07), // f=1, s= 0, t=  82
                    new astro_top_term_t(      267, -2.5074175872274202e-06, -1.7701423128546709e-07), // f=1, s= 0, t=  83
                    new astro_top_term_t(     3162,  2.0698330022272070e-06, -1.3093627746428469e-06), // f=1, s= 0, t=  84
                    new astro_top_term_t(      832, -1.9285373422971391e-06,  1.4160461431375550e-06), // f=1, s= 0, t=  85
                    new astro_top_term_t(     1190,  2.2341390259691221e-06, -7.2482790590962755e-07), // f=1, s= 0, t=  86
                    new astro_top_term_t(      526,  1.7171371907181510e-06,  1.4529625626913769e-06), // f=1, s= 0, t=  87
                    new astro_top_term_t(     2733,  1.0269773239622380e-06, -2.0007533912028411e-06), // f=1, s= 0, t=  88
                    new astro_top_term_t(      574, -5.1232176836951890e-07, -2.1792410718247190e-06)  // f=1, s= 0, t=  89
                },
                new astro_top_term_t[]  // f=1, s= 1
                {
                    new astro_top_term_t(        0,  2.5335660204370001e+01,  0.0000000000000000e+00), // f=1, s= 1, t=   0
                    new astro_top_term_t(        4, -2.1897916824442529e-04,  1.7741955864290151e-03), // f=1, s= 1, t=   1
                    new astro_top_term_t(     1402, -1.3029527067277161e-03, -5.8987171410210541e-04), // f=1, s= 1, t=   2
                    new astro_top_term_t(        8, -5.7857466108113462e-05,  6.2766572590063064e-04), // f=1, s= 1, t=   3
                    new astro_top_term_t(      522, -3.0338218770459020e-04, -4.6787236333551149e-05), // f=1, s= 1, t=   4
                    new astro_top_term_t(     1331, -1.6234312756262061e-04, -1.3226650376835280e-04), // f=1, s= 1, t=   5
                    new astro_top_term_t(       35, -1.6846578574212679e-04,  5.1474337733700072e-05), // f=1, s= 1, t=   6
                    new astro_top_term_t(     1473,  1.3716672333829441e-04,  2.7605165443775319e-05), // f=1, s= 1, t=   7
                    new astro_top_term_t(       12, -1.1426985570984979e-04,  1.7062621376287921e-05), // f=1, s= 1, t=   8
                    new astro_top_term_t(     2875, -8.2654240464229338e-05, -1.4230432684507880e-05), // f=1, s= 1, t=   9
                    new astro_top_term_t(     2945,  3.6084370654255798e-05, -3.8286603564628976e-06), // f=1, s= 1, t=  10
                    new astro_top_term_t(      593,  3.1050555447062818e-05, -2.3120488185252789e-06), // f=1, s= 1, t=  11
                    new astro_top_term_t(       16, -8.7245446550032187e-06, -2.2703299980076952e-05), // f=1, s= 1, t=  12
                    new astro_top_term_t(      137,  5.5783261240786444e-06, -2.0642718216052691e-05), // f=1, s= 1, t=  13
                    new astro_top_term_t(     1115, -1.3513503169043030e-05, -1.1380785261757990e-05), // f=1, s= 1, t=  14
                    new astro_top_term_t(       71,  1.0574660674287450e-05,  1.3281184472910570e-05), // f=1, s= 1, t=  15
                    new astro_top_term_t(     1261, -8.2709801392281214e-06, -1.1688074277507450e-05), // f=1, s= 1, t=  16
                    new astro_top_term_t(     2804, -1.1505534938221521e-05, -5.1883886955282330e-06), // f=1, s= 1, t=  17
                    new astro_top_term_t(       27, -9.0174735475683120e-06, -5.2285448475511598e-06), // f=1, s= 1, t=  18
                    new astro_top_term_t(      452, -9.3099342968225354e-06, -3.8372164516008766e-06), // f=1, s= 1, t=  19
                    new astro_top_term_t(     1543,  2.4828575108302229e-06, -9.6169483680059550e-06), // f=1, s= 1, t=  20
                    new astro_top_term_t(       63, -6.5055787623425069e-06, -7.4385376578198243e-06), // f=1, s= 1, t=  21
                    new astro_top_term_t(      106, -4.5512156458035857e-06,  8.6658259433823040e-06), // f=1, s= 1, t=  22
                    new astro_top_term_t(      133, -6.3072609842624613e-06, -6.9148116074105151e-06), // f=1, s= 1, t=  23
                    new astro_top_term_t(    28337, -4.4530735415054223e-06, -6.7124377273408592e-06), // f=1, s= 1, t=  24
                    new astro_top_term_t(      141,  6.9466051844107804e-06,  4.0449257224296616e-06), // f=1, s= 1, t=  25
                    new astro_top_term_t(     1398,  1.0522658062293921e-06, -7.0563446922276261e-06), // f=1, s= 1, t=  26
                    new astro_top_term_t(       59, -5.7183266993529056e-06,  1.6273101086020201e-06), // f=1, s= 1, t=  27
                    new astro_top_term_t(     1383, -5.6807273131313079e-06, -9.5798703870467336e-07), // f=1, s= 1, t=  28
                    new astro_top_term_t(       20,  4.8130148431204634e-06, -2.5299767388566989e-06), // f=1, s= 1, t=  29
                    new astro_top_term_t(     4348, -5.3269871823774314e-06,  6.0216192993201565e-07), // f=1, s= 1, t=  30
                    new astro_top_term_t(      129, -4.2816786763917080e-06,  1.2524207693334740e-06)  // f=1, s= 1, t=  31
                },
                new astro_top_term_t[]  // f=1, s= 2
                {
                    new astro_top_term_t(        0, -1.8272218839163919e-02,  0.0000000000000000e+00), // f=1, s= 2, t=   0
                    new astro_top_term_t(        4, -4.2382205514535369e-04,  6.0951225139627217e-05), // f=1, s= 2, t=   1
                    new astro_top_term_t(        8, -2.4214950161520300e-04,  1.3555498866999700e-04), // f=1, s= 2, t=   2
                    new astro_top_term_t(     1402, -1.6731953542354990e-04,  1.4877590360992571e-04), // f=1, s= 2, t=   3
                    new astro_top_term_t(       12, -4.0514580902466651e-05, -4.6859429953625132e-05), // f=1, s= 2, t=   4
                    new astro_top_term_t(     1331, -4.4255479363996823e-05,  2.8523966290637259e-05), // f=1, s= 2, t=   5
                    new astro_top_term_t(      522, -2.6626163791203940e-05,  4.0414506964273448e-05), // f=1, s= 2, t=   6
                    new astro_top_term_t(       35, -2.3510812229573040e-06,  1.7218835328103482e-05), // f=1, s= 2, t=   7
                    new astro_top_term_t(     2875, -7.9375576768093987e-06,  1.3976740463566499e-05), // f=1, s= 2, t=   8
                    new astro_top_term_t(       16,  9.7305978044426803e-06, -1.0841117038814680e-05), // f=1, s= 2, t=   9
                    new astro_top_term_t(     2945, -4.5741854272088092e-07, -9.4120776625816716e-06)  // f=1, s= 2, t=  10
                },
                new astro_top_term_t[]  // f=1, s= 3
                {
                    new astro_top_term_t(        0,  1.9409931667071581e-03,  0.0000000000000000e+00), // f=1, s= 3, t=   0
                    new astro_top_term_t(        8, -5.2528177259165953e-05, -6.1055439661395280e-05), // f=1, s= 3, t=   1
                    new astro_top_term_t(        4, -6.0738316603187738e-05,  3.5387155556321078e-05), // f=1, s= 3, t=   2
                    new astro_top_term_t(     1402,  1.5967276563962451e-05,  3.5538891371840628e-05), // f=1, s= 3, t=   3
                    new astro_top_term_t(       12,  1.5332065032938759e-05, -2.1295985440213221e-05)  // f=1, s= 3, t=   4
                },
                new astro_top_term_t[]  // f=1, s= 4
                {
                    new astro_top_term_t(        0,  8.6099959150566781e-05,  0.0000000000000000e+00), // f=1, s= 4, t=   0
                    new astro_top_term_t(        4, -5.3544872037241911e-05, -3.2079251781364850e-05)  // f=1, s= 4, t=   1
                }
            },
            new astro_top_term_t[][]  // f=2
            {
                new astro_top_term_t[]  // f=2, s= 0
                {
                    new astro_top_term_t(        0, -1.7873895940349999e-01,  0.0000000000000000e+00), // f=2, s= 0, t=   0
                    new astro_top_term_t(     1473,  3.1629832749992988e-03, -2.0985870294942081e-03), // f=2, s= 0, t=   1
                    new astro_top_term_t(       71, -6.8180357663860398e-04,  1.1113519163940930e-03), // f=2, s= 0, t=   2
                    new astro_top_term_t(     1331,  1.6911781266743710e-04,  1.1885125258969610e-03), // f=2, s= 0, t=   3
                    new astro_top_term_t(      593,  5.4774166542017916e-04, -6.3844158559171749e-04), // f=2, s= 0, t=   4
                    new astro_top_term_t(     1402, -2.1438473609329679e-05, -5.2026929830140383e-04), // f=2, s= 0, t=   5
                    new astro_top_term_t(     1261, -5.7368849879358177e-05,  4.6888445010058260e-04), // f=2, s= 0, t=   6
                    new astro_top_term_t(      141, -9.6040157452198532e-05,  3.0415286563684890e-04), // f=2, s= 0, t=   7
                    new astro_top_term_t(      452,  1.2124121714193030e-04,  2.7470351829330812e-04), // f=2, s= 0, t=   8
                    new astro_top_term_t(     2945,  1.1035693239374779e-04, -1.4867545478165761e-04), // f=2, s= 0, t=   9
                    new astro_top_term_t(     1190, -5.9767419503565077e-05,  1.5052450817429230e-04), // f=2, s= 0, t=  10
                    new astro_top_term_t(     1543, -1.0081841648710751e-04,  1.1874253775582299e-04), // f=2, s= 0, t=  11
                    new astro_top_term_t(      522, -3.7998691236472437e-05, -1.1710923038633279e-04), // f=2, s= 0, t=  12
                    new astro_top_term_t(      381,  1.7939534333569889e-05,  1.2165766508553840e-04), // f=2, s= 0, t=  13
                    new astro_top_term_t(        8, -5.1065826531599063e-05, -1.0380178618772450e-04), // f=2, s= 0, t=  14
                    new astro_top_term_t(        4, -9.6101012457693182e-05, -3.7409967137145340e-05), // f=2, s= 0, t=  15
                    new astro_top_term_t(      212, -1.0011374095546651e-05,  9.1861092027193070e-05), // f=2, s= 0, t=  16
                    new astro_top_term_t(      208,  6.3118700287276614e-05,  6.6551813665925998e-05), // f=2, s= 0, t=  17
                    new astro_top_term_t(      106,  4.6362839297978442e-05,  7.1895265631125496e-05), // f=2, s= 0, t=  18
                    new astro_top_term_t(     2804,  2.5951137645879960e-05,  4.8810749391069221e-05), // f=2, s= 0, t=  19
                    new astro_top_term_t(     1119, -3.1798344550251363e-05,  4.3442443246298743e-05), // f=2, s= 0, t=  20
                    new astro_top_term_t(       35, -2.6365309213148771e-05, -3.7318748312868848e-05), // f=2, s= 0, t=  21
                    new astro_top_term_t(     1186,  4.5373467873211630e-05, -3.9315678706168707e-06), // f=2, s= 0, t=  22
                    new astro_top_term_t(      310, -5.7440916735318752e-06,  4.2608008327900133e-05), // f=2, s= 0, t=  23
                    new astro_top_term_t(       67, -3.2628206231552743e-05,  1.4188485967101389e-05), // f=2, s= 0, t=  24
                    new astro_top_term_t(      664, -1.3138986439777871e-05,  2.8745825031870391e-05), // f=2, s= 0, t=  25
                    new astro_top_term_t(    17476, -4.5717070878285386e-06, -2.7229104695341480e-05), // f=2, s= 0, t=  26
                    new astro_top_term_t(      283,  4.1649905603620409e-06,  2.5778764278582031e-05), // f=2, s= 0, t=  27
                    new astro_top_term_t(    28407, -2.6086891392275169e-05,  6.5019643042548361e-07), // f=2, s= 0, t=  28
                    new astro_top_term_t(     2875, -1.2604181155312930e-05, -2.2155626126791911e-05), // f=2, s= 0, t=  29
                    new astro_top_term_t(     2733,  4.6350251894966598e-06,  2.1642396781861640e-05), // f=2, s= 0, t=  30
                    new astro_top_term_t(       12,  1.3787761564064851e-05, -1.6882625493196851e-05), // f=2, s= 0, t=  31
                    new astro_top_term_t(       63,  7.3243482301750183e-06, -2.0430415756143061e-05), // f=2, s= 0, t=  32
                    new astro_top_term_t(      137,  1.7688914790024399e-05, -2.6466558968384412e-06), // f=2, s= 0, t=  33
                    new astro_top_term_t(     1048, -1.3743478634746570e-05,  1.1178051438197709e-05), // f=2, s= 0, t=  34
                    new astro_top_term_t(     1614, -1.4527685432017941e-05, -2.5021418914365550e-06), // f=2, s= 0, t=  35
                    new astro_top_term_t(     1044, -5.5188522930647246e-06,  1.3429997819426629e-05), // f=2, s= 0, t=  36
                    new astro_top_term_t(      239, -6.3224380303283023e-06,  1.2677314304305631e-05), // f=2, s= 0, t=  37
                    new astro_top_term_t(      133,  2.9376054111895790e-06, -1.3077911011252820e-05), // f=2, s= 0, t=  38
                    new astro_top_term_t(     1492, -9.7624713292315362e-06,  4.8639192065138997e-06), // f=2, s= 0, t=  39
                    new astro_top_term_t(     1454,  8.1915641060518482e-06, -7.1484181790528786e-06), // f=2, s= 0, t=  40
                    new astro_top_term_t(     4418,  2.9126133839920300e-06, -9.7148097661500440e-06), // f=2, s= 0, t=  41
                    new astro_top_term_t(      177, -8.1108828396444682e-06, -5.8479936691072763e-06), // f=2, s= 0, t=  42
                    new astro_top_term_t(     2663, -7.4109426380852051e-07,  8.1320485449777430e-06), // f=2, s= 0, t=  43
                    new astro_top_term_t(    17334,  7.7062782105104396e-06,  2.1832991677727231e-06), // f=2, s= 0, t=  44
                    new astro_top_term_t(     3016, -2.7097054305923721e-06,  7.2985152622098958e-06), // f=2, s= 0, t=  45
                    new astro_top_term_t(    28266,  3.0978305562474201e-06, -6.9523478329398256e-06), // f=2, s= 0, t=  46
                    new astro_top_term_t(       75, -7.5005963705916821e-06,  5.8552558907898478e-07), // f=2, s= 0, t=  47
                    new astro_top_term_t(      354,  2.2767150385710288e-06,  7.1310390742790356e-06), // f=2, s= 0, t=  48
                    new astro_top_term_t(      974, -2.8814074492595782e-06,  6.7739175425797296e-06), // f=2, s= 0, t=  49
                    new astro_top_term_t(       59, -4.8255490447201409e-06, -5.3633858025303309e-06), // f=2, s= 0, t=  50
                    new astro_top_term_t(     1115,  2.4165363833927252e-06, -5.9760260615249186e-06), // f=2, s= 0, t=  51
                    new astro_top_term_t(      204, -1.5420460797033710e-07, -6.3645330720578637e-06), // f=2, s= 0, t=  52
                    new astro_top_term_t(      612,  4.7791204064686441e-06, -4.1989180399775663e-06), // f=2, s= 0, t=  53
                    new astro_top_term_t(      574, -3.2533619740383079e-06,  4.8942036839749464e-06), // f=2, s= 0, t=  54
                    new astro_top_term_t(      978, -5.2816161970547299e-06,  2.4194318535884879e-06), // f=2, s= 0, t=  55
                    new astro_top_term_t(      129, -4.1704589160563732e-06, -3.8074741128128591e-06), // f=2, s= 0, t=  56
                    new astro_top_term_t(      200, -4.3860533690857434e-06, -2.8426817721421870e-06), // f=2, s= 0, t=  57
                    new astro_top_term_t(     1685, -4.4937857977759297e-06, -2.6664132359831861e-06), // f=2, s= 0, t=  58
                    new astro_top_term_t(     1327, -4.0394363815867892e-06,  3.0479150446235980e-06), // f=2, s= 0, t=  59
                    new astro_top_term_t(     1257, -4.4570496522011644e-06,  1.8417244601604851e-06)  // f=2, s= 0, t=  60
                },
                new astro_top_term_t[]  // f=2, s= 1
                {
                    new astro_top_term_t(        0, -6.1339663802007860e-04,  0.0000000000000000e+00), // f=2, s= 1, t=   0
                    new astro_top_term_t(     1331,  5.6664110172313892e-04, -8.1133128701713292e-05), // f=2, s= 1, t=   1
                    new astro_top_term_t(     1473, -1.9632805654171590e-04, -2.9843120844069362e-04), // f=2, s= 1, t=   2
                    new astro_top_term_t(       71, -1.9721831884265259e-04, -1.2980246221969450e-04), // f=2, s= 1, t=   3
                    new astro_top_term_t(     1402, -1.4909525632022481e-04,  5.5470899732068926e-06), // f=2, s= 1, t=   4
                    new astro_top_term_t(     1261,  1.4379398981882231e-04,  1.8301247697520619e-05), // f=2, s= 1, t=   5
                    new astro_top_term_t(     2945, -7.2143717725975376e-05, -6.1374574108146777e-05), // f=2, s= 1, t=   6
                    new astro_top_term_t(     1190,  7.4611359353024315e-05,  3.0289539777815409e-05), // f=2, s= 1, t=   7
                    new astro_top_term_t(      593, -5.9939495914062063e-05, -5.1836225465555612e-05), // f=2, s= 1, t=   8
                    new astro_top_term_t(        8,  3.9984073707170713e-05, -3.0799116570219068e-05), // f=2, s= 1, t=   9
                    new astro_top_term_t(     1543,  3.1094179263405468e-05,  2.6861784078602320e-05), // f=2, s= 1, t=  10
                    new astro_top_term_t(      381,  3.7436713283265682e-05, -5.1576484771147526e-06), // f=2, s= 1, t=  11
                    new astro_top_term_t(     1119,  2.9759191459225299e-05,  2.2175266912892650e-05), // f=2, s= 1, t=  12
                    new astro_top_term_t(      452,  3.2286698454013622e-05, -1.4330028054859710e-05), // f=2, s= 1, t=  13
                    new astro_top_term_t(      522, -3.3610639556594409e-05,  1.0732560764842151e-05), // f=2, s= 1, t=  14
                    new astro_top_term_t(     2804,  2.6892727962532569e-05, -1.2202088159044841e-05), // f=2, s= 1, t=  15
                    new astro_top_term_t(        4,  1.6145362682689730e-05, -2.1535180644982459e-05), // f=2, s= 1, t=  16
                    new astro_top_term_t(       35,  1.1619842301292181e-05,  2.1259793782462840e-05), // f=2, s= 1, t=  17
                    new astro_top_term_t(      310,  2.1159419772930708e-05,  3.2365387780938209e-06), // f=2, s= 1, t=  18
                    new astro_top_term_t(      212, -1.7383881837993001e-05, -1.2570139978872199e-07), // f=2, s= 1, t=  19
                    new astro_top_term_t(     2733,  1.5857546421435610e-05, -2.6232455173985150e-06), // f=2, s= 1, t=  20
                    new astro_top_term_t(     1048,  9.8372722970562282e-06,  1.2182701032807050e-05), // f=2, s= 1, t=  21
                    new astro_top_term_t(       12,  1.1952587407675961e-05,  7.5377841861106607e-06), // f=2, s= 1, t=  22
                    new astro_top_term_t(      283, -1.0224581925761621e-05,  1.2303488561812810e-06), // f=2, s= 1, t=  23
                    new astro_top_term_t(      239,  8.6293198115122691e-06,  4.6174857696641109e-06)  // f=2, s= 1, t=  24
                },
                new astro_top_term_t[]  // f=2, s= 2
                {
                    new astro_top_term_t(     1331,  2.3820657956875651e-05, -1.4117217046583029e-04), // f=2, s= 2, t=   0
                    new astro_top_term_t(        0,  5.7169784317623151e-05,  0.0000000000000000e+00), // f=2, s= 2, t=   1
                    new astro_top_term_t(     1261,  2.8584658843704731e-05, -1.9153125455819128e-05), // f=2, s= 2, t=   2
                    new astro_top_term_t(       71, -6.3863669075967202e-06, -3.0182516098252429e-05), // f=2, s= 2, t=   3
                    new astro_top_term_t(     2945, -1.6831473476249729e-05,  1.7778057659448761e-05), // f=2, s= 2, t=   4
                    new astro_top_term_t(     1402, -8.7282659443619542e-06,  2.1969730180690399e-05), // f=2, s= 2, t=   5
                    new astro_top_term_t(     1190,  1.8655664239774939e-05, -1.4275349333837280e-05), // f=2, s= 2, t=   6
                    new astro_top_term_t(        8,  2.0572423386833478e-05,  2.9970650695923720e-06)  // f=2, s= 2, t=   7
                }
            },
            new astro_top_term_t[][]  // f=3
            {
                new astro_top_term_t[]  // f=3, s= 0
                {
                    new astro_top_term_t(        0, -1.7340471864230000e-01,  0.0000000000000000e+00), // f=3, s= 0, t=   0
                    new astro_top_term_t(     1473,  2.1234813115076170e-03,  3.0130058621335031e-03), // f=3, s= 0, t=   1
                    new astro_top_term_t(       71, -1.0505499200359431e-03, -7.0849649152681434e-04), // f=3, s= 0, t=   2
                    new astro_top_term_t(     1331,  1.1926139685020671e-03, -1.2467514555563171e-04), // f=3, s= 0, t=   3
                    new astro_top_term_t(      593,  6.3328818822552336e-04,  5.1838245948530805e-04), // f=3, s= 0, t=   4
                    new astro_top_term_t(     1402, -4.0198902705438299e-04,  3.4879411820834029e-04), // f=3, s= 0, t=   5
                    new astro_top_term_t(     1261,  4.6754614871725769e-04,  6.5799333440031267e-05), // f=3, s= 0, t=   6
                    new astro_top_term_t(      141, -3.1574545029112658e-04, -1.1057917696717649e-04), // f=3, s= 0, t=   7
                    new astro_top_term_t(      452,  2.7842572840314938e-04, -1.1076531434037220e-04), // f=3, s= 0, t=   8
                    new astro_top_term_t(     2945,  1.4727652777389929e-04,  1.0316600256389621e-04), // f=3, s= 0, t=   9
                    new astro_top_term_t(     1190,  1.4974726446281261e-04,  6.1544751468406387e-05), // f=3, s= 0, t=  10
                    new astro_top_term_t(     1543, -1.1418008898386521e-04, -9.8159556344747838e-05), // f=3, s= 0, t=  11
                    new astro_top_term_t(      522, -6.9432011255429296e-05,  1.0500648997119170e-04), // f=3, s= 0, t=  12
                    new astro_top_term_t(      381,  1.2168346615540400e-04, -1.5724958971846889e-05), // f=3, s= 0, t=  13
                    new astro_top_term_t(        4,  3.5694868279606838e-05, -1.1563400718698000e-04), // f=3, s= 0, t=  14
                    new astro_top_term_t(        8,  1.0601495232676501e-04, -5.3792770194159239e-05), // f=3, s= 0, t=  15
                    new astro_top_term_t(      212, -8.8206829660081605e-05, -5.6186044809908770e-06), // f=3, s= 0, t=  16
                    new astro_top_term_t(      208, -6.2116275994912254e-05,  6.2681468828265711e-05), // f=3, s= 0, t=  17
                    new astro_top_term_t(      106, -5.9881700950469839e-05,  1.8173312312951261e-05), // f=3, s= 0, t=  18
                    new astro_top_term_t(     2804,  4.9734881336053070e-05, -2.4663551556821301e-05), // f=3, s= 0, t=  19
                    new astro_top_term_t(     1119,  4.3142176920002318e-05,  3.2166118196149167e-05), // f=3, s= 0, t=  20
                    new astro_top_term_t(     1186,  5.3160131153765831e-06,  4.6681999520445390e-05), // f=3, s= 0, t=  21
                    new astro_top_term_t(      310,  4.2381917128956993e-05,  6.2114355060802278e-06), // f=3, s= 0, t=  22
                    new astro_top_term_t(       67,  9.7801450369382322e-06,  3.8008739374555199e-05), // f=3, s= 0, t=  23
                    new astro_top_term_t(       35,  3.4777523130477253e-05,  8.0807437621346330e-06), // f=3, s= 0, t=  24
                    new astro_top_term_t(      664, -2.7772391338347379e-05, -1.2933056019872690e-05), // f=3, s= 0, t=  25
                    new astro_top_term_t(    17476,  2.6116611382398249e-05, -5.2512707918252350e-06), // f=3, s= 0, t=  26
                    new astro_top_term_t(      283, -2.5901894192097461e-05,  7.6613799535167943e-07), // f=3, s= 0, t=  27
                    new astro_top_term_t(    28407, -1.3727149833568929e-06, -2.5571899427495732e-05), // f=3, s= 0, t=  28
                    new astro_top_term_t(       12,  1.7840323650284740e-05,  1.4017491541647661e-05), // f=3, s= 0, t=  29
                    new astro_top_term_t(     2875, -1.2752702589739361e-05,  1.8674367648851481e-05), // f=3, s= 0, t=  30
                    new astro_top_term_t(       63, -2.1141833621778831e-05, -6.5305966080253417e-06), // f=3, s= 0, t=  31
                    new astro_top_term_t(     2733,  2.1699001266656308e-05, -4.2325482084098446e-06), // f=3, s= 0, t=  32
                    new astro_top_term_t(     1048,  1.1071845430997370e-05,  1.3823004095404599e-05), // f=3, s= 0, t=  33
                    new astro_top_term_t(      137, -1.3154331920868999e-05, -9.3564722901572934e-06), // f=3, s= 0, t=  34
                    new astro_top_term_t(     1044,  1.3602483536936790e-05,  6.6866532885834768e-06), // f=3, s= 0, t=  35
                    new astro_top_term_t(     1614,  3.0521160394305390e-06, -1.4295846223173739e-05), // f=3, s= 0, t=  36
                    new astro_top_term_t(      239,  1.2491201410693211e-05,  6.3469005652658392e-06), // f=3, s= 0, t=  37
                    new astro_top_term_t(      133, -1.3409179117794041e-05, -2.6113349355778870e-06), // f=3, s= 0, t=  38
                    new astro_top_term_t(     1492, -4.9926812258289137e-06, -9.3442347365693383e-06), // f=3, s= 0, t=  39
                    new astro_top_term_t(     1454,  7.1182846121023137e-06,  7.7835699509488784e-06), // f=3, s= 0, t=  40
                    new astro_top_term_t(     4418,  9.5085270250543767e-06,  2.5720788995216492e-06), // f=3, s= 0, t=  41
                    new astro_top_term_t(      354, -7.7105865781446408e-06,  3.5195735852526461e-06), // f=3, s= 0, t=  42
                    new astro_top_term_t(     2663,  8.0823599962064017e-06,  8.8289336639692517e-07), // f=3, s= 0, t=  43
                    new astro_top_term_t(       75, -6.1089854320482681e-07, -8.0284473434343767e-06), // f=3, s= 0, t=  44
                    new astro_top_term_t(    17334,  2.4666413292489791e-06, -7.6105265820617346e-06), // f=3, s= 0, t=  45
                    new astro_top_term_t(      974,  6.9130746927758178e-06,  3.5048585089334658e-06), // f=3, s= 0, t=  46
                    new astro_top_term_t(    28266, -6.8235793370187026e-06, -3.3487790585021339e-06), // f=3, s= 0, t=  47
                    new astro_top_term_t(     3016, -7.0427593206860361e-06, -2.6580682576526168e-06), // f=3, s= 0, t=  48
                    new astro_top_term_t(       59, -5.5456494164143571e-06,  4.9667283828556126e-06), // f=3, s= 0, t=  49
                    new astro_top_term_t(      204, -7.0017722378178446e-06, -1.6817751033546979e-06), // f=3, s= 0, t=  50
                    new astro_top_term_t(     1115, -6.4024937184914472e-06,  2.1192315231268140e-06), // f=3, s= 0, t=  51
                    new astro_top_term_t(      978,  2.3623585675585379e-06,  5.3260537432133881e-06), // f=3, s= 0, t=  52
                    new astro_top_term_t(      129, -3.9867297079366098e-06,  4.2341276287979132e-06), // f=3, s= 0, t=  53
                    new astro_top_term_t(      574, -4.8314043929159330e-06, -2.9544972252560498e-06), // f=3, s= 0, t=  54
                    new astro_top_term_t(      612,  4.0663883671886490e-06,  3.8452439907704313e-06), // f=3, s= 0, t=  55
                    new astro_top_term_t(      200, -3.2782990554685469e-06,  4.3746397303989681e-06), // f=3, s= 0, t=  56
                    new astro_top_term_t(      247,  4.4086088444075100e-06, -3.0262449322858620e-06), // f=3, s= 0, t=  57
                    new astro_top_term_t(     1685,  2.7665561507113492e-06, -4.4616466628081997e-06), // f=3, s= 0, t=  58
                    new astro_top_term_t(     1335, -2.3134382826680641e-06,  4.5013914052228881e-06), // f=3, s= 0, t=  59
                    new astro_top_term_t(     1327,  3.1175433103076271e-06,  3.7924860549293981e-06), // f=3, s= 0, t=  60
                    new astro_top_term_t(      903,  4.4741994730303496e-06,  1.7822080022821150e-06), // f=3, s= 0, t=  61
                    new astro_top_term_t(       16, -1.4927342901310379e-06,  4.4884086825602869e-06), // f=3, s= 0, t=  62
                    new astro_top_term_t(      169,  2.5354018835911621e-06,  3.5730000800789471e-06), // f=3, s= 0, t=  63
                    new astro_top_term_t(      271, -2.1688629242033941e-06,  3.7820832945120280e-06), // f=3, s= 0, t=  64
                    new astro_top_term_t(       79,  3.8979055894631164e-06, -1.8923407502835230e-06), // f=3, s= 0, t=  65
                    new astro_top_term_t(      416,  3.9376438761867171e-06, -5.2538703904847715e-07), // f=3, s= 0, t=  66
                    new astro_top_term_t(    17405,  1.0904344065136349e-06,  3.5409195078133210e-06), // f=3, s= 0, t=  67
                    new astro_top_term_t(    28337,  3.4736718287656798e-06, -5.2407682497120616e-07), // f=3, s= 0, t=  68
                    new astro_top_term_t(     1351, -3.4436468577006500e-06, -7.4140986208801818e-08), // f=3, s= 0, t=  69
                    new astro_top_term_t(     1312,  3.3433376977429141e-06, -7.9547674738811048e-07), // f=3, s= 0, t=  70
                    new astro_top_term_t(      832,  3.0647904819735712e-06,  1.2969722855848811e-06), // f=3, s= 0, t=  71
                    new astro_top_term_t(      275,  1.5394331893610239e-06, -2.7751033895625021e-06), // f=3, s= 0, t=  72
                    new astro_top_term_t(     2592,  2.8273974520595848e-06,  1.3333677990888271e-06), // f=3, s= 0, t=  73
                    new astro_top_term_t(    17264,  1.6047868553221060e-06, -2.5731990436490451e-06)  // f=3, s= 0, t=  74
                },
                new astro_top_term_t[]  // f=3, s= 1
                {
                    new astro_top_term_t(     1331, -6.0044174932694893e-05, -5.6841814115066798e-04), // f=3, s= 1, t=   0
                    new astro_top_term_t(     1473,  2.8299001111328808e-04, -1.9899683117934169e-04), // f=3, s= 1, t=   1
                    new astro_top_term_t(       71,  1.3109847892346901e-04, -1.9645547029607899e-04), // f=3, s= 1, t=   2
                    new astro_top_term_t(     1402,  9.8753503746441094e-05,  1.1624141830784459e-04), // f=3, s= 1, t=   3
                    new astro_top_term_t(     1261,  2.0857393493693379e-05, -1.4333407898431410e-04), // f=3, s= 1, t=   4
                    new astro_top_term_t(     2945,  5.7600265965337831e-05, -7.1668848804472786e-05), // f=3, s= 1, t=   5
                    new astro_top_term_t(     1190,  3.1158128562030883e-05, -7.4221508541133025e-05), // f=3, s= 1, t=   6
                    new astro_top_term_t(      593,  4.9568824746840162e-05, -5.9443745710499677e-05), // f=3, s= 1, t=   7
                    new astro_top_term_t(        0, -6.3916210233836622e-05,  0.0000000000000000e+00), // f=3, s= 1, t=   8
                    new astro_top_term_t(        8,  3.2326459763508557e-05,  4.1097067153751401e-05), // f=3, s= 1, t=   9
                    new astro_top_term_t(        4,  3.8505374914800768e-05, -2.4690750688454580e-05), // f=3, s= 1, t=  10
                    new astro_top_term_t(     1543, -2.6127774810849349e-05,  2.9873942333388641e-05), // f=3, s= 1, t=  11
                    new astro_top_term_t(      381, -4.4772545431527864e-06, -3.7396119251182150e-05), // f=3, s= 1, t=  12
                    new astro_top_term_t(     1119,  2.2430310753365549e-05, -2.9557761386193919e-05), // f=3, s= 1, t=  13
                    new astro_top_term_t(      522,  2.9886224978272480e-05,  2.0259791365064171e-05), // f=3, s= 1, t=  14
                    new astro_top_term_t(      452, -1.3106498259692980e-05, -3.2637289766253343e-05), // f=3, s= 1, t=  15
                    new astro_top_term_t(     2804, -1.1508149132666279e-05, -2.7316790371599908e-05), // f=3, s= 1, t=  16
                    new astro_top_term_t(      310,  3.4738350608151140e-06, -2.1017437982180379e-05), // f=3, s= 1, t=  17
                    new astro_top_term_t(      212,  2.1402716283134239e-06, -1.9178262093378069e-05), // f=3, s= 1, t=  18
                    new astro_top_term_t(     2733, -2.3283762042050202e-06, -1.5882216399278110e-05), // f=3, s= 1, t=  19
                    new astro_top_term_t(     1048,  1.2246829111767810e-05, -9.7512980427413643e-06), // f=3, s= 1, t=  20
                    new astro_top_term_t(       12, -7.6135340345472092e-06,  1.2831304341856990e-05), // f=3, s= 1, t=  21
                    new astro_top_term_t(       35, -5.7732549173295510e-06,  1.0016273100035810e-05), // f=3, s= 1, t=  22
                    new astro_top_term_t(      283, -2.2363814820487809e-06, -9.6218108123675468e-06), // f=3, s= 1, t=  23
                    new astro_top_term_t(      239,  4.6765179736328208e-06, -8.5111533123445213e-06), // f=3, s= 1, t=  24
                    new astro_top_term_t(      106, -1.9836445213635969e-06,  7.8016311549488939e-06), // f=3, s= 1, t=  25
                    new astro_top_term_t(     2875,  5.8334802150818700e-06,  5.1997278670085084e-06), // f=3, s= 1, t=  26
                    new astro_top_term_t(     1044,  2.4978011849958161e-06, -6.6276491176842909e-06)  // f=3, s= 1, t=  27
                },
                new astro_top_term_t[]  // f=3, s= 2
                {
                    new astro_top_term_t(     1331, -1.3993900805601181e-04, -2.9006879879181548e-05), // f=3, s= 2, t=   0
                    new astro_top_term_t(        0,  6.7829783611580237e-05,  0.0000000000000000e+00), // f=3, s= 2, t=   1
                    new astro_top_term_t(     1261, -1.8612445466711130e-05, -2.8908877766737650e-05), // f=3, s= 2, t=   2
                    new astro_top_term_t(       71,  2.9853088727300880e-05, -7.1012915175069879e-06), // f=3, s= 2, t=   3
                    new astro_top_term_t(        4,  2.6321188906063949e-05,  1.3464234343799510e-05), // f=3, s= 2, t=   4
                    new astro_top_term_t(     1402,  2.3284325058403039e-05, -6.8144038337101704e-06), // f=3, s= 2, t=   5
                    new astro_top_term_t(     2945, -1.7707706192590389e-05, -1.5833741011362669e-05), // f=3, s= 2, t=   6
                    new astro_top_term_t(     1190, -1.4044910483808361e-05, -1.8816964596903641e-05), // f=3, s= 2, t=   7
                    new astro_top_term_t(        8, -3.2757789019448882e-06,  2.1268129579065079e-05), // f=3, s= 2, t=   8
                    new astro_top_term_t(     1473, -9.8933125150553904e-06, -1.3083286541079480e-05), // f=3, s= 2, t=   9
                    new astro_top_term_t(     1119, -7.2140288335202669e-06, -1.1729524860644351e-05)  // f=3, s= 2, t=  10
                },
                new astro_top_term_t[]  // f=3, s= 3
                {
                    new astro_top_term_t(     1331, -1.7737816039767501e-05,  2.8055420335966871e-05)  // f=3, s= 3, t=   0
                }
            },
            new astro_top_term_t[][]  // f=4
            {
                new astro_top_term_t[]  // f=4, s= 0
                {
                    new astro_top_term_t(        0, -5.1702307822780000e-02,  0.0000000000000000e+00), // f=4, s= 0, t=   0
                    new astro_top_term_t(     1402, -1.3296787157385281e-04,  1.3523784099823399e-04), // f=4, s= 0, t=   1
                    new astro_top_term_t(     1543,  1.6300899640352639e-04,  5.5927755421839437e-05), // f=4, s= 0, t=   2
                    new astro_top_term_t(     1473, -2.3179975765177140e-05, -9.8184622098344091e-05), // f=4, s= 0, t=   3
                    new astro_top_term_t(      522, -1.9097390976957099e-05,  3.6607666661423369e-05), // f=4, s= 0, t=   4
                    new astro_top_term_t(      664,  3.2543564829276270e-05,  1.0672339397019160e-06), // f=4, s= 0, t=   5
                    new astro_top_term_t(     1331, -2.0841818514170151e-05,  1.2410223144740100e-05), // f=4, s= 0, t=   6
                    new astro_top_term_t(        4, -8.0356073429830634e-06,  1.9857752253856971e-05), // f=4, s= 0, t=   7
                    new astro_top_term_t(      593, -1.0459691293338759e-05, -1.7857031017673810e-05), // f=4, s= 0, t=   8
                    new astro_top_term_t(     1614,  2.0002336741252269e-05,  1.4700260033157791e-06), // f=4, s= 0, t=   9
                    new astro_top_term_t(     2875, -3.5880881863000781e-06,  8.0256130260543578e-06), // f=4, s= 0, t=  10
                    new astro_top_term_t(     3016,  8.5254320868063325e-06, -1.8479247218425171e-07), // f=4, s= 0, t=  11
                    new astro_top_term_t(        8, -5.4541260810660720e-06,  3.9805913226076404e-06), // f=4, s= 0, t=  12
                    new astro_top_term_t(      452, -3.7355080445397231e-06,  4.0025013123176398e-06), // f=4, s= 0, t=  13
                    new astro_top_term_t(       35, -2.1277230177098351e-06, -4.9607027782605714e-06), // f=4, s= 0, t=  14
                    new astro_top_term_t(      137, -3.8233051639183096e-06, -3.1849844866680141e-06)  // f=4, s= 0, t=  15
                },
                new astro_top_term_t[]  // f=4, s= 1
                {
                    new astro_top_term_t(        0,  1.9166844047183211e-04,  0.0000000000000000e+00), // f=4, s= 1, t=   0
                    new astro_top_term_t(     1402,  3.9019701387765488e-05,  3.8671378812923788e-05), // f=4, s= 1, t=   1
                    new astro_top_term_t(     1543,  1.5119207079965939e-05, -4.3332491992581899e-05), // f=4, s= 1, t=   2
                    new astro_top_term_t(     1331,  5.8922122899993431e-06,  1.0029385340576861e-05), // f=4, s= 1, t=   3
                    new astro_top_term_t(      522,  1.0280806202111021e-05,  5.2916455231996232e-06)  // f=4, s= 1, t=   4
                }
            },
            new astro_top_term_t[][]  // f=5
            {
                new astro_top_term_t[]  // f=5, s= 0
                {
                    new astro_top_term_t(        0,  1.3977992515640000e-01,  0.0000000000000000e+00), // f=5, s= 0, t=   0
                    new astro_top_term_t(     1402,  1.2883200572372361e-04,  1.3471156433694759e-04), // f=5, s= 0, t=   1
                    new astro_top_term_t(     1543, -4.9979310098879121e-05,  1.6180896042292819e-04), // f=5, s= 0, t=   2
                    new astro_top_term_t(     1473, -2.2060691195269469e-05, -9.3668360843677731e-05), // f=5, s= 0, t=   3
                    new astro_top_term_t(      522,  3.5537595493916502e-05,  1.9797390613745271e-05), // f=5, s= 0, t=   4
                    new astro_top_term_t(      664, -7.4709042862279015e-08,  3.1983280221550278e-05), // f=5, s= 0, t=   5
                    new astro_top_term_t(     1331,  1.1695376261855580e-05,  2.0798562232570440e-05), // f=5, s= 0, t=   6
                    new astro_top_term_t(      593, -9.9523808317349038e-06, -1.7074224315439731e-05), // f=5, s= 0, t=   7
                    new astro_top_term_t(     1614, -9.4118113483089596e-07,  1.9684674172624949e-05), // f=5, s= 0, t=   8
                    new astro_top_term_t(        4, -1.5932093431444279e-05,  9.4977011164352140e-06), // f=5, s= 0, t=   9
                    new astro_top_term_t(     2875,  7.8417800319575014e-06,  3.7228551547920001e-06), // f=5, s= 0, t=  10
                    new astro_top_term_t(     3016,  4.3444197029476981e-07,  8.3694483911743903e-06), // f=5, s= 0, t=  11
                    new astro_top_term_t(        8, -5.9902893836260643e-06, -4.3498982251067908e-06), // f=5, s= 0, t=  12
                    new astro_top_term_t(      137, -3.0046395968887028e-06,  4.7462972937701063e-06), // f=5, s= 0, t=  13
                    new astro_top_term_t(      452,  3.9819715036126689e-06,  3.7540995120279268e-06), // f=5, s= 0, t=  14
                    new astro_top_term_t(       35, -4.0365965056395762e-08,  5.3159364396026402e-06), // f=5, s= 0, t=  15
                    new astro_top_term_t(      141,  1.2649492629014499e-07,  4.8536627302669076e-06), // f=5, s= 0, t=  16
                    new astro_top_term_t(     1261,  1.1901325673451860e-06,  4.3994722705730970e-06), // f=5, s= 0, t=  17
                    new astro_top_term_t(     2945, -2.7264749792923272e-06, -3.6454310686215340e-06), // f=5, s= 0, t=  18
                    new astro_top_term_t(     1685,  6.9410338646847810e-07,  3.3566050467992739e-06), // f=5, s= 0, t=  19
                    new astro_top_term_t(      734,  8.2952840171332113e-07,  3.2356209378486682e-06), // f=5, s= 0, t=  20
                    new astro_top_term_t(      208,  2.4256846662067760e-06, -9.5181484992359002e-07), // f=5, s= 0, t=  21
                    new astro_top_term_t(      279, -2.5271926303381422e-06, -3.7922031983877281e-07), // f=5, s= 0, t=  22
                    new astro_top_term_t(     1115,  5.6515925473070844e-07,  2.1078436053637770e-06), // f=5, s= 0, t=  23
                    new astro_top_term_t(     1257, -1.3558265037946960e-06,  1.4215822110290420e-06), // f=5, s= 0, t=  24
                    new astro_top_term_t(       12,  1.3204389683748120e-07, -1.8835278008457151e-06), // f=5, s= 0, t=  25
                    new astro_top_term_t(       71, -1.9042145379339269e-07, -1.8630215116559800e-06), // f=5, s= 0, t=  26
                    new astro_top_term_t(      212,  1.6053042001407380e-06,  8.4671630560881122e-07), // f=5, s= 0, t=  27
                    new astro_top_term_t(    17405,  1.3113174080732559e-06, -4.9311123653462177e-07), // f=5, s= 0, t=  28
                    new astro_top_term_t(    17547,  1.0602949344425471e-06,  8.9833952832208403e-07), // f=5, s= 0, t=  29
                    new astro_top_term_t(     2804,  9.2561678972169579e-07,  9.7236819416524451e-07), // f=5, s= 0, t=  30
                    new astro_top_term_t(    28337, -9.4657252626238933e-08, -1.1160255453211769e-06), // f=5, s= 0, t=  31
                    new astro_top_term_t(    28478,  9.2430781502463316e-07, -6.2265716348343233e-07), // f=5, s= 0, t=  32
                    new astro_top_term_t(       75, -9.2130387409465058e-07,  5.9403283989994324e-07), // f=5, s= 0, t=  33
                    new astro_top_term_t(      106, -2.0304770919600219e-07, -1.0388792453459521e-06), // f=5, s= 0, t=  34
                    new astro_top_term_t(      381,  5.9915368756519242e-07,  8.3180916441107633e-07), // f=5, s= 0, t=  35
                    new astro_top_term_t(     1190,  2.7201162062442490e-08,  1.0190741100380801e-06), // f=5, s= 0, t=  36
                    new astro_top_term_t(     3087,  2.9873098731161120e-07,  9.6483981362737291e-07), // f=5, s= 0, t=  37
                    new astro_top_term_t(       63,  5.8956418822020129e-07, -7.7745585971202028e-07), // f=5, s= 0, t=  38
                    new astro_top_term_t(     1186,  2.4117256641652101e-07, -9.2801537744071498e-07), // f=5, s= 0, t=  39
                    new astro_top_term_t(       67, -3.2941946279216570e-07,  7.9963658829793748e-07), // f=5, s= 0, t=  40
                    new astro_top_term_t(    17476, -7.4361213537121858e-07, -1.2891468267067801e-07), // f=5, s= 0, t=  41
                    new astro_top_term_t(      283,  2.8076344536095831e-08,  7.3853891759583075e-07), // f=5, s= 0, t=  42
                    new astro_top_term_t(     1756,  3.1161152193772379e-07,  6.3104215205486464e-07), // f=5, s= 0, t=  43
                    new astro_top_term_t(      247, -3.0658465478480348e-07, -6.2524043128321752e-07), // f=5, s= 0, t=  44
                    new astro_top_term_t(    28407, -2.6107491970896922e-07,  5.4499817140536954e-07), // f=5, s= 0, t=  45
                    new astro_top_term_t(      354,  4.9460726403829815e-07,  2.9447867230664199e-07), // f=5, s= 0, t=  46
                    new astro_top_term_t(     1421, -3.1957781759109811e-07, -4.3033184225872391e-07), // f=5, s= 0, t=  47
                    new astro_top_term_t(     1383,  4.0796989384481828e-07,  3.2647659230997358e-07), // f=5, s= 0, t=  48
                    new astro_top_term_t(      133,  1.1815490710962760e-07, -5.0524398677283852e-07), // f=5, s= 0, t=  49
                    new astro_top_term_t(      177, -4.6755777925983522e-07,  2.0641392675487111e-07), // f=5, s= 0, t=  50
                    new astro_top_term_t(      805,  2.4510909183821792e-07,  4.3308098936875489e-07), // f=5, s= 0, t=  51
                    new astro_top_term_t(     1563,  2.0041064699998999e-07, -4.4445832345851582e-07), // f=5, s= 0, t=  52
                    new astro_top_term_t(       59, -1.2658294965005670e-07, -4.6671870382330612e-07), // f=5, s= 0, t=  53
                    new astro_top_term_t(     1524, -7.8614887212231827e-08,  4.6760439016103602e-07), // f=5, s= 0, t=  54
                    new astro_top_term_t(     4348,  4.6767825285878958e-07,  3.5290700804681861e-08), // f=5, s= 0, t=  55
                    new astro_top_term_t(       16,  4.2127049178471667e-07, -1.8684252468076170e-07), // f=5, s= 0, t=  56
                    new astro_top_term_t(     4489,  1.8055735116211371e-07,  4.2107641051871891e-07), // f=5, s= 0, t=  57
                    new astro_top_term_t(       79, -4.0548063730912672e-07, -2.1065225409216059e-07), // f=5, s= 0, t=  58
                    new astro_top_term_t(     1398, -1.2435454796683390e-07,  3.9919922981210330e-07), // f=5, s= 0, t=  59
                    new astro_top_term_t(      145, -3.0487589061275071e-07,  2.8316787420643652e-07), // f=5, s= 0, t=  60
                    new astro_top_term_t(      345,  3.4373900983948560e-07, -2.2193498580565619e-07), // f=5, s= 0, t=  61
                    new astro_top_term_t(      275, -7.6720036516022037e-09, -4.0650745658300768e-07), // f=5, s= 0, t=  62
                    new astro_top_term_t(     1406, -3.7816142635595027e-07,  9.1996049592574050e-08), // f=5, s= 0, t=  63
                    new astro_top_term_t(     1044,  7.4671401050977241e-08,  3.7193907363226299e-07), // f=5, s= 0, t=  64
                    new astro_top_term_t(      318,  7.1927751245935557e-08,  3.6928457890682560e-07), // f=5, s= 0, t=  65
                    new astro_top_term_t(      416,  4.1798174498739052e-08, -3.6796074151811958e-07), // f=5, s= 0, t=  66
                    new astro_top_term_t(      204, -1.5624747021437499e-08, -3.6071195814331651e-07), // f=5, s= 0, t=  67
                    new astro_top_term_t(     1539,  3.3199296973509958e-07, -1.0438655407426570e-07), // f=5, s= 0, t=  68
                    new astro_top_term_t(     1547,  2.1465474690281490e-07,  2.6396972342688539e-07), // f=5, s= 0, t=  69
                    new astro_top_term_t(     1327, -1.5719732197123969e-07,  2.8693854391332578e-07), // f=5, s= 0, t=  70
                    new astro_top_term_t(      541,  2.6198971229951878e-07,  1.7511268865594321e-07)  // f=5, s= 0, t=  71
                },
                new astro_top_term_t[]  // f=5, s= 1
                {
                    new astro_top_term_t(        0,  7.7993297015907634e-05,  0.0000000000000000e+00), // f=5, s= 1, t=   0
                    new astro_top_term_t(     1402,  3.9152687150310170e-05, -3.7172486736475712e-05), // f=5, s= 1, t=   1
                    new astro_top_term_t(     1543,  4.3031136375758243e-05,  1.3517906141260201e-05), // f=5, s= 1, t=   2
                    new astro_top_term_t(     1331,  1.0005173914454110e-05, -5.5510992766576859e-06), // f=5, s= 1, t=   3
                    new astro_top_term_t(      522,  5.4799934835759310e-06, -9.9803123122233115e-06), // f=5, s= 1, t=   4
                    new astro_top_term_t(     1473, -9.2379228749259814e-06,  2.0049455605464851e-06), // f=5, s= 1, t=   5
                    new astro_top_term_t(      664, -3.2895419619690672e-06, -4.5974891215164072e-08), // f=5, s= 1, t=   6
                    new astro_top_term_t(        8,  1.3908106265782320e-06, -2.9292856582721550e-06), // f=5, s= 1, t=   7
                    new astro_top_term_t(     2875,  1.5562630279974180e-06, -2.5756858183805818e-06), // f=5, s= 1, t=   8
                    new astro_top_term_t(        4, -2.4245388783587939e-06, -1.5068482062721909e-06), // f=5, s= 1, t=   9
                    new astro_top_term_t(     3016,  2.7144173403334680e-06,  1.5168858882407361e-07), // f=5, s= 1, t=  10
                    new astro_top_term_t(     2945, -1.9556025171954529e-06,  1.2707620314125440e-06), // f=5, s= 1, t=  11
                    new astro_top_term_t(       35,  1.3344840689621260e-06,  1.4001396639490341e-06), // f=5, s= 1, t=  12
                    new astro_top_term_t(      593, -1.5045227426363480e-06,  9.0293863845346919e-07), // f=5, s= 1, t=  13
                    new astro_top_term_t(     1614,  1.4931074066475380e-06,  1.3901882598918031e-07), // f=5, s= 1, t=  14
                    new astro_top_term_t(     1261,  1.3653087138211129e-06, -3.4669288959066708e-07), // f=5, s= 1, t=  15
                    new astro_top_term_t(       12,  1.2315835838874919e-06, -5.9561839689487617e-08), // f=5, s= 1, t=  16
                    new astro_top_term_t(      137,  9.0428610682743405e-07,  5.7865365381345505e-07), // f=5, s= 1, t=  17
                    new astro_top_term_t(     2804,  5.6172637460901586e-07, -4.7331838000806132e-07), // f=5, s= 1, t=  18
                    new astro_top_term_t(      212, -6.0813020695312143e-07,  1.1354550736721950e-07)  // f=5, s= 1, t=  19
                },
                new astro_top_term_t[]  // f=5, s= 2
                {
                    new astro_top_term_t(     1402, -2.9412500484509639e-06, -8.0598049136548207e-06), // f=5, s= 2, t=   0
                    new astro_top_term_t(     1543, -1.1369975198335059e-06, -6.6800431913703534e-06), // f=5, s= 2, t=   1
                    new astro_top_term_t(     1331, -5.6564967655275266e-07, -2.8456118566588939e-06), // f=5, s= 2, t=   2
                    new astro_top_term_t(        0,  2.4849999999999999e-06,  0.0000000000000000e+00), // f=5, s= 2, t=   3
                    new astro_top_term_t(      522, -1.0188268325841671e-06, -1.4049711115470319e-06), // f=5, s= 2, t=   4
                    new astro_top_term_t(        4, -1.2789041308060199e-06, -7.7265243668823065e-07), // f=5, s= 2, t=   5
                    new astro_top_term_t(        8,  1.2656119872224440e-06, -4.3448350220923409e-07)  // f=5, s= 2, t=   6
                }
            }
        };



        private static readonly double[] top_freq = new double[]
        {
            0.5296909622785881e+03,
            0.2132990811942489e+03,
            0.7478166163181234e+02,
            0.3813297236217556e+02,
            0.2533566020437000e+02
        };

        private static double calc_elliptical_coord(astro_top_term_t[][] formula, double dmu, int f, double t1)
        {
            double el = 0.0;
            double tpower = 1.0;
            for (int s=0; s < formula.Length; ++s)
            {
                astro_top_term_t[] series = formula[s];
                for (int t=0; t < series.Length; ++t)
                {
                    astro_top_term_t term = series[t];
                    if (f==1 && s==1 && term.k==0)
                        continue;
                    double arg = term.k * dmu * t1;
                    el += tpower * (term.c*Math.Cos(arg) + term.s*Math.Sin(arg));
                }
                tpower *= t1;
            }
            return el;
        }

        private static top_elliptical_t TopCalcElliptical(int planet, astro_top_term_t[][][] model, double tt)
        {
            /* Translated from: TOP2013.f */
            /* See: https://github.com/cosinekitty/ephemeris/tree/master/top2013 */
            /* Copied from: ftp://ftp.imcce.fr/pub/ephem/planets/top2013 */
            double t1 = tt / 365250.0;
            double dmu = (top_freq[0] - top_freq[1]) / 880.0;

            var ellip = new top_elliptical_t();
            ellip.a = calc_elliptical_coord(model[0], dmu, 0, t1);
            ellip.lambda = calc_elliptical_coord(model[1], dmu, 1, t1);
            ellip.k = calc_elliptical_coord(model[2], dmu, 2, t1);
            ellip.h = calc_elliptical_coord(model[3], dmu, 3, t1);
            ellip.q = calc_elliptical_coord(model[4], dmu, 4, t1);
            ellip.p = calc_elliptical_coord(model[5], dmu, 5, t1);

            double xl = ellip.lambda + top_freq[planet - 5] * t1;
            xl %= PI2;
            if (xl < 0.0)
                xl += PI2;
            ellip.lambda = xl;
            return ellip;
        }

        static AstroVector TopEcliptic(top_elliptical_t ellip, AstroTime time)
        {
            double xa = ellip.a;
            double xl = ellip.lambda;
            double xk = ellip.k;
            double xh = ellip.h;
            double xq = ellip.q;
            double xp = ellip.p;

            double xfi = Math.Sqrt(1.0 - xk*xk - xh*xh);
            double xki = Math.Sqrt(1.0 - xq*xq - xp*xp);
            double zr = xk;
            double zi = xh;
            double u = 1.0 / (1.0 + xfi);
            double ex2 = zr*zr + zi*zi;
            double ex = Math.Sqrt(ex2);
            double ex3 = ex * ex2;
            double z1r = zr;
            double z1i = -zi;

            double gl = xl % PI2;
            double gm = gl - Math.Atan2(xh, xk);
            double e = gl + (ex - 0.125*ex3)*Math.Sin(gm) + 0.5*ex2*Math.Sin(2.0*gm) + 0.375*ex3*Math.Sin(3.0*gm);

            double dl, z2r, z2i, zteta_r, zteta_i, z3r, z3i, rsa;
            do
            {
                z2r = 0.0;
                z2i = e;
                zteta_r = Math.Cos(z2i);
                zteta_i = Math.Sin(z2i);
                z3r = z1r*zteta_r - z1i*zteta_i;
                z3i = z1r*zteta_i + z1i*zteta_r;
                dl = gl - e + z3i;
                rsa = 1.0 - z3r;
                e += dl/rsa;
            } while (Math.Abs(dl) >= 1.0e-15);

            z1r = z3i * u * zr;
            z1i = z3i * u * zi;
            z2r = +z1i;
            z2i = -z1r;
            double zto_r = (-zr + zteta_r + z2r) / rsa;
            double zto_i = (-zi + zteta_i + z2i) / rsa;
            double xcw = zto_r;
            double xsw = zto_i;
            double xm = xp*xcw - xq*xsw;
            double xr = xa*rsa;

            return new AstroVector(
                xr*(xcw - 2.0*xp*xm),
                xr*(xsw + 2.0*xq*xm),
                -2.0*xr*xki*xm,
                time
            );
        }

        private static double[,] top_eq_rot;

        private static AstroVector TopEquatorial(AstroVector ecl)
        {
            if (top_eq_rot == null)
            {
                const double sdrad = DEG2RAD / 3600.0;
                const double eps = (23.0 + 26.0/60.0 + 21.41136/3600.0)*DEG2RAD;
                const double phi = -0.05188 * sdrad;
                double ceps = Math.Cos(eps);
                double seps = Math.Sin(eps);
                double cphi = Math.Cos(phi);
                double sphi = Math.Sin(phi);

                top_eq_rot = new double[3,3];

                top_eq_rot[0,0] =  cphi;
                top_eq_rot[0,1] = -sphi*ceps;
                top_eq_rot[0,2] =  sphi*seps;
                top_eq_rot[1,0] =  sphi;
                top_eq_rot[1,1] =  cphi*ceps;
                top_eq_rot[1,2] = -cphi*seps;
                top_eq_rot[2,0] =  0.0;
                top_eq_rot[2,1] =  seps;
                top_eq_rot[2,2] =  ceps;
            }

            return new AstroVector(
                (top_eq_rot[0,0] * ecl.x) + (top_eq_rot[0,1] * ecl.y) + (top_eq_rot[0,2] * ecl.z),
                (top_eq_rot[1,0] * ecl.x) + (top_eq_rot[1,1] * ecl.y) + (top_eq_rot[1,2] * ecl.z),
                (top_eq_rot[2,0] * ecl.x) + (top_eq_rot[2,1] * ecl.y) + (top_eq_rot[2,2] * ecl.z),
                ecl.t
            );
        }

        private static AstroVector TopPosition(astro_top_term_t[][][] model, int planet, AstroTime time)
        {
            top_elliptical_t ellip = TopCalcElliptical(planet, model, time.tt);
            AstroVector ecl = TopEcliptic(ellip, time);
            return TopEquatorial(ecl);
        }

        private static RotationMatrix precession_rot(double tt1, double tt2)
        {
            double xx, yx, zx, xy, yy, zy, xz, yz, zz;
            double t, psia, omegaa, chia, sa, ca, sb, cb, sc, cc, sd, cd;
            double eps0 = 84381.406;

            if ((tt1 != 0.0) && (tt2 != 0.0))
                throw new ArgumentException("precession_rot: one of (tt1, tt2) must be zero.");

            t = (tt2 - tt1) / 36525;
            if (tt2 == 0)
                t = -t;

            psia   = (((((-    0.0000000951  * t
                        +    0.000132851 ) * t
                        -    0.00114045  ) * t
                        -    1.0790069   ) * t
                        + 5038.481507    ) * t);

            omegaa = (((((+    0.0000003337  * t
                        -    0.000000467 ) * t
                        -    0.00772503  ) * t
                        +    0.0512623   ) * t
                        -    0.025754    ) * t + eps0);

            chia   = (((((-    0.0000000560  * t
                        +    0.000170663 ) * t
                        -    0.00121197  ) * t
                        -    2.3814292   ) * t
                        +   10.556403    ) * t);

            eps0 = eps0 * ASEC2RAD;
            psia = psia * ASEC2RAD;
            omegaa = omegaa * ASEC2RAD;
            chia = chia * ASEC2RAD;

            sa = Math.Sin(eps0);
            ca = Math.Cos(eps0);
            sb = Math.Sin(-psia);
            cb = Math.Cos(-psia);
            sc = Math.Sin(-omegaa);
            cc = Math.Cos(-omegaa);
            sd = Math.Sin(chia);
            cd = Math.Cos(chia);

            xx =  cd * cb - sb * sd * cc;
            yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc;
            zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc;
            xy = -sd * cb - sb * cd * cc;
            yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc;
            zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc;
            xz =  sb * sc;
            yz = -sc * cb * ca - sa * cc;
            zz = -sc * cb * sa + cc * ca;

            var rot = new double[3,3];
            if (tt2 == 0.0)
            {
                /* Perform rotation from other epoch to J2000.0. */
                rot[0, 0] = xx;
                rot[0, 1] = yx;
                rot[0, 2] = zx;
                rot[1, 0] = xy;
                rot[1, 1] = yy;
                rot[1, 2] = zy;
                rot[2, 0] = xz;
                rot[2, 1] = yz;
                rot[2, 2] = zz;
            }
            else
            {
                /* Perform rotation from J2000.0 to other epoch. */
                rot[0, 0] = xx;
                rot[0, 1] = xy;
                rot[0, 2] = xz;
                rot[1, 0] = yx;
                rot[1, 1] = yy;
                rot[1, 2] = yz;
                rot[2, 0] = zx;
                rot[2, 1] = zy;
                rot[2, 2] = zz;
            }

            return new RotationMatrix(rot);
        }

        private static AstroVector precession(double tt1, AstroVector pos, double tt2)
        {
            RotationMatrix r = precession_rot(tt1, tt2);
            return new AstroVector(
                r.rot[0, 0]*pos.x + r.rot[1, 0]*pos.y + r.rot[2, 0]*pos.z,
                r.rot[0, 1]*pos.x + r.rot[1, 1]*pos.y + r.rot[2, 1]*pos.z,
                r.rot[0, 2]*pos.x + r.rot[1, 2]*pos.y + r.rot[2, 2]*pos.z,
                null
            );
        }

        private struct earth_tilt_t
        {
            public double tt;
            public double dpsi;
            public double deps;
            public double ee;
            public double mobl;
            public double tobl;

            public earth_tilt_t(double tt, double dpsi, double deps, double ee, double mobl, double tobl)
            {
                this.tt = tt;
                this.dpsi = dpsi;
                this.deps = deps;
                this.ee = ee;
                this.mobl = mobl;
                this.tobl = tobl;
            }
        }

        private struct iau_row_t
        {
            public int nals0;
            public int nals1;
            public int nals2;
            public int nals3;
            public int nals4;

            public double cls0;
            public double cls1;
            public double cls2;
            public double cls3;
            public double cls4;
            public double cls5;
        }

        private static readonly iau_row_t[] iau_row = new iau_row_t[]
        {

        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  1 , cls0 =   -172064161, cls1 =      -174666, cls2 =        33386, cls3 =     92052331, cls4 =         9086, cls5 =        15377 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  2, nals3 = -2, nals4 =  2 , cls0 =    -13170906, cls1 =        -1675, cls2 =       -13696, cls3 =      5730336, cls4 =        -3015, cls5 =        -4587 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =     -2276413, cls1 =         -234, cls2 =         2796, cls3 =       978459, cls4 =         -485, cls5 =         1374 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  2 , cls0 =      2074554, cls1 =          207, cls2 =         -698, cls3 =      -897492, cls4 =          470, cls5 =         -291 },
        new iau_row_t { nals0 =  0, nals1 =  1, nals2 =  0, nals3 =  0, nals4 =  0 , cls0 =      1475877, cls1 =        -3633, cls2 =        11817, cls3 =        73871, cls4 =         -184, cls5 =        -1924 },
        new iau_row_t { nals0 =  0, nals1 =  1, nals2 =  2, nals3 = -2, nals4 =  2 , cls0 =      -516821, cls1 =         1226, cls2 =         -524, cls3 =       224386, cls4 =         -677, cls5 =         -174 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  0 , cls0 =       711159, cls1 =           73, cls2 =         -872, cls3 =        -6750, cls4 =            0, cls5 =          358 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  1 , cls0 =      -387298, cls1 =         -367, cls2 =          380, cls3 =       200728, cls4 =           18, cls5 =          318 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =      -301461, cls1 =          -36, cls2 =          816, cls3 =       129025, cls4 =          -63, cls5 =          367 },
        new iau_row_t { nals0 =  0, nals1 = -1, nals2 =  2, nals3 = -2, nals4 =  2 , cls0 =       215829, cls1 =         -494, cls2 =          111, cls3 =       -95929, cls4 =          299, cls5 =          132 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  2, nals3 = -2, nals4 =  1 , cls0 =       128227, cls1 =          137, cls2 =          181, cls3 =       -68982, cls4 =           -9, cls5 =           39 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =       123457, cls1 =           11, cls2 =           19, cls3 =       -53311, cls4 =           32, cls5 =           -4 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  0, nals3 =  2, nals4 =  0 , cls0 =       156994, cls1 =           10, cls2 =         -168, cls3 =        -1235, cls4 =            0, cls5 =           82 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  1 , cls0 =        63110, cls1 =           63, cls2 =           27, cls3 =       -33228, cls4 =            0, cls5 =           -9 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  1 , cls0 =       -57976, cls1 =          -63, cls2 =         -189, cls3 =        31429, cls4 =            0, cls5 =          -75 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  2, nals3 =  2, nals4 =  2 , cls0 =       -59641, cls1 =          -11, cls2 =          149, cls3 =        25543, cls4 =          -11, cls5 =           66 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  1 , cls0 =       -51613, cls1 =          -42, cls2 =          129, cls3 =        26366, cls4 =            0, cls5 =           78 },
        new iau_row_t { nals0 = -2, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  1 , cls0 =        45893, cls1 =           50, cls2 =           31, cls3 =       -24236, cls4 =          -10, cls5 =           20 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  0, nals3 =  2, nals4 =  0 , cls0 =        63384, cls1 =           11, cls2 =         -150, cls3 =        -1220, cls4 =            0, cls5 =           29 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  2, nals3 =  2, nals4 =  2 , cls0 =       -38571, cls1 =           -1, cls2 =          158, cls3 =        16452, cls4 =          -11, cls5 =           68 },
        new iau_row_t { nals0 =  0, nals1 = -2, nals2 =  2, nals3 = -2, nals4 =  2 , cls0 =        32481, cls1 =            0, cls2 =            0, cls3 =       -13870, cls4 =            0, cls5 =            0 },
        new iau_row_t { nals0 = -2, nals1 =  0, nals2 =  0, nals3 =  2, nals4 =  0 , cls0 =       -47722, cls1 =            0, cls2 =          -18, cls3 =          477, cls4 =            0, cls5 =          -25 },
        new iau_row_t { nals0 =  2, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =       -31046, cls1 =           -1, cls2 =          131, cls3 =        13238, cls4 =          -11, cls5 =           59 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  2, nals3 = -2, nals4 =  2 , cls0 =        28593, cls1 =            0, cls2 =           -1, cls3 =       -12338, cls4 =           10, cls5 =           -3 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  1 , cls0 =        20441, cls1 =           21, cls2 =           10, cls3 =       -10758, cls4 =            0, cls5 =           -3 },
        new iau_row_t { nals0 =  2, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  0 , cls0 =        29243, cls1 =            0, cls2 =          -74, cls3 =         -609, cls4 =            0, cls5 =           13 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  0 , cls0 =        25887, cls1 =            0, cls2 =          -66, cls3 =         -550, cls4 =            0, cls5 =           11 },
        new iau_row_t { nals0 =  0, nals1 =  1, nals2 =  0, nals3 =  0, nals4 =  1 , cls0 =       -14053, cls1 =          -25, cls2 =           79, cls3 =         8551, cls4 =           -2, cls5 =          -45 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  0, nals3 =  2, nals4 =  1 , cls0 =        15164, cls1 =           10, cls2 =           11, cls3 =        -8001, cls4 =            0, cls5 =           -1 },
        new iau_row_t { nals0 =  0, nals1 =  2, nals2 =  2, nals3 = -2, nals4 =  2 , cls0 =       -15794, cls1 =           72, cls2 =          -16, cls3 =         6850, cls4 =          -42, cls5 =           -5 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 = -2, nals3 =  2, nals4 =  0 , cls0 =        21783, cls1 =            0, cls2 =           13, cls3 =         -167, cls4 =            0, cls5 =           13 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  0, nals3 = -2, nals4 =  1 , cls0 =       -12873, cls1 =          -10, cls2 =          -37, cls3 =         6953, cls4 =            0, cls5 =          -14 },
        new iau_row_t { nals0 =  0, nals1 = -1, nals2 =  0, nals3 =  0, nals4 =  1 , cls0 =       -12654, cls1 =           11, cls2 =           63, cls3 =         6415, cls4 =            0, cls5 =           26 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  2, nals3 =  2, nals4 =  1 , cls0 =       -10204, cls1 =            0, cls2 =           25, cls3 =         5222, cls4 =            0, cls5 =           15 },
        new iau_row_t { nals0 =  0, nals1 =  2, nals2 =  0, nals3 =  0, nals4 =  0 , cls0 =        16707, cls1 =          -85, cls2 =          -10, cls3 =          168, cls4 =           -1, cls5 =           10 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  2, nals3 =  2, nals4 =  2 , cls0 =        -7691, cls1 =            0, cls2 =           44, cls3 =         3268, cls4 =            0, cls5 =           19 },
        new iau_row_t { nals0 = -2, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  0 , cls0 =       -11024, cls1 =            0, cls2 =          -14, cls3 =          104, cls4 =            0, cls5 =            2 },
        new iau_row_t { nals0 =  0, nals1 =  1, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =         7566, cls1 =          -21, cls2 =          -11, cls3 =        -3250, cls4 =            0, cls5 =           -5 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  2, nals3 =  2, nals4 =  1 , cls0 =        -6637, cls1 =          -11, cls2 =           25, cls3 =         3353, cls4 =            0, cls5 =           14 },
        new iau_row_t { nals0 =  0, nals1 = -1, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =        -7141, cls1 =           21, cls2 =            8, cls3 =         3070, cls4 =            0, cls5 =            4 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  0, nals3 =  2, nals4 =  1 , cls0 =        -6302, cls1 =          -11, cls2 =            2, cls3 =         3272, cls4 =            0, cls5 =            4 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  2, nals3 = -2, nals4 =  1 , cls0 =         5800, cls1 =           10, cls2 =            2, cls3 =        -3045, cls4 =            0, cls5 =           -1 },
        new iau_row_t { nals0 =  2, nals1 =  0, nals2 =  2, nals3 = -2, nals4 =  2 , cls0 =         6443, cls1 =            0, cls2 =           -7, cls3 =        -2768, cls4 =            0, cls5 =           -4 },
        new iau_row_t { nals0 = -2, nals1 =  0, nals2 =  0, nals3 =  2, nals4 =  1 , cls0 =        -5774, cls1 =          -11, cls2 =          -15, cls3 =         3041, cls4 =            0, cls5 =           -5 },
        new iau_row_t { nals0 =  2, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  1 , cls0 =        -5350, cls1 =            0, cls2 =           21, cls3 =         2695, cls4 =            0, cls5 =           12 },
        new iau_row_t { nals0 =  0, nals1 = -1, nals2 =  2, nals3 = -2, nals4 =  1 , cls0 =        -4752, cls1 =          -11, cls2 =           -3, cls3 =         2719, cls4 =            0, cls5 =           -3 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  0, nals3 = -2, nals4 =  1 , cls0 =        -4940, cls1 =          -11, cls2 =          -21, cls3 =         2720, cls4 =            0, cls5 =           -9 },
        new iau_row_t { nals0 = -1, nals1 = -1, nals2 =  0, nals3 =  2, nals4 =  0 , cls0 =         7350, cls1 =            0, cls2 =           -8, cls3 =          -51, cls4 =            0, cls5 =            4 },
        new iau_row_t { nals0 =  2, nals1 =  0, nals2 =  0, nals3 = -2, nals4 =  1 , cls0 =         4065, cls1 =            0, cls2 =            6, cls3 =        -2206, cls4 =            0, cls5 =            1 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  0, nals3 =  2, nals4 =  0 , cls0 =         6579, cls1 =            0, cls2 =          -24, cls3 =         -199, cls4 =            0, cls5 =            2 },
        new iau_row_t { nals0 =  0, nals1 =  1, nals2 =  2, nals3 = -2, nals4 =  1 , cls0 =         3579, cls1 =            0, cls2 =            5, cls3 =        -1900, cls4 =            0, cls5 =            1 },
        new iau_row_t { nals0 =  1, nals1 = -1, nals2 =  0, nals3 =  0, nals4 =  0 , cls0 =         4725, cls1 =            0, cls2 =           -6, cls3 =          -41, cls4 =            0, cls5 =            3 },
        new iau_row_t { nals0 = -2, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =        -3075, cls1 =            0, cls2 =           -2, cls3 =         1313, cls4 =            0, cls5 =           -1 },
        new iau_row_t { nals0 =  3, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =        -2904, cls1 =            0, cls2 =           15, cls3 =         1233, cls4 =            0, cls5 =            7 },
        new iau_row_t { nals0 =  0, nals1 = -1, nals2 =  0, nals3 =  2, nals4 =  0 , cls0 =         4348, cls1 =            0, cls2 =          -10, cls3 =          -81, cls4 =            0, cls5 =            2 },
        new iau_row_t { nals0 =  1, nals1 = -1, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =        -2878, cls1 =            0, cls2 =            8, cls3 =         1232, cls4 =            0, cls5 =            4 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  0, nals3 =  1, nals4 =  0 , cls0 =        -4230, cls1 =            0, cls2 =            5, cls3 =          -20, cls4 =            0, cls5 =           -2 },
        new iau_row_t { nals0 = -1, nals1 = -1, nals2 =  2, nals3 =  2, nals4 =  2 , cls0 =        -2819, cls1 =            0, cls2 =            7, cls3 =         1207, cls4 =            0, cls5 =            3 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  0 , cls0 =        -4056, cls1 =            0, cls2 =            5, cls3 =           40, cls4 =            0, cls5 =           -2 },
        new iau_row_t { nals0 =  0, nals1 = -1, nals2 =  2, nals3 =  2, nals4 =  2 , cls0 =        -2647, cls1 =            0, cls2 =           11, cls3 =         1129, cls4 =            0, cls5 =            5 },
        new iau_row_t { nals0 = -2, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  1 , cls0 =        -2294, cls1 =            0, cls2 =          -10, cls3 =         1266, cls4 =            0, cls5 =           -4 },
        new iau_row_t { nals0 =  1, nals1 =  1, nals2 =  2, nals3 =  0, nals4 =  2 , cls0 =         2481, cls1 =            0, cls2 =           -7, cls3 =        -1062, cls4 =            0, cls5 =           -3 },
        new iau_row_t { nals0 =  2, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  1 , cls0 =         2179, cls1 =            0, cls2 =           -2, cls3 =        -1129, cls4 =            0, cls5 =           -2 },
        new iau_row_t { nals0 = -1, nals1 =  1, nals2 =  0, nals3 =  1, nals4 =  0 , cls0 =         3276, cls1 =            0, cls2 =            1, cls3 =           -9, cls4 =            0, cls5 =            0 },
        new iau_row_t { nals0 =  1, nals1 =  1, nals2 =  0, nals3 =  0, nals4 =  0 , cls0 =        -3389, cls1 =            0, cls2 =            5, cls3 =           35, cls4 =            0, cls5 =           -2 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  2, nals3 =  0, nals4 =  0 , cls0 =         3339, cls1 =            0, cls2 =          -13, cls3 =         -107, cls4 =            0, cls5 =            1 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  2, nals3 = -2, nals4 =  1 , cls0 =        -1987, cls1 =            0, cls2 =           -6, cls3 =         1073, cls4 =            0, cls5 =           -2 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  2 , cls0 =        -1981, cls1 =            0, cls2 =            0, cls3 =          854, cls4 =            0, cls5 =            0 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  0, nals3 =  1, nals4 =  0 , cls0 =         4026, cls1 =            0, cls2 =         -353, cls3 =         -553, cls4 =            0, cls5 =         -139 },
        new iau_row_t { nals0 =  0, nals1 =  0, nals2 =  2, nals3 =  1, nals4 =  2 , cls0 =         1660, cls1 =            0, cls2 =           -5, cls3 =         -710, cls4 =            0, cls5 =           -2 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  2, nals3 =  4, nals4 =  2 , cls0 =        -1521, cls1 =            0, cls2 =            9, cls3 =          647, cls4 =            0, cls5 =            4 },
        new iau_row_t { nals0 = -1, nals1 =  1, nals2 =  0, nals3 =  1, nals4 =  1 , cls0 =         1314, cls1 =            0, cls2 =            0, cls3 =         -700, cls4 =            0, cls5 =            0 },
        new iau_row_t { nals0 =  0, nals1 = -2, nals2 =  2, nals3 = -2, nals4 =  1 , cls0 =        -1283, cls1 =            0, cls2 =            0, cls3 =          672, cls4 =            0, cls5 =            0 },
        new iau_row_t { nals0 =  1, nals1 =  0, nals2 =  2, nals3 =  2, nals4 =  1 , cls0 =        -1331, cls1 =            0, cls2 =            8, cls3 =          663, cls4 =            0, cls5 =            4 },
        new iau_row_t { nals0 = -2, nals1 =  0, nals2 =  2, nals3 =  2, nals4 =  2 , cls0 =         1383, cls1 =            0, cls2 =           -2, cls3 =         -594, cls4 =            0, cls5 =           -2 },
        new iau_row_t { nals0 = -1, nals1 =  0, nals2 =  0, nals3 =  0, nals4 =  2 , cls0 =         1405, cls1 =            0, cls2 =            4, cls3 =         -610, cls4 =            0, cls5 =            2 },
        new iau_row_t { nals0 =  1, nals1 =  1, nals2 =  2, nals3 = -2, nals4 =  2 , cls0 =         1290, cls1 =            0, cls2 =            0, cls3 =         -556, cls4 =            0, cls5 =            0 }

        };

        private static void iau2000b(AstroTime time)
        {
            /* Adapted from the NOVAS C 3.1 function of the same name. */

            double t, el, elp, f, d, om, arg, dp, de, sarg, carg;
            int i;

            if (double.IsNaN(time.psi))
            {
                t = time.tt / 36525.0;
                el  = ((485868.249036 + t * 1717915923.2178) % ASEC360) * ASEC2RAD;
                elp = ((1287104.79305 + t * 129596581.0481)  % ASEC360) * ASEC2RAD;
                f   = ((335779.526232 + t * 1739527262.8478) % ASEC360) * ASEC2RAD;
                d   = ((1072260.70369 + t * 1602961601.2090) % ASEC360) * ASEC2RAD;
                om  = ((450160.398036 - t * 6962890.5431)    % ASEC360) * ASEC2RAD;
                dp = 0;
                de = 0;
                for (i=76; i >= 0; --i)
                {
                    arg = (iau_row[i].nals0*el + iau_row[i].nals1*elp + iau_row[i].nals2*f + iau_row[i].nals3*d + iau_row[i].nals4*om) % PI2;
                    sarg = Math.Sin(arg);
                    carg = Math.Cos(arg);
                    dp += (iau_row[i].cls0 + iau_row[i].cls1*t) * sarg + iau_row[i].cls2*carg;
                    de += (iau_row[i].cls3 + iau_row[i].cls4*t) * carg + iau_row[i].cls5*sarg;
                }

                time.psi = -0.000135 + (dp * 1.0e-7);
                time.eps = +0.000388 + (de * 1.0e-7);
            }
        }

        private static double mean_obliq(double tt)
        {
            double t = tt / 36525.0;
            double asec =
                (((( -  0.0000000434   * t
                    -  0.000000576  ) * t
                    +  0.00200340   ) * t
                    -  0.0001831    ) * t
                    - 46.836769     ) * t + 84381.406;

            return asec / 3600.0;
        }

        private static earth_tilt_t e_tilt(AstroTime time)
        {
            iau2000b(time);

            double mobl = mean_obliq(time.tt);
            double tobl = mobl + (time.eps / 3600.0);
            double ee = time.psi * Math.Cos(mobl * DEG2RAD) / 15.0;
            return new earth_tilt_t(time.tt, time.psi, time.eps, ee, mobl, tobl);
        }

        private static double era(double ut)        /* Earth Rotation Angle */
        {
            double thet1 = 0.7790572732640 + 0.00273781191135448 * ut;
            double thet3 = ut % 1.0;
            double theta = 360.0 *((thet1 + thet3) % 1.0);
            if (theta < 0.0)
                theta += 360.0;

            return theta;
        }

        private static double sidereal_time(AstroTime time)
        {
            double t = time.tt / 36525.0;
            double eqeq = 15.0 * e_tilt(time).ee;    /* Replace with eqeq=0 to get GMST instead of GAST (if we ever need it) */
            double theta = era(time.ut);
            double st = (eqeq + 0.014506 +
                (((( -    0.0000000368   * t
                    -    0.000029956  ) * t
                    -    0.00000044   ) * t
                    +    1.3915817    ) * t
                    + 4612.156534     ) * t);

            double gst = ((st/3600.0 + theta) % 360.0) / 15.0;
            if (gst < 0.0)
                gst += 24.0;

            return gst;
        }

        private static AstroVector terra(Observer observer, double st)
        {
            double df = 1.0 - 0.003352819697896;    /* flattening of the Earth */
            double df2 = df * df;
            double phi = observer.latitude * DEG2RAD;
            double sinphi = Math.Sin(phi);
            double cosphi = Math.Cos(phi);
            double c = 1.0 / Math.Sqrt(cosphi*cosphi + df2*sinphi*sinphi);
            double s = df2 * c;
            double ht_km = observer.height / 1000.0;
            double ach = EARTH_EQUATORIAL_RADIUS_KM*c + ht_km;
            double ash = EARTH_EQUATORIAL_RADIUS_KM*s + ht_km;
            double stlocl = (15.0*st + observer.longitude) * DEG2RAD;
            double sinst = Math.Sin(stlocl);
            double cosst = Math.Cos(stlocl);

            return new AstroVector(
                ach * cosphi * cosst / KM_PER_AU,
                ach * cosphi * sinst / KM_PER_AU,
                ash * sinphi / KM_PER_AU,
                null
            );
        }

        private static RotationMatrix nutation_rot(AstroTime time, int direction)
        {
            earth_tilt_t tilt = e_tilt(time);
            double oblm = tilt.mobl * DEG2RAD;
            double oblt = tilt.tobl * DEG2RAD;
            double psi = tilt.dpsi * ASEC2RAD;
            double cobm = Math.Cos(oblm);
            double sobm = Math.Sin(oblm);
            double cobt = Math.Cos(oblt);
            double sobt = Math.Sin(oblt);
            double cpsi = Math.Cos(psi);
            double spsi = Math.Sin(psi);

            double xx = cpsi;
            double yx = -spsi * cobm;
            double zx = -spsi * sobm;
            double xy = spsi * cobt;
            double yy = cpsi * cobm * cobt + sobm * sobt;
            double zy = cpsi * sobm * cobt - cobm * sobt;
            double xz = spsi * sobt;
            double yz = cpsi * cobm * sobt - sobm * cobt;
            double zz = cpsi * sobm * sobt + cobm * cobt;

            var rot = new double[3,3];

            if (direction == 0)
            {
                /* forward rotation */
                rot[0, 0] = xx;
                rot[0, 1] = xy;
                rot[0, 2] = xz;
                rot[1, 0] = yx;
                rot[1, 1] = yy;
                rot[1, 2] = yz;
                rot[2, 0] = zx;
                rot[2, 1] = zy;
                rot[2, 2] = zz;
            }
            else
            {
                /* inverse rotation */
                rot[0, 0] = xx;
                rot[0, 1] = yx;
                rot[0, 2] = zx;
                rot[1, 0] = xy;
                rot[1, 1] = yy;
                rot[1, 2] = zy;
                rot[2, 0] = xz;
                rot[2, 1] = yz;
                rot[2, 2] = zz;
            }

            return new RotationMatrix(rot);
        }


        private static AstroVector nutation(AstroTime time, int direction, AstroVector pos)
        {
            RotationMatrix r = nutation_rot(time, direction);
            return new AstroVector(
                r.rot[0, 0]*pos.x + r.rot[1, 0]*pos.y + r.rot[2, 0]*pos.z,
                r.rot[0, 1]*pos.x + r.rot[1, 1]*pos.y + r.rot[2, 1]*pos.z,
                r.rot[0, 2]*pos.x + r.rot[1, 2]*pos.y + r.rot[2, 2]*pos.z,
                time
            );
        }


        private static Equatorial vector2radec(AstroVector pos)
        {
            double ra, dec, dist;
            double xyproj;

            xyproj = pos.x*pos.x + pos.y*pos.y;
            dist = Math.Sqrt(xyproj + pos.z*pos.z);
            if (xyproj == 0.0)
            {
                if (pos.z == 0.0)
                {
                    /* Indeterminate coordinates; pos vector has zero length. */
                    throw new ArgumentException("Bad vector");
                }

                if (pos.z < 0)
                {
                    ra = 0.0;
                    dec = -90.0;
                }
                else
                {
                    ra = 0.0;
                    dec = +90.0;
                }
            }
            else
            {
                ra = Math.Atan2(pos.y, pos.x) / (DEG2RAD * 15.0);
                if (ra < 0)
                    ra += 24.0;

                dec = RAD2DEG * Math.Atan2(pos.z, Math.Sqrt(xyproj));
            }

            return new Equatorial(ra, dec, dist);
        }

        private static AstroVector geo_pos(AstroTime time, Observer observer)
        {
            double gast = sidereal_time(time);
            AstroVector pos1 = terra(observer, gast);
            AstroVector pos2 = nutation(time, -1, pos1);
            return precession(time.tt, pos2, 0.0);
        }

        private static AstroVector spin(double angle, AstroVector pos)
        {
            double angr = angle * DEG2RAD;
            double cosang = Math.Cos(angr);
            double sinang = Math.Sin(angr);
            return new AstroVector(
                +cosang*pos.x + sinang*pos.y,
                -sinang*pos.x + cosang*pos.y,
                pos.z,
                null
            );
        }

        private static AstroVector ecl2equ_vec(AstroTime time, AstroVector ecl)
        {
            double obl = mean_obliq(time.tt) * DEG2RAD;
            double cos_obl = Math.Cos(obl);
            double sin_obl = Math.Sin(obl);

            return new AstroVector(
                ecl.x,
                ecl.y*cos_obl - ecl.z*sin_obl,
                ecl.y*sin_obl + ecl.z*cos_obl,
                time
            );
        }

        private static AstroVector GeoMoon(AstroTime time)
        {
            var context = new MoonContext(time.tt / 36525.0);
            MoonResult moon = context.CalcMoon();

            /* Convert geocentric ecliptic spherical coordinates to Cartesian coordinates. */
            double dist_cos_lat = moon.distance_au * Math.Cos(moon.geo_eclip_lat);

            var gepos = new AstroVector(
                dist_cos_lat * Math.Cos(moon.geo_eclip_lon),
                dist_cos_lat * Math.Sin(moon.geo_eclip_lon),
                moon.distance_au * Math.Sin(moon.geo_eclip_lat),
                null
            );

            /* Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date. */
            AstroVector mpos1 = ecl2equ_vec(time, gepos);

            /* Convert from mean equinox of date to J2000. */
            AstroVector mpos2 = precession(time.tt, mpos1, 0);

            /* Patch in the correct time value into the returned vector. */
            return new AstroVector(mpos2.x, mpos2.y, mpos2.z, time);
        }

        private static AstroVector BarycenterContrib(AstroTime time, Body body, double pmass)
        {
            double shift = pmass / (pmass + SUN_MASS);
            AstroVector p = CalcVsop(vsop[(int)body], time);
            return new AstroVector(
                shift * p.x,
                shift * p.y,
                shift * p.z,
                time
            );
        }

        private static AstroVector CalcSolarSystemBarycenter(AstroTime time)
        {
            AstroVector j = BarycenterContrib(time, Body.Jupiter, JUPITER_MASS);
            AstroVector s = BarycenterContrib(time, Body.Saturn,  SATURN_MASS);
            AstroVector u = BarycenterContrib(time, Body.Uranus,  URANUS_MASS);
            AstroVector n = BarycenterContrib(time, Body.Neptune, NEPTUNE_MASS);
            return new AstroVector(
                j.x + s.x + u.x + n.x,
                j.y + s.y + u.y + n.y,
                j.z + s.z + u.z + n.z,
                time
            );
        }

        /// <summary>
        /// Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.
        /// </summary>
        /// <remarks>
        /// This function calculates the position of the given celestial body as a vector,
        /// using the center of the Sun as the origin.  The result is expressed as a Cartesian
        /// vector in the J2000 equatorial system: the coordinates are based on the mean equator
        /// of the Earth at noon UTC on 1 January 2000.
        ///
        /// The position is not corrected for light travel time or aberration.
        /// This is different from the behavior of #Astronomy.GeoVector.
        ///
        /// If given an invalid value for `body`, or the body is `Body.Pluto` and the `time` is outside
        /// the year range 1700..2200, this function will throw an `ArgumentException`.
        /// </remarks>
        /// <param name="body">A body for which to calculate a heliocentric position: the Sun, Moon, EMB, SSB, or any of the planets.</param>
        /// <param name="time">The date and time for which to calculate the position.</param>
        /// <returns>A heliocentric position vector of the center of the given body.</returns>
        public static AstroVector HelioVector(Body body, AstroTime time)
        {
            AstroVector earth, geomoon;

            switch (body)
            {
                case Body.Sun:
                    return new AstroVector(0.0, 0.0, 0.0, time);

                case Body.Mercury:
                case Body.Venus:
                case Body.Earth:
                case Body.Mars:
                case Body.Jupiter:
                case Body.Saturn:
                case Body.Uranus:
                case Body.Neptune:
                    return CalcVsop(vsop[(int)body], time);

                case Body.Pluto:
                    return TopPosition(top_model_8, 9, time);

                case Body.Moon:
                    geomoon = GeoMoon(time);
                    earth = CalcEarth(time);
                    return new AstroVector(
                        earth.x + geomoon.x,
                        earth.y + geomoon.y,
                        earth.z + geomoon.z,
                        time
                    );

                case Body.EMB:
                    geomoon = GeoMoon(time);
                    earth = CalcEarth(time);
                    double denom = 1.0 + EARTH_MOON_MASS_RATIO;
                    return new AstroVector(
                        earth.x + (geomoon.x / denom),
                        earth.y + (geomoon.y / denom),
                        earth.z + (geomoon.z / denom),
                        time
                    );

                case Body.SSB:
                    return CalcSolarSystemBarycenter(time);

                default:
                    throw new InvalidBodyException(body);
            }
        }

        /// <summary>
        /// Calculates the distance between a body and the Sun at a given time.
        /// </summary>
        /// <remarks>
        /// Given a date and time, this function calculates the distance between
        /// the center of `body` and the center of the Sun.
        /// For the planets Mercury through Neptune, this function is significantly
        /// more efficient than calling #Astronomy.HelioVector followed by taking the length
        /// of the resulting vector.
        /// </remarks>
        /// <param name="body">
        /// A body for which to calculate a heliocentric distance:
        /// the Sun, Moon, or any of the planets.
        /// </param>
        /// <param name="time">
        /// The date and time for which to calculate the heliocentric distance.
        /// </param>
        /// <returns>
        /// The heliocentric distance in AU.
        /// </returns>
        public static double HelioDistance(Body body, AstroTime time)
        {
            switch (body)
            {
                case Body.Sun:
                    return 0.0;

                case Body.Mercury:
                case Body.Venus:
                case Body.Earth:
                case Body.Mars:
                case Body.Jupiter:
                case Body.Saturn:
                case Body.Uranus:
                case Body.Neptune:
                    return VsopHelioDistance(vsop[(int)body], time);

                default:
                    /* For non-VSOP objects, fall back to taking the length of the heliocentric vector. */
                    return HelioVector(body, time).Length();
            }
        }

        private static AstroVector CalcEarth(AstroTime time)
        {
            return CalcVsop(vsop[(int)Body.Earth], time);
        }

        ///
        /// <summary>
        /// Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.
        /// </summary>
        /// <remarks>
        /// This function calculates the position of the given celestial body as a vector,
        /// using the center of the Earth as the origin.  The result is expressed as a Cartesian
        /// vector in the J2000 equatorial system: the coordinates are based on the mean equator
        /// of the Earth at noon UTC on 1 January 2000.
        ///
        /// If given an invalid value for `body`, or the body is `Body.Pluto` and the `time` is outside
        /// the year range 1700..2200, this function will throw an exception.
        ///
        /// Unlike #Astronomy.HelioVector, this function always corrects for light travel time.
        /// This means the position of the body is "back-dated" by the amount of time it takes
        /// light to travel from that body to an observer on the Earth.
        ///
        /// Also, the position can optionally be corrected for
        /// [aberration](https://en.wikipedia.org/wiki/Aberration_of_light), an effect
        /// causing the apparent direction of the body to be shifted due to transverse
        /// movement of the Earth with respect to the rays of light coming from that body.
        /// </remarks>
        /// <param name="body">A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.</param>
        /// <param name="time">The date and time for which to calculate the position.</param>
        /// <param name="aberration">`Aberration.Corrected` to correct for aberration, or `Aberration.None` to leave uncorrected.</param>
        /// <returns>A geocentric position vector of the center of the given body.</returns>
        public static AstroVector GeoVector(
            Body body,
            AstroTime time,
            Aberration aberration)
        {
            AstroVector vector;
            AstroVector earth = new AstroVector(0.0, 0.0, 0.0, null);
            AstroTime ltime;
            AstroTime ltime2;
            double dt;
            int iter;

            if (aberration != Aberration.Corrected && aberration != Aberration.None)
                throw new ArgumentException(string.Format("Unsupported aberration option {0}", aberration));

            switch (body)
            {
            case Body.Earth:
                /* The Earth's geocentric coordinates are always (0,0,0). */
                return new AstroVector(0.0, 0.0, 0.0, time);

            case Body.Moon:
                return GeoMoon(time);

            default:
                /* For all other bodies, apply light travel time correction. */

                if (aberration == Aberration.None)
                {
                    /* No aberration, so calculate Earth's position once, at the time of observation. */
                    earth = CalcEarth(time);
                }

                ltime = time;
                for (iter=0; iter < 10; ++iter)
                {
                    vector = HelioVector(body, ltime);
                    if (aberration == Aberration.Corrected)
                    {
                        /*
                            Include aberration, so make a good first-order approximation
                            by backdating the Earth's position also.
                            This is confusing, but it works for objects within the Solar System
                            because the distance the Earth moves in that small amount of light
                            travel time (a few minutes to a few hours) is well approximated
                            by a line segment that substends the angle seen from the remote
                            body viewing Earth. That angle is pretty close to the aberration
                            angle of the moving Earth viewing the remote body.
                            In other words, both of the following approximate the aberration angle:
                                (transverse distance Earth moves) / (distance to body)
                                (transverse speed of Earth) / (speed of light).
                        */
                        earth = CalcEarth(ltime);
                    }

                    /* Convert heliocentric vector to geocentric vector. */
                    vector = new AstroVector(vector.x - earth.x, vector.y - earth.y, vector.z - earth.z, time);
                    ltime2 = time.AddDays(-vector.Length() / C_AUDAY);
                    dt = Math.Abs(ltime2.tt - ltime.tt);
                    if (dt < 1.0e-9)
                        return vector;

                    ltime = ltime2;
                }
                throw new Exception("Light travel time correction did not converge");
            }
        }

        /// <summary>
        /// Calculates equatorial coordinates of a celestial body as seen by an observer on the Earth's surface.
        /// </summary>
        /// <remarks>
        /// Calculates topocentric equatorial coordinates in one of two different systems:
        /// J2000 or true-equator-of-date, depending on the value of the `equdate` parameter.
        /// Equatorial coordinates include right ascension, declination, and distance in astronomical units.
        ///
        /// This function corrects for light travel time: it adjusts the apparent location
        /// of the observed body based on how long it takes for light to travel from the body to the Earth.
        ///
        /// This function corrects for *topocentric parallax*, meaning that it adjusts for the
        /// angular shift depending on where the observer is located on the Earth. This is most
        /// significant for the Moon, because it is so close to the Earth. However, parallax corection
        /// has a small effect on the apparent positions of other bodies.
        ///
        /// Correction for aberration is optional, using the `aberration` parameter.
        /// </remarks>
        /// <param name="body">The celestial body to be observed. Not allowed to be `Body.Earth`.</param>
        /// <param name="time">The date and time at which the observation takes place.</param>
        /// <param name="observer">A location on or near the surface of the Earth.</param>
        /// <param name="equdate">Selects the date of the Earth's equator in which to express the equatorial coordinates.</param>
        /// <param name="aberration">Selects whether or not to correct for aberration.</param>
        public static Equatorial Equator(
            Body body,
            AstroTime time,
            Observer observer,
            EquatorEpoch equdate,
            Aberration aberration)
        {
            AstroVector gc_observer = geo_pos(time, observer);
            AstroVector gc = GeoVector(body, time, aberration);
            AstroVector j2000 = new AstroVector(gc.x - gc_observer.x, gc.y - gc_observer.y, gc.z - gc_observer.z, time);

            switch (equdate)
            {
                case EquatorEpoch.OfDate:
                    AstroVector temp = precession(0.0, j2000, time.tt);
                    AstroVector datevect = nutation(time, 0, temp);
                    return vector2radec(datevect);

                case EquatorEpoch.J2000:
                    return vector2radec(j2000);

                default:
                    throw new ArgumentException(string.Format("Unsupported equator epoch {0}", equdate));
            }
        }

        /// <summary>
        /// Calculates the apparent location of a body relative to the local horizon of an observer on the Earth.
        /// </summary>
        /// <remarks>
        /// Given a date and time, the geographic location of an observer on the Earth, and
        /// equatorial coordinates (right ascension and declination) of a celestial body,
        /// this function returns horizontal coordinates (azimuth and altitude angles) for the body
        /// relative to the horizon at the geographic location.
        ///
        /// The right ascension `ra` and declination `dec` passed in must be *equator of date*
        /// coordinates, based on the Earth's true equator at the date and time of the observation.
        /// Otherwise the resulting horizontal coordinates will be inaccurate.
        /// Equator of date coordinates can be obtained by calling #Astronomy.Equator, passing in
        /// `EquatorEpoch.OfDate` as its `equdate` parameter. It is also recommended to enable
        /// aberration correction by passing in `Aberration.Corrected` as the `aberration` parameter.
        ///
        /// This function optionally corrects for atmospheric refraction.
        /// For most uses, it is recommended to pass `Refraction.Normal` in the `refraction` parameter to
        /// correct for optical lensing of the Earth's atmosphere that causes objects
        /// to appear somewhat higher above the horizon than they actually are.
        /// However, callers may choose to avoid this correction by passing in `Refraction.None`.
        /// If refraction correction is enabled, the azimuth, altitude, right ascension, and declination
        /// in the #Topocentric structure returned by this function will all be corrected for refraction.
        /// If refraction is disabled, none of these four coordinates will be corrected; in that case,
        /// the right ascension and declination in the returned structure will be numerically identical
        /// to the respective `ra` and `dec` values passed in.
        /// </remarks>
        /// <param name="time">The date and time of the observation.</param>
        /// <param name="observer">The geographic location of the observer.</param>
        /// <param name="ra">The right ascension of the body in sidereal hours. See remarks above for more details.</param>
        /// <param name="dec">The declination of the body in degrees. See remarks above for more details.</param>
        /// <param name="refraction">
        /// Selects whether to correct for atmospheric refraction, and if so, which model to use.
        /// The recommended value for most uses is `Refraction.Normal`.
        /// See remarks above for more details.
        /// </param>
        /// <returns>
        /// The body's apparent horizontal coordinates and equatorial coordinates, both optionally corrected for refraction.
        /// </returns>
        public static Topocentric Horizon(
            AstroTime time,
            Observer observer,
            double ra,
            double dec,
            Refraction refraction)
        {
            double sinlat = Math.Sin(observer.latitude * DEG2RAD);
            double coslat = Math.Cos(observer.latitude * DEG2RAD);
            double sinlon = Math.Sin(observer.longitude * DEG2RAD);
            double coslon = Math.Cos(observer.longitude * DEG2RAD);
            double sindc = Math.Sin(dec * DEG2RAD);
            double cosdc = Math.Cos(dec * DEG2RAD);
            double sinra = Math.Sin(ra * 15 * DEG2RAD);
            double cosra = Math.Cos(ra * 15 * DEG2RAD);

            var uze = new AstroVector(coslat * coslon, coslat * sinlon, sinlat, null);
            var une = new AstroVector(-sinlat * coslon, -sinlat * sinlon, coslat, null);
            var uwe = new AstroVector(sinlon, -coslon, 0.0, null);

            double spin_angle = -15.0 * sidereal_time(time);
            AstroVector uz = spin(spin_angle, uze);
            AstroVector un = spin(spin_angle, une);
            AstroVector uw = spin(spin_angle, uwe);

            var p = new AstroVector(cosdc * cosra, cosdc * sinra, sindc, null);
            double pz = p.x*uz.x + p.y*uz.y + p.z*uz.z;
            double pn = p.x*un.x + p.y*un.y + p.z*un.z;
            double pw = p.x*uw.x + p.y*uw.y + p.z*uw.z;

            double proj = Math.Sqrt(pn*pn + pw*pw);
            double az = 0.0;
            if (proj > 0.0)
            {
                az = -Math.Atan2(pw, pn) * RAD2DEG;
                if (az < 0.0)
                    az += 360.0;
                else if (az >= 360.0)
                    az -= 360.0;
            }
            double zd = Math.Atan2(proj, pz) * RAD2DEG;
            double hor_ra = ra;
            double hor_dec = dec;

            if (refraction == Refraction.Normal || refraction == Refraction.JplHor)
            {
                double zd0 = zd;
                double refr = RefractionAngle(refraction, 90.0 - zd);
                zd -= refr;

                if (refr > 0.0 && zd > 3.0e-4)
                {
                    double sinzd = Math.Sin(zd * DEG2RAD);
                    double coszd = Math.Cos(zd * DEG2RAD);
                    double sinzd0 = Math.Sin(zd0 * DEG2RAD);
                    double coszd0 = Math.Cos(zd0 * DEG2RAD);

                    double prx = ((p.x - coszd0 * uz.x) / sinzd0)*sinzd + uz.x*coszd;
                    double pry = ((p.y - coszd0 * uz.y) / sinzd0)*sinzd + uz.y*coszd;
                    double prz = ((p.z - coszd0 * uz.z) / sinzd0)*sinzd + uz.z*coszd;

                    proj = Math.Sqrt(prx*prx + pry*pry);
                    if (proj > 0.0)
                    {
                        hor_ra = Math.Atan2(pry, prx) * (RAD2DEG / 15.0);
                        if (hor_ra < 0.0)
                            hor_ra += 24.0;
                        else if (hor_ra >= 24.0)
                            hor_ra -= 24.0;
                    }
                    else
                    {
                        hor_ra = 0.0;
                    }
                    hor_dec = Math.Atan2(prz, proj) * RAD2DEG;
                }
            }
            else if (refraction != Refraction.None)
                throw new ArgumentException(string.Format("Unsupported refraction option {0}", refraction));

            return new Topocentric(az, 90.0 - zd, hor_ra, hor_dec);
        }

        /// <summary>
        /// Calculates geocentric ecliptic coordinates for the Sun.
        /// </summary>
        /// <remarks>
        /// This function calculates the position of the Sun as seen from the Earth.
        /// The returned value includes both Cartesian and spherical coordinates.
        /// The x-coordinate and longitude values in the returned structure are based
        /// on the *true equinox of date*: one of two points in the sky where the instantaneous
        /// plane of the Earth's equator at the given date and time (the *equatorial plane*)
        /// intersects with the plane of the Earth's orbit around the Sun (the *ecliptic plane*).
        /// By convention, the apparent location of the Sun at the March equinox is chosen
        /// as the longitude origin and x-axis direction, instead of the one for September.
        ///
        /// `SunPosition` corrects for precession and nutation of the Earth's axis
        /// in order to obtain the exact equatorial plane at the given time.
        ///
        /// This function can be used for calculating changes of seasons: equinoxes and solstices.
        /// In fact, the function #Astronomy.Seasons does use this function for that purpose.
        /// </remarks>
        /// <param name="time">
        /// The date and time for which to calculate the Sun's position.
        /// </param>
        /// <returns>
        /// The ecliptic coordinates of the Sun using the Earth's true equator of date.
        /// </returns>
        public static Ecliptic SunPosition(AstroTime time)
        {
            /* Correct for light travel time from the Sun. */
            /* Otherwise season calculations (equinox, solstice) will all be early by about 8 minutes! */
            AstroTime adjusted_time = time.AddDays(-1.0 / C_AUDAY);

            AstroVector earth2000 = CalcEarth(adjusted_time);

            /* Convert heliocentric location of Earth to geocentric location of Sun. */
            AstroVector sun2000 = new AstroVector(-earth2000.x, -earth2000.y, -earth2000.z, adjusted_time);

            /* Convert to equatorial Cartesian coordinates of date. */
            AstroVector stemp = precession(0.0, sun2000, adjusted_time.tt);
            AstroVector sun_ofdate = nutation(adjusted_time, 0, stemp);

            /* Convert equatorial coordinates to ecliptic coordinates. */
            double true_obliq = DEG2RAD * e_tilt(adjusted_time).tobl;
            return RotateEquatorialToEcliptic(sun_ofdate, true_obliq);
        }

        private static Ecliptic RotateEquatorialToEcliptic(AstroVector pos, double obliq_radians)
        {
            double cos_ob = Math.Cos(obliq_radians);
            double sin_ob = Math.Sin(obliq_radians);

            double ex = +pos.x;
            double ey = +pos.y*cos_ob + pos.z*sin_ob;
            double ez = -pos.y*sin_ob + pos.z*cos_ob;

            double xyproj = Math.Sqrt(ex*ex + ey*ey);
            double elon = 0.0;
            if (xyproj > 0.0)
            {
                elon = RAD2DEG * Math.Atan2(ey, ex);
                if (elon < 0.0)
                    elon += 360.0;
            }

            double elat = RAD2DEG * Math.Atan2(ez, xyproj);

            return new Ecliptic(ex, ey, ez, elat, elon);
        }

        /// <summary>
        /// Converts J2000 equatorial Cartesian coordinates to J2000 ecliptic coordinates.
        /// </summary>
        /// <remarks>
        /// Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
        /// on 1 January 2000), this function converts those coordinates to J2000 ecliptic coordinates,
        /// which are relative to the plane of the Earth's orbit around the Sun.
        /// </remarks>
        /// <param name="equ">
        /// Equatorial coordinates in the J2000 frame of reference.
        /// You can call #Astronomy.GeoVector to obtain suitable equatorial coordinates.
        /// </param>
        /// <returns>Ecliptic coordinates in the J2000 frame of reference.</returns>
        public static Ecliptic EquatorialToEcliptic(AstroVector equ)
        {
            /* Based on NOVAS functions equ2ecl() and equ2ecl_vec(). */
            const double ob2000 = 0.40909260059599012;   /* mean obliquity of the J2000 ecliptic in radians */
            return RotateEquatorialToEcliptic(equ, ob2000);
        }

        /// <summary>
        /// Finds both equinoxes and both solstices for a given calendar year.
        /// </summary>
        /// <remarks>
        /// The changes of seasons are defined by solstices and equinoxes.
        /// Given a calendar year number, this function calculates the
        /// March and September equinoxes and the June and December solstices.
        ///
        /// The equinoxes are the moments twice each year when the plane of the
        /// Earth's equator passes through the center of the Sun. In other words,
        /// the Sun's declination is zero at both equinoxes.
        /// The March equinox defines the beginning of spring in the northern hemisphere
        /// and the beginning of autumn in the southern hemisphere.
        /// The September equinox defines the beginning of autumn in the northern hemisphere
        /// and the beginning of spring in the southern hemisphere.
        ///
        /// The solstices are the moments twice each year when one of the Earth's poles
        /// is most tilted toward the Sun. More precisely, the Sun's declination reaches
        /// its minimum value at the December solstice, which defines the beginning of
        /// winter in the northern hemisphere and the beginning of summer in the southern
        /// hemisphere. The Sun's declination reaches its maximum value at the June solstice,
        /// which defines the beginning of summer in the northern hemisphere and the beginning
        /// of winter in the southern hemisphere.
        /// </remarks>
        /// <param name="year">
        /// The calendar year number for which to calculate equinoxes and solstices.
        /// The value may be any integer, but only the years 1800 through 2100 have been
        /// validated for accuracy: unit testing against data from the
        /// United States Naval Observatory confirms that all equinoxes and solstices
        /// for that range of years are within 2 minutes of the correct time.
        /// </param>
        /// <returns>
        /// A #SeasonsInfo structure that contains four #AstroTime values:
        /// the March and September equinoxes and the June and December solstices.
        /// </returns>
        public static SeasonsInfo Seasons(int year)
        {
            return new SeasonsInfo(
                FindSeasonChange(  0, year,  3, 19),
                FindSeasonChange( 90, year,  6, 19),
                FindSeasonChange(180, year,  9, 21),
                FindSeasonChange(270, year, 12, 20)
            );
        }

        private static AstroTime FindSeasonChange(double targetLon, int year, int month, int day)
        {
            var startTime = new AstroTime(year, month, day, 0, 0, 0);
            return SearchSunLongitude(targetLon, startTime, 4.0);
        }

        /// <summary>
        /// Searches for the time when the Sun reaches an apparent ecliptic longitude as seen from the Earth.
        /// </summary>
        /// <remarks>
        /// This function finds the moment in time, if any exists in the given time window,
        /// that the center of the Sun reaches a specific ecliptic longitude as seen from the center of the Earth.
        ///
        /// This function can be used to determine equinoxes and solstices.
        /// However, it is usually more convenient and efficient to call #Astronomy.Seasons
        /// to calculate all equinoxes and solstices for a given calendar year.
        ///
        /// The function searches the window of time specified by `startTime` and `startTime+limitDays`.
        /// The search will return an error if the Sun never reaches the longitude `targetLon` or
        /// if the window is so large that the longitude ranges more than 180 degrees within it.
        /// It is recommended to keep the window smaller than 10 days when possible.
        /// </remarks>
        /// <param name="targetLon">
        /// The desired ecliptic longitude in degrees, relative to the true equinox of date.
        /// This may be any value in the range [0, 360), although certain values have
        /// conventional meanings:
        /// 0 = March equinox, 90 = June solstice, 180 = September equinox, 270 = December solstice.
        /// </param>
        /// <param name="startTime">
        /// The date and time for starting the search for the desired longitude event.
        /// </param>
        /// <param name="limitDays">
        /// The real-valued number of days, which when added to `startTime`, limits the
        /// range of time over which the search looks.
        /// It is recommended to keep this value between 1 and 10 days.
        /// See remarks above for more details.
        /// </param>
        /// <returns>
        /// The date and time when the Sun reaches the specified apparent ecliptic longitude.
        /// </returns>
        public static AstroTime SearchSunLongitude(double targetLon, AstroTime startTime, double limitDays)
        {
            var sun_offset = new SearchContext_SunOffset(targetLon);
            AstroTime t2 = startTime.AddDays(limitDays);
            return Search(sun_offset, startTime, t2, 1.0);
        }

        /// <summary>
        /// Searches for a time at which a function's value increases through zero.
        /// </summary>
        /// <remarks>
        /// Certain astronomy calculations involve finding a time when an event occurs.
        /// Often such events can be defined as the root of a function:
        /// the time at which the function's value becomes zero.
        ///
        /// `Search` finds the *ascending root* of a function: the time at which
        /// the function's value becomes zero while having a positive slope. That is, as time increases,
        /// the function transitions from a negative value, through zero at a specific moment,
        /// to a positive value later. The goal of the search is to find that specific moment.
        ///
        /// The `func` parameter is an instance of the abstract class #SearchContext.
        /// As an example, a caller may wish to find the moment a celestial body reaches a certain
        /// ecliptic longitude. In that case, the caller might derive a class that contains
        /// a #Body member to specify the body and a `double` to hold the target longitude.
        /// It could subtract the target longitude from the actual longitude at a given time;
        /// thus the difference would equal zero at the moment in time the planet reaches the
        /// desired longitude.
        ///
        /// Every call to `func.Eval` must either return a valid #AstroTime or throw an exception.
        ///
        /// The search calls `func.Eval` repeatedly to rapidly narrow in on any ascending
        /// root within the time window specified by `t1` and `t2`. The search never
        /// reports a solution outside this time window.
        ///
        /// `Search` uses a combination of bisection and quadratic interpolation
        /// to minimize the number of function calls. However, it is critical that the
        /// supplied time window be small enough that there cannot be more than one root
        /// (ascedning or descending) within it; otherwise the search can fail.
        /// Beyond that, it helps to make the time window as small as possible, ideally
        /// such that the function itself resembles a smooth parabolic curve within that window.
        ///
        /// If an ascending root is not found, or more than one root
        /// (ascending and/or descending) exists within the window `t1`..`t2`,
        /// the search will return `null`.
        ///
        /// If the search does not converge within 20 iterations, it will throw an exception.
        /// </remarks>
        /// <param name="func">
        /// The function for which to find the time of an ascending root.
        /// See remarks above for more details.
        /// </param>
        /// <param name="t1">
        /// The lower time bound of the search window.
        /// See remarks above for more details.
        /// </param>
        /// <param name="t2">
        /// The upper time bound of the search window.
        /// See remarks above for more details.
        /// </param>
        /// <param name="dt_tolerance_seconds">
        /// Specifies an amount of time in seconds within which a bounded ascending root
        /// is considered accurate enough to stop. A typical value is 1 second.
        /// </param>
        /// <returns>
        /// If successful, returns an #AstroTime value indicating a date and time
        /// that is within `dt_tolerance_seconds` of an ascending root.
        /// If no ascending root is found, or more than one root exists in the time
        /// window `t1`..`t2`, the function returns `null`.
        /// If the search does not converge within 20 iterations, an exception is thrown.
        /// </returns>
        public static AstroTime Search(
            SearchContext func,
            AstroTime t1,
            AstroTime t2,
            double dt_tolerance_seconds)
        {
            const int iter_limit = 20;
            double dt_days = Math.Abs(dt_tolerance_seconds / SECONDS_PER_DAY);
            double f1 = func.Eval(t1);
            double f2 = func.Eval(t2);
            int iter = 0;
            bool calc_fmid = true;
            double fmid = 0.0;
            for(;;)
            {
                if (++iter > iter_limit)
                    throw new Exception(string.Format("Search did not converge within {0} iterations.", iter_limit));

                double dt = (t2.tt - t1.tt) / 2.0;
                AstroTime tmid = t1.AddDays(dt);
                if (Math.Abs(dt) < dt_days)
                {
                    /* We are close enough to the event to stop the search. */
                    return tmid;
                }

                if (calc_fmid)
                    fmid = func.Eval(tmid);
                else
                    calc_fmid = true;   /* we already have the correct value of fmid from the previous loop */

                /* Quadratic interpolation: */
                /* Try to find a parabola that passes through the 3 points we have sampled: */
                /* (t1,f1), (tmid,fmid), (t2,f2) */

                double q_x, q_ut, q_df_dt;
                if (QuadInterp(tmid.ut, t2.ut - tmid.ut, f1, fmid, f2, out q_x, out q_ut, out q_df_dt))
                {
                    var tq = new AstroTime(q_ut);
                    double fq = func.Eval(tq);
                    if (q_df_dt != 0.0)
                    {
                        double dt_guess = Math.Abs(fq / q_df_dt);
                        if (dt_guess < dt_days)
                        {
                            /* The estimated time error is small enough that we can quit now. */
                            return tq;
                        }

                        /* Try guessing a tighter boundary with the interpolated root at the center. */
                        dt_guess *= 1.2;
                        if (dt_guess < dt/10.0)
                        {
                            AstroTime tleft = tq.AddDays(-dt_guess);
                            AstroTime tright = tq.AddDays(+dt_guess);
                            if ((tleft.ut - t1.ut)*(tleft.ut - t2.ut) < 0.0)
                            {
                                if ((tright.ut - t1.ut)*(tright.ut - t2.ut) < 0.0)
                                {
                                    double fleft, fright;
                                    fleft = func.Eval(tleft);
                                    fright = func.Eval(tright);
                                    if (fleft<0.0 && fright>=0.0)
                                    {
                                        f1 = fleft;
                                        f2 = fright;
                                        t1 = tleft;
                                        t2 = tright;
                                        fmid = fq;
                                        calc_fmid = false;  /* save a little work -- no need to re-calculate fmid next time around the loop */
                                        continue;
                                    }
                                }
                            }
                        }
                    }
                }

                /* After quadratic interpolation attempt. */
                /* Now just divide the region in two parts and pick whichever one appears to contain a root. */
                if (f1 < 0.0 && fmid >= 0.0)
                {
                    t2 = tmid;
                    f2 = fmid;
                    continue;
                }

                if (fmid < 0.0 && f2 >= 0.0)
                {
                    t1 = tmid;
                    f1 = fmid;
                    continue;
                }

                /* Either there is no ascending zero-crossing in this range */
                /* or the search window is too wide (more than one zero-crossing). */
                return null;
            }
        }

        private static bool QuadInterp(
            double tm, double dt, double fa, double fm, double fb,
            out double out_x, out double out_t, out double out_df_dt)
        {
            double Q, R, S;
            double u, ru, x1, x2;

            out_x = out_t = out_df_dt = 0.0;

            Q = (fb + fa)/2.0 - fm;
            R = (fb - fa)/2.0;
            S = fm;

            if (Q == 0.0)
            {
                /* This is a line, not a parabola. */
                if (R == 0.0)
                    return false;       /* This is a HORIZONTAL line... can't make progress! */
                out_x = -S / R;
                if (out_x < -1.0 || out_x > +1.0)
                    return false;   /* out of bounds */
            }
            else
            {
                /* This really is a parabola. Find roots x1, x2. */
                u = R*R - 4*Q*S;
                if (u <= 0.0)
                    return false;   /* can't solve if imaginary, or if vertex of parabola is tangent. */

                ru = Math.Sqrt(u);
                x1 = (-R + ru) / (2.0 * Q);
                x2 = (-R - ru) / (2.0 * Q);
                if (-1.0 <= x1 && x1 <= +1.0)
                {
                    if (-1.0 <= x2 && x2 <= +1.0)
                        return false;   /* two roots are within bounds; we require a unique zero-crossing. */
                    out_x = x1;
                }
                else if (-1.0 <= x2 && x2 <= +1.0)
                    out_x = x2;
                else
                    return false;   /* neither root is within bounds */
            }

            out_t = tm + out_x*dt;
            out_df_dt = (2*Q*out_x + R) / dt;
            return true;   /* success */
        }

        ///
        /// <summary>
        /// Returns a body's ecliptic longitude with respect to the Sun, as seen from the Earth.
        /// </summary>
        /// <remarks>
        /// This function can be used to determine where a planet appears around the ecliptic plane
        /// (the plane of the Earth's orbit around the Sun) as seen from the Earth,
        /// relative to the Sun's apparent position.
        ///
        /// The angle starts at 0 when the body and the Sun are at the same ecliptic longitude
        /// as seen from the Earth. The angle increases in the prograde direction
        /// (the direction that the planets orbit the Sun and the Moon orbits the Earth).
        ///
        /// When the angle is 180 degrees, it means the Sun and the body appear on opposite sides
        /// of the sky for an Earthly observer. When `body` is a planet whose orbit around the
        /// Sun is farther than the Earth's, 180 degrees indicates opposition. For the Moon,
        /// it indicates a full moon.
        ///
        /// The angle keeps increasing up to 360 degrees as the body's apparent prograde
        /// motion continues relative to the Sun. When the angle reaches 360 degrees, it starts
        /// over at 0 degrees.
        ///
        /// Values between 0 and 180 degrees indicate that the body is visible in the evening sky
        /// after sunset.  Values between 180 degrees and 360 degrees indicate that the body
        /// is visible in the morning sky before sunrise.
        /// </remarks>
        /// <param name="body">The celestial body for which to find longitude from the Sun.</param>
        /// <param name="time">The date and time of the observation.</param>
        /// <returns>
        /// A value in the range [0, 360), expressed in degrees.
        /// </returns>
        public static double LongitudeFromSun(Body body, AstroTime time)
        {
            if (body == Body.Earth)
                throw new EarthNotAllowedException();

            AstroVector sv = GeoVector(Body.Sun, time, Aberration.None);
            Ecliptic se = EquatorialToEcliptic(sv);

            AstroVector bv = GeoVector(body, time, Aberration.None);
            Ecliptic be = EquatorialToEcliptic(bv);

            return NormalizeLongitude(be.elon - se.elon);
        }

        /// <summary>
        /// Returns the Moon's phase as an angle from 0 to 360 degrees.
        /// </summary>
        /// <remarks>
        /// This function determines the phase of the Moon using its apparent
        /// ecliptic longitude relative to the Sun, as seen from the center of the Earth.
        /// Certain values of the angle have conventional definitions:
        ///
        /// - 0 = new moon
        /// - 90 = first quarter
        /// - 180 = full moon
        /// - 270 = third quarter
        /// </remarks>
        /// <param name="time">The date and time of the observation.</param>
        /// <returns>The angle as described above, a value in the range 0..360 degrees.</returns>
        public static double MoonPhase(AstroTime time)
        {
            return LongitudeFromSun(Body.Moon, time);
        }

        /// <summary>
        /// Finds the first lunar quarter after the specified date and time.
        /// </summary>
        /// <remarks>
        /// A lunar quarter is one of the following four lunar phase events:
        /// new moon, first quarter, full moon, third quarter.
        /// This function finds the lunar quarter that happens soonest
        /// after the specified date and time.
        ///
        /// To continue iterating through consecutive lunar quarters, call this function once,
        /// followed by calls to #Astronomy.NextMoonQuarter as many times as desired.
        /// </remarks>
        /// <param name="startTime">The date and time at which to start the search.</param>
        /// <returns>
        /// A #MoonQuarterInfo structure reporting the next quarter phase and the time it will occur.
        /// </returns>
        public static MoonQuarterInfo SearchMoonQuarter(AstroTime startTime)
        {
            double angres = MoonPhase(startTime);
            int quarter = (1 + (int)Math.Floor(angres / 90.0)) % 4;
            AstroTime qtime = SearchMoonPhase(90.0 * quarter, startTime, 10.0);
            return new MoonQuarterInfo(quarter, qtime);
        }

        /// <summary>
        /// Continues searching for lunar quarters from a previous search.
        /// </summary>
        /// <remarks>
        /// After calling #Astronomy.SearchMoonQuarter, this function can be called
        /// one or more times to continue finding consecutive lunar quarters.
        /// This function finds the next consecutive moon quarter event after
        /// the one passed in as the parameter `mq`.
        /// </remarks>
        /// <param name="mq">The previous moon quarter found by a call to #Astronomy.SearchMoonQuarter or `Astronomy.NextMoonQuarter`.</param>
        /// <returns>The moon quarter that occurs next in time after the one passed in `mq`.</returns>
        public static MoonQuarterInfo NextMoonQuarter(MoonQuarterInfo mq)
        {
            /* Skip 6 days past the previous found moon quarter to find the next one. */
            /* This is less than the minimum possible increment. */
            /* So far I have seen the interval well contained by the range (6.5, 8.3) days. */

            AstroTime time = mq.time.AddDays(6.0);
            MoonQuarterInfo next_mq = SearchMoonQuarter(time);
            /* Verify that we found the expected moon quarter. */
            if (next_mq.quarter != (1 + mq.quarter) % 4)
                throw new Exception("Internal error: found the wrong moon quarter.");
            return next_mq;
        }

        ///
        /// <summary>Searches for the time that the Moon reaches a specified phase.</summary>
        /// <remarks>
        /// Lunar phases are conventionally defined in terms of the Moon's geocentric ecliptic
        /// longitude with respect to the Sun's geocentric ecliptic longitude.
        /// When the Moon and the Sun have the same longitude, that is defined as a new moon.
        /// When their longitudes are 180 degrees apart, that is defined as a full moon.
        ///
        /// This function searches for any value of the lunar phase expressed as an
        /// angle in degrees in the range [0, 360).
        ///
        /// If you want to iterate through lunar quarters (new moon, first quarter, full moon, third quarter)
        /// it is much easier to call the functions #Astronomy.SearchMoonQuarter and #Astronomy.NextMoonQuarter.
        /// This function is useful for finding general phase angles outside those four quarters.
        /// </remarks>
        /// <param name="targetLon">
        /// The difference in geocentric longitude between the Sun and Moon
        /// that specifies the lunar phase being sought. This can be any value
        /// in the range [0, 360).  Certain values have conventional names:
        /// 0 = new moon, 90 = first quarter, 180 = full moon, 270 = third quarter.
        /// </param>
        /// <param name="startTime">
        /// The beginning of the time window in which to search for the Moon reaching the specified phase.
        /// </param>
        /// <param name="limitDays">
        /// The number of days after `startTime` that limits the time window for the search.
        /// </param>
        /// <returns>
        /// If successful, returns the date and time the moon reaches the phase specified by
        /// `targetlon`. This function will return throw an exception if the phase does not
        /// occur within `limitDays` of `startTime`; that is, if the search window is too small.
        /// </returns>
        public static AstroTime SearchMoonPhase(double targetLon, AstroTime startTime, double limitDays)
        {
            /*
                To avoid discontinuities in the moon_offset function causing problems,
                we need to approximate when that function will next return 0.
                We probe it with the start time and take advantage of the fact
                that every lunar phase repeats roughly every 29.5 days.
                There is a surprising uncertainty in the quarter timing,
                due to the eccentricity of the moon's orbit.
                I have seen up to 0.826 days away from the simple prediction.
                To be safe, we take the predicted time of the event and search
                +/-0.9 days around it (a 1.8-day wide window).
                Return null if the final result goes beyond limitDays after startTime.
            */

            const double uncertainty = 0.9;
            var moon_offset = new SearchContext_MoonOffset(targetLon);

            double ya = moon_offset.Eval(startTime);
            if (ya > 0.0) ya -= 360.0;  /* force searching forward in time, not backward */
            double est_dt = -(MEAN_SYNODIC_MONTH * ya) / 360.0;
            double dt1 = est_dt - uncertainty;
            if (dt1 > limitDays)
                return null;    /* not possible for moon phase to occur within specified window (too short) */
            double dt2 = est_dt + uncertainty;
            if (limitDays < dt2)
                dt2 = limitDays;
            AstroTime t1 = startTime.AddDays(dt1);
            AstroTime t2 = startTime.AddDays(dt2);
            AstroTime time = Search(moon_offset, t1, t2, 1.0);
            if (time == null)
                throw new Exception(string.Format("Could not find moon longitude {0} within {1} days of {2}", targetLon, limitDays, startTime));
            return time;
        }

        /// <summary>
        /// Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.
        /// </summary>
        /// <remarks>
        /// This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
        /// Rise time is when the body first starts to be visible above the horizon.
        /// For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
        /// Set time is the moment when the body appears to vanish below the horizon.
        ///
        /// This function corrects for typical atmospheric refraction, which causes celestial
        /// bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
        /// It also adjusts for the apparent angular radius of the observed body (significant only for the Sun and Moon).
        ///
        /// Note that rise or set may not occur in every 24 hour period.
        /// For example, near the Earth's poles, there are long periods of time where
        /// the Sun stays below the horizon, never rising.
        /// Also, it is possible for the Moon to rise just before midnight but not set during the subsequent 24-hour day.
        /// This is because the Moon sets nearly an hour later each day due to orbiting the Earth a
        /// significant amount during each rotation of the Earth.
        /// Therefore callers must not assume that the function will always succeed.
        /// </remarks>
        ///
        /// <param name="body">The Sun, Moon, or any planet other than the Earth.</param>
        ///
        /// <param name="observer">The location where observation takes place.</param>
        ///
        /// <param name="direction">
        ///      Either `Direction.Rise` to find a rise time or `Direction.Set` to find a set time.
        /// </param>
        ///
        /// <param name="startTime">The date and time at which to start the search.</param>
        ///
        /// <param name="limitDays">
        /// Limits how many days to search for a rise or set time.
        /// To limit a rise or set time to the same day, you can use a value of 1 day.
        /// In cases where you want to find the next rise or set time no matter how far
        /// in the future (for example, for an observer near the south pole), you can
        /// pass in a larger value like 365.
        /// </param>
        ///
        /// <returns>
        /// On success, returns the date and time of the rise or set time as requested.
        /// If the function returns `null`, it means the rise or set event does not occur
        /// within `limitDays` days of `startTime`. This is a normal condition,
        /// not an error.
        /// </returns>
        public static AstroTime SearchRiseSet(
            Body body,
            Observer observer,
            Direction direction,
            AstroTime startTime,
            double limitDays)
        {
            if (body == Body.Earth)
                throw new EarthNotAllowedException();

            double ha_before, ha_after;
            switch (direction)
            {
                case Direction.Rise:
                    ha_before = 12.0;   /* minimum altitude (bottom) happens BEFORE the body rises. */
                    ha_after = 0.0;     /* maximum altitude (culmination) happens AFTER the body rises. */
                    break;

                case Direction.Set:
                    ha_before = 0.0;    /* culmination happens BEFORE the body sets. */
                    ha_after = 12.0;    /* bottom happens AFTER the body sets. */
                    break;

                default:
                    throw new ArgumentException(string.Format("Unsupported direction value {0}", direction));
            }

            var peak_altitude = new SearchContext_PeakAltitude(body, direction, observer);
            /*
                See if the body is currently above/below the horizon.
                If we are looking for next rise time and the body is below the horizon,
                we use the current time as the lower time bound and the next culmination
                as the upper bound.
                If the body is above the horizon, we search for the next bottom and use it
                as the lower bound and the next culmination after that bottom as the upper bound.
                The same logic applies for finding set times, only we swap the hour angles.
            */

            HourAngleInfo evt_before, evt_after;
            AstroTime time_start = startTime;
            double alt_before = peak_altitude.Eval(time_start);
            AstroTime time_before;
            if (alt_before > 0.0)
            {
                /* We are past the sought event, so we have to wait for the next "before" event (culm/bottom). */
                evt_before = SearchHourAngle(body, observer, ha_before, time_start);
                time_before = evt_before.time;
                alt_before = peak_altitude.Eval(time_before);
            }
            else
            {
                /* We are before or at the sought event, so we find the next "after" event (bottom/culm), */
                /* and use the current time as the "before" event. */
                time_before = time_start;
            }

            evt_after = SearchHourAngle(body, observer, ha_after, time_before);
            double alt_after = peak_altitude.Eval(evt_after.time);

            for(;;)
            {
                if (alt_before <= 0.0 && alt_after > 0.0)
                {
                    /* Search between evt_before and evt_after for the desired event. */
                    AstroTime result = Search(peak_altitude, time_before, evt_after.time, 1.0);
                    if (result != null)
                        return result;
                }

                /* If we didn't find the desired event, use evt_after.time to find the next before-event. */
                evt_before = SearchHourAngle(body, observer, ha_before, evt_after.time);
                evt_after = SearchHourAngle(body, observer, ha_after, evt_before.time);

                if (evt_before.time.ut >= time_start.ut + limitDays)
                    return null;

                time_before = evt_before.time;

                alt_before = peak_altitude.Eval(evt_before.time);
                alt_after = peak_altitude.Eval(evt_after.time);
            }
        }

        /// <summary>
        /// Searches for the time when a celestial body reaches a specified hour angle as seen by an observer on the Earth.
        /// </summary>
        ///
        /// <remarks>
        /// The *hour angle* of a celestial body indicates its position in the sky with respect
        /// to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
        /// The hour angle is 0 when the body reaches its highest angle above the horizon in a given day.
        /// The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
        /// to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
        /// the number of hours that have passed since the most recent time that the body has culminated,
        /// or reached its highest point.
        ///
        /// This function searches for the next time a celestial body reaches the given hour angle
        /// after the date and time specified by `startTime`.
        /// To find when a body culminates, pass 0 for `hourAngle`.
        /// To find when a body reaches its lowest point in the sky, pass 12 for `hourAngle`.
        ///
        /// Note that, especially close to the Earth's poles, a body as seen on a given day
        /// may always be above the horizon or always below the horizon, so the caller cannot
        /// assume that a culminating object is visible nor that an object is below the horizon
        /// at its minimum altitude.
        ///
        /// On success, the function reports the date and time, along with the horizontal coordinates
        /// of the body at that time, as seen by the given observer.
        /// </remarks>
        ///
        /// <param name="body">
        /// The celestial body, which can the Sun, the Moon, or any planet other than the Earth.
        /// </param>
        ///
        /// <param name="observer">
        /// Indicates a location on or near the surface of the Earth where the observer is located.
        /// </param>
        ///
        /// <param name="hourAngle">
        /// An hour angle value in the range [0, 24) indicating the number of sidereal hours after the
        /// body's most recent culmination.
        /// </param>
        ///
        /// <param name="startTime">
        /// The date and time at which to start the search.
        /// </param>
        ///
        /// <returns>
        /// This function returns a valid #HourAngleInfo object on success.
        /// If any error occurs, it throws an exception.
        /// It never returns a null value.
        /// </returns>
        public static HourAngleInfo SearchHourAngle(
            Body body,
            Observer observer,
            double hourAngle,
            AstroTime startTime)
        {
            int iter = 0;

            if (body == Body.Earth)
                throw new EarthNotAllowedException();

            if (hourAngle < 0.0 || hourAngle >= 24.0)
                throw new ArgumentException("hourAngle is out of the allowed range [0, 24).");

            AstroTime time = startTime;
            for(;;)
            {
                ++iter;

                /* Calculate Greenwich Apparent Sidereal Time (GAST) at the given time. */
                double gast = sidereal_time(time);

                /* Obtain equatorial coordinates of date for the body. */
                Equatorial ofdate = Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);

                /* Calculate the adjustment needed in sidereal time */
                /* to bring the hour angle to the desired value. */

                double delta_sidereal_hours = ((hourAngle + ofdate.ra - observer.longitude/15.0) - gast) % 24.0;
                if (iter == 1)
                {
                    /* On the first iteration, always search forward in time. */
                    if (delta_sidereal_hours < 0.0)
                        delta_sidereal_hours += 24.0;
                }
                else
                {
                    /* On subsequent iterations, we make the smallest possible adjustment, */
                    /* either forward or backward in time. */
                    if (delta_sidereal_hours < -12.0)
                        delta_sidereal_hours += 24.0;
                    else if (delta_sidereal_hours > +12.0)
                        delta_sidereal_hours -= 24.0;
                }

                /* If the error is tolerable (less than 0.1 seconds), the search has succeeded. */
                if (Math.Abs(delta_sidereal_hours) * 3600.0 < 0.1)
                {
                    Topocentric hor = Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.Normal);
                    return new HourAngleInfo(time, hor);
                }

                /* We need to loop another time to get more accuracy. */
                /* Update the terrestrial time (in solar days) adjusting by sidereal time (sidereal hours). */
                time = time.AddDays((delta_sidereal_hours / 24.0) * SOLAR_DAYS_PER_SIDEREAL_DAY);
            }
        }

        /// <summary>
        ///      Searches for the time when the Earth and another planet are separated by a specified angle
        ///      in ecliptic longitude, as seen from the Sun.
        /// </summary>
        ///
        /// <remarks>
        /// A relative longitude is the angle between two bodies measured in the plane of the Earth's orbit
        /// (the ecliptic plane). The distance of the bodies above or below the ecliptic plane is ignored.
        /// If you imagine the shadow of the body cast onto the ecliptic plane, and the angle measured around
        /// that plane from one body to the other in the direction the planets orbit the Sun, you will get an
        /// angle somewhere between 0 and 360 degrees. This is the relative longitude.
        ///
        /// Given a planet other than the Earth in `body` and a time to start the search in `startTime`,
        /// this function searches for the next time that the relative longitude measured from the planet
        /// to the Earth is `targetRelLon`.
        ///
        /// Certain astronomical events are defined in terms of relative longitude between the Earth and another planet:
        ///
        /// - When the relative longitude is 0 degrees, it means both planets are in the same direction from the Sun.
        ///   For planets that orbit closer to the Sun (Mercury and Venus), this is known as *inferior conjunction*,
        ///   a time when the other planet becomes very difficult to see because of being lost in the Sun's glare.
        ///   (The only exception is in the rare event of a transit, when we see the silhouette of the planet passing
        ///   between the Earth and the Sun.)
        ///
        /// - When the relative longitude is 0 degrees and the other planet orbits farther from the Sun,
        ///   this is known as *opposition*.  Opposition is when the planet is closest to the Earth, and
        ///   also when it is visible for most of the night, so it is considered the best time to observe the planet.
        ///
        /// - When the relative longitude is 180 degrees, it means the other planet is on the opposite side of the Sun
        ///   from the Earth. This is called *superior conjunction*. Like inferior conjunction, the planet is
        ///   very difficult to see from the Earth. Superior conjunction is possible for any planet other than the Earth.
        /// </remarks>
        ///
        /// <param name="body">
        ///      A planet other than the Earth.
        ///      If `body` is `Body.Earth`, `Body.Sun`, or `Body.Moon`, this function throws an exception.
        /// </param>
        ///
        /// <param name="targetRelLon">
        ///      The desired relative longitude, expressed in degrees. Must be in the range [0, 360).
        /// </param>
        ///
        /// <param name="startTime">
        ///      The date and time at which to begin the search.
        /// </param>
        ///
        /// <returns>
        ///      If successful, returns the date and time of the relative longitude event.
        ///      Otherwise this function returns null.
        /// </returns>
        public static AstroTime SearchRelativeLongitude(Body body, double targetRelLon, AstroTime startTime)
        {
            if (body == Body.Earth || body == Body.Sun || body == Body.Moon)
                throw new InvalidBodyException(body);

            double syn = SynodicPeriod(body);
            int direction = IsSuperiorPlanet(body) ? +1 : -1;

            /* Iterate until we converge on the desired event. */
            /* Calculate the error angle, which will be a negative number of degrees, */
            /* meaning we are "behind" the target relative longitude. */

            double error_angle = rlon_offset(body, startTime, direction, targetRelLon);
            if (error_angle > 0.0)
                error_angle -= 360.0;    /* force searching forward in time */

            AstroTime time = startTime;
            for (int iter = 0; iter < 100; ++iter)
            {
                /* Estimate how many days in the future (positive) or past (negative) */
                /* we have to go to get closer to the target relative longitude. */
                double day_adjust = (-error_angle/360.0) * syn;
                time = time.AddDays(day_adjust);
                if (Math.Abs(day_adjust) * SECONDS_PER_DAY < 1.0)
                    return time;

                double prev_angle = error_angle;
                error_angle = rlon_offset(body, time, direction, targetRelLon);
                if (Math.Abs(prev_angle) < 30.0 && (prev_angle != error_angle))
                {
                    /* Improve convergence for Mercury/Mars (eccentric orbits) */
                    /* by adjusting the synodic period to more closely match the */
                    /* variable speed of both planets in this part of their respective orbits. */
                    double ratio = prev_angle / (prev_angle - error_angle);
                    if (ratio > 0.5 && ratio < 2.0)
                        syn *= ratio;
                }
            }

            throw new Exception("Relative longitude search failed to converge.");
        }

        private static double rlon_offset(Body body, AstroTime time, int direction, double targetRelLon)
        {
            double plon = EclipticLongitude(body, time);
            double elon = EclipticLongitude(Body.Earth, time);
            double diff = direction * (elon - plon);
            return LongitudeOffset(diff - targetRelLon);
        }

        private static double SynodicPeriod(Body body)
        {
            /* The Earth does not have a synodic period as seen from itself. */
            if (body == Body.Earth)
                throw new EarthNotAllowedException();

            if (body == Body.Moon)
                return MEAN_SYNODIC_MONTH;

            double Tp = PlanetOrbitalPeriod(body);
            return Math.Abs(EARTH_ORBITAL_PERIOD / (EARTH_ORBITAL_PERIOD/Tp - 1.0));
        }

        /// <summary>Calculates heliocentric ecliptic longitude of a body based on the J2000 equinox.</summary>
        /// <remarks>
        /// This function calculates the angle around the plane of the Earth's orbit
        /// of a celestial body, as seen from the center of the Sun.
        /// The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
        /// in degrees from the J2000 equinox. The ecliptic longitude is always in the range [0, 360).
        /// </remarks>
        ///
        /// <param name="body">A body other than the Sun.</param>
        ///
        /// <param name="time">The date and time at which the body's ecliptic longitude is to be calculated.</param>
        ///
        /// <returns>
        ///      Returns the ecliptic longitude in degrees of the given body at the given time.
        /// </returns>
        public static double EclipticLongitude(Body body, AstroTime time)
        {
            if (body == Body.Sun)
                throw new ArgumentException("Cannot calculate heliocentric longitude of the Sun.");

            AstroVector hv = HelioVector(body, time);
            Ecliptic eclip = EquatorialToEcliptic(hv);
            return eclip.elon;
        }

        private static double PlanetOrbitalPeriod(Body body)
        {
            /* Returns the number of days it takes for a planet to orbit the Sun. */
            switch (body)
            {
                case Body.Mercury:  return     87.969;
                case Body.Venus:    return    224.701;
                case Body.Earth:    return    EARTH_ORBITAL_PERIOD;
                case Body.Mars:     return    686.980;
                case Body.Jupiter:  return   4332.589;
                case Body.Saturn:   return  10759.22;
                case Body.Uranus:   return  30685.4;
                case Body.Neptune:  return  NEPTUNE_ORBITAL_PERIOD;
                case Body.Pluto:    return  90560.0;
                default:
                    throw new InvalidBodyException(body);
            }
        }

        private static bool IsSuperiorPlanet(Body body)
        {
            switch (body)
            {
                case Body.Mars:
                case Body.Jupiter:
                case Body.Saturn:
                case Body.Uranus:
                case Body.Neptune:
                case Body.Pluto:
                    return true;

                default:
                    return false;
            }
        }

        /// <summary>
        /// Determines visibility of a celestial body relative to the Sun, as seen from the Earth.
        /// </summary>
        ///
        /// <remarks>
        /// This function returns an #ElongationInfo structure, which provides the following
        /// information about the given celestial body at the given time:
        ///
        /// - `visibility` is an enumerated type that specifies whether the body is more easily seen
        ///    in the morning before sunrise, or in the evening after sunset.
        ///
        /// - `elongation` is the angle in degrees between two vectors: one from the center of the Earth to the
        ///    center of the Sun, the other from the center of the Earth to the center of the specified body.
        ///    This angle indicates how far away the body is from the glare of the Sun.
        ///    The elongation angle is always in the range [0, 180].
        ///
        /// - `ecliptic_separation` is the absolute value of the difference between the body's ecliptic longitude
        ///   and the Sun's ecliptic longitude, both as seen from the center of the Earth. This angle measures
        ///   around the plane of the Earth's orbit, and ignores how far above or below that plane the body is.
        ///   The ecliptic separation is measured in degrees and is always in the range [0, 180].
        /// </remarks>
        ///
        /// <param name="body">
        ///      The celestial body whose visibility is to be calculated.
        /// </param>
        ///
        /// <param name="time">
        ///      The date and time of the observation.
        /// </param>
        ///
        /// <returns>
        /// Returns a valid #ElongationInfo structure, or throws an exception if there is an error.
        /// </returns>
        public static ElongationInfo Elongation(Body body, AstroTime time)
        {
            Visibility visibility;
            double ecliptic_separation = LongitudeFromSun(body, time);
            if (ecliptic_separation > 180.0)
            {
                visibility = Visibility.Morning;
                ecliptic_separation = 360.0 - ecliptic_separation;
            }
            else
            {
                visibility = Visibility.Evening;
            }

            double elongation = AngleFromSun(body, time);
            return new ElongationInfo(time, visibility, elongation, ecliptic_separation);
        }

        /// <summary>
        /// Finds a date and time when Mercury or Venus reaches its maximum angle from the Sun as seen from the Earth.
        /// </summary>
        ///
        /// <remarks>
        /// Mercury and Venus are are often difficult to observe because they are closer to the Sun than the Earth is.
        /// Mercury especially is almost always impossible to see because it gets lost in the Sun's glare.
        /// The best opportunities for spotting Mercury, and the best opportunities for viewing Venus through
        /// a telescope without atmospheric interference, are when these planets reach maximum elongation.
        /// These are events where the planets reach the maximum angle from the Sun as seen from the Earth.
        ///
        /// This function solves for those times, reporting the next maximum elongation event's date and time,
        /// the elongation value itself, the relative longitude with the Sun, and whether the planet is best
        /// observed in the morning or evening. See #Astronomy.Elongation for more details about the returned structure.
        /// </remarks>
        ///
        /// <param name="body">
        /// Either `Body.Mercury` or `Body.Venus`. Any other value will result in an exception.
        /// To find the best viewing opportunites for planets farther from the Sun than the Earth is (Mars through Pluto)
        /// use #Astronomy.SearchRelativeLongitude to find the next opposition event.
        /// </param>
        ///
        /// <param name="startTime">
        /// The date and time at which to begin the search. The maximum elongation event found will always
        /// be the first one that occurs after this date and time.
        /// </param>
        ///
        /// <returns>
        /// Either an exception will be thrown, or the function will return a valid value.
        /// </returns>
        public static ElongationInfo SearchMaxElongation(Body body, AstroTime startTime)
        {
            double s1, s2;
            switch (body)
            {
                case Body.Mercury:
                    s1 = 50.0;
                    s2 = 85.0;
                    break;

                case Body.Venus:
                    s1 = 40.0;
                    s2 = 50.0;
                    break;

                default:
                    throw new InvalidBodyException(body);
            }

            double syn = SynodicPeriod(body);
            var neg_elong_slope = new SearchContext_NegElongSlope(body);

            for (int iter=0; ++iter <= 2;)
            {
                double plon = EclipticLongitude(body, startTime);
                double elon = EclipticLongitude(Body.Earth, startTime);
                double rlon = LongitudeOffset(plon - elon);     /* clamp to (-180, +180] */

                /* The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees */
                /* because there is a cusp there that causes a discontinuity in the derivative. */
                /* So we need to guard against searching near such times. */
                double adjust_days, rlon_lo, rlon_hi;
                if (rlon >= -s1 && rlon < +s1)
                {
                    /* Seek to the window [+s1, +s2]. */
                    adjust_days = 0.0;
                    /* Search forward for the time t1 when rel lon = +s1. */
                    rlon_lo = +s1;
                    /* Search forward for the time t2 when rel lon = +s2. */
                    rlon_hi = +s2;
                }
                else if (rlon > +s2 || rlon < -s2)
                {
                    /* Seek to the next search window at [-s2, -s1]. */
                    adjust_days = 0.0;
                    /* Search forward for the time t1 when rel lon = -s2. */
                    rlon_lo = -s2;
                    /* Search forward for the time t2 when rel lon = -s1. */
                    rlon_hi = -s1;
                }
                else if (rlon >= 0.0)
                {
                    /* rlon must be in the middle of the window [+s1, +s2]. */
                    /* Search BACKWARD for the time t1 when rel lon = +s1. */
                    adjust_days = -syn / 4.0;
                    rlon_lo = +s1;
                    rlon_hi = +s2;
                    /* Search forward from t1 to find t2 such that rel lon = +s2. */
                }
                else
                {
                    /* rlon must be in the middle of the window [-s2, -s1]. */
                    /* Search BACKWARD for the time t1 when rel lon = -s2. */
                    adjust_days = -syn / 4.0;
                    rlon_lo = -s2;
                    /* Search forward from t1 to find t2 such that rel lon = -s1. */
                    rlon_hi = -s1;
                }

                AstroTime t_start = startTime.AddDays(adjust_days);

                AstroTime t1 = SearchRelativeLongitude(body, rlon_lo, t_start);
                AstroTime t2 = SearchRelativeLongitude(body, rlon_hi, t1);

                /* Now we have a time range [t1,t2] that brackets a maximum elongation event. */
                /* Confirm the bracketing. */
                double m1 = neg_elong_slope.Eval(t1);
                if (m1 >= 0.0)
                    throw new Exception("There is a bug in the bracketing algorithm! m1 = " + m1);

                double m2 = neg_elong_slope.Eval(t2);
                if (m2 <= 0.0)
                    throw new Exception("There is a bug in the bracketing algorithm! m2 = " + m2);

                /* Use the generic search algorithm to home in on where the slope crosses from negative to positive. */
                AstroTime searchx = Search(neg_elong_slope, t1, t2, 10.0);
                if (searchx == null)
                    throw new Exception("Maximum elongation search failed.");

                if (searchx.tt >= startTime.tt)
                    return Elongation(body, searchx);

                /* This event is in the past (earlier than startTime). */
                /* We need to search forward from t2 to find the next possible window. */
                /* We never need to search more than twice. */
                startTime = t2.AddDays(1.0);
            }

            throw new Exception("Maximum elongation search iterated too many times.");
        }

        ///
        /// <summary>Returns the angle between the given body and the Sun, as seen from the Earth.</summary>
        ///
        /// <remarks>
        /// This function calculates the angular separation between the given body and the Sun,
        /// as seen from the center of the Earth. This angle is helpful for determining how
        /// easy it is to see the body away from the glare of the Sun.
        /// </remarks>
        ///
        /// <param name="body">
        /// The celestial body whose angle from the Sun is to be measured.
        /// Not allowed to be `Body.Earth`.
        /// </param>
        ///
        /// <param name="time">
        /// The time at which the observation is made.
        /// </param>
        ///
        /// <returns>
        /// Returns the angle in degrees between the Sun and the specified body as
        /// seen from the center of the Earth.
        /// </returns>
        public static double AngleFromSun(Body body, AstroTime time)
        {
            if (body == Body.Earth)
                throw new EarthNotAllowedException();

            AstroVector sv = GeoVector(Body.Sun, time, Aberration.Corrected);
            AstroVector bv = GeoVector(body, time, Aberration.Corrected);
            return AngleBetween(sv, bv);
        }

        private static double AngleBetween(AstroVector a, AstroVector b)
        {
            double r = a.Length() * b.Length();
            if (r < 1.0e-8)
                throw new Exception("Cannot find angle between vectors because they are too short.");

            double dot = (a.x*b.x + a.y*b.y + a.z*b.z) / r;

            if (dot <= -1.0)
                return 180.0;

            if (dot >= +1.0)
                return 0.0;

            return RAD2DEG * Math.Acos(dot);
        }

        /// <summary>
        ///      Finds the date and time of the Moon's closest distance (perigee)
        ///      or farthest distance (apogee) with respect to the Earth.
        /// </summary>
        /// <remarks>
        /// Given a date and time to start the search in `startTime`, this function finds the
        /// next date and time that the center of the Moon reaches the closest or farthest point
        /// in its orbit with respect to the center of the Earth, whichever comes first
        /// after `startTime`.
        ///
        /// The closest point is called *perigee* and the farthest point is called *apogee*.
        /// The word *apsis* refers to either event.
        ///
        /// To iterate through consecutive alternating perigee and apogee events, call `Astronomy.SearchLunarApsis`
        /// once, then use the return value to call #Astronomy.NextLunarApsis. After that,
        /// keep feeding the previous return value from `Astronomy.NextLunarApsis` into another
        /// call of `Astronomy.NextLunarApsis` as many times as desired.
        /// </remarks>
        /// <param name="startTime">
        ///      The date and time at which to start searching for the next perigee or apogee.
        /// </param>
        /// <returns>
        /// Returns an #ApsisInfo structure containing information about the next lunar apsis.
        /// </returns>
        public static ApsisInfo SearchLunarApsis(AstroTime startTime)
        {
            const double increment = 5.0;   /* number of days to skip in each iteration */
            var positive_slope = new SearchContext_MoonDistanceSlope(+1);
            var negative_slope = new SearchContext_MoonDistanceSlope(-1);

            /*
                Check the rate of change of the distance dr/dt at the start time.
                If it is positive, the Moon is currently getting farther away,
                so start looking for apogee.
                Conversely, if dr/dt < 0, start looking for perigee.
                Either way, the polarity of the slope will change, so the product will be negative.
                Handle the crazy corner case of exactly touching zero by checking for m1*m2 <= 0.
            */
            AstroTime t1 = startTime;
            double m1 = positive_slope.Eval(t1);
            for (int iter=0; iter * increment < 2.0 * Astronomy.MEAN_SYNODIC_MONTH; ++iter)
            {
                AstroTime t2 = t1.AddDays(increment);
                double m2 = positive_slope.Eval(t2);
                if (m1 * m2 <= 0.0)
                {
                    /* There is a change of slope polarity within the time range [t1, t2]. */
                    /* Therefore this time range contains an apsis. */
                    /* Figure out whether it is perigee or apogee. */

                    AstroTime search;
                    ApsisKind kind;
                    if (m1 < 0.0 || m2 > 0.0)
                    {
                        /* We found a minimum-distance event: perigee. */
                        /* Search the time range for the time when the slope goes from negative to positive. */
                        search = Search(positive_slope, t1, t2, 1.0);
                        kind = ApsisKind.Pericenter;
                    }
                    else if (m1 > 0.0 || m2 < 0.0)
                    {
                        /* We found a maximum-distance event: apogee. */
                        /* Search the time range for the time when the slope goes from positive to negative. */
                        search = Search(negative_slope, t1, t2, 1.0);
                        kind = ApsisKind.Apocenter;
                    }
                    else
                    {
                        /* This should never happen. It should not be possible for both slopes to be zero. */
                        throw new Exception("Internal error with slopes in SearchLunarApsis");
                    }

                    if (search == null)
                        throw new Exception("Failed to find slope transition in lunar apsis search.");

                    double dist_au = SearchContext_MoonDistanceSlope.MoonDistance(search);
                    return new ApsisInfo(search, kind, dist_au);
                }
                /* We have not yet found a slope polarity change. Keep searching. */
                t1 = t2;
                m1 = m2;
            }

            /* It should not be possible to fail to find an apsis within 2 synodic months. */
            throw new Exception("Internal error: should have found lunar apsis within 2 synodic months.");
        }

        /// <summary>
        /// Finds the next lunar perigee or apogee event in a series.
        /// </summary>
        /// <remarks>
        /// This function requires an #ApsisInfo value obtained from a call
        /// to #Astronomy.SearchLunarApsis or `Astronomy.NextLunarApsis`. Given
        /// an apogee event, this function finds the next perigee event, and vice versa.
        ///
        /// See #Astronomy.SearchLunarApsis for more details.
        /// </remarks>
        /// <param name="apsis">
        /// An apsis event obtained from a call to #Astronomy.SearchLunarApsis or `Astronomy.NextLunarApsis`.
        /// See #Astronomy.SearchLunarApsis for more details.
        /// </param>
        /// <returns>
        /// Same as the return value for #Astronomy.SearchLunarApsis.
        /// </returns>
        public static ApsisInfo NextLunarApsis(ApsisInfo apsis)
        {
            const double skip = 11.0;   // number of days to skip to start looking for next apsis event

            if (apsis.kind != ApsisKind.Pericenter && apsis.kind != ApsisKind.Apocenter)
                throw new ArgumentException("Invalid apsis kind");

            AstroTime time = apsis.time.AddDays(skip);
            ApsisInfo next =  SearchLunarApsis(time);
            if ((int)next.kind + (int)apsis.kind != 1)
                throw new Exception(string.Format("Internal error: previous apsis was {0}, but found {1} for next apsis.", apsis.kind, next.kind));
            return next;
        }

        private static ApsisInfo PlanetExtreme(Body body, ApsisKind kind, AstroTime start_time, double dayspan)
        {
            double direction = (kind == ApsisKind.Apocenter) ? +1.0 : -1.0;
            const int npoints = 10;

            for(;;)
            {
                double interval = dayspan / (npoints - 1);

                if (interval < 1.0 / 1440.0)    /* iterate until uncertainty is less than one minute */
                {
                    AstroTime apsis_time = start_time.AddDays(interval / 2.0);
                    double dist_au = HelioDistance(body, apsis_time);
                    return new ApsisInfo(apsis_time, kind, dist_au);
                }

                int best_i = -1;
                double best_dist = 0.0;
                for (int i=0; i < npoints; ++i)
                {
                    AstroTime time = start_time.AddDays(i * interval);
                    double dist = direction * HelioDistance(body, time);
                    if (i==0 || dist > best_dist)
                    {
                        best_i = i;
                        best_dist = dist;
                    }
                }

                /* Narrow in on the extreme point. */
                start_time = start_time.AddDays((best_i - 1) * interval);
                dayspan = 2.0 * interval;
            }
        }

        private static ApsisInfo BruteSearchPlanetApsis(Body body, AstroTime startTime)
        {
            const int npoints = 100;
            int i;
            var perihelion = new ApsisInfo();
            var aphelion = new ApsisInfo();

            /*
                Neptune is a special case for two reasons:
                1. Its orbit is nearly circular (low orbital eccentricity).
                2. It is so distant from the Sun that the orbital period is very long.
                Put together, this causes wobbling of the Sun around the Solar System Barycenter (SSB)
                to be so significant that there are 3 local minima in the distance-vs-time curve
                near each apsis. Therefore, unlike for other planets, we can't use an optimized
                algorithm for finding dr/dt = 0.
                Instead, we use a dumb, brute-force algorithm of sampling and finding min/max
                heliocentric distance.

                There is a similar problem in the TOP2013 model for Pluto:
                Its position vector has high-frequency oscillations that confuse the
                slope-based determination of apsides.
            */

            /*
                Rewind approximately 30 degrees in the orbit,
                then search forward for 270 degrees.
                This is a very cautious way to prevent missing an apsis.
                Typically we will find two apsides, and we pick whichever
                apsis is ealier, but after startTime.
                Sample points around this orbital arc and find when the distance
                is greatest and smallest.
            */
            double period = PlanetOrbitalPeriod(body);
            AstroTime t1 = startTime.AddDays(period * ( -30.0 / 360.0));
            AstroTime t2 = startTime.AddDays(period * (+270.0 / 360.0));
            AstroTime t_min = t1;
            AstroTime t_max = t1;
            double min_dist = -1.0;
            double max_dist = -1.0;
            double interval = (t2.ut - t1.ut) / (npoints - 1.0);

            for (i=0; i < npoints; ++i)
            {
                AstroTime time = t1.AddDays(i * interval);
                double dist = HelioDistance(body, time);
                if (i == 0)
                {
                    max_dist = min_dist = dist;
                }
                else
                {
                    if (dist > max_dist)
                    {
                        max_dist = dist;
                        t_max = time;
                    }
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        t_min = time;
                    }
                }
            }

            t1 = t_min.AddDays(-2 * interval);
            perihelion = PlanetExtreme(body, ApsisKind.Pericenter, t1, 4 * interval);

            t1 = t_max.AddDays(-2 * interval);
            aphelion = PlanetExtreme(body, ApsisKind.Apocenter, t1, 4 * interval);

            if (perihelion.time.tt >= startTime.tt)
            {
                if (aphelion.time.tt >= startTime.tt)
                {
                    /* Perihelion and aphelion are both valid. Pick the one that comes first. */
                    if (aphelion.time.tt < perihelion.time.tt)
                        return aphelion;
                }
                return perihelion;
            }

            if (aphelion.time.tt >= startTime.tt)
                return aphelion;

            throw new Exception("Internal error: failed to find planet apsis.");
        }


        /// <summary>
        /// Finds the date and time of a planet's perihelion (closest approach to the Sun)
        /// or aphelion (farthest distance from the Sun) after a given time.
        /// </summary>
        /// <remarks>
        /// Given a date and time to start the search in `startTime`, this function finds the
        /// next date and time that the center of the specified planet reaches the closest or farthest point
        /// in its orbit with respect to the center of the Sun, whichever comes first
        /// after `startTime`.
        ///
        /// The closest point is called *perihelion* and the farthest point is called *aphelion*.
        /// The word *apsis* refers to either event.
        ///
        /// To iterate through consecutive alternating perihelion and aphelion events,
        /// call `Astronomy.SearchPlanetApsis` once, then use the return value to call
        /// #Astronomy.NextPlanetApsis. After that, keep feeding the previous return value
        /// from `Astronomy.NextPlanetApsis` into another call of `Astronomy.NextPlanetApsis`
        /// as many times as desired.
        /// </remarks>
        /// <param name="body">
        /// The planet for which to find the next perihelion/aphelion event.
        /// Not allowed to be `Body.Sun` or `Body.Moon`.
        /// </param>
        /// <param name="startTime">
        /// The date and time at which to start searching for the next perihelion or aphelion.
        /// </param>
        /// <returns>
        /// Returns a structure in which `time` holds the date and time of the next planetary apsis,
        /// `kind` holds either `ApsisKind.Pericenter` for perihelion or `ApsisKind.Apocenter` for aphelion.
        /// and distance values `dist_au` (astronomical units) and `dist_km` (kilometers).
        /// </returns>
        public static ApsisInfo SearchPlanetApsis(Body body, AstroTime startTime)
        {
            if (body == Body.Neptune || body == Body.Pluto)
                return BruteSearchPlanetApsis(body, startTime);

            var positive_slope = new SearchContext_PlanetDistanceSlope(+1.0, body);
            var negative_slope = new SearchContext_PlanetDistanceSlope(-1.0, body);
            double orbit_period_days = PlanetOrbitalPeriod(body);
            double increment = orbit_period_days / 6.0;
            AstroTime t1 = startTime;
            double m1 = positive_slope.Eval(t1);
            for (int iter = 0; iter * increment < 2.0 * orbit_period_days; ++iter)
            {
                AstroTime t2 = t1.AddDays(increment);
                double m2 = positive_slope.Eval(t2);
                if (m1 * m2 <= 0.0)
                {
                    /* There is a change of slope polarity within the time range [t1, t2]. */
                    /* Therefore this time range contains an apsis. */
                    /* Figure out whether it is perihelion or aphelion. */

                    SearchContext_PlanetDistanceSlope slope_func;
                    ApsisKind kind;
                    if (m1 < 0.0 || m2 > 0.0)
                    {
                        /* We found a minimum-distance event: perihelion. */
                        /* Search the time range for the time when the slope goes from negative to positive. */
                        slope_func = positive_slope;
                        kind = ApsisKind.Pericenter;
                    }
                    else if (m1 > 0.0 || m2 < 0.0)
                    {
                        /* We found a maximum-distance event: aphelion. */
                        /* Search the time range for the time when the slope goes from positive to negative. */
                        slope_func = negative_slope;
                        kind = ApsisKind.Apocenter;
                    }
                    else
                    {
                        /* This should never happen. It should not be possible for both slopes to be zero. */
                        throw new Exception("Internal error with slopes in SearchPlanetApsis");
                    }

                    AstroTime search = Search(slope_func, t1, t2, 1.0);
                    if (search == null)
                        throw new Exception("Failed to find slope transition in planetary apsis search.");

                    double dist = HelioDistance(body, search);
                    return new ApsisInfo(search, kind, dist);
                }
                /* We have not yet found a slope polarity change. Keep searching. */
                t1 = t2;
                m1 = m2;
            }
            /* It should not be possible to fail to find an apsis within 2 planet orbits. */
            throw new Exception("Internal error: should have found planetary apsis within 2 orbital periods.");
        }

        /// <summary>
        /// Finds the next planetary perihelion or aphelion event in a series.
        /// </summary>
        /// <remarks>
        /// This function requires an #ApsisInfo value obtained from a call
        /// to #Astronomy.SearchPlanetApsis or `Astronomy.NextPlanetApsis`.
        /// Given an aphelion event, this function finds the next perihelion event, and vice versa.
        /// See #Astronomy.SearchPlanetApsis for more details.
        /// </remarks>
        /// <param name="body">
        /// The planet for which to find the next perihelion/aphelion event.
        /// Not allowed to be `Body.Sun` or `Body.Moon`.
        /// Must match the body passed into the call that produced the `apsis` parameter.
        /// </param>
        /// <param name="apsis">
        /// An apsis event obtained from a call to #Astronomy.SearchPlanetApsis or `Astronomy.NextPlanetApsis`.
        /// </param>
        /// <returns>
        /// Same as the return value for #Astronomy.SearchPlanetApsis.
        /// </returns>
        public static ApsisInfo NextPlanetApsis(Body body, ApsisInfo apsis)
        {
            if (apsis.kind != ApsisKind.Apocenter && apsis.kind != ApsisKind.Pericenter)
                throw new ArgumentException("Invalid apsis kind");

            /* skip 1/4 of an orbit before starting search again */
            double skip = 0.25 * PlanetOrbitalPeriod(body);
            if (skip <= 0.0)
                throw new InvalidBodyException(body);

            AstroTime time = apsis.time.AddDays(skip);
            ApsisInfo next = SearchPlanetApsis(body, time);

            /* Verify that we found the opposite apsis from the previous one. */
            if ((int)next.kind + (int)apsis.kind != 1)
                throw new Exception(string.Format("Internal error: previous apsis was {0}, but found {1} for next apsis.", apsis.kind, next.kind));

            return next;
        }


        // We can get away with creating a single EarthShadowSlope context
        // because it contains no state and it has no side-effects.
        // This reduces memory allocation overhead.
        private static readonly SearchContext_EarthShadowSlope earthShadowSlopeContext = new SearchContext_EarthShadowSlope();

        private static ShadowInfo PeakEarthShadow(AstroTime search_center_time)
        {
            const double window = 0.03;        /* initial search window, in days, before/after given time */
            AstroTime t1 = search_center_time.AddDays(-window);
            AstroTime t2 = search_center_time.AddDays(+window);
            AstroTime tx = Search(earthShadowSlopeContext, t1, t2, 1.0);
            return EarthShadow(tx);
        }


        /// <summary>Searches for a lunar eclipse.</summary>
        /// <remarks>
        /// This function finds the first lunar eclipse that occurs after `startTime`.
        /// A lunar eclipse may be penumbral, partial, or total.
        /// See #LunarEclipseInfo for more information.
        /// To find a series of lunar eclipses, call this function once,
        /// then keep calling #Astronomy.NextLunarEclipse as many times as desired,
        /// passing in the `center` value returned from the previous call.
        /// </remarks>
        /// <param name="startTime">
        ///      The date and time for starting the search for a lunar eclipse.
        /// </param>
        /// <returns>
        ///      A #LunarEclipseInfo structure containing information about the lunar eclipse.
        /// </returns>
        public static LunarEclipseInfo SearchLunarEclipse(AstroTime startTime)
        {
            const double PruneLatitude = 1.8;   /* full Moon's ecliptic latitude above which eclipse is impossible */
            // Iterate through consecutive full moons until we find any kind of lunar eclipse.
            AstroTime fmtime = startTime;
            for (int fmcount=0; fmcount < 12; ++fmcount)
            {
                // Search for the next full moon. Any eclipse will be near it.
                AstroTime fullmoon = SearchMoonPhase(180.0, fmtime, 40.0);

                /*
                    Pruning: if the full Moon's ecliptic latitude is too large,
                    a lunar eclipse is not possible. Avoid needless work searching for
                    the minimum moon distance.
                */
                var mc = new MoonContext(fullmoon.tt / 36525.0);
                MoonResult mr = mc.CalcMoon();
                if (RAD2DEG * Math.Abs(mr.geo_eclip_lat) < PruneLatitude)
                {
                    // Search near the full moon for the time when the center of the Moon
                    // is closest to the line passing through the centers of the Sun and Earth.
                    ShadowInfo shadow = PeakEarthShadow(fullmoon);

                    if (shadow.r < shadow.p + MOON_MEAN_RADIUS_KM)
                    {
                        // This is at least a penumbral eclipse. We will return a result.
                        EclipseKind kind = EclipseKind.Penumbral;
                        double sd_total = 0.0;
                        double sd_partial = 0.0;
                        double sd_penum = ShadowSemiDurationMinutes(shadow.time, shadow.p + MOON_MEAN_RADIUS_KM, 200.0);

                        if (shadow.r < shadow.k + MOON_MEAN_RADIUS_KM)
                        {
                            // This is at least a partial eclipse.
                            kind = EclipseKind.Partial;
                            sd_partial = ShadowSemiDurationMinutes(shadow.time, shadow.k + MOON_MEAN_RADIUS_KM, sd_penum);

                            if (shadow.r + MOON_MEAN_RADIUS_KM < shadow.k)
                            {
                                // This is a total eclipse.
                                kind = EclipseKind.Total;
                                sd_total = ShadowSemiDurationMinutes(shadow.time, shadow.k - MOON_MEAN_RADIUS_KM, sd_partial);
                            }
                        }
                        return new LunarEclipseInfo(kind, shadow.time, sd_penum, sd_partial, sd_total);
                    }
                }

                // We didn't find an eclipse on this full moon, so search for the next one.
                fmtime = fullmoon.AddDays(10.0);
            }

            // This should never happen, because there should be at least 2 lunar eclipses per year.
            throw new Exception("Internal error: failed to find lunar eclipse within 12 full moons.");
        }


        /// <summary>Searches for the next lunar eclipse in a series.</summary>
        /// <remarks>
        /// After using #Astronomy.SearchLunarEclipse to find the first lunar eclipse
        /// in a series, you can call this function to find the next consecutive lunar eclipse.
        /// Pass in the `center` value from the #LunarEclipseInfo returned by the
        /// previous call to `Astronomy.SearchLunarEclipse` or `Astronomy.NextLunarEclipse`
        /// to find the next lunar eclipse.
        /// </remarks>
        ///
        /// <param name="prevEclipseTime">
        /// A date and time near a full moon. Lunar eclipse search will start at the next full moon.
        /// </param>
        ///
        /// <returns>
        /// A #LunarEclipseInfo structure containing information about the lunar eclipse.
        /// </returns>
        public static LunarEclipseInfo NextLunarEclipse(AstroTime prevEclipseTime)
        {
            AstroTime startTime = prevEclipseTime.AddDays(10.0);
            return SearchLunarEclipse(startTime);
        }


        private static double ShadowSemiDurationMinutes(AstroTime center_time, double radius_limit, double window_minutes)
        {
            // Search backwards and forwards from the center time until shadow axis distance crosses radius limit.
            double window = window_minutes / (24.0 * 60.0);
            AstroTime before = center_time.AddDays(-window);
            AstroTime after  = center_time.AddDays(+window);
            AstroTime t1 = Search(new SearchContext_EarthShadow(radius_limit, -1.0), before, center_time, 1.0);
            AstroTime t2 = Search(new SearchContext_EarthShadow(radius_limit, +1.0), center_time, after, 1.0);
            return (t2.ut - t1.ut) * ((24.0 * 60.0) / 2.0);    // convert days to minutes and average the semi-durations.
        }


        /// <summary>
        /// Searches for a solar eclipse visible anywhere on the Earth's surface.
        /// </summary>
        /// <remarks>
        /// This function finds the first solar eclipse that occurs after `startTime`.
        /// A solar eclipse may be partial, annular, or total.
        /// See #GlobalSolarEclipseInfo for more information.
        /// To find a series of solar eclipses, call this function once,
        /// then keep calling #Astronomy.NextGlobalSolarEclipse as many times as desired,
        /// passing in the `peak` value returned from the previous call.
        /// </remarks>
        /// <param name="startTime">The date and time for starting the search for a solar eclipse.</param>
        public static GlobalSolarEclipseInfo SearchGlobalSolarEclipse(AstroTime startTime)
        {
            const double PruneLatitude = 1.8;   /* Moon's ecliptic latitude beyond which eclipse is impossible */

            /* Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth. */
            AstroTime nmtime = startTime;
            for (int nmcount=0; nmcount < 12; ++nmcount)
            {
                /* Search for the next new moon. Any eclipse will be near it. */
                AstroTime newmoon = SearchMoonPhase(0.0, nmtime, 40.0);

                /* Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible. */
                double eclip_lat = MoonEclipticLatitudeDegrees(newmoon);
                if (Math.Abs(eclip_lat) < PruneLatitude)
                {
                    /* Search near the new moon for the time when the center of the Earth */
                    /* is closest to the line passing through the centers of the Sun and Moon. */
                    ShadowInfo shadow = PeakMoonShadow(newmoon);
                    if (shadow.r < shadow.p + EARTH_MEAN_RADIUS_KM)
                    {
                        /* This is at least a partial solar eclipse visible somewhere on Earth. */
                        /* Try to find an intersection between the shadow axis and the Earth's oblate geoid. */
                        return GeoidIntersect(shadow);
                    }
                }

                /* We didn't find an eclipse on this new moon, so search for the next one. */
                nmtime = newmoon.AddDays(10.0);
            }

            /* Safety valve to prevent infinite loop. */
            /* This should never happen, because at least 2 solar eclipses happen per year. */
            throw new Exception("Failure to find global solar eclipse.");
        }


        /// <summary>
        /// Searches for the next global solar eclipse in a series.
        /// </summary>
        /// <remarks>
        /// After using #Astronomy.SearchGlobalSolarEclipse to find the first solar eclipse
        /// in a series, you can call this function to find the next consecutive solar eclipse.
        /// Pass in the `peak` value from the #GlobalSolarEclipseInfo returned by the
        /// previous call to `Astronomy.SearchGlobalSolarEclipse` or `Astronomy.NextGlobalSolarEclipse`
        /// to find the next solar eclipse.
        /// </remarks>
        /// <param name="prevEclipseTime">
        /// A date and time near a new moon. Solar eclipse search will start at the next new moon.
        /// </param>
        public static GlobalSolarEclipseInfo NextGlobalSolarEclipse(AstroTime prevEclipseTime)
        {
            AstroTime startTime = prevEclipseTime.AddDays(10.0);
            return SearchGlobalSolarEclipse(startTime);
        }


        private static GlobalSolarEclipseInfo GeoidIntersect(ShadowInfo shadow)
        {
            var eclipse = new GlobalSolarEclipseInfo();
            eclipse.kind = EclipseKind.Partial;
            eclipse.peak = shadow.time;
            eclipse.distance = shadow.r;
            eclipse.latitude = eclipse.longitude = double.NaN;

            /*
                We want to calculate the intersection of the shadow axis with the Earth's geoid.
                First we must convert EQJ (equator of J2000) coordinates to EQD (equator of date)
                coordinates that are perfectly aligned with the Earth's equator at this
                moment in time.
            */
            RotationMatrix rot = Rotation_EQJ_EQD(shadow.time);

            AstroVector v = RotateVector(rot, shadow.dir);        /* shadow-axis vector in equator-of-date coordinates */
            AstroVector e = RotateVector(rot, shadow.target);     /* lunacentric Earth in equator-of-date coordinates */

            /*
                Convert all distances from AU to km.
                But dilate the z-coordinates so that the Earth becomes a perfect sphere.
                Then find the intersection of the vector with the sphere.
                See p 184 in Montenbruck & Pfleger's "Astronomy on the Personal Computer", second edition.
            */
            v.x *= KM_PER_AU;
            v.y *= KM_PER_AU;
            v.z *= KM_PER_AU / EARTH_FLATTENING;

            e.x *= KM_PER_AU;
            e.y *= KM_PER_AU;
            e.z *= KM_PER_AU / EARTH_FLATTENING;

            /*
                Solve the quadratic equation that finds whether and where
                the shadow axis intersects with the Earth in the dilated coordinate system.
            */
            double R = EARTH_EQUATORIAL_RADIUS_KM;
            double A = v.x*v.x + v.y*v.y + v.z*v.z;
            double B = -2.0 * (v.x*e.x + v.y*e.y + v.z*e.z);
            double C = (e.x*e.x + e.y*e.y + e.z*e.z) - R*R;
            double radic = B*B - 4*A*C;

            if (radic > 0.0)
            {
                /* Calculate the closer of the two intersection points. */
                /* This will be on the day side of the Earth. */
                double u = (-B - Math.Sqrt(radic)) / (2 * A);

                /* Convert lunacentric dilated coordinates to geocentric coordinates. */
                double px = u*v.x - e.x;
                double py = u*v.y - e.y;
                double pz = (u*v.z - e.z) * EARTH_FLATTENING;

                /* Convert cartesian coordinates into geodetic latitude/longitude. */
                double proj = Math.Sqrt(px*px + py*py) * (EARTH_FLATTENING * EARTH_FLATTENING);
                if (proj == 0.0)
                    eclipse.latitude = (pz > 0.0) ? +90.0 : -90.0;
                else
                    eclipse.latitude = RAD2DEG * Math.Atan(pz / proj);

                /* Adjust longitude for Earth's rotation at the given UT. */
                double gast = sidereal_time(eclipse.peak);
                eclipse.longitude = ((RAD2DEG*Math.Atan2(py, px)) - (15*gast)) % 360.0;
                if (eclipse.longitude <= -180.0)
                    eclipse.longitude += 360.0;
                else if (eclipse.longitude > +180.0)
                    eclipse.longitude -= 360.0;

                /* We want to determine whether the observer sees a total eclipse or an annular eclipse. */
                /* We need to perform a series of vector calculations... */
                /* Calculate the inverse rotation matrix, so we can convert EQD to EQJ. */
                RotationMatrix inv = InverseRotation(rot);

                /* Put the EQD geocentric coordinates of the observer into the vector 'o'. */
                /* Also convert back from kilometers to astronomical units. */
                var o = new AstroVector(px / KM_PER_AU, py / KM_PER_AU, pz / KM_PER_AU, shadow.time);

                /* Rotate the observer's geocentric EQD back to the EQJ system. */
                o = RotateVector(inv, o);

                /* Convert geocentric vector to lunacentric vector. */
                o.x += shadow.target.x;
                o.y += shadow.target.y;
                o.z += shadow.target.z;

                /* Recalculate the shadow using a vector from the Moon's center toward the observer. */
                ShadowInfo surface = CalcShadow(MOON_POLAR_RADIUS_KM, shadow.time, o, shadow.dir);

                /* If we did everything right, the shadow distance should be very close to zero. */
                /* That's because we already determined the observer 'o' is on the shadow axis! */
                if (surface.r > 1.0e-9 || surface.r < 0.0)
                    throw new Exception("Invalid surface distance from intersection.");

                eclipse.kind = EclipseKindFromUmbra(surface.k);
            }

            return eclipse;
        }


        private static EclipseKind EclipseKindFromUmbra(double k)
        {
            // The umbra radius tells us what kind of eclipse the observer sees.
            // If the umbra radius is positive, this is a total eclipse. Otherwise, it's annular.
            // HACK: I added a tiny bias (14 meters) to match Espenak test data.
            return (k > 0.014) ? EclipseKind.Total : EclipseKind.Annular;
        }


        private static readonly SearchContext_MoonShadowSlope moonShadowSlopeContext = new SearchContext_MoonShadowSlope();

        private static ShadowInfo PeakMoonShadow(AstroTime search_center_time)
        {
            /* Search for when the Moon's shadow axis is closest to the center of the Earth. */

            const double window = 0.03;     /* days before/after new moon to search for minimum shadow distance */
            AstroTime t1 = search_center_time.AddDays(-window);
            AstroTime t2 = search_center_time.AddDays(+window);
            AstroTime time = Search(moonShadowSlopeContext, t1, t2, 1.0);
            return MoonShadow(time);
        }

        private static ShadowInfo PeakLocalMoonShadow(AstroTime search_center_time, Observer observer)
        {
            /*
                Search for the time near search_center_time that the Moon's shadow comes
                closest to the given observer.
            */
            const double window = 0.2;
            AstroTime t1 = search_center_time.AddDays(-window);
            AstroTime t2 = search_center_time.AddDays(+window);
            var context = new SearchContext_LocalMoonShadowSlope(observer);
            AstroTime time = Search(context, t1, t2, 1.0);
            return LocalMoonShadow(time, observer);
        }

        private static ShadowInfo PeakPlanetShadow(Body body, double planet_radius_km, AstroTime search_center_time)
        {
            /* Search for when the body's shadow is closest to the center of the Earth. */
            const double window = 1.0;     /* days before/after inferior conjunction to search for minimum shadow distance */
            AstroTime t1 = search_center_time.AddDays(-window);
            AstroTime t2 = search_center_time.AddDays(+window);
            var context = new SearchContext_PlanetShadowSlope(body, planet_radius_km);
            AstroTime time = Search(context, t1, t2, 1.0);
            return PlanetShadow(body, planet_radius_km, time);
        }

        private static ShadowInfo CalcShadow(
            double body_radius_km,
            AstroTime time,
            AstroVector target,
            AstroVector dir)
        {
            double u = (dir.x*target.x + dir.y*target.y + dir.z*target.z) / (dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
            double dx = (u * dir.x) - target.x;
            double dy = (u * dir.y) - target.y;
            double dz = (u * dir.z) - target.z;
            double r = KM_PER_AU * Math.Sqrt(dx*dx + dy*dy + dz*dz);
            double k = +SUN_RADIUS_KM - (1.0 + u)*(SUN_RADIUS_KM - body_radius_km);
            double p = -SUN_RADIUS_KM + (1.0 + u)*(SUN_RADIUS_KM + body_radius_km);
            return new ShadowInfo(time, u, r, k, p, target, dir);
        }


        internal static ShadowInfo EarthShadow(AstroTime time)
        {
            /* This function helps find when the Earth's shadow falls upon the Moon. */
            AstroVector e = CalcEarth(time);            /* This function never fails; no need to check return value */
            AstroVector m = GeoMoon(time);    /* This function never fails; no need to check return value */

            return CalcShadow(EARTH_ECLIPSE_RADIUS_KM, time, m, e);
        }


        internal static ShadowInfo MoonShadow(AstroTime time)
        {
            /* This function helps find when the Moon's shadow falls upon the Earth. */

            /*
                This is a variation on the logic in EarthShadow().
                Instead of a heliocentric Earth and a geocentric Moon,
                we want a heliocentric Moon and a lunacentric Earth.
            */

            AstroVector h = CalcEarth(time);    /* heliocentric Earth */
            AstroVector m = GeoMoon(time);      /* geocentric Moon */

            /* Calculate lunacentric Earth. */
            var e = new AstroVector(-m.x, -m.y, -m.z, m.t);

            /* Convert geocentric moon to heliocentric Moon. */
            m.x += h.x;
            m.y += h.y;
            m.z += h.z;

            return CalcShadow(MOON_MEAN_RADIUS_KM, time, e, m);
        }


        internal static ShadowInfo LocalMoonShadow(AstroTime time, Observer observer)
        {
            /* Calculate observer's geocentric position. */
            /* For efficiency, do this first, to populate the earth rotation parameters in 'time'. */
            /* That way they can be recycled instead of recalculated. */
            AstroVector pos = geo_pos(time, observer);

            AstroVector h = CalcEarth(time);    /* heliocentric Earth */
            AstroVector m = GeoMoon(time);      /* geocentric Moon */

            /* Calculate lunacentric location of an observer on the Earth's surface. */
            var o = new AstroVector(pos.x - m.x, pos.y - m.y, pos.z - m.z, time);

            /* Convert geocentric moon to heliocentric Moon. */
            m.x += h.x;
            m.y += h.y;
            m.z += h.z;

            return CalcShadow(MOON_MEAN_RADIUS_KM, time, o, m);
        }


        internal static ShadowInfo PlanetShadow(Body body, double planet_radius_km, AstroTime time)
        {
            /* Calculate light-travel-corrected vector from Earth to planet. */
            AstroVector g = GeoVector(body, time, Aberration.None);

            /* Calculate light-travel-corrected vector from Earth to Sun. */
            AstroVector e = GeoVector(Body.Sun, time, Aberration.None);

            /* Deduce light-travel-corrected vector from Sun to planet. */
            var p = new AstroVector(g.x - e.x, g.y - e.y, g.z - e.z, time);

            /* Calcluate Earth's position from the planet's point of view. */
            e.x = -g.x;
            e.y = -g.y;
            e.z = -g.z;

            return CalcShadow(planet_radius_km, time, e, p);
        }


        private static double MoonEclipticLatitudeDegrees(AstroTime time)
        {
            var context = new MoonContext(time.tt / 36525.0);
            MoonResult moon = context.CalcMoon();
            return RAD2DEG * moon.geo_eclip_lat;
        }

        /// <summary>
        /// Searches for a solar eclipse visible at a specific location on the Earth's surface.
        /// </summary>
        /// <remarks>
        /// This function finds the first solar eclipse that occurs after `startTime`.
        /// A solar eclipse may be partial, annular, or total.
        /// See #LocalSolarEclipseInfo for more information.
        ///
        /// To find a series of solar eclipses, call this function once,
        /// then keep calling #Astronomy.NextLocalSolarEclipse as many times as desired,
        /// passing in the `peak` value returned from the previous call.
        ///
        /// IMPORTANT: An eclipse reported by this function might be partly or
        /// completely invisible to the observer due to the time of day.
        /// See #LocalSolarEclipseInfo for more information about this topic.
        /// </remarks>
        ///
        /// <param name="startTime">The date and time for starting the search for a solar eclipse.</param>
        /// <param name="observer">The geographic location of the observer.</param>
        public static LocalSolarEclipseInfo SearchLocalSolarEclipse(AstroTime startTime, Observer observer)
        {
            const double PruneLatitude = 1.8;   /* Moon's ecliptic latitude beyond which eclipse is impossible */

            /* Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth. */
            AstroTime nmtime = startTime;
            for(;;)
            {
                /* Search for the next new moon. Any eclipse will be near it. */
                AstroTime newmoon = SearchMoonPhase(0.0, nmtime, 40.0);

                /* Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible. */
                double eclip_lat = MoonEclipticLatitudeDegrees(newmoon);
                if (Math.Abs(eclip_lat) < PruneLatitude)
                {
                    /* Search near the new moon for the time when the observer */
                    /* is closest to the line passing through the centers of the Sun and Moon. */
                    ShadowInfo shadow = PeakLocalMoonShadow(newmoon, observer);
                    if (shadow.r < shadow.p)
                    {
                        /* This is at least a partial solar eclipse for the observer. */
                        LocalSolarEclipseInfo eclipse = LocalEclipse(shadow, observer);

                        /* Ignore any eclipse that happens completely at night. */
                        /* More precisely, the center of the Sun must be above the horizon */
                        /* at the beginning or the end of the eclipse, or we skip the event. */
                        if (eclipse.partial_begin.altitude > 0.0 || eclipse.partial_end.altitude > 0.0)
                            return eclipse;
                    }
                }

                /* We didn't find an eclipse on this new moon, so search for the next one. */
                nmtime = newmoon.AddDays(10.0);
            }
        }


        /// <summary>
        /// Searches for the next local solar eclipse in a series.
        /// </summary>
        ///
        /// <remarks>
        /// After using #Astronomy.SearchLocalSolarEclipse to find the first solar eclipse
        /// in a series, you can call this function to find the next consecutive solar eclipse.
        /// Pass in the `peak` value from the #LocalSolarEclipseInfo returned by the
        /// previous call to `Astronomy.SearchLocalSolarEclipse` or `Astronomy.NextLocalSolarEclipse`
        /// to find the next solar eclipse.
        /// </remarks>
        ///
        /// <param name="prevEclipseTime">
        ///      A date and time near a new moon. Solar eclipse search will start at the next new moon.
        /// </param>
        ///
        /// <param name="observer">
        ///      The geographic location of the observer.
        /// </param>
        public static LocalSolarEclipseInfo NextLocalSolarEclipse(AstroTime prevEclipseTime, Observer observer)
        {
            AstroTime startTime = prevEclipseTime.AddDays(10.0);
            return SearchLocalSolarEclipse(startTime, observer);
        }


        private static double local_partial_distance(ShadowInfo shadow)
        {
            return shadow.p - shadow.r;
        }

        private static double local_total_distance(ShadowInfo shadow)
        {
            /* Must take the absolute value of the umbra radius 'k' */
            /* because it can be negative for an annular eclipse. */
            return Math.Abs(shadow.k) - shadow.r;
        }

        private static LocalSolarEclipseInfo LocalEclipse(ShadowInfo shadow, Observer observer)
        {
            const double PARTIAL_WINDOW = 0.2;
            const double TOTAL_WINDOW = 0.01;

            var eclipse = new LocalSolarEclipseInfo();
            eclipse.peak = CalcEvent(observer, shadow.time);
            AstroTime t1 = shadow.time.AddDays(-PARTIAL_WINDOW);
            AstroTime t2 = shadow.time.AddDays(+PARTIAL_WINDOW);
            eclipse.partial_begin = LocalEclipseTransition(observer, +1.0, local_partial_distance, t1, shadow.time);
            eclipse.partial_end   = LocalEclipseTransition(observer, -1.0, local_partial_distance, shadow.time, t2);

            if (shadow.r < Math.Abs(shadow.k))      /* take absolute value of 'k' to handle annular eclipses too. */
            {
                t1 = shadow.time.AddDays(-TOTAL_WINDOW);
                t2 = shadow.time.AddDays(+TOTAL_WINDOW);
                eclipse.total_begin = LocalEclipseTransition(observer, +1.0, local_total_distance, t1, shadow.time);
                eclipse.total_end = LocalEclipseTransition(observer, -1.0, local_total_distance, shadow.time, t2);
                eclipse.kind = EclipseKindFromUmbra(shadow.k);
            }
            else
            {
                eclipse.kind = EclipseKind.Partial;
            }

            return eclipse;
        }

        private static EclipseEvent LocalEclipseTransition(
            Observer observer,
            double direction,
            Func<ShadowInfo,double> func,
            AstroTime t1,
            AstroTime t2)
        {
            var context = new SearchContext_LocalEclipseTransition(func, direction, observer);
            AstroTime search = Search(context, t1, t2, 1.0);
            if (search == null)
                throw new Exception("Local eclipse transition search failed.");
            return CalcEvent(observer, search);
        }

        private static EclipseEvent CalcEvent(Observer observer, AstroTime time)
        {
            var evt = new EclipseEvent();
            evt.time = time;
            evt.altitude = SunAltitude(time, observer);
            return evt;
        }

        private static double SunAltitude(AstroTime time, Observer observer)
        {
            Equatorial equ = Equator(Body.Sun, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
            Topocentric hor = Horizon(time, observer, equ.ra, equ.dec, Refraction.Normal);
            return hor.altitude;
        }


        private static AstroTime PlanetTransitBoundary(
            Body body,
            double planet_radius_km,
            AstroTime t1,
            AstroTime t2,
            double direction)
        {
            /* Search for the time the planet's penumbra begins/ends making contact with the center of the Earth. */
            var context = new SearchContext_PlanetShadowBoundary(body, planet_radius_km, direction);
            AstroTime time = Search(context, t1, t2, 1.0);
            if (time == null)
                throw new Exception("Planet transit boundary search failed");
            return time;
        }


        /// <summary>
        /// Searches for the first transit of Mercury or Venus after a given date.
        /// </summary>
        /// <remarks>
        /// Finds the first transit of Mercury or Venus after a specified date.
        /// A transit is when an inferior planet passes between the Sun and the Earth
        /// so that the silhouette of the planet is visible against the Sun in the background.
        /// To continue the search, pass the `finish` time in the returned structure to
        /// #Astronomy.NextTransit.
        /// </remarks>
        /// <param name="body">
        /// The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`.
        /// </param>
        /// <param name="startTime">
        /// The date and time for starting the search for a transit.
        /// </param>
        public static TransitInfo SearchTransit(Body body, AstroTime startTime)
        {
            const double threshold_angle = 0.4;     /* maximum angular separation to attempt transit calculation */
            const double dt_days = 1.0;

            // Validate the planet and find its mean radius.
            double planet_radius_km;
            switch (body)
            {
                case Body.Mercury:
                    planet_radius_km = 2439.7;
                    break;

                case Body.Venus:
                    planet_radius_km = 6051.8;
                    break;

                default:
                    throw new InvalidBodyException(body);
            }

            AstroTime search_time = startTime;
            for(;;)
            {
                /*
                    Search for the next inferior conjunction of the given planet.
                    This is the next time the Earth and the other planet have the same
                    ecliptic longitude as seen from the Sun.
                */
                AstroTime conj = SearchRelativeLongitude(body, 0.0, search_time);

                /* Calculate the angular separation between the body and the Sun at this time. */
                double separation = AngleFromSun(body, conj);

                if (separation < threshold_angle)
                {
                    /*
                        The planet's angular separation from the Sun is small enough
                        to consider it a transit candidate.
                        Search for the moment when the line passing through the Sun
                        and planet are closest to the Earth's center.
                    */
                    ShadowInfo shadow = PeakPlanetShadow(body, planet_radius_km, conj);

                    if (shadow.r < shadow.p)        /* does the planet's penumbra touch the Earth's center? */
                    {
                        var transit = new TransitInfo();

                        /* Find the beginning and end of the penumbral contact. */
                        AstroTime tx = shadow.time.AddDays(-dt_days);
                        transit.start = PlanetTransitBoundary(body, planet_radius_km, tx, shadow.time, -1.0);

                        tx = shadow.time.AddDays(+dt_days);
                        transit.finish = PlanetTransitBoundary(body, planet_radius_km, shadow.time, tx, +1.0);

                        transit.peak = shadow.time;
                        transit.separation = 60.0 * AngleFromSun(body, shadow.time);
                        return transit;
                    }
                }

                /* This inferior conjunction was not a transit. Try the next inferior conjunction. */
                search_time = conj.AddDays(10.0);
            }
        }


        /// <summary>
        /// Searches for another transit of Mercury or Venus.
        /// </summary>
        /// <remarks>
        /// After calling #Astronomy.SearchTransit to find a transit of Mercury or Venus,
        /// this function finds the next transit after that.
        /// Keep calling this function as many times as you want to keep finding more transits.
        /// </remarks>
        /// <param name="body">
        /// The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`.
        /// </param>
        /// <param name="prevTransitTime">
        /// A date and time near the previous transit.
        /// </param>
        public static TransitInfo NextTransit(Body body, AstroTime prevTransitTime)
        {
            AstroTime startTime = prevTransitTime.AddDays(100.0);
            return SearchTransit(body, startTime);
        }

        /// <summary>
        /// Finds visual magnitude, phase angle, and other illumination information about a celestial body.
        /// </summary>
        /// <remarks>
        /// This function calculates information about how bright a celestial body appears from the Earth,
        /// reported as visual magnitude, which is a smaller (or even negative) number for brighter objects
        /// and a larger number for dimmer objects.
        ///
        /// For bodies other than the Sun, it reports a phase angle, which is the angle in degrees between
        /// the Sun and the Earth, as seen from the center of the body. Phase angle indicates what fraction
        /// of the body appears illuminated as seen from the Earth. For example, when the phase angle is
        /// near zero, it means the body appears "full" as seen from the Earth.  A phase angle approaching
        /// 180 degrees means the body appears as a thin crescent as seen from the Earth.  A phase angle
        /// of 90 degrees means the body appears "half full".
        /// For the Sun, the phase angle is always reported as 0; the Sun emits light rather than reflecting it,
        /// so it doesn't have a phase angle.
        ///
        /// When the body is Saturn, the returned structure contains a field `ring_tilt` that holds
        /// the tilt angle in degrees of Saturn's rings as seen from the Earth. A value of 0 means
        /// the rings appear edge-on, and are thus nearly invisible from the Earth. The `ring_tilt` holds
        /// 0 for all bodies other than Saturn.
        /// </remarks>
        /// <param name="body">The Sun, Moon, or any planet other than the Earth.</param>
        /// <param name="time">The date and time of the observation.</param>
        /// <returns>An #IllumInfo structure with fields as documented above.</returns>
        public static IllumInfo Illumination(Body body, AstroTime time)
        {
            if (body == Body.Earth)
                throw new EarthNotAllowedException();

            AstroVector earth = CalcEarth(time);

            AstroVector gc;
            AstroVector hc;
            double phase_angle;
            if (body == Body.Sun)
            {
                gc = new AstroVector(-earth.x, -earth.y, -earth.z, time);
                hc = new AstroVector(0.0, 0.0, 0.0, time);
                // The Sun emits light instead of reflecting it,
                // so we report a placeholder phase angle of 0.
                phase_angle = 0.0;
            }
            else
            {
                if (body == Body.Moon)
                {
                    // For extra numeric precision, use geocentric Moon formula directly.
                    gc = GeoMoon(time);
                    hc = new AstroVector(earth.x + gc.x, earth.y + gc.y, earth.z + gc.z, time);
                }
                else
                {
                    // For planets, the heliocentric vector is more direct to calculate.
                    hc = HelioVector(body, time);
                    gc = new AstroVector(hc.x - earth.x, hc.y - earth.y, hc.z - earth.z, time);
                }

                phase_angle = AngleBetween(gc, hc);
            }

            double geo_dist = gc.Length();
            double helio_dist = hc.Length();
            double ring_tilt = 0.0;

            double mag;
            switch (body)
            {
                case Body.Sun:
                    mag = -0.17 + 5.0*Math.Log10(geo_dist / AU_PER_PARSEC);
                    break;

                case Body.Moon:
                    mag = MoonMagnitude(phase_angle, helio_dist, geo_dist);
                    break;

                case Body.Saturn:
                    mag = SaturnMagnitude(phase_angle, helio_dist, geo_dist, gc, time, out ring_tilt);
                    break;

                default:
                    mag = VisualMagnitude(body, phase_angle, helio_dist, geo_dist);
                    break;
            }

            return new IllumInfo(time, mag, phase_angle, helio_dist, ring_tilt);
        }

        private static double MoonMagnitude(double phase, double helio_dist, double geo_dist)
        {
            /* https://astronomy.stackexchange.com/questions/10246/is-there-a-simple-analytical-formula-for-the-lunar-phase-brightness-curve */
            double rad = phase * DEG2RAD;
            double rad2 = rad * rad;
            double rad4 = rad2 * rad2;
            double mag = -12.717 + 1.49*Math.Abs(rad) + 0.0431*rad4;
            double moon_mean_distance_au = 385000.6 / KM_PER_AU;
            double geo_au = geo_dist / moon_mean_distance_au;
            mag += 5.0 * Math.Log10(helio_dist * geo_au);
            return mag;
        }

        private static double VisualMagnitude(
            Body body,
            double phase,
            double helio_dist,
            double geo_dist)
        {
            /* For Mercury and Venus, see:  https://iopscience.iop.org/article/10.1086/430212 */
            double c0, c1=0, c2=0, c3=0;
            switch (body)
            {
                case Body.Mercury:
                    c0 = -0.60; c1 = +4.98; c2 = -4.88; c3 = +3.02; break;
                case Body.Venus:
                    if (phase < 163.6)
                    {
                        c0 = -4.47; c1 = +1.03; c2 = +0.57; c3 = +0.13;
                    }
                    else
                    {
                        c0 = 0.98; c1 = -1.02;
                    }
                    break;
                case Body.Mars:        c0 = -1.52; c1 = +1.60;   break;
                case Body.Jupiter:     c0 = -9.40; c1 = +0.50;   break;
                case Body.Uranus:      c0 = -7.19; c1 = +0.25;   break;
                case Body.Neptune:     c0 = -6.87;               break;
                case Body.Pluto:       c0 = -1.00; c1 = +4.00;   break;
                default:
                    throw new InvalidBodyException(body);
            }

            double x = phase / 100;
            double mag = c0 + x*(c1 + x*(c2 + x*c3));
            mag += 5.0 * Math.Log10(helio_dist * geo_dist);
            return mag;
        }

        private static double SaturnMagnitude(
            double phase,
            double helio_dist,
            double geo_dist,
            AstroVector gc,
            AstroTime time,
            out double ring_tilt)
        {
            /* Based on formulas by Paul Schlyter found here: */
            /* http://www.stjarnhimlen.se/comp/ppcomp.html#15 */

            /* We must handle Saturn's rings as a major component of its visual magnitude. */
            /* Find geocentric ecliptic coordinates of Saturn. */
            Ecliptic eclip = EquatorialToEcliptic(gc);

            double ir = DEG2RAD * 28.06;   /* tilt of Saturn's rings to the ecliptic, in radians */
            double Nr = DEG2RAD * (169.51 + (3.82e-5 * time.tt));    /* ascending node of Saturn's rings, in radians */

            /* Find tilt of Saturn's rings, as seen from Earth. */
            double lat = DEG2RAD * eclip.elat;
            double lon = DEG2RAD * eclip.elon;
            double tilt = Math.Asin(Math.Sin(lat)*Math.Cos(ir) - Math.Cos(lat)*Math.Sin(ir)*Math.Sin(lon-Nr));
            double sin_tilt = Math.Sin(Math.Abs(tilt));

            double mag = -9.0 + 0.044*phase;
            mag += sin_tilt*(-2.6 + 1.2*sin_tilt);
            mag += 5.0 * Math.Log10(helio_dist * geo_dist);

            ring_tilt = RAD2DEG * tilt;

            return mag;
        }

        /// <summary>Searches for the date and time Venus will next appear brightest as seen from the Earth.</summary>
        /// <remarks>
        /// This function searches for the date and time Venus appears brightest as seen from the Earth.
        /// Currently only Venus is supported for the `body` parameter, though this could change in the future.
        /// Mercury's peak magnitude occurs at superior conjunction, when it is virtually impossible to see from the Earth,
        /// so peak magnitude events have little practical value for that planet.
        /// Planets other than Venus and Mercury reach peak magnitude at opposition, which can
        /// be found using #Astronomy.SearchRelativeLongitude.
        /// The Moon reaches peak magnitude at full moon, which can be found using
        /// #Astronomy.SearchMoonQuarter or #Astronomy.SearchMoonPhase.
        /// The Sun reaches peak magnitude at perihelion, which occurs each year in January.
        /// However, the difference is minor and has little practical value.
        /// </remarks>
        ///
        /// <param name="body">
        ///      Currently only `Body.Venus` is allowed. Any other value causes an exception.
        ///      See remarks above for more details.
        /// </param>
        /// <param name="startTime">
        ///     The date and time to start searching for the next peak magnitude event.
        /// </param>
        /// <returns>
        ///      See documentation about the return value from #Astronomy.Illumination.
        /// </returns>
        public static IllumInfo SearchPeakMagnitude(Body body, AstroTime startTime)
        {
            /* s1 and s2 are relative longitudes within which peak magnitude of Venus can occur. */
            const double s1 = 10.0;
            const double s2 = 30.0;

            if (body != Body.Venus)
                throw new ArgumentException("Peak magnitude currently is supported for Venus only.");

            var mag_slope = new SearchContext_MagnitudeSlope(body);

            int iter = 0;
            while (++iter <= 2)
            {
                /* Find current heliocentric relative longitude between the */
                /* inferior planet and the Earth. */
                double plon = EclipticLongitude(body, startTime);
                double elon = EclipticLongitude(Body.Earth, startTime);
                double rlon = LongitudeOffset(plon - elon);     // clamp to (-180, +180].

                /* The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees */
                /* because there is a cusp there that causes a discontinuity in the derivative. */
                /* So we need to guard against searching near such times. */

                double rlon_lo, rlon_hi, adjust_days, syn;
                if (rlon >= -s1 && rlon < +s1)
                {
                    /* Seek to the window [+s1, +s2]. */
                    adjust_days = 0.0;
                    /* Search forward for the time t1 when rel lon = +s1. */
                    rlon_lo = +s1;
                    /* Search forward for the time t2 when rel lon = +s2. */
                    rlon_hi = +s2;
                }
                else if (rlon >= +s2 || rlon < -s2)
                {
                    /* Seek to the next search window at [-s2, -s1]. */
                    adjust_days = 0.0;
                    /* Search forward for the time t1 when rel lon = -s2. */
                    rlon_lo = -s2;
                    /* Search forward for the time t2 when rel lon = -s1. */
                    rlon_hi = -s1;
                }
                else if (rlon >= 0)
                {
                    /* rlon must be in the middle of the window [+s1, +s2]. */
                    /* Search BACKWARD for the time t1 when rel lon = +s1. */
                    syn = SynodicPeriod(body);
                    adjust_days = -syn / 4;
                    rlon_lo = +s1;
                    /* Search forward from t1 to find t2 such that rel lon = +s2. */
                    rlon_hi = +s2;
                }
                else
                {
                    /* rlon must be in the middle of the window [-s2, -s1]. */
                    /* Search BACKWARD for the time t1 when rel lon = -s2. */
                    syn = SynodicPeriod(body);
                    adjust_days = -syn / 4;
                    rlon_lo = -s2;
                    /* Search forward from t1 to find t2 such that rel lon = -s1. */
                    rlon_hi = -s1;
                }
                AstroTime t_start = startTime.AddDays(adjust_days);
                AstroTime t1 = SearchRelativeLongitude(body, rlon_lo, t_start);
                AstroTime t2 = SearchRelativeLongitude(body, rlon_hi, t1);

                /* Now we have a time range [t1,t2] that brackets a maximum magnitude event. */
                /* Confirm the bracketing. */
                double m1 = mag_slope.Eval(t1);
                if (m1 >= 0.0)
                    throw new Exception("Internal error: m1 >= 0");    /* should never happen! */

                double m2 = mag_slope.Eval(t2);
                if (m2 <= 0.0)
                    throw new Exception("Internal error: m2 <= 0");    /* should never happen! */

                /* Use the generic search algorithm to home in on where the slope crosses from negative to positive. */
                AstroTime tx = Search(mag_slope, t1, t2, 10.0);
                if (tx == null)
                    throw new Exception("Failed to find magnitude slope transition.");

                if (tx.tt >= startTime.tt)
                    return Illumination(body, tx);

                /* This event is in the past (earlier than startTime). */
                /* We need to search forward from t2 to find the next possible window. */
                /* We never need to search more than twice. */
                startTime = t2.AddDays(1.0);
            }
            // This should never happen. If it does, please report as a bug in Astronomy Engine.
            throw new Exception("Peak magnitude search failed.");
        }

        /// <summary>Calculates the inverse of a rotation matrix.</summary>
        /// <remarks>
        /// Given a rotation matrix that performs some coordinate transform,
        /// this function returns the matrix that reverses that trasnform.
        /// </remarks>
        /// <param name="rotation">The rotation matrix to be inverted.</param>
        /// <returns>A rotation matrix that performs the opposite transformation.</returns>
        public static RotationMatrix InverseRotation(RotationMatrix rotation)
        {
            var inverse = new RotationMatrix(new double[3,3]);

            inverse.rot[0, 0] = rotation.rot[0, 0];
            inverse.rot[0, 1] = rotation.rot[1, 0];
            inverse.rot[0, 2] = rotation.rot[2, 0];
            inverse.rot[1, 0] = rotation.rot[0, 1];
            inverse.rot[1, 1] = rotation.rot[1, 1];
            inverse.rot[1, 2] = rotation.rot[2, 1];
            inverse.rot[2, 0] = rotation.rot[0, 2];
            inverse.rot[2, 1] = rotation.rot[1, 2];
            inverse.rot[2, 2] = rotation.rot[2, 2];

            return inverse;
        }

        /// <summary>Creates a rotation based on applying one rotation followed by another.</summary>
        /// <remarks>
        /// Given two rotation matrices, returns a combined rotation matrix that is
        /// equivalent to rotating based on the first matrix, followed by the second.
        /// </remarks>
        /// <param name="a">The first rotation to apply.</param>
        /// <param name="b">The second rotation to apply.</param>
        /// <returns>The combined rotation matrix.</returns>
        public static RotationMatrix CombineRotation(RotationMatrix a, RotationMatrix b)
        {
            var rot = new double[3,3];

            // Use matrix multiplication: c = b*a.
            // We put 'b' on the left and 'a' on the right because,
            // just like when you use a matrix M to rotate a vector V,
            // you put the M on the left in the product M*V.
            // We can think of this as 'b' rotating all the 3 column vectors in 'a'.

            rot[0, 0] = b.rot[0, 0]*a.rot[0, 0] + b.rot[1, 0]*a.rot[0, 1] + b.rot[2, 0]*a.rot[0, 2];
            rot[1, 0] = b.rot[0, 0]*a.rot[1, 0] + b.rot[1, 0]*a.rot[1, 1] + b.rot[2, 0]*a.rot[1, 2];
            rot[2, 0] = b.rot[0, 0]*a.rot[2, 0] + b.rot[1, 0]*a.rot[2, 1] + b.rot[2, 0]*a.rot[2, 2];
            rot[0, 1] = b.rot[0, 1]*a.rot[0, 0] + b.rot[1, 1]*a.rot[0, 1] + b.rot[2, 1]*a.rot[0, 2];
            rot[1, 1] = b.rot[0, 1]*a.rot[1, 0] + b.rot[1, 1]*a.rot[1, 1] + b.rot[2, 1]*a.rot[1, 2];
            rot[2, 1] = b.rot[0, 1]*a.rot[2, 0] + b.rot[1, 1]*a.rot[2, 1] + b.rot[2, 1]*a.rot[2, 2];
            rot[0, 2] = b.rot[0, 2]*a.rot[0, 0] + b.rot[1, 2]*a.rot[0, 1] + b.rot[2, 2]*a.rot[0, 2];
            rot[1, 2] = b.rot[0, 2]*a.rot[1, 0] + b.rot[1, 2]*a.rot[1, 1] + b.rot[2, 2]*a.rot[1, 2];
            rot[2, 2] = b.rot[0, 2]*a.rot[2, 0] + b.rot[1, 2]*a.rot[2, 1] + b.rot[2, 2]*a.rot[2, 2];

            return new RotationMatrix(rot);
        }

        /// <summary>Applies a rotation to a vector, yielding a rotated vector.</summary>
        /// <remarks>
        /// This function transforms a vector in one orientation to a vector
        /// in another orientation.
        /// </remarks>
        /// <param name="rotation">A rotation matrix that specifies how the orientation of the vector is to be changed.</param>
        /// <param name="vector">The vector whose orientation is to be changed.</param>
        /// <returns>A vector in the orientation specified by `rotation`.</returns>
        public static AstroVector RotateVector(RotationMatrix rotation, AstroVector vector)
        {
            return new AstroVector(
                rotation.rot[0, 0]*vector.x + rotation.rot[1, 0]*vector.y + rotation.rot[2, 0]*vector.z,
                rotation.rot[0, 1]*vector.x + rotation.rot[1, 1]*vector.y + rotation.rot[2, 1]*vector.z,
                rotation.rot[0, 2]*vector.x + rotation.rot[1, 2]*vector.y + rotation.rot[2, 2]*vector.z,
                vector.t
            );
        }

        /// <summary>Converts spherical coordinates to Cartesian coordinates.</summary>
        /// <remarks>
        /// Given spherical coordinates and a time at which they are valid,
        /// returns a vector of Cartesian coordinates. The returned value
        /// includes the time, as required by the type #AstroVector.
        /// </remarks>
        /// <param name="sphere">Spherical coordinates to be converted.</param>
        /// <param name="time">The time that should be included in the return value.</param>
        /// <returns>The vector form of the supplied spherical coordinates.</returns>
        public static AstroVector VectorFromSphere(Spherical sphere, AstroTime time)
        {
            double radlat = sphere.lat * DEG2RAD;
            double radlon = sphere.lon * DEG2RAD;
            double rcoslat = sphere.dist * Math.Cos(radlat);
            return new AstroVector(
                rcoslat * Math.Cos(radlon),
                rcoslat * Math.Sin(radlon),
                sphere.dist * Math.Sin(radlat),
                time
            );
        }

        /// <summary>Converts Cartesian coordinates to spherical coordinates.</summary>
        /// <remarks>
        /// Given a Cartesian vector, returns latitude, longitude, and distance.
        /// </remarks>
        /// <param name="vector">Cartesian vector to be converted to spherical coordinates.</param>
        /// <returns>Spherical coordinates that are equivalent to the given vector.</returns>
        public static Spherical SphereFromVector(AstroVector vector)
        {
            double xyproj = vector.x*vector.x + vector.y*vector.y;
            double dist = Math.Sqrt(xyproj + vector.z*vector.z);
            double lat, lon;
            if (xyproj == 0.0)
            {
                if (vector.z == 0.0)
                {
                    /* Indeterminate coordinates; pos vector has zero length. */
                    throw new ArgumentException("Cannot convert zero-length vector to spherical coordinates.");
                }

                lon = 0.0;
                lat = (vector.z < 0.0) ? -90.0 : +90.0;
            }
            else
            {
                lon = RAD2DEG * Math.Atan2(vector.y, vector.x);
                if (lon < 0.0)
                    lon += 360.0;

                lat = RAD2DEG * Math.Atan2(vector.z, Math.Sqrt(xyproj));
            }

            return new Spherical(lat, lon, dist);
        }


        /// <summary>Given angular equatorial coordinates in `equ`, calculates equatorial vector.</summary>
        /// <param name="equ">Angular equatorial coordinates to be converted to a vector.</param>
        /// <param name="time">
        /// The date and time of the observation. This is needed because the returned
        /// vector requires a valid time value when passed to certain other functions.
        /// </param>
        /// <returns>A vector in the equatorial system.</returns>
        public static AstroVector VectorFromEquator(Equatorial equ, AstroTime time)
        {
            var sphere = new Spherical(equ.dec, 15.0 * equ.ra, equ.dist);
            return VectorFromSphere(sphere, time);
        }


        /// <summary>Given an equatorial vector, calculates equatorial angular coordinates.</summary>
        /// <param name="vector">A vector in an equatorial coordinate system.</param>
        /// <returns>Angular coordinates expressed in the same equatorial system as `vector`.</returns>
        public static Equatorial EquatorFromVector(AstroVector vector)
        {
            Spherical sphere = SphereFromVector(vector);
            return new Equatorial(sphere.lon / 15.0, sphere.lat, sphere.dist);
        }


        private static double ToggleAzimuthDirection(double az)
        {
            az = 360.0 - az;
            if (az >= 360.0)
                az -= 360.0;
            else if (az < 0.0)
                az += 360.0;
            return az;
        }


        /// <summary>
        /// Converts Cartesian coordinates to horizontal coordinates.
        /// </summary>
        /// <remarks>
        /// Given a horizontal Cartesian vector, returns horizontal azimuth and altitude.
        ///
        /// *IMPORTANT:* This function differs from #Astronomy.SphereFromVector in two ways:
        /// - `Astronomy.SphereFromVector` returns a `lon` value that represents azimuth defined counterclockwise
        ///   from north (e.g., west = +90), but this function represents a clockwise rotation
        ///   (e.g., east = +90). The difference is because `Astronomy.SphereFromVector` is intended
        ///   to preserve the vector "right-hand rule", while this function defines azimuth in a more
        ///   traditional way as used in navigation and cartography.
        /// - This function optionally corrects for atmospheric refraction, while `Astronomy.SphereFromVector`
        ///   does not.
        ///
        /// The returned structure contains the azimuth in `lon`.
        /// It is measured in degrees clockwise from north: east = +90 degrees, west = +270 degrees.
        ///
        /// The altitude is stored in `lat`.
        ///
        /// The distance to the observed object is stored in `dist`,
        /// and is expressed in astronomical units (AU).
        /// </remarks>
        /// <param name="vector">Cartesian vector to be converted to horizontal coordinates.</param>
        /// <param name="refraction">
        /// `Refraction.Normal`: correct altitude for atmospheric refraction (recommended).
        /// `Refraction.None`: no atmospheric refraction correction is performed.
        /// `Refraction.JplHor`: for JPL Horizons compatibility testing only; not recommended for normal use.
        /// </param>
        /// <returns>
        /// Horizontal spherical coordinates as described above.
        /// </returns>
        public static Spherical HorizonFromVector(AstroVector vector, Refraction refraction)
        {
            Spherical sphere = SphereFromVector(vector);
            return new Spherical(
                sphere.lat + RefractionAngle(refraction, sphere.lat),
                ToggleAzimuthDirection(sphere.lon),
                sphere.dist
            );
        }


        /// <summary>
        /// Given apparent angular horizontal coordinates in `sphere`, calculate horizontal vector.
        /// </summary>
        /// <param name="sphere">
        /// A structure that contains apparent horizontal coordinates:
        /// `lat` holds the refracted azimuth angle,
        /// `lon` holds the azimuth in degrees clockwise from north,
        /// and `dist` holds the distance from the observer to the object in AU.
        /// </param>
        /// <param name="time">
        /// The date and time of the observation. This is needed because the returned
        /// #AstroVector requires a valid time value when passed to certain other functions.
        /// </param>
        /// <param name="refraction">
        /// The refraction option used to model atmospheric lensing. See #Astronomy.RefractionAngle.
        /// This specifies how refraction is to be removed from the altitude stored in `sphere.lat`.
        /// </param>
        /// <returns>
        /// A vector in the horizontal system: `x` = north, `y` = west, and `z` = zenith (up).
        /// </returns>
        public static AstroVector VectorFromHorizon(Spherical sphere, AstroTime time, Refraction refraction)
        {
            return VectorFromSphere(
                new Spherical(
                    sphere.lat + InverseRefractionAngle(refraction, sphere.lat),
                    ToggleAzimuthDirection(sphere.lon),
                    sphere.dist
                ),
                time
            );
        }


        /// <summary>
        /// Calculates the amount of "lift" to an altitude angle caused by atmospheric refraction.
        /// </summary>
        /// <remarks>
        /// Given an altitude angle and a refraction option, calculates
        /// the amount of "lift" caused by atmospheric refraction.
        /// This is the number of degrees higher in the sky an object appears
        /// due to the lensing of the Earth's atmosphere.
        /// </remarks>
        /// <param name="refraction">
        /// The option selecting which refraction correction to use.
        /// If `Refraction.Normal`, uses a well-behaved refraction model that works well for
        /// all valid values (-90 to +90) of `altitude`.
        /// If `Refraction.JplHor`, this function returns a compatible value with the JPL Horizons tool.
        /// If any other value (including `Refraction.None`), this function returns 0.
        /// </param>
        /// <param name="altitude">
        /// An altitude angle in a horizontal coordinate system. Must be a value between -90 and +90.
        /// </param>
        /// <returns>
        /// The angular adjustment in degrees to be added to the altitude angle to correct for atmospheric lensing.
        /// </returns>
        public static double RefractionAngle(Refraction refraction, double altitude)
        {
            if (altitude < -90.0 || altitude > +90.0)
                return 0.0;     /* no attempt to correct an invalid altitude */

            double refr;
            if (refraction == Refraction.Normal || refraction == Refraction.JplHor)
            {
                // http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
                // JPL Horizons says it uses refraction algorithm from
                // Meeus "Astronomical Algorithms", 1991, p. 101-102.
                // I found the following Go implementation:
                // https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
                // This is a translation from the function "Saemundsson" there.
                // I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
                // This is important because the 'refr' formula below goes crazy near hd = -5.11.
                double hd = altitude;
                if (hd < -1.0)
                    hd = -1.0;

                refr = (1.02 / Math.Tan((hd+10.3/(hd+5.11))*DEG2RAD)) / 60.0;

                if (refraction == Refraction.Normal && altitude < -1.0)
                {
                    // In "normal" mode we gradually reduce refraction toward the nadir
                    // so that we never get an altitude angle less than -90 degrees.
                    // When horizon angle is -1 degrees, the factor is exactly 1.
                    // As altitude approaches -90 (the nadir), the fraction approaches 0 linearly.
                    refr *= (altitude + 90.0) / 89.0;
                }
            }
            else
            {
                /* No refraction, or the refraction option is invalid. */
                refr = 0.0;
            }

            return refr;
        }


        /// <summary>
        /// Calculates the inverse of an atmospheric refraction angle.
        /// </summary>
        /// <remarks>
        /// Given an observed altitude angle that includes atmospheric refraction,
        /// calculate the negative angular correction to obtain the unrefracted
        /// altitude. This is useful for cases where observed horizontal
        /// coordinates are to be converted to another orientation system,
        /// but refraction first must be removed from the observed position.
        /// </remarks>
        /// <param name="refraction">
        /// The option selecting which refraction correction to use.
        /// See #Astronomy.RefractionAngle.
        /// </param>
        /// <param name="bent_altitude">
        /// The apparent altitude that includes atmospheric refraction.
        /// </param>
        /// <returns>
        /// The angular adjustment in degrees to be added to the
        /// altitude angle to correct for atmospheric lensing.
        /// This will be less than or equal to zero.
        /// </returns>
        public static double InverseRefractionAngle(Refraction refraction, double bent_altitude)
        {
            if (bent_altitude < -90.0 || bent_altitude > +90.0)
                return 0.0;     /* no attempt to correct an invalid altitude */

            /* Find the pre-adjusted altitude whose refraction correction leads to 'altitude'. */
            double altitude = bent_altitude - RefractionAngle(refraction, bent_altitude);
            for(;;)
            {
                /* See how close we got. */
                double diff = (altitude + RefractionAngle(refraction, altitude)) - bent_altitude;
                if (Math.Abs(diff) < 1.0e-14)
                    return altitude - bent_altitude;

                altitude -= diff;
            }
        }


        /// <summary>Calculates a rotation matrix from equatorial J2000 (EQJ) to ecliptic J2000 (ECL).</summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQJ = equatorial system, using equator at J2000 epoch.
        /// Target: ECL = ecliptic system, using equator at J2000 epoch.
        /// </remarks>
        /// <returns>A rotation matrix that converts EQJ to ECL.</returns>
        public static RotationMatrix Rotation_EQJ_ECL()
        {
            /* ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians. */
            const double c = 0.9174821430670688;    /* cos(ob) */
            const double s = 0.3977769691083922;    /* sin(ob) */
            var r = new RotationMatrix(new double[3,3]);

            r.rot[0, 0] = 1.0;  r.rot[1, 0] = 0.0;  r.rot[2, 0] = 0.0;
            r.rot[0, 1] = 0.0;  r.rot[1, 1] = +c;   r.rot[2, 1] = +s;
            r.rot[0, 2] = 0.0;  r.rot[1, 2] = -s;   r.rot[2, 2] = +c;

            return r;
        }


        /// <summary>Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial J2000 (EQJ).</summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: ECL = ecliptic system, using equator at J2000 epoch.
        /// Target: EQJ = equatorial system, using equator at J2000 epoch.
        /// </remarks>
        /// <returns>A rotation matrix that converts ECL to EQJ.</returns>
        public static RotationMatrix Rotation_ECL_EQJ()
        {
            /* ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians. */
            const double c = 0.9174821430670688;    /* cos(ob) */
            const double s = 0.3977769691083922;    /* sin(ob) */
            var r = new RotationMatrix(new double[3,3]);

            r.rot[0, 0] = 1.0;  r.rot[1, 0] = 0.0;  r.rot[2, 0] = 0.0;
            r.rot[0, 1] = 0.0;  r.rot[1, 1] = +c;   r.rot[2, 1] = -s;
            r.rot[0, 2] = 0.0;  r.rot[1, 2] = +s;   r.rot[2, 2] = +c;

            return r;
        }


        /// <summary>
        /// Calculates a rotation matrix from equatorial J2000 (EQJ) to equatorial of-date (EQD).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQJ = equatorial system, using equator at J2000 epoch.
        /// Target: EQD = equatorial system, using equator of the specified date/time.
        /// </remarks>
        /// <param name="time">
        /// The date and time at which the Earth's equator defines the target orientation.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts EQJ to EQD at `time`.
        /// </returns>
        public static RotationMatrix Rotation_EQJ_EQD(AstroTime time)
        {
            RotationMatrix prec = precession_rot(0.0, time.tt);
            RotationMatrix nut = nutation_rot(time, 0);
            return CombineRotation(prec, nut);
        }


        /// <summary>
        /// Calculates a rotation matrix from equatorial of-date (EQD) to equatorial J2000 (EQJ).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQD = equatorial system, using equator of the specified date/time.
        /// Target: EQJ = equatorial system, using equator at J2000 epoch.
        /// </remarks>
        /// <param name="time">
        /// The date and time at which the Earth's equator defines the source orientation.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts EQD at `time` to EQJ.
        /// </returns>
        public static RotationMatrix Rotation_EQD_EQJ(AstroTime time)
        {
            RotationMatrix nut = nutation_rot(time, 1);
            RotationMatrix prec = precession_rot(time.tt, 0.0);
            return CombineRotation(nut, prec);
        }


        /// <summary>
        /// Calculates a rotation matrix from equatorial of-date (EQD) to horizontal (HOR).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQD = equatorial system, using equator of the specified date/time.
        /// Target: HOR = horizontal system.
        /// </remarks>
        /// <param name="time">
        /// The date and time at which the Earth's equator applies.
        /// </param>
        /// <param name="observer">
        /// A location near the Earth's mean sea level that defines the observer's horizon.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts EQD to HOR at `time` and for `observer`.
        /// The components of the horizontal vector are:
        /// x = north, y = west, z = zenith (straight up from the observer).
        /// These components are chosen so that the "right-hand rule" works for the vector
        /// and so that north represents the direction where azimuth = 0.
        /// </returns>
        public static RotationMatrix Rotation_EQD_HOR(AstroTime time, Observer observer)
        {
            double sinlat = Math.Sin(observer.latitude * DEG2RAD);
            double coslat = Math.Cos(observer.latitude * DEG2RAD);
            double sinlon = Math.Sin(observer.longitude * DEG2RAD);
            double coslon = Math.Cos(observer.longitude * DEG2RAD);

            var uze = new AstroVector(coslat * coslon, coslat * sinlon, sinlat, null);
            var une = new AstroVector(-sinlat * coslon, -sinlat * sinlon, coslat, null);
            var uwe = new AstroVector(sinlon, -coslon, 0.0, null);

            double spin_angle = -15.0 * sidereal_time(time);
            AstroVector uz = spin(spin_angle, uze);
            AstroVector un = spin(spin_angle, une);
            AstroVector uw = spin(spin_angle, uwe);

            var rot = new double[3,3];
            rot[0, 0] = un.x; rot[1, 0] = un.y; rot[2, 0] = un.z;
            rot[0, 1] = uw.x; rot[1, 1] = uw.y; rot[2, 1] = uw.z;
            rot[0, 2] = uz.x; rot[1, 2] = uz.y; rot[2, 2] = uz.z;

            return new RotationMatrix(rot);
        }


        /// <summary>
        /// Calculates a rotation matrix from horizontal (HOR) to equatorial of-date (EQD).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: HOR = horizontal system (x=North, y=West, z=Zenith).
        /// Target: EQD = equatorial system, using equator of the specified date/time.
        /// </remarks>
        /// <param name="time">
        /// The date and time at which the Earth's equator applies.
        /// </param>
        /// <param name="observer">
        /// A location near the Earth's mean sea level that defines the observer's horizon.
        /// </param>
        /// <returns>
        ///  A rotation matrix that converts HOR to EQD at `time` and for `observer`.
        /// </returns>
        public static RotationMatrix Rotation_HOR_EQD(AstroTime time, Observer observer)
        {
            RotationMatrix rot = Rotation_EQD_HOR(time, observer);
            return InverseRotation(rot);
        }


        /// <summary>
        /// Calculates a rotation matrix from horizontal (HOR) to J2000 equatorial (EQJ).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: HOR = horizontal system (x=North, y=West, z=Zenith).
        /// Target: EQJ = equatorial system, using equator at the J2000 epoch.
        /// </remarks>
        /// <param name="time">
        /// The date and time of the observation.
        /// </param>
        /// <param name="observer">
        /// A location near the Earth's mean sea level that defines the observer's horizon.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts HOR to EQD at `time` and for `observer`.
        /// </returns>
        public static RotationMatrix Rotation_HOR_EQJ(AstroTime time, Observer observer)
        {
            RotationMatrix hor_eqd = Rotation_HOR_EQD(time, observer);
            RotationMatrix eqd_eqj = Rotation_EQD_EQJ(time);
            return CombineRotation(hor_eqd, eqd_eqj);
        }


        /// <summary>
        /// Calculates a rotation matrix from equatorial J2000 (EQJ) to horizontal (HOR).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQJ = equatorial system, using the equator at the J2000 epoch.
        /// Target: HOR = horizontal system.
        /// </remarks>
        /// <param name="time">
        /// The date and time of the desired horizontal orientation.
        /// </param>
        /// <param name="observer">
        /// A location near the Earth's mean sea level that defines the observer's horizon.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts EQJ to HOR at `time` and for `observer`.
        /// The components of the horizontal vector are:
        /// x = north, y = west, z = zenith (straight up from the observer).
        /// These components are chosen so that the "right-hand rule" works for the vector
        /// and so that north represents the direction where azimuth = 0.
        /// </returns>
        public static RotationMatrix Rotation_EQJ_HOR(AstroTime time, Observer observer)
        {
            RotationMatrix rot = Rotation_HOR_EQJ(time, observer);
            return InverseRotation(rot);
        }


        /// <summary>
        /// Calculates a rotation matrix from equatorial of-date (EQD) to ecliptic J2000 (ECL).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQD = equatorial system, using equator of date.
        /// Target: ECL = ecliptic system, using equator at J2000 epoch.
        /// </remarks>
        /// <param name="time">
        /// The date and time of the source equator.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts EQD to ECL.
        /// </returns>
        public static RotationMatrix Rotation_EQD_ECL(AstroTime time)
        {
            RotationMatrix eqd_eqj = Rotation_EQD_EQJ(time);
            RotationMatrix eqj_ecl = Rotation_EQJ_ECL();
            return CombineRotation(eqd_eqj, eqj_ecl);
        }


        /// <summary>
        /// Calculates a rotation matrix from ecliptic J2000 (ECL) to equatorial of-date (EQD).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: ECL = ecliptic system, using equator at J2000 epoch.
        /// Target: EQD = equatorial system, using equator of date.
        /// </remarks>
        /// <param name="time">
        /// The date and time of the desired equator.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts ECL to EQD.
        /// </returns>
        public static RotationMatrix Rotation_ECL_EQD(AstroTime time)
        {
            RotationMatrix rot = Rotation_EQD_ECL(time);
            return InverseRotation(rot);
        }


        /// <summary>
        /// Calculates a rotation matrix from ecliptic J2000 (ECL) to horizontal (HOR).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: ECL = ecliptic system, using equator at J2000 epoch.
        /// Target: HOR = horizontal system.
        /// </remarks>
        /// <param name="time">
        /// The date and time of the desired horizontal orientation.
        /// </param>
        /// <param name="observer">
        /// A location near the Earth's mean sea level that defines the observer's horizon.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts ECL to HOR at `time` and for `observer`.
        /// The components of the horizontal vector are:
        /// x = north, y = west, z = zenith (straight up from the observer).
        /// These components are chosen so that the "right-hand rule" works for the vector
        /// and so that north represents the direction where azimuth = 0.
        /// </returns>
        public static RotationMatrix Rotation_ECL_HOR(AstroTime time, Observer observer)
        {
            RotationMatrix ecl_eqd = Rotation_ECL_EQD(time);
            RotationMatrix eqd_hor = Rotation_EQD_HOR(time, observer);
            return CombineRotation(ecl_eqd, eqd_hor);
        }

        /// <summary>
        /// Calculates a rotation matrix from horizontal (HOR) to ecliptic J2000 (ECL).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: HOR = horizontal system.
        /// Target: ECL = ecliptic system, using equator at J2000 epoch.
        /// </remarks>
        /// <param name="time">
        /// The date and time of the horizontal observation.
        /// </param>
        /// <param name="observer">
        /// The location of the horizontal observer.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts HOR to ECL.
        /// </returns>
        public static RotationMatrix Rotation_HOR_ECL(AstroTime time, Observer observer)
        {
            RotationMatrix rot = Rotation_ECL_HOR(time, observer);
            return InverseRotation(rot);
        }

        private struct constel_info_t
        {
            public readonly string symbol;
            public readonly string name;

            public constel_info_t(string symbol, string name)
            {
                this.symbol = symbol;
                this.name = name;
            }
        }

        private struct constel_boundary_t
        {
            public readonly int index;
            public readonly double ra_lo;
            public readonly double ra_hi;
            public readonly double dec_lo;

            public constel_boundary_t(int index, double ra_lo, double ra_hi, double dec_lo)
            {
                this.index = index;
                this.ra_lo = ra_lo;
                this.ra_hi = ra_hi;
                this.dec_lo = dec_lo;
            }
        }

        private static readonly object ConstelLock = new object();
        private static RotationMatrix ConstelRot;
        private static AstroTime Epoch2000;

        /// <summary>
        /// Determines the constellation that contains the given point in the sky.
        /// </summary>
        /// <remarks>
        /// Given J2000 equatorial (EQJ) coordinates of a point in the sky, determines the
        /// constellation that contains that point.
        /// </remarks>
        /// <param name="ra">
        /// The right ascension (RA) of a point in the sky, using the J2000 equatorial system.
        /// </param>
        /// <param name="dec">
        /// The declination (DEC) of a point in the sky, using the J2000 equatorial system.
        /// </param>
        /// <returns>
        /// A structure that contains the 3-letter abbreviation and full name
        /// of the constellation that contains the given (ra,dec), along with
        /// the converted B1875 (ra,dec) for that point.
        /// </returns>
        public static ConstellationInfo Constellation(double ra, double dec)
        {
            if (dec < -90.0 || dec > +90.0)
                throw new ArgumentException("Invalid declination angle. Must be -90..+90.");

            // Allow right ascension to "wrap around". Clamp to [0, 24) sidereal hours.
            ra %= 24.0;
            if (ra < 0.0)
                ra += 24.0;

            lock (ConstelLock)
            {
                if (ConstelRot.rot == null)
                {
                    // Lazy-initialize the rotation matrix for converting J2000 to B1875.
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
                    var time = new AstroTime(-45655.74141261017);
                    ConstelRot = Rotation_EQJ_EQD(time);
                    Epoch2000 = new AstroTime(0.0);
                }
            }

            // Convert coordinates from J2000 to B1875.
            var equ2000 = new Equatorial(ra, dec, 1.0);
            AstroVector vec2000 = VectorFromEquator(equ2000, Epoch2000);
            AstroVector vec1875 = RotateVector(ConstelRot, vec2000);
            Equatorial equ1875 = EquatorFromVector(vec1875);

            // Search for the constellation using the B1875 coordinates.
            foreach (constel_boundary_t b in ConstelBounds)
                if ((b.dec_lo <= equ1875.dec) && (b.ra_hi > equ1875.ra) && (b.ra_lo <= equ1875.ra))
                    return new ConstellationInfo(ConstelNames[b.index].symbol, ConstelNames[b.index].name, equ1875.ra, equ1875.dec);

            // This should never happen!
            throw new Exception("Unable to find constellation for given coordinates.");
        }

        private static readonly constel_info_t[] ConstelNames = new constel_info_t[]
        {
            new constel_info_t("And", "Andromeda"           )  //  0
        ,   new constel_info_t("Ant", "Antila"              )  //  1
        ,   new constel_info_t("Aps", "Apus"                )  //  2
        ,   new constel_info_t("Aql", "Aquila"              )  //  3
        ,   new constel_info_t("Aqr", "Aquarius"            )  //  4
        ,   new constel_info_t("Ara", "Ara"                 )  //  5
        ,   new constel_info_t("Ari", "Aries"               )  //  6
        ,   new constel_info_t("Aur", "Auriga"              )  //  7
        ,   new constel_info_t("Boo", "Bootes"              )  //  8
        ,   new constel_info_t("Cae", "Caelum"              )  //  9
        ,   new constel_info_t("Cam", "Camelopardis"        )  // 10
        ,   new constel_info_t("Cap", "Capricornus"         )  // 11
        ,   new constel_info_t("Car", "Carina"              )  // 12
        ,   new constel_info_t("Cas", "Cassiopeia"          )  // 13
        ,   new constel_info_t("Cen", "Centaurus"           )  // 14
        ,   new constel_info_t("Cep", "Cepheus"             )  // 15
        ,   new constel_info_t("Cet", "Cetus"               )  // 16
        ,   new constel_info_t("Cha", "Chamaeleon"          )  // 17
        ,   new constel_info_t("Cir", "Circinus"            )  // 18
        ,   new constel_info_t("CMa", "Canis Major"         )  // 19
        ,   new constel_info_t("CMi", "Canis Minor"         )  // 20
        ,   new constel_info_t("Cnc", "Cancer"              )  // 21
        ,   new constel_info_t("Col", "Columba"             )  // 22
        ,   new constel_info_t("Com", "Coma Berenices"      )  // 23
        ,   new constel_info_t("CrA", "Corona Australis"    )  // 24
        ,   new constel_info_t("CrB", "Corona Borealis"     )  // 25
        ,   new constel_info_t("Crt", "Crater"              )  // 26
        ,   new constel_info_t("Cru", "Crux"                )  // 27
        ,   new constel_info_t("Crv", "Corvus"              )  // 28
        ,   new constel_info_t("CVn", "Canes Venatici"      )  // 29
        ,   new constel_info_t("Cyg", "Cygnus"              )  // 30
        ,   new constel_info_t("Del", "Delphinus"           )  // 31
        ,   new constel_info_t("Dor", "Dorado"              )  // 32
        ,   new constel_info_t("Dra", "Draco"               )  // 33
        ,   new constel_info_t("Equ", "Equuleus"            )  // 34
        ,   new constel_info_t("Eri", "Eridanus"            )  // 35
        ,   new constel_info_t("For", "Fornax"              )  // 36
        ,   new constel_info_t("Gem", "Gemini"              )  // 37
        ,   new constel_info_t("Gru", "Grus"                )  // 38
        ,   new constel_info_t("Her", "Hercules"            )  // 39
        ,   new constel_info_t("Hor", "Horologium"          )  // 40
        ,   new constel_info_t("Hya", "Hydra"               )  // 41
        ,   new constel_info_t("Hyi", "Hydrus"              )  // 42
        ,   new constel_info_t("Ind", "Indus"               )  // 43
        ,   new constel_info_t("Lac", "Lacerta"             )  // 44
        ,   new constel_info_t("Leo", "Leo"                 )  // 45
        ,   new constel_info_t("Lep", "Lepus"               )  // 46
        ,   new constel_info_t("Lib", "Libra"               )  // 47
        ,   new constel_info_t("LMi", "Leo Minor"           )  // 48
        ,   new constel_info_t("Lup", "Lupus"               )  // 49
        ,   new constel_info_t("Lyn", "Lynx"                )  // 50
        ,   new constel_info_t("Lyr", "Lyra"                )  // 51
        ,   new constel_info_t("Men", "Mensa"               )  // 52
        ,   new constel_info_t("Mic", "Microscopium"        )  // 53
        ,   new constel_info_t("Mon", "Monoceros"           )  // 54
        ,   new constel_info_t("Mus", "Musca"               )  // 55
        ,   new constel_info_t("Nor", "Norma"               )  // 56
        ,   new constel_info_t("Oct", "Octans"              )  // 57
        ,   new constel_info_t("Oph", "Ophiuchus"           )  // 58
        ,   new constel_info_t("Ori", "Orion"               )  // 59
        ,   new constel_info_t("Pav", "Pavo"                )  // 60
        ,   new constel_info_t("Peg", "Pegasus"             )  // 61
        ,   new constel_info_t("Per", "Perseus"             )  // 62
        ,   new constel_info_t("Phe", "Phoenix"             )  // 63
        ,   new constel_info_t("Pic", "Pictor"              )  // 64
        ,   new constel_info_t("PsA", "Pisces Austrinus"    )  // 65
        ,   new constel_info_t("Psc", "Pisces"              )  // 66
        ,   new constel_info_t("Pup", "Puppis"              )  // 67
        ,   new constel_info_t("Pyx", "Pyxis"               )  // 68
        ,   new constel_info_t("Ret", "Reticulum"           )  // 69
        ,   new constel_info_t("Scl", "Sculptor"            )  // 70
        ,   new constel_info_t("Sco", "Scorpius"            )  // 71
        ,   new constel_info_t("Sct", "Scutum"              )  // 72
        ,   new constel_info_t("Ser", "Serpens"             )  // 73
        ,   new constel_info_t("Sex", "Sextans"             )  // 74
        ,   new constel_info_t("Sge", "Sagitta"             )  // 75
        ,   new constel_info_t("Sgr", "Sagittarius"         )  // 76
        ,   new constel_info_t("Tau", "Taurus"              )  // 77
        ,   new constel_info_t("Tel", "Telescopium"         )  // 78
        ,   new constel_info_t("TrA", "Triangulum Australe" )  // 79
        ,   new constel_info_t("Tri", "Triangulum"          )  // 80
        ,   new constel_info_t("Tuc", "Tucana"              )  // 81
        ,   new constel_info_t("UMa", "Ursa Major"          )  // 82
        ,   new constel_info_t("UMi", "Ursa Minor"          )  // 83
        ,   new constel_info_t("Vel", "Vela"                )  // 84
        ,   new constel_info_t("Vir", "Virgo"               )  // 85
        ,   new constel_info_t("Vol", "Volans"              )  // 86
        ,   new constel_info_t("Vul", "Vulpecula"           )  // 87
        };

        private static readonly constel_boundary_t[] ConstelBounds = new constel_boundary_t[]
        {
            new constel_boundary_t(83,  0.00000000000000, 24.00000000000000, 88.00000000000000)    // UMi
        ,   new constel_boundary_t(83,  8.00000000000000, 14.50000000000000, 86.50000000000000)    // UMi
        ,   new constel_boundary_t(83, 21.00000000000000, 23.00000000000000, 86.16666666666667)    // UMi
        ,   new constel_boundary_t(83, 18.00000000000000, 21.00000000000000, 86.00000000000000)    // UMi
        ,   new constel_boundary_t(15,  0.00000000000000,  8.00000000000000, 85.00000000000000)    // Cep
        ,   new constel_boundary_t(10,  9.16666666666667, 10.66666666666667, 82.00000000000000)    // Cam
        ,   new constel_boundary_t(15,  0.00000000000000,  5.00000000000000, 80.00000000000000)    // Cep
        ,   new constel_boundary_t(10, 10.66666666666667, 14.50000000000000, 80.00000000000000)    // Cam
        ,   new constel_boundary_t(83, 17.50000000000000, 18.00000000000000, 80.00000000000000)    // UMi
        ,   new constel_boundary_t(33, 20.16666666666667, 21.00000000000000, 80.00000000000000)    // Dra
        ,   new constel_boundary_t(15,  0.00000000000000,  3.50833333333333, 77.00000000000000)    // Cep
        ,   new constel_boundary_t(10, 11.50000000000000, 13.58333333333333, 77.00000000000000)    // Cam
        ,   new constel_boundary_t(83, 16.53333333333333, 17.50000000000000, 75.00000000000000)    // UMi
        ,   new constel_boundary_t(15, 20.16666666666667, 20.66666666666667, 75.00000000000000)    // Cep
        ,   new constel_boundary_t(10,  7.96666666666667,  9.16666666666667, 73.50000000000000)    // Cam
        ,   new constel_boundary_t(33,  9.16666666666667, 11.33333333333333, 73.50000000000000)    // Dra
        ,   new constel_boundary_t(83, 13.00000000000000, 16.53333333333333, 70.00000000000000)    // UMi
        ,   new constel_boundary_t(13,  3.10000000000000,  3.41666666666667, 68.00000000000000)    // Cas
        ,   new constel_boundary_t(33, 20.41666666666667, 20.66666666666667, 67.00000000000000)    // Dra
        ,   new constel_boundary_t(33, 11.33333333333333, 12.00000000000000, 66.50000000000000)    // Dra
        ,   new constel_boundary_t(15,  0.00000000000000,  0.33333333333333, 66.00000000000000)    // Cep
        ,   new constel_boundary_t(83, 14.00000000000000, 15.66666666666667, 66.00000000000000)    // UMi
        ,   new constel_boundary_t(15, 23.58333333333333, 24.00000000000000, 66.00000000000000)    // Cep
        ,   new constel_boundary_t(33, 12.00000000000000, 13.50000000000000, 64.00000000000000)    // Dra
        ,   new constel_boundary_t(33, 13.50000000000000, 14.41666666666667, 63.00000000000000)    // Dra
        ,   new constel_boundary_t(15, 23.16666666666667, 23.58333333333333, 63.00000000000000)    // Cep
        ,   new constel_boundary_t(10,  6.10000000000000,  7.00000000000000, 62.00000000000000)    // Cam
        ,   new constel_boundary_t(33, 20.00000000000000, 20.41666666666667, 61.50000000000000)    // Dra
        ,   new constel_boundary_t(15, 20.53666666666667, 20.60000000000000, 60.91666666666666)    // Cep
        ,   new constel_boundary_t(10,  7.00000000000000,  7.96666666666667, 60.00000000000000)    // Cam
        ,   new constel_boundary_t(82,  7.96666666666667,  8.41666666666667, 60.00000000000000)    // UMa
        ,   new constel_boundary_t(33, 19.76666666666667, 20.00000000000000, 59.50000000000000)    // Dra
        ,   new constel_boundary_t(15, 20.00000000000000, 20.53666666666667, 59.50000000000000)    // Cep
        ,   new constel_boundary_t(15, 22.86666666666667, 23.16666666666667, 59.08333333333334)    // Cep
        ,   new constel_boundary_t(13,  0.00000000000000,  2.43333333333333, 58.50000000000000)    // Cas
        ,   new constel_boundary_t(33, 19.41666666666667, 19.76666666666667, 58.00000000000000)    // Dra
        ,   new constel_boundary_t(13,  1.70000000000000,  1.90833333333333, 57.50000000000000)    // Cas
        ,   new constel_boundary_t(13,  2.43333333333333,  3.10000000000000, 57.00000000000000)    // Cas
        ,   new constel_boundary_t(10,  3.10000000000000,  3.16666666666667, 57.00000000000000)    // Cam
        ,   new constel_boundary_t(15, 22.31666666666667, 22.86666666666667, 56.25000000000000)    // Cep
        ,   new constel_boundary_t(10,  5.00000000000000,  6.10000000000000, 56.00000000000000)    // Cam
        ,   new constel_boundary_t(82, 14.03333333333333, 14.41666666666667, 55.50000000000000)    // UMa
        ,   new constel_boundary_t(33, 14.41666666666667, 19.41666666666667, 55.50000000000000)    // Dra
        ,   new constel_boundary_t(10,  3.16666666666667,  3.33333333333333, 55.00000000000000)    // Cam
        ,   new constel_boundary_t(15, 22.13333333333333, 22.31666666666667, 55.00000000000000)    // Cep
        ,   new constel_boundary_t(15, 20.60000000000000, 21.96666666666667, 54.83333333333334)    // Cep
        ,   new constel_boundary_t(13,  0.00000000000000,  1.70000000000000, 54.00000000000000)    // Cas
        ,   new constel_boundary_t(50,  6.10000000000000,  6.50000000000000, 54.00000000000000)    // Lyn
        ,   new constel_boundary_t(82, 12.08333333333333, 13.50000000000000, 53.00000000000000)    // UMa
        ,   new constel_boundary_t(33, 15.25000000000000, 15.75000000000000, 53.00000000000000)    // Dra
        ,   new constel_boundary_t(15, 21.96666666666667, 22.13333333333333, 52.75000000000000)    // Cep
        ,   new constel_boundary_t(10,  3.33333333333333,  5.00000000000000, 52.50000000000000)    // Cam
        ,   new constel_boundary_t(13, 22.86666666666667, 23.33333333333333, 52.50000000000000)    // Cas
        ,   new constel_boundary_t(33, 15.75000000000000, 17.00000000000000, 51.50000000000000)    // Dra
        ,   new constel_boundary_t(62,  2.04166666666667,  2.51666666666667, 50.50000000000000)    // Per
        ,   new constel_boundary_t(33, 17.00000000000000, 18.23333333333333, 50.50000000000000)    // Dra
        ,   new constel_boundary_t(13,  0.00000000000000,  1.36666666666667, 50.00000000000000)    // Cas
        ,   new constel_boundary_t(62,  1.36666666666667,  1.66666666666667, 50.00000000000000)    // Per
        ,   new constel_boundary_t(50,  6.50000000000000,  6.80000000000000, 50.00000000000000)    // Lyn
        ,   new constel_boundary_t(13, 23.33333333333333, 24.00000000000000, 50.00000000000000)    // Cas
        ,   new constel_boundary_t(82, 13.50000000000000, 14.03333333333333, 48.50000000000000)    // UMa
        ,   new constel_boundary_t(13,  0.00000000000000,  1.11666666666667, 48.00000000000000)    // Cas
        ,   new constel_boundary_t(13, 23.58333333333333, 24.00000000000000, 48.00000000000000)    // Cas
        ,   new constel_boundary_t(39, 18.17500000000000, 18.23333333333333, 47.50000000000000)    // Her
        ,   new constel_boundary_t(33, 18.23333333333333, 19.08333333333333, 47.50000000000000)    // Dra
        ,   new constel_boundary_t(30, 19.08333333333333, 19.16666666666667, 47.50000000000000)    // Cyg
        ,   new constel_boundary_t(62,  1.66666666666667,  2.04166666666667, 47.00000000000000)    // Per
        ,   new constel_boundary_t(82,  8.41666666666667,  9.16666666666667, 47.00000000000000)    // UMa
        ,   new constel_boundary_t(13,  0.16666666666667,  0.86666666666667, 46.00000000000000)    // Cas
        ,   new constel_boundary_t(82, 12.00000000000000, 12.08333333333333, 45.00000000000000)    // UMa
        ,   new constel_boundary_t(50,  6.80000000000000,  7.36666666666667, 44.50000000000000)    // Lyn
        ,   new constel_boundary_t(30, 21.90833333333333, 21.96666666666667, 44.00000000000000)    // Cyg
        ,   new constel_boundary_t(30, 21.87500000000000, 21.90833333333333, 43.75000000000000)    // Cyg
        ,   new constel_boundary_t(30, 19.16666666666667, 19.40000000000000, 43.50000000000000)    // Cyg
        ,   new constel_boundary_t(82,  9.16666666666667, 10.16666666666667, 42.00000000000000)    // UMa
        ,   new constel_boundary_t(82, 10.16666666666667, 10.78333333333333, 40.00000000000000)    // UMa
        ,   new constel_boundary_t( 8, 15.43333333333333, 15.75000000000000, 40.00000000000000)    // Boo
        ,   new constel_boundary_t(39, 15.75000000000000, 16.33333333333333, 40.00000000000000)    // Her
        ,   new constel_boundary_t(50,  9.25000000000000,  9.58333333333333, 39.75000000000000)    // Lyn
        ,   new constel_boundary_t( 0,  0.00000000000000,  2.51666666666667, 36.75000000000000)    // And
        ,   new constel_boundary_t(62,  2.51666666666667,  2.56666666666667, 36.75000000000000)    // Per
        ,   new constel_boundary_t(51, 19.35833333333333, 19.40000000000000, 36.50000000000000)    // Lyr
        ,   new constel_boundary_t(62,  4.50000000000000,  4.69166666666667, 36.00000000000000)    // Per
        ,   new constel_boundary_t(30, 21.73333333333333, 21.87500000000000, 36.00000000000000)    // Cyg
        ,   new constel_boundary_t(44, 21.87500000000000, 22.00000000000000, 36.00000000000000)    // Lac
        ,   new constel_boundary_t( 7,  6.53333333333333,  7.36666666666667, 35.50000000000000)    // Aur
        ,   new constel_boundary_t(50,  7.36666666666667,  7.75000000000000, 35.50000000000000)    // Lyn
        ,   new constel_boundary_t( 0,  0.00000000000000,  2.00000000000000, 35.00000000000000)    // And
        ,   new constel_boundary_t(44, 22.00000000000000, 22.81666666666667, 35.00000000000000)    // Lac
        ,   new constel_boundary_t(44, 22.81666666666667, 22.86666666666667, 34.50000000000000)    // Lac
        ,   new constel_boundary_t( 0, 22.86666666666667, 23.50000000000000, 34.50000000000000)    // And
        ,   new constel_boundary_t(62,  2.56666666666667,  2.71666666666667, 34.00000000000000)    // Per
        ,   new constel_boundary_t(82, 10.78333333333333, 11.00000000000000, 34.00000000000000)    // UMa
        ,   new constel_boundary_t(29, 12.00000000000000, 12.33333333333333, 34.00000000000000)    // CVn
        ,   new constel_boundary_t(50,  7.75000000000000,  9.25000000000000, 33.50000000000000)    // Lyn
        ,   new constel_boundary_t(48,  9.25000000000000,  9.88333333333333, 33.50000000000000)    // LMi
        ,   new constel_boundary_t( 0,  0.71666666666667,  1.40833333333333, 33.00000000000000)    // And
        ,   new constel_boundary_t( 8, 15.18333333333333, 15.43333333333333, 33.00000000000000)    // Boo
        ,   new constel_boundary_t( 0, 23.50000000000000, 23.75000000000000, 32.08333333333334)    // And
        ,   new constel_boundary_t(29, 12.33333333333333, 13.25000000000000, 32.00000000000000)    // CVn
        ,   new constel_boundary_t( 0, 23.75000000000000, 24.00000000000000, 31.33333333333333)    // And
        ,   new constel_boundary_t(29, 13.95833333333333, 14.03333333333333, 30.75000000000000)    // CVn
        ,   new constel_boundary_t(80,  2.41666666666667,  2.71666666666667, 30.66666666666667)    // Tri
        ,   new constel_boundary_t(62,  2.71666666666667,  4.50000000000000, 30.66666666666667)    // Per
        ,   new constel_boundary_t( 7,  4.50000000000000,  4.75000000000000, 30.00000000000000)    // Aur
        ,   new constel_boundary_t(51, 18.17500000000000, 19.35833333333333, 30.00000000000000)    // Lyr
        ,   new constel_boundary_t(82, 11.00000000000000, 12.00000000000000, 29.00000000000000)    // UMa
        ,   new constel_boundary_t(30, 19.66666666666667, 20.91666666666667, 29.00000000000000)    // Cyg
        ,   new constel_boundary_t( 7,  4.75000000000000,  5.88333333333333, 28.50000000000000)    // Aur
        ,   new constel_boundary_t(48,  9.88333333333333, 10.50000000000000, 28.50000000000000)    // LMi
        ,   new constel_boundary_t(29, 13.25000000000000, 13.95833333333333, 28.50000000000000)    // CVn
        ,   new constel_boundary_t( 0,  0.00000000000000,  0.06666666666667, 28.00000000000000)    // And
        ,   new constel_boundary_t(80,  1.40833333333333,  1.66666666666667, 28.00000000000000)    // Tri
        ,   new constel_boundary_t( 7,  5.88333333333333,  6.53333333333333, 28.00000000000000)    // Aur
        ,   new constel_boundary_t(37,  7.88333333333333,  8.00000000000000, 28.00000000000000)    // Gem
        ,   new constel_boundary_t(30, 20.91666666666667, 21.73333333333333, 28.00000000000000)    // Cyg
        ,   new constel_boundary_t(30, 19.25833333333333, 19.66666666666667, 27.50000000000000)    // Cyg
        ,   new constel_boundary_t(80,  1.91666666666667,  2.41666666666667, 27.25000000000000)    // Tri
        ,   new constel_boundary_t(25, 16.16666666666667, 16.33333333333333, 27.00000000000000)    // CrB
        ,   new constel_boundary_t( 8, 15.08333333333333, 15.18333333333333, 26.00000000000000)    // Boo
        ,   new constel_boundary_t(25, 15.18333333333333, 16.16666666666667, 26.00000000000000)    // CrB
        ,   new constel_boundary_t(51, 18.36666666666667, 18.86666666666667, 26.00000000000000)    // Lyr
        ,   new constel_boundary_t(48, 10.75000000000000, 11.00000000000000, 25.50000000000000)    // LMi
        ,   new constel_boundary_t(51, 18.86666666666667, 19.25833333333333, 25.50000000000000)    // Lyr
        ,   new constel_boundary_t(80,  1.66666666666667,  1.91666666666667, 25.00000000000000)    // Tri
        ,   new constel_boundary_t(66,  0.71666666666667,  0.85000000000000, 23.75000000000000)    // Psc
        ,   new constel_boundary_t(48, 10.50000000000000, 10.75000000000000, 23.50000000000000)    // LMi
        ,   new constel_boundary_t(87, 21.25000000000000, 21.41666666666667, 23.50000000000000)    // Vul
        ,   new constel_boundary_t(77,  5.70000000000000,  5.88333333333333, 22.83333333333333)    // Tau
        ,   new constel_boundary_t( 0,  0.06666666666667,  0.14166666666667, 22.00000000000000)    // And
        ,   new constel_boundary_t(73, 15.91666666666667, 16.03333333333333, 22.00000000000000)    // Ser
        ,   new constel_boundary_t(37,  5.88333333333333,  6.21666666666667, 21.50000000000000)    // Gem
        ,   new constel_boundary_t(87, 19.83333333333333, 20.25000000000000, 21.25000000000000)    // Vul
        ,   new constel_boundary_t(87, 18.86666666666667, 19.25000000000000, 21.08333333333333)    // Vul
        ,   new constel_boundary_t( 0,  0.14166666666667,  0.85000000000000, 21.00000000000000)    // And
        ,   new constel_boundary_t(87, 20.25000000000000, 20.56666666666667, 20.50000000000000)    // Vul
        ,   new constel_boundary_t(37,  7.80833333333333,  7.88333333333333, 20.00000000000000)    // Gem
        ,   new constel_boundary_t(87, 20.56666666666667, 21.25000000000000, 19.50000000000000)    // Vul
        ,   new constel_boundary_t(87, 19.25000000000000, 19.83333333333333, 19.16666666666667)    // Vul
        ,   new constel_boundary_t( 6,  3.28333333333333,  3.36666666666667, 19.00000000000000)    // Ari
        ,   new constel_boundary_t(75, 18.86666666666667, 19.00000000000000, 18.50000000000000)    // Sge
        ,   new constel_boundary_t(59,  5.70000000000000,  5.76666666666667, 18.00000000000000)    // Ori
        ,   new constel_boundary_t(37,  6.21666666666667,  6.30833333333333, 17.50000000000000)    // Gem
        ,   new constel_boundary_t(75, 19.00000000000000, 19.83333333333333, 16.16666666666667)    // Sge
        ,   new constel_boundary_t(77,  4.96666666666667,  5.33333333333333, 16.00000000000000)    // Tau
        ,   new constel_boundary_t(39, 15.91666666666667, 16.08333333333333, 16.00000000000000)    // Her
        ,   new constel_boundary_t(75, 19.83333333333333, 20.25000000000000, 15.75000000000000)    // Sge
        ,   new constel_boundary_t(77,  4.61666666666667,  4.96666666666667, 15.50000000000000)    // Tau
        ,   new constel_boundary_t(77,  5.33333333333333,  5.60000000000000, 15.50000000000000)    // Tau
        ,   new constel_boundary_t(23, 12.83333333333333, 13.50000000000000, 15.00000000000000)    // Com
        ,   new constel_boundary_t(39, 17.25000000000000, 18.25000000000000, 14.33333333333333)    // Her
        ,   new constel_boundary_t(23, 11.86666666666667, 12.83333333333333, 14.00000000000000)    // Com
        ,   new constel_boundary_t(37,  7.50000000000000,  7.80833333333333, 13.50000000000000)    // Gem
        ,   new constel_boundary_t(39, 16.75000000000000, 17.25000000000000, 12.83333333333333)    // Her
        ,   new constel_boundary_t(61,  0.00000000000000,  0.14166666666667, 12.50000000000000)    // Peg
        ,   new constel_boundary_t(77,  5.60000000000000,  5.76666666666667, 12.50000000000000)    // Tau
        ,   new constel_boundary_t(37,  7.00000000000000,  7.50000000000000, 12.50000000000000)    // Gem
        ,   new constel_boundary_t(61, 21.11666666666667, 21.33333333333333, 12.50000000000000)    // Peg
        ,   new constel_boundary_t(37,  6.30833333333333,  6.93333333333333, 12.00000000000000)    // Gem
        ,   new constel_boundary_t(39, 18.25000000000000, 18.86666666666667, 12.00000000000000)    // Her
        ,   new constel_boundary_t(31, 20.87500000000000, 21.05000000000000, 11.83333333333333)    // Del
        ,   new constel_boundary_t(61, 21.05000000000000, 21.11666666666667, 11.83333333333333)    // Peg
        ,   new constel_boundary_t(45, 11.51666666666667, 11.86666666666667, 11.00000000000000)    // Leo
        ,   new constel_boundary_t(59,  6.24166666666667,  6.30833333333333, 10.00000000000000)    // Ori
        ,   new constel_boundary_t(37,  6.93333333333333,  7.00000000000000, 10.00000000000000)    // Gem
        ,   new constel_boundary_t(21,  7.80833333333333,  7.92500000000000, 10.00000000000000)    // Cnc
        ,   new constel_boundary_t(61, 23.83333333333333, 24.00000000000000, 10.00000000000000)    // Peg
        ,   new constel_boundary_t( 6,  1.66666666666667,  3.28333333333333,  9.91666666666667)    // Ari
        ,   new constel_boundary_t(31, 20.14166666666667, 20.30000000000000,  8.50000000000000)    // Del
        ,   new constel_boundary_t( 8, 13.50000000000000, 15.08333333333333,  8.00000000000000)    // Boo
        ,   new constel_boundary_t(61, 22.75000000000000, 23.83333333333333,  7.50000000000000)    // Peg
        ,   new constel_boundary_t(21,  7.92500000000000,  9.25000000000000,  7.00000000000000)    // Cnc
        ,   new constel_boundary_t(45,  9.25000000000000, 10.75000000000000,  7.00000000000000)    // Leo
        ,   new constel_boundary_t(58, 18.25000000000000, 18.66222222222222,  6.25000000000000)    // Oph
        ,   new constel_boundary_t( 3, 18.66222222222222, 18.86666666666667,  6.25000000000000)    // Aql
        ,   new constel_boundary_t(31, 20.83333333333333, 20.87500000000000,  6.00000000000000)    // Del
        ,   new constel_boundary_t(20,  7.00000000000000,  7.01666666666667,  5.50000000000000)    // CMi
        ,   new constel_boundary_t(73, 18.25000000000000, 18.42500000000000,  4.50000000000000)    // Ser
        ,   new constel_boundary_t(39, 16.08333333333333, 16.75000000000000,  4.00000000000000)    // Her
        ,   new constel_boundary_t(58, 18.25000000000000, 18.42500000000000,  3.00000000000000)    // Oph
        ,   new constel_boundary_t(61, 21.46666666666667, 21.66666666666667,  2.75000000000000)    // Peg
        ,   new constel_boundary_t(66,  0.00000000000000,  2.00000000000000,  2.00000000000000)    // Psc
        ,   new constel_boundary_t(73, 18.58333333333333, 18.86666666666667,  2.00000000000000)    // Ser
        ,   new constel_boundary_t(31, 20.30000000000000, 20.83333333333333,  2.00000000000000)    // Del
        ,   new constel_boundary_t(34, 20.83333333333333, 21.33333333333333,  2.00000000000000)    // Equ
        ,   new constel_boundary_t(61, 21.33333333333333, 21.46666666666667,  2.00000000000000)    // Peg
        ,   new constel_boundary_t(61, 22.00000000000000, 22.75000000000000,  2.00000000000000)    // Peg
        ,   new constel_boundary_t(61, 21.66666666666667, 22.00000000000000,  1.75000000000000)    // Peg
        ,   new constel_boundary_t(20,  7.01666666666667,  7.20000000000000,  1.50000000000000)    // CMi
        ,   new constel_boundary_t(77,  3.58333333333333,  4.61666666666667,  0.00000000000000)    // Tau
        ,   new constel_boundary_t(59,  4.61666666666667,  4.66666666666667,  0.00000000000000)    // Ori
        ,   new constel_boundary_t(20,  7.20000000000000,  8.08333333333333,  0.00000000000000)    // CMi
        ,   new constel_boundary_t(85, 14.66666666666667, 15.08333333333333,  0.00000000000000)    // Vir
        ,   new constel_boundary_t(58, 17.83333333333333, 18.25000000000000,  0.00000000000000)    // Oph
        ,   new constel_boundary_t(16,  2.65000000000000,  3.28333333333333, -1.75000000000000)    // Cet
        ,   new constel_boundary_t(77,  3.28333333333333,  3.58333333333333, -1.75000000000000)    // Tau
        ,   new constel_boundary_t(73, 15.08333333333333, 16.26666666666667, -3.25000000000000)    // Ser
        ,   new constel_boundary_t(59,  4.66666666666667,  5.08333333333333, -4.00000000000000)    // Ori
        ,   new constel_boundary_t(59,  5.83333333333333,  6.24166666666667, -4.00000000000000)    // Ori
        ,   new constel_boundary_t(73, 17.83333333333333, 17.96666666666667, -4.00000000000000)    // Ser
        ,   new constel_boundary_t(73, 18.25000000000000, 18.58333333333333, -4.00000000000000)    // Ser
        ,   new constel_boundary_t( 3, 18.58333333333333, 18.86666666666667, -4.00000000000000)    // Aql
        ,   new constel_boundary_t(66, 22.75000000000000, 23.83333333333333, -4.00000000000000)    // Psc
        ,   new constel_boundary_t(45, 10.75000000000000, 11.51666666666667, -6.00000000000000)    // Leo
        ,   new constel_boundary_t(85, 11.51666666666667, 11.83333333333333, -6.00000000000000)    // Vir
        ,   new constel_boundary_t(66,  0.00000000000000,  0.33333333333333, -7.00000000000000)    // Psc
        ,   new constel_boundary_t(66, 23.83333333333333, 24.00000000000000, -7.00000000000000)    // Psc
        ,   new constel_boundary_t(85, 14.25000000000000, 14.66666666666667, -8.00000000000000)    // Vir
        ,   new constel_boundary_t(58, 15.91666666666667, 16.26666666666667, -8.00000000000000)    // Oph
        ,   new constel_boundary_t( 3, 20.00000000000000, 20.53333333333333, -9.00000000000000)    // Aql
        ,   new constel_boundary_t( 4, 21.33333333333333, 21.86666666666667, -9.00000000000000)    // Aqr
        ,   new constel_boundary_t(58, 17.16666666666667, 17.96666666666667, -10.00000000000000)    // Oph
        ,   new constel_boundary_t(54,  5.83333333333333,  8.08333333333333, -11.00000000000000)    // Mon
        ,   new constel_boundary_t(35,  4.91666666666667,  5.08333333333333, -11.00000000000000)    // Eri
        ,   new constel_boundary_t(59,  5.08333333333333,  5.83333333333333, -11.00000000000000)    // Ori
        ,   new constel_boundary_t(41,  8.08333333333333,  8.36666666666667, -11.00000000000000)    // Hya
        ,   new constel_boundary_t(74,  9.58333333333333, 10.75000000000000, -11.00000000000000)    // Sex
        ,   new constel_boundary_t(85, 11.83333333333333, 12.83333333333333, -11.00000000000000)    // Vir
        ,   new constel_boundary_t(58, 17.58333333333333, 17.66666666666667, -11.66666666666667)    // Oph
        ,   new constel_boundary_t( 3, 18.86666666666667, 20.00000000000000, -12.03333333333333)    // Aql
        ,   new constel_boundary_t(35,  4.83333333333333,  4.91666666666667, -14.50000000000000)    // Eri
        ,   new constel_boundary_t( 4, 20.53333333333333, 21.33333333333333, -15.00000000000000)    // Aqr
        ,   new constel_boundary_t(73, 17.16666666666667, 18.25000000000000, -16.00000000000000)    // Ser
        ,   new constel_boundary_t(72, 18.25000000000000, 18.86666666666667, -16.00000000000000)    // Sct
        ,   new constel_boundary_t(41,  8.36666666666667,  8.58333333333333, -17.00000000000000)    // Hya
        ,   new constel_boundary_t(58, 16.26666666666667, 16.37500000000000, -18.25000000000000)    // Oph
        ,   new constel_boundary_t(41,  8.58333333333333,  9.08333333333333, -19.00000000000000)    // Hya
        ,   new constel_boundary_t(26, 10.75000000000000, 10.83333333333333, -19.00000000000000)    // Crt
        ,   new constel_boundary_t(71, 16.26666666666667, 16.37500000000000, -19.25000000000000)    // Sco
        ,   new constel_boundary_t(47, 15.66666666666667, 15.91666666666667, -20.00000000000000)    // Lib
        ,   new constel_boundary_t(28, 12.58333333333333, 12.83333333333333, -22.00000000000000)    // Crv
        ,   new constel_boundary_t(85, 12.83333333333333, 14.25000000000000, -22.00000000000000)    // Vir
        ,   new constel_boundary_t(41,  9.08333333333333,  9.75000000000000, -24.00000000000000)    // Hya
        ,   new constel_boundary_t(16,  1.66666666666667,  2.65000000000000, -24.38333333333333)    // Cet
        ,   new constel_boundary_t(35,  2.65000000000000,  3.75000000000000, -24.38333333333333)    // Eri
        ,   new constel_boundary_t(26, 10.83333333333333, 11.83333333333333, -24.50000000000000)    // Crt
        ,   new constel_boundary_t(28, 11.83333333333333, 12.58333333333333, -24.50000000000000)    // Crv
        ,   new constel_boundary_t(47, 14.25000000000000, 14.91666666666667, -24.50000000000000)    // Lib
        ,   new constel_boundary_t(58, 16.26666666666667, 16.75000000000000, -24.58333333333333)    // Oph
        ,   new constel_boundary_t(16,  0.00000000000000,  1.66666666666667, -25.50000000000000)    // Cet
        ,   new constel_boundary_t(11, 21.33333333333333, 21.86666666666667, -25.50000000000000)    // Cap
        ,   new constel_boundary_t( 4, 21.86666666666667, 23.83333333333333, -25.50000000000000)    // Aqr
        ,   new constel_boundary_t(16, 23.83333333333333, 24.00000000000000, -25.50000000000000)    // Cet
        ,   new constel_boundary_t(41,  9.75000000000000, 10.25000000000000, -26.50000000000000)    // Hya
        ,   new constel_boundary_t(35,  4.70000000000000,  4.83333333333333, -27.25000000000000)    // Eri
        ,   new constel_boundary_t(46,  4.83333333333333,  6.11666666666667, -27.25000000000000)    // Lep
        ,   new constel_boundary_t(11, 20.00000000000000, 21.33333333333333, -28.00000000000000)    // Cap
        ,   new constel_boundary_t(41, 10.25000000000000, 10.58333333333333, -29.16666666666667)    // Hya
        ,   new constel_boundary_t(41, 12.58333333333333, 14.91666666666667, -29.50000000000000)    // Hya
        ,   new constel_boundary_t(47, 14.91666666666667, 15.66666666666667, -29.50000000000000)    // Lib
        ,   new constel_boundary_t(71, 15.66666666666667, 16.00000000000000, -29.50000000000000)    // Sco
        ,   new constel_boundary_t(35,  4.58333333333333,  4.70000000000000, -30.00000000000000)    // Eri
        ,   new constel_boundary_t(58, 16.75000000000000, 17.60000000000000, -30.00000000000000)    // Oph
        ,   new constel_boundary_t(76, 17.60000000000000, 17.83333333333333, -30.00000000000000)    // Sgr
        ,   new constel_boundary_t(41, 10.58333333333333, 10.83333333333333, -31.16666666666667)    // Hya
        ,   new constel_boundary_t(19,  6.11666666666667,  7.36666666666667, -33.00000000000000)    // CMa
        ,   new constel_boundary_t(41, 12.25000000000000, 12.58333333333333, -33.00000000000000)    // Hya
        ,   new constel_boundary_t(41, 10.83333333333333, 12.25000000000000, -35.00000000000000)    // Hya
        ,   new constel_boundary_t(36,  3.50000000000000,  3.75000000000000, -36.00000000000000)    // For
        ,   new constel_boundary_t(68,  8.36666666666667,  9.36666666666667, -36.75000000000000)    // Pyx
        ,   new constel_boundary_t(35,  4.26666666666667,  4.58333333333333, -37.00000000000000)    // Eri
        ,   new constel_boundary_t(76, 17.83333333333333, 19.16666666666667, -37.00000000000000)    // Sgr
        ,   new constel_boundary_t(65, 21.33333333333333, 23.00000000000000, -37.00000000000000)    // PsA
        ,   new constel_boundary_t(70, 23.00000000000000, 23.33333333333333, -37.00000000000000)    // Scl
        ,   new constel_boundary_t(36,  3.00000000000000,  3.50000000000000, -39.58333333333334)    // For
        ,   new constel_boundary_t( 1,  9.36666666666667, 11.00000000000000, -39.75000000000000)    // Ant
        ,   new constel_boundary_t(70,  0.00000000000000,  1.66666666666667, -40.00000000000000)    // Scl
        ,   new constel_boundary_t(36,  1.66666666666667,  3.00000000000000, -40.00000000000000)    // For
        ,   new constel_boundary_t(35,  3.86666666666667,  4.26666666666667, -40.00000000000000)    // Eri
        ,   new constel_boundary_t(70, 23.33333333333333, 24.00000000000000, -40.00000000000000)    // Scl
        ,   new constel_boundary_t(14, 14.16666666666667, 14.91666666666667, -42.00000000000000)    // Cen
        ,   new constel_boundary_t(49, 15.66666666666667, 16.00000000000000, -42.00000000000000)    // Lup
        ,   new constel_boundary_t(71, 16.00000000000000, 16.42083333333333, -42.00000000000000)    // Sco
        ,   new constel_boundary_t( 9,  4.83333333333333,  5.00000000000000, -43.00000000000000)    // Cae
        ,   new constel_boundary_t(22,  5.00000000000000,  6.58333333333333, -43.00000000000000)    // Col
        ,   new constel_boundary_t(67,  8.00000000000000,  8.36666666666667, -43.00000000000000)    // Pup
        ,   new constel_boundary_t(35,  3.41666666666667,  3.86666666666667, -44.00000000000000)    // Eri
        ,   new constel_boundary_t(71, 16.42083333333333, 17.83333333333333, -45.50000000000000)    // Sco
        ,   new constel_boundary_t(24, 17.83333333333333, 19.16666666666667, -45.50000000000000)    // CrA
        ,   new constel_boundary_t(76, 19.16666666666667, 20.33333333333333, -45.50000000000000)    // Sgr
        ,   new constel_boundary_t(53, 20.33333333333333, 21.33333333333333, -45.50000000000000)    // Mic
        ,   new constel_boundary_t(35,  3.00000000000000,  3.41666666666667, -46.00000000000000)    // Eri
        ,   new constel_boundary_t( 9,  4.50000000000000,  4.83333333333333, -46.50000000000000)    // Cae
        ,   new constel_boundary_t(49, 15.33333333333333, 15.66666666666667, -48.00000000000000)    // Lup
        ,   new constel_boundary_t(63,  0.00000000000000,  2.33333333333333, -48.16666666666666)    // Phe
        ,   new constel_boundary_t(35,  2.66666666666667,  3.00000000000000, -49.00000000000000)    // Eri
        ,   new constel_boundary_t(40,  4.08333333333333,  4.26666666666667, -49.00000000000000)    // Hor
        ,   new constel_boundary_t( 9,  4.26666666666667,  4.50000000000000, -49.00000000000000)    // Cae
        ,   new constel_boundary_t(38, 21.33333333333333, 22.00000000000000, -50.00000000000000)    // Gru
        ,   new constel_boundary_t(67,  6.00000000000000,  8.00000000000000, -50.75000000000000)    // Pup
        ,   new constel_boundary_t(84,  8.00000000000000,  8.16666666666667, -50.75000000000000)    // Vel
        ,   new constel_boundary_t(35,  2.41666666666667,  2.66666666666667, -51.00000000000000)    // Eri
        ,   new constel_boundary_t(40,  3.83333333333333,  4.08333333333333, -51.00000000000000)    // Hor
        ,   new constel_boundary_t(63,  0.00000000000000,  1.83333333333333, -51.50000000000000)    // Phe
        ,   new constel_boundary_t(12,  6.00000000000000,  6.16666666666667, -52.50000000000000)    // Car
        ,   new constel_boundary_t(84,  8.16666666666667,  8.45000000000000, -53.00000000000000)    // Vel
        ,   new constel_boundary_t(40,  3.50000000000000,  3.83333333333333, -53.16666666666666)    // Hor
        ,   new constel_boundary_t(32,  3.83333333333333,  4.00000000000000, -53.16666666666666)    // Dor
        ,   new constel_boundary_t(63,  0.00000000000000,  1.58333333333333, -53.50000000000000)    // Phe
        ,   new constel_boundary_t(35,  2.16666666666667,  2.41666666666667, -54.00000000000000)    // Eri
        ,   new constel_boundary_t(64,  4.50000000000000,  5.00000000000000, -54.00000000000000)    // Pic
        ,   new constel_boundary_t(49, 15.05000000000000, 15.33333333333333, -54.00000000000000)    // Lup
        ,   new constel_boundary_t(84,  8.45000000000000,  8.83333333333333, -54.50000000000000)    // Vel
        ,   new constel_boundary_t(12,  6.16666666666667,  6.50000000000000, -55.00000000000000)    // Car
        ,   new constel_boundary_t(14, 11.83333333333333, 12.83333333333333, -55.00000000000000)    // Cen
        ,   new constel_boundary_t(49, 14.16666666666667, 15.05000000000000, -55.00000000000000)    // Lup
        ,   new constel_boundary_t(56, 15.05000000000000, 15.33333333333333, -55.00000000000000)    // Nor
        ,   new constel_boundary_t(32,  4.00000000000000,  4.33333333333333, -56.50000000000000)    // Dor
        ,   new constel_boundary_t(84,  8.83333333333333, 11.00000000000000, -56.50000000000000)    // Vel
        ,   new constel_boundary_t(14, 11.00000000000000, 11.25000000000000, -56.50000000000000)    // Cen
        ,   new constel_boundary_t( 5, 17.50000000000000, 18.00000000000000, -57.00000000000000)    // Ara
        ,   new constel_boundary_t(78, 18.00000000000000, 20.33333333333333, -57.00000000000000)    // Tel
        ,   new constel_boundary_t(38, 22.00000000000000, 23.33333333333333, -57.00000000000000)    // Gru
        ,   new constel_boundary_t(40,  3.20000000000000,  3.50000000000000, -57.50000000000000)    // Hor
        ,   new constel_boundary_t(64,  5.00000000000000,  5.50000000000000, -57.50000000000000)    // Pic
        ,   new constel_boundary_t(12,  6.50000000000000,  6.83333333333333, -58.00000000000000)    // Car
        ,   new constel_boundary_t(63,  0.00000000000000,  1.33333333333333, -58.50000000000000)    // Phe
        ,   new constel_boundary_t(35,  1.33333333333333,  2.16666666666667, -58.50000000000000)    // Eri
        ,   new constel_boundary_t(63, 23.33333333333333, 24.00000000000000, -58.50000000000000)    // Phe
        ,   new constel_boundary_t(32,  4.33333333333333,  4.58333333333333, -59.00000000000000)    // Dor
        ,   new constel_boundary_t(56, 15.33333333333333, 16.42083333333333, -60.00000000000000)    // Nor
        ,   new constel_boundary_t(43, 20.33333333333333, 21.33333333333333, -60.00000000000000)    // Ind
        ,   new constel_boundary_t(64,  5.50000000000000,  6.00000000000000, -61.00000000000000)    // Pic
        ,   new constel_boundary_t(18, 15.16666666666667, 15.33333333333333, -61.00000000000000)    // Cir
        ,   new constel_boundary_t( 5, 16.42083333333333, 16.58333333333333, -61.00000000000000)    // Ara
        ,   new constel_boundary_t(18, 14.91666666666667, 15.16666666666667, -63.58333333333334)    // Cir
        ,   new constel_boundary_t( 5, 16.58333333333333, 16.75000000000000, -63.58333333333334)    // Ara
        ,   new constel_boundary_t(64,  6.00000000000000,  6.83333333333333, -64.00000000000000)    // Pic
        ,   new constel_boundary_t(12,  6.83333333333333,  9.03333333333333, -64.00000000000000)    // Car
        ,   new constel_boundary_t(14, 11.25000000000000, 11.83333333333333, -64.00000000000000)    // Cen
        ,   new constel_boundary_t(27, 11.83333333333333, 12.83333333333333, -64.00000000000000)    // Cru
        ,   new constel_boundary_t(14, 12.83333333333333, 14.53333333333333, -64.00000000000000)    // Cen
        ,   new constel_boundary_t(18, 13.50000000000000, 13.66666666666667, -65.00000000000000)    // Cir
        ,   new constel_boundary_t( 5, 16.75000000000000, 16.83333333333333, -65.00000000000000)    // Ara
        ,   new constel_boundary_t(40,  2.16666666666667,  3.20000000000000, -67.50000000000000)    // Hor
        ,   new constel_boundary_t(69,  3.20000000000000,  4.58333333333333, -67.50000000000000)    // Ret
        ,   new constel_boundary_t(18, 14.75000000000000, 14.91666666666667, -67.50000000000000)    // Cir
        ,   new constel_boundary_t( 5, 16.83333333333333, 17.50000000000000, -67.50000000000000)    // Ara
        ,   new constel_boundary_t(60, 17.50000000000000, 18.00000000000000, -67.50000000000000)    // Pav
        ,   new constel_boundary_t(81, 22.00000000000000, 23.33333333333333, -67.50000000000000)    // Tuc
        ,   new constel_boundary_t(32,  4.58333333333333,  6.58333333333333, -70.00000000000000)    // Dor
        ,   new constel_boundary_t(18, 13.66666666666667, 14.75000000000000, -70.00000000000000)    // Cir
        ,   new constel_boundary_t(79, 14.75000000000000, 17.00000000000000, -70.00000000000000)    // TrA
        ,   new constel_boundary_t(81,  0.00000000000000,  1.33333333333333, -75.00000000000000)    // Tuc
        ,   new constel_boundary_t(42,  3.50000000000000,  4.58333333333333, -75.00000000000000)    // Hyi
        ,   new constel_boundary_t(86,  6.58333333333333,  9.03333333333333, -75.00000000000000)    // Vol
        ,   new constel_boundary_t(12,  9.03333333333333, 11.25000000000000, -75.00000000000000)    // Car
        ,   new constel_boundary_t(55, 11.25000000000000, 13.66666666666667, -75.00000000000000)    // Mus
        ,   new constel_boundary_t(60, 18.00000000000000, 21.33333333333333, -75.00000000000000)    // Pav
        ,   new constel_boundary_t(43, 21.33333333333333, 23.33333333333333, -75.00000000000000)    // Ind
        ,   new constel_boundary_t(81, 23.33333333333333, 24.00000000000000, -75.00000000000000)    // Tuc
        ,   new constel_boundary_t(81,  0.75000000000000,  1.33333333333333, -76.00000000000000)    // Tuc
        ,   new constel_boundary_t(42,  0.00000000000000,  3.50000000000000, -82.50000000000000)    // Hyi
        ,   new constel_boundary_t(17,  7.66666666666667, 13.66666666666667, -82.50000000000000)    // Cha
        ,   new constel_boundary_t( 2, 13.66666666666667, 18.00000000000000, -82.50000000000000)    // Aps
        ,   new constel_boundary_t(52,  3.50000000000000,  7.66666666666667, -85.00000000000000)    // Men
        ,   new constel_boundary_t(57,  0.00000000000000, 24.00000000000000, -90.00000000000000)    // Oct
        };



    }
}
