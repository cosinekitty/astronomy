/*
    Astronomy Engine for C# / .NET.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019-2023 Don Cross <cosinekitty@gmail.com>

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
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

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
        internal EarthNotAllowedException():
            base("The Earth is not allowed as the body parameter.")
            {}
    }

    /// <summary>
    /// This exception is thrown by certain Astronomy Engine functions
    /// when a body is specified that is not appropriate for the given operation.
    /// </summary>
    public class InvalidBodyException: ArgumentException
    {
        internal InvalidBodyException(Body body):
            base("Invalid body: " + body)
            {}
    }

    /// <summary>
    /// This exception indicates an unexpected error occurred inside Astronomy Engine.
    /// Please report any such errors by creating an issue at:
    /// https://github.com/cosinekitty/astronomy/issues
    /// </summary>
    public class InternalError: Exception
    {
        internal InternalError(string message):
            base("Internal error. Please report an issue at: https://github.com/cosinekitty/astronomy/issues. Diagnostic: " + message)
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

        /// <summary>
        /// User-defined star #1.
        /// </summary>
        Star1 = 101,

        /// <summary>
        /// User-defined star #2.
        /// </summary>
        Star2,

        /// <summary>
        /// User-defined star #3.
        /// </summary>
        Star3,

        /// <summary>
        /// User-defined star #4.
        /// </summary>
        Star4,

        /// <summary>
        /// User-defined star #5.
        /// </summary>
        Star5,

        /// <summary>
        /// User-defined star #6.
        /// </summary>
        Star6,

        /// <summary>
        /// User-defined star #7.
        /// </summary>
        Star7,

        /// <summary>
        /// User-defined star #8.
        /// </summary>
        Star8,
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

        internal double psi = double.NaN;    // For internal use only. Used to optimize Earth tilt calculations.
        internal double eps = double.NaN;    // For internal use only. Used to optimize Earth tilt calculations.
        internal double st = double.NaN;     // For internal use only.  Lazy-caches sidereal time (Earth rotation).

        private AstroTime(double ut, double tt)
        {
            this.ut = ut;
            this.tt = tt;
        }

        /// <summary>
        /// Creates an `AstroTime` object from a Universal Time day value.
        /// </summary>
        /// <param name="ut">The number of days after the J2000 epoch.</param>
        public AstroTime(double ut)
            : this(ut, Astronomy.TerrestrialTime(ut))
        {
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
        /// <param name="second">The UTC second value [0, 60).</param>
        public AstroTime(int year, int month, int day, int hour, int minute, double second)
            : this(UniversalTimeFromCalendar(year, month, day, hour, minute, second))
        {
        }

        /// <summary>
        /// Creates an `AstroTime` object from a Terrestrial Time day value.
        /// </summary>
        /// <remarks>
        /// This function can be used in rare cases where a time must be based
        /// on Terrestrial Time (TT) rather than Universal Time (UT).
        /// Most developers will want to invoke `new AstroTime(ut)` with a universal time
        /// instead of this function, because usually time is based on civil time adjusted
        /// by leap seconds to match the Earth's rotation, rather than the uniformly
        /// flowing TT used to calculate solar system dynamics. In rare cases
        /// where the caller already knows TT, this function is provided to create
        /// an `AstroTime` value that can be passed to Astronomy Engine functions.
        /// </remarks>
        /// <param name="tt">The number of days after the J2000 epoch.</param>
        public static AstroTime FromTerrestrialTime(double tt)
        {
            return new AstroTime(Astronomy.UniversalTime(tt), tt);
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
        /// Converts this object to our custom type #CalendarDateTime.
        /// </summary>
        /// <remarks>
        /// The .NET type `DateTime` can only represent years in the range 0000..9999.
        /// However, the Astronomy Engine type #CalendarDateTime can represent
        /// years in the range -999999..+999999. This is a time span of nearly 2 million years.
        /// This function converts this `AstroTime` object to an equivalent Gregorian calendar representation.
        /// </remarks>
        public CalendarDateTime ToCalendarDateTime() => new CalendarDateTime(ut);

        /// <summary>
        /// Converts this `AstroTime` to ISO 8601 format, expressed in UTC with millisecond resolution.
        /// </summary>
        /// <returns>Example: "2019-08-30T17:45:22.763Z".</returns>
        public override string ToString() => ToCalendarDateTime().ToString();

        private static Regex re = new Regex(
            @"^
            ([\+\-]?[0-9]{4,6})     # 1 : year  : should be 4-digit, or 6-digit with +/- prefix, but be flexible
            -(0[0-9]|1[012])        # 2 : month
            -([012][0-9]|3[01])     # 3 : day
            T([01][0-9]|2[0-3])     # 4 : hour
            :([0-5][0-9])           # 5 : minute
            (                       # 6
                :
                (                   # 7
                    [0-5][0-9]      # optional seconds
                    (\.[0-9]+)?     # optional fraction of a second
                )
            )?
            Z$                      # terminator",
            RegexOptions.Compiled | RegexOptions.CultureInvariant | RegexOptions.IgnorePatternWhitespace
        );

        /// <summary>
        /// Converts a string of the format returned by #AstroTime.ToString back into an `AstroTime`.
        /// </summary>
        /// <remarks>
        /// This function attempts to parse an ISO 8601 formatted date and time string
        /// into an `AstroTime` object.
        /// If the string is valid, sets `time` to a new object and returns `true`.
        /// If the string is not valid, sets `time` to `null` and returns `false`.
        /// </remarks>
        /// <param name="text">The string from which to parse a date and time.</param>
        /// <param name="time">On success, receives the date and time value. On failure, receives `null`.</param>
        public static bool TryParse(string text, out AstroTime time)
        {
            time = null;

            if (text == null)
                return false;

            Match m = re.Match(text);
            if (!m.Success)
                return false;

            if (!int.TryParse(m.Groups[1].Value, out int year))
                return false;

            if (!int.TryParse(m.Groups[2].Value, out int month))
                return false;

            if (!int.TryParse(m.Groups[3].Value, out int day))
                return false;

            if (!int.TryParse(m.Groups[4].Value, out int hour))
                return false;

            if (!int.TryParse(m.Groups[5].Value, out int minute))
                return false;

            double second = 0.0;
            string stext = m.Groups[7].Value;
            if (!string.IsNullOrEmpty(stext))
            {
                var styles = NumberStyles.AllowLeadingSign | NumberStyles.AllowDecimalPoint | NumberStyles.AllowExponent;
                if (!double.TryParse(stext, styles, CultureInfo.InvariantCulture, out second))
                    return false;
            }

            time = new AstroTime(year, month, day, hour, minute, second);
            return true;
        }

        /// <summary>
        /// Calculates the sum or difference of an #AstroTime with a specified floating point number of days.
        /// </summary>
        /// <remarks>
        /// Sometimes we need to adjust a given #AstroTime value by a certain amount of time.
        /// This function adds the given real number of days in `days` to the date and time in this object.
        ///
        /// More precisely, the result's Universal Time field `ut` is exactly adjusted by `days` and
        /// the Terrestrial Time field `tt` is adjusted for the resulting UTC date and time,
        /// using a best-fit piecewise polynomial model devised by
        /// [Espenak and Meeus](https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html).
        /// </remarks>
        /// <param name="days">A floating point number of days by which to adjust `time`. May be negative, 0, or positive.</param>
        /// <returns>A date and time that is conceptually equal to `time + days`.</returns>
        public AstroTime AddDays(double days)
        {
            return new AstroTime(this.ut + days);
        }

        /// <summary>
        /// Nutation angle `psi`. Intended for unit testing only.
        /// </summary>
        public double Psi => psi;

        /// <summary>
        /// Nutation angle `eps`. Intended for unit testing only.
        /// </summary>
        public double Eps => eps;

        private static double UniversalTimeFromCalendar(int year, int month, int day, int hour, int minute, double second)
        {
            // This formula is adapted from NOVAS C 3.1 function julian_date(),
            // which in turn comes from Henry F. Fliegel & Thomas C. Van Flendern:
            // Communications of the ACM, Vol 11, No 10, October 1968, p. 657.
            // See: https://dl.acm.org/doi/pdf/10.1145/364096.364097
            //
            // [Don Cross - 2023-02-25] I modified the formula so that it will
            // work correctly with years as far back as -999999.

            long y = (long)year;
            long m = (long)month;
            long d = (long)day;
            long f = (14 - m) / 12;

            long y2000 = (
                (d - 365972956)
                + (1461*(y + 1000000 - f))/4
                + (367*(m - 2 + 12*f))/12
                - (3*((y + 1000100 - f) / 100))/4
            );

            double ut = (y2000 - 0.5) + (hour / 24.0) + (minute / 1440.0) + (second / 86400.0);
            return ut;
        }
    }

    /// <summary>
    /// Represents a Gregorian calendar date and time within plus or minus 1 million years from the year 0.
    /// </summary>
    /// <remarks>
    /// The C# standard type `System.DateTime` only allows years from 0001 to 9999.
    /// However, the #AstroTime class can represent years in the range -999999 to +999999.
    /// In order to support formatting dates with extreme year values in an extrapolated
    /// Gregorian calendar, the `CalendarDateTime` class breaks out the components of
    /// a date into separate fields.
    /// </remarks>
    public struct CalendarDateTime
    {
        /// <summary>The year value in the range -999999 to +999999.</summary>
        public int year;

        /// <summary>The calendar month in the range 1..12.</summary>
        public int month;

        /// <summary>The day of the month in the reange 1..31.</summary>
        public int day;

        /// <summary>The hour in the range 0..23.</summary>
        public int hour;

        /// <summary>The minute in the range 0..59.</summary>
        public int minute;

        /// <summary>The real-valued second in the half-open range [0, 60).</summary>
        public double second;

        /// <summary>Convert a J2000 day value to a Gregorian calendar date.</summary>
        /// <param name="ut">The real-valued number of days since the J2000 epoch.</param>
        public CalendarDateTime(double ut)
        {
            // Adapted from the NOVAS C 3.1 function cal_date().
            // Convert fractional days since J2000 into Gregorian calendar date/time.

            double djd = ut + 2451545.5;
            long jd = (long)Math.Floor(djd);
            double x = 24.0 * (djd % 1.0);
            if (x < 0.0)
                x += 24.0;
            hour = (int)x;
            x = 60.0 * (x % 1.0);
            minute = (int)x;
            second = 60.0 * (x % 1.0);

            // This is my own adjustment to the NOVAS cal_date logic
            // so that it can handle dates much farther back in the past.
            // I add c*400 years worth of days at the front,
            // then subtract c*400 years at the back,
            // which avoids negative values in the formulas that mess up
            // the calendar date calculations.
            // Any multiple of 400 years has the same number of days,
            // because it eliminates all the special cases for leap years.
            const long c = 2500;

            long k = jd + (68569 + c*146097);
            long n = (4 * k) / 146097;
            k = k - (146097*n + 3) / 4;
            long m = (4000 * (k+1)) / 1461001;
            k = k - (1461 * m)/4 + 31;

            month = (int) ((80 * k) / 2447);
            day = (int) (k - (2447 * month)/80);
            k = month / 11;

            month = (int) (month + 2 - 12*k);
            year = (int) (100 * (n - 49) + m + k - 400*c);

            if (year < -999999 || year > +999999)
                throw new ArgumentOutOfRangeException("The supplied time is too far from the year 2000 to be represented.");

            if (month < 1 || month > 12 || day < 1 || day > 31)
                throw new InternalError($"Invalid calendar date calculated: month={month}, day={day}.");
        }

        /// <summary>
        /// Converts this `CalendarDateTime` to ISO 8601 format, expressed in UTC with millisecond resolution.
        /// </summary>
        /// <returns>Example: "2019-08-30T17:45:22.763Z".</returns>
        public override string ToString()
        {
            int millis = Math.Max(0, Math.Min(59999, (int)Math.Round(second * 1000.0)));
            string y;
            if (year < 0)
                y = "-" + (-year).ToString("000000");
            else if (year <= 9999)
                y = year.ToString("0000");
            else
                y = "+" + year.ToString("000000");
            return $"{y}-{month:00}-{day:00}T{hour:00}:{minute:00}:{millis/1000:00}.{millis%1000:000}Z";
        }
    }

    internal struct TerseVector
    {
        public double x;
        public double y;
        public double z;

        public TerseVector(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public static readonly TerseVector Zero = new TerseVector(0.0, 0.0, 0.0);

        public AstroVector ToAstroVector(AstroTime time)
        {
            return new AstroVector(x, y, z, time);
        }

        public static TerseVector operator +(TerseVector a, TerseVector b)
        {
            return new TerseVector(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        public static TerseVector operator -(TerseVector a, TerseVector b)
        {
            return new TerseVector(a.x - b.x, a.y - b.y, a.z - b.z);
        }

        public static TerseVector operator -(TerseVector a)
        {
            return new TerseVector(-a.x, -a.y, -a.z);
        }

        public static TerseVector operator *(double s, TerseVector v)
        {
            return new TerseVector(s*v.x, s*v.y, s*v.z);
        }

        public static TerseVector operator /(TerseVector v, double s)
        {
            return new TerseVector(v.x/s, v.y/s, v.z/s);
        }

        public double Quadrature()
        {
            return x*x + y*y + z*z;
        }

        public double Magnitude()
        {
            return Math.Sqrt(Quadrature());
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
            if (t == null)
                throw new NullReferenceException("AstroTime parameter is not allowed to be null.");

            this.x = x;
            this.y = y;
            this.z = z;
            this.t = t;
        }

        /// <summary>
        /// Converts the vector to a string of the format (x, y, z, t).
        /// </summary>
        public override string ToString()
        {
            return $"({x:G16}, {y:G16}, {z:G16}, {t})";
        }

        // (0.1428571428571428, 1.333333333333333, 3.846153846153846E-07, 2023-02-14T09:45:30.000Z)
        private static Regex re = new Regex(
            @"^\s*\(\s*                 # (
            ([^\s,]+) \s* , \s*         # x ,
            ([^\s,]+) \s* , \s*         # y ,
            ([^\s,]+) \s* , \s*         # z ,
            ([^\s\)]+) \s* \) \s* $     # t )",
            RegexOptions.Compiled | RegexOptions.CultureInvariant | RegexOptions.IgnorePatternWhitespace
        );

        /// <summary>
        /// Parses a vector from a string as formatted by #AstroVector.ToString.
        /// On success, `vector` receives the vector and the function returns `true`.
        /// Otherwise, `vector` receives the value (0, 0, 0, null) and the function returns `false`.
        /// </summary>
        /// <param name="text">A string of the form "(x, y, z, t)".</param>
        /// <param name="vector">Receives the output vector.</param>
        public static bool TryParse(string text, out AstroVector vector)
        {
            vector = new AstroVector();
            if (text != null)
            {
                Match m = re.Match(text);
                if (m.Success)
                {
                    var styles = NumberStyles.AllowLeadingSign | NumberStyles.AllowDecimalPoint | NumberStyles.AllowExponent;
                    return (
                        double.TryParse(m.Groups[1].Value, styles, CultureInfo.InvariantCulture, out vector.x) &&
                        double.TryParse(m.Groups[2].Value, styles, CultureInfo.InvariantCulture, out vector.y) &&
                        double.TryParse(m.Groups[3].Value, styles, CultureInfo.InvariantCulture, out vector.z) &&
                        AstroTime.TryParse(m.Groups[4].Value, out vector.t)
                    );
                }
            }
            return false;
        }

        /// <summary>
        /// Calculates the total distance in AU represented by this vector.
        /// </summary>
        /// <returns>The nonnegative length of the Cartisian vector in AU.</returns>
        public double Length()
        {
            return Astronomy.hypot(x, y, z);
        }

#pragma warning disable 1591        // we don't need XML documentation for these operator overloads
        public static AstroVector operator - (AstroVector a)
        {
            return new AstroVector(-a.x, -a.y, -a.z, a.t);
        }

        public static AstroVector operator - (AstroVector a, AstroVector b)
        {
            return new AstroVector (
                a.x - b.x,
                a.y - b.y,
                a.z - b.z,
                VerifyIdenticalTimes(a.t, b.t)
            );
        }

        public static AstroVector operator + (AstroVector a, AstroVector b)
        {
            return new AstroVector (
                a.x + b.x,
                a.y + b.y,
                a.z + b.z,
                VerifyIdenticalTimes(a.t, b.t)
            );
        }

        public static double operator * (AstroVector a, AstroVector b)
        {
            // the scalar dot product of two vectors
            VerifyIdenticalTimes(a.t, b.t);
            return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
        }

        public static AstroVector operator * (double factor, AstroVector a)
        {
            return new AstroVector(
                factor * a.x,
                factor * a.y,
                factor * a.z,
                a.t
            );
        }

        public static AstroVector operator / (AstroVector a, double denom)
        {
            if (denom == 0.0)
                throw new ArgumentException("Attempt to divide a vector by zero.");

            return new AstroVector(
                a.x / denom,
                a.y / denom,
                a.z / denom,
                a.t
            );
        }
#pragma warning restore 1591

        private static AstroTime VerifyIdenticalTimes(AstroTime a, AstroTime b)
        {
            if (a.tt != b.tt)
                throw new ArgumentException("Attempt to operate on two vectors from different times.");

            // If either time has already had its nutation calculated, retain that work.
            return !double.IsNaN(a.psi) ? a : b;
        }
    }

    /// <summary>
    /// A combination of a position vector and a velocity vector at a given moment in time.
    /// </summary>
    /// <remarks>
    /// A state vector represents the dynamic state of a point at a given moment.
    /// It includes the position vector of the point, expressed in Astronomical Units (AU)
    /// along with the velocity vector of the point, expressed in AU/day.
    /// </remarks>
    public struct StateVector
    {
        /// <summary>
        /// The position x-coordinate in AU.
        /// </summary>
        public double x;

        /// <summary>
        /// The position y-coordinate in AU.
        /// </summary>
        public double y;

        /// <summary>
        /// The position z-coordinate in AU.
        /// </summary>
        public double z;

        /// <summary>
        /// The velocity x-component in AU/day.
        /// </summary>
        public double vx;

        /// <summary>
        /// The velocity y-component in AU/day.
        /// </summary>
        public double vy;

        /// <summary>
        /// The velocity z-component in AU/day.
        /// </summary>
        public double vz;

        /// <summary>
        /// The date and time at which this vector is valid.
        /// </summary>
        public AstroTime t;

        /// <summary>
        /// Creates an AstroVector.
        /// </summary>
        /// <param name="x">A position x-coordinate expressed in AU.</param>
        /// <param name="y">A position y-coordinate expressed in AU.</param>
        /// <param name="z">A position z-coordinate expressed in AU.</param>
        /// <param name="vx">A velocity x-component expressed in AU/day.</param>
        /// <param name="vy">A velocity y-component expressed in AU/day.</param>
        /// <param name="vz">A velocity z-component expressed in AU/day.</param>
        /// <param name="t">The date and time at which this state vector is valid.</param>
        public StateVector(double x, double y, double z, double vx, double vy, double vz, AstroTime t)
        {
            if (t == null)
                throw new NullReferenceException("AstroTime parameter is not allowed to be null.");

            this.x = x;
            this.y = y;
            this.z = z;
            this.vx = vx;
            this.vy = vy;
            this.vz = vz;
            this.t = t;
        }

        /// <summary>
        /// Combines a position vector and a velocity vector into a single state vector.
        /// </summary>
        /// <param name="pos">A position vector.</param>
        /// <param name="vel">A velocity vector.</param>
        /// <param name="time">The common time that represents the given position and velocity.</param>
        public StateVector(AstroVector pos, AstroVector vel, AstroTime time)
        {
            if (time == null)
                throw new NullReferenceException("AstroTime parameter is not allowed to be null.");

            this.x = pos.x;
            this.y = pos.y;
            this.z = pos.z;
            this.vx = vel.x;
            this.vy = vel.y;
            this.vz = vel.z;
            this.t = time;
        }

        /// <summary>
        /// Converts the state vector to a string of the format (x, y, z, vx, vy, vz, t).
        /// </summary>
        public override string ToString()
        {
            return $"({x:G16}, {y:G16}, {z:G16}, {vx:G16}, {vy:G16}, {vz:G16}, {t})";
        }

        /// <summary>
        /// Returns the position vector associated with this state vector.
        /// </summary>
        public AstroVector Position()
        {
            return new AstroVector(x, y, z, t);
        }

        /// <summary>
        /// Returns the velocity vector associated with this state vector.
        /// </summary>
        public AstroVector Velocity()
        {
            return new AstroVector(vx, vy, vz, t);
        }
    }

    /// <summary>
    /// Holds the positions and velocities of Jupiter's major 4 moons.
    /// </summary>
    /// <remarks>
    /// The #Astronomy.JupiterMoons function returns an object of this type
    /// to report position and velocity vectors for Jupiter's largest 4 moons
    /// Io, Europa, Ganymede, and Callisto. Each position vector is relative
    /// to the center of Jupiter. Both position and velocity are oriented in
    /// the EQJ system (that is, using Earth's equator at the J2000 epoch).
    /// The positions are expressed in astronomical units (AU),
    /// and the velocities in AU/day.
    /// </remarks>
    public struct JupiterMoonsInfo
    {
        /// <summary>The position and velocity of Jupiter's moon Io.</summary>
        public StateVector io;

        /// <summary>The position and velocity of Jupiter's moon Europa.</summary>
        public StateVector europa;

        /// <summary>The position and velocity of Jupiter's moon Ganymede.</summary>
        public StateVector ganymede;

        /// <summary>The position and velocity of Jupiter's moon Callisto.</summary>
        public StateVector callisto;
    }

    /// <summary>
    /// A rotation matrix that can be used to transform one coordinate system to another.
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

        /// <summary>
        /// Converts an `Observer` to a string representation like `(N 26.728965, W 093.157562, 1234.567 m)`.
        /// </summary>
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("(");
            sb.Append(latitude < 0.0 ? "S " : "N ");
            sb.Append(Math.Abs(latitude).ToString("00.000000"));
            sb.Append(", ");
            sb.Append(longitude < 0.0 ? "W " : "E ");
            sb.Append(Math.Abs(longitude).ToString("000.000000"));
            sb.Append(", ");
            sb.Append(height.ToString("0.000"));
            sb.Append(" m)");
            return sb.ToString();
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
    /// Equatorial angular and cartesian coordinates.
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
        /// Equatorial coordinates in cartesian vector form: x = March equinox, y = June solstice, z = north.
        /// </summary>
        public readonly AstroVector vec;

        internal Equatorial(double ra, double dec, double dist, AstroVector vec)
        {
            this.ra = ra;
            this.dec = dec;
            this.dist = dist;
            this.vec = vec;
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
        /// Cartesian ecliptic vector, with components as follows:
        /// x: the direction of the equinox along the ecliptic plane.
        /// y: in the ecliptic plane 90 degrees prograde from the equinox.
        /// z: perpendicular to the ecliptic plane. Positive is north.
        /// </summary>
        public readonly AstroVector vec;

        /// <summary>
        /// Latitude in degrees north (positive) or south (negative) of the ecliptic plane.
        /// </summary>
        public readonly double elat;

        /// <summary>
        /// Longitude in degrees around the ecliptic plane prograde from the equinox.
        /// </summary>
        public readonly double elon;

        internal Ecliptic(AstroVector vec, double elat, double elon)
        {
            this.vec = vec;
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

        internal Topocentric(double azimuth, double altitude, double ra, double dec)
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
    /// Information about idealized atmospheric variables at a given elevation.
    /// </summary>
    public struct AtmosphereInfo
    {
        /// <summary>
        /// Atmospheric pressure in pascals.
        /// </summary>
        public double pressure;

        /// <summary>
        /// Atmospheric temperature in kelvins.
        /// </summary>
        public double temperature;

        /// <summary>
        /// Atmospheric density relative to sea level.
        /// </summary>
        public double density;
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
    /// Lunar libration angles, returned by #Astronomy.Libration.
    /// </summary>
    public struct LibrationInfo
    {
        /// <summary>Sub-Earth libration ecliptic latitude angle, in degrees.</summary>
        public double elat;

        /// <summary>Sub-Earth libration ecliptic longitude angle, in degrees.</summary>
        public double elon;

        /// <summary>Moon's geocentric ecliptic latitude in degrees.</summary>
        public double mlat;

        /// <summary>Moon's geocentric ecliptic longitude in degrees.</summary>
        public double mlon;

        /// <summary>Distance between the centers of the Earth and Moon in kilometers.</summary>
        public double dist_km;

        /// <summary>The apparent angular diameter of the Moon, in degrees, as seen from the center of the Earth.</summary>
        public double diam_deg;
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

        internal HourAngleInfo(AstroTime time, Topocentric hor)
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

        internal ElongationInfo(AstroTime time, Visibility visibility, double elongation, double ecliptic_separation)
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

    /// <summary>The different kinds of lunar/solar eclipses.</summary>
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
    /// Penumbral eclipses are difficult to observe, because the Moon is only slightly dimmed
    /// by the Earth's penumbra; no part of the Moon touches the Earth's umbra.
    /// Partial eclipses occur when part, but not all, of the Moon touches the Earth's umbra.
    /// Total eclipses occur when the entire Moon passes into the Earth's umbra.
    ///
    /// The `kind` field thus holds `EclipseKind.Penumbral`, `EclipseKind.Partial`,
    /// or `EclipseKind.Total`, depending on the kind of lunar eclipse found.
    ///
    /// The `obscuration` field holds a value in the range [0, 1] that indicates what fraction
    /// of the Moon's apparent disc area is covered by the Earth's umbra at the eclipse's peak.
    /// This indicates how dark the peak eclipse appears. For penumbral eclipses, the obscuration
    /// is 0, because the Moon does not pass through the Earth's umbra. For partial eclipses,
    /// the obscuration is somewhere between 0 and 1. For total lunar eclipses, the obscuration is 1.
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

        /// <summary>The peak fraction of the Moon's apparent disc that is covered by the Earth's umbra.</summary>
        public double obscuration;

        /// <summary>The time of the eclipse at its peak.</summary>
        public AstroTime peak;

        /// <summary>The semi-duration of the penumbral phase in minutes.</summary>
        public double sd_penum;

        /// <summary>The semi-duration of the partial phase in minutes, or 0.0 if none.</summary>
        public double sd_partial;

        /// <summary>The semi-duration of the total phase in minutes, or 0.0 if none.</summary>
        public double sd_total;

        internal LunarEclipseInfo(EclipseKind kind, double obscuration, AstroTime peak, double sd_penum, double sd_partial, double sd_total)
        {
            this.kind = kind;
            this.obscuration = obscuration;
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
    ///
    /// For total or annular eclipses, the `obscuration` field holds the fraction (0, 1]
    /// of the Sun's apparent disc area that is blocked from view by the Moon's silhouette,
    /// as seen by an observer located at the geographic coordinates `latitude`, `longitude`
    /// at the darkest time `peak`. The value will always be 1 for total eclipses, and less than
    /// 1 for annular eclipses.
    /// For partial eclipses, `obscuration` is undefined and should not be used.
    /// This is because there is little practical use for an obscuration value of
    /// a partial eclipse without supplying a particular observation location.
    /// Developers who wish to find an obscuration value for partial solar eclipses should therefore use
    /// #Astronomy.SearchLocalSolarEclipse and provide the geographic coordinates of an observer.
    /// </remarks>
    public struct GlobalSolarEclipseInfo
    {
        /// <summary>The type of solar eclipse: `EclipseKind.Partial`, `EclipseKind.Annular`, or `EclipseKind.Total`.</summary>
        public EclipseKind kind;

        /// <summary>The peak fraction of the Sun's apparent disc area obscured by the Moon (total and annular eclipses only).</summary>
        public double obscuration;

        /// <summary>
        /// The date and time when the solar eclipse is at its darkest.
        /// This is the instant when the axis of the Moon's shadow cone passes closest to the Earth's center.
        /// </summary>
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
    /// the horizon.
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
    /// The `obscuration` field reports what fraction of the Sun's disc appears blocked
    /// by the Moon when viewed by the observer at the peak eclipse time.
    /// This is a value that ranges from 0 (no blockage) to 1 (total eclipse).
    /// The obscuration value will be between 0 and 1 for partial eclipses and annular eclipses.
    /// The value will be exactly 1 for total eclipses. Obscuration gives an indication
    /// of how dark the eclipse appears.
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

        /// <summary>The fraction of the Sun's apparent disc area obscured by the Moon at the eclipse peak.</summary>
        public double obscuration;

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

        /// <summary>A value in the range [0.0, 1.0] indicating what fraction of the body's apparent disc is illuminated, as seen from the Earth.</summary>
        public readonly double phase_fraction;

        /// <summary>The distance between the Sun and the body at the observation time.</summary>
        public readonly double helio_dist;

        /// <summary>For Saturn, the tilt angle in degrees of its rings as seen from Earth. For all other bodies, 0.</summary>
        public readonly double ring_tilt;

        internal IllumInfo(AstroTime time, double mag, double phase_angle, double helio_dist, double ring_tilt)
        {
            this.time = time;
            this.mag = mag;
            this.phase_angle = phase_angle;
            this.phase_fraction = (1.0 + Math.Cos(Astronomy.DEG2RAD * phase_angle)) / 2.0;
            this.helio_dist = helio_dist;
            this.ring_tilt = ring_tilt;
        }
    }

    /// <summary>
    /// Information about a body's rotation axis at a given time.
    /// </summary>
    /// <remarks>
    /// This structure is returned by #Astronomy.RotationAxis to report
    /// the orientation of a body's rotation axis at a given moment in time.
    /// The axis is specified by the direction in space that the body's north pole
    /// points, using angular equatorial coordinates in the J2000 system (EQJ).
    ///
    /// Thus `ra` is the right ascension, and `dec` is the declination, of the
    /// body's north pole vector at the given moment in time. The north pole
    /// of a body is defined as the pole that lies on the north side of the
    /// [Solar System's invariable plane](https://en.wikipedia.org/wiki/Invariable_plane),
    /// regardless of the body's direction of rotation.
    ///
    /// The `spin` field indicates the angular position of a prime meridian
    /// arbitrarily recommended for the body by the International Astronomical
    /// Union (IAU).
    ///
    /// The fields `ra`, `dec`, and `spin` correspond to the variables
    /// α0, δ0, and W, respectively, from
    /// [Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015](https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf).
    ///
    /// The field `north` is a unit vector pointing in the direction of the body's north pole.
    /// It is expressed in the J2000 mean equator system (EQJ).
    /// </remarks>
    public struct AxisInfo
    {
        /// <summary>The J2000 right ascension of the body's north pole direction, in sidereal hours.</summary>
        public double ra;

        /// <summary>The J2000 declination of the body's north pole direction, in degrees.</summary>
        public double dec;

        /// <summary>Rotation angle of the body's prime meridian, in degrees.</summary>
        public double spin;

        /// <summary>A J2000 dimensionless unit vector pointing in the direction of the body's north pole.</summary>
        public AstroVector north;
    }

    /// <summary>
    /// Indicates whether a crossing through the ecliptic plane is ascending or descending.
    /// </summary>
    public enum NodeEventKind
    {
        /// <summary>Placeholder value for a missing or invalid node.</summary>
        Invalid = 0,

        /// <summary>The body passes through the ecliptic plane from south to north.</summary>
        Ascending = +1,

        /// <summary>The body passes through the ecliptic plane from north to south.</summary>
        Descending = -1,
    }

    /// <summary>
    /// Information about an ascending or descending node of a body.
    /// </summary>
    /// <remarks>
    /// This structure is returned by #Astronomy.SearchMoonNode and #Astronomy.NextMoonNode
    /// to report information about the center of the Moon passing through the ecliptic plane.
    /// </remarks>
    public struct NodeEventInfo
    {
        /// <summary>The time when the body passes through the ecliptic plane.</summary>
        public AstroTime        time;

        /// <summary>Whether the node is ascending (south to north) or descending (north to south).</summary>
        public NodeEventKind    kind;
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
            // The Search() function finds a transition from negative to positive values.
            // The derivative of magnitude y with respect to time t (dy/dt)
            // is negative as an object gets brighter, because the magnitude numbers
            // get smaller. At peak magnitude dy/dt = 0, then as the object gets dimmer,
            // dy/dt > 0.
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

    internal class SearchContext_MoonNode: SearchContext
    {
        public double Direction;

        public override double Eval(AstroTime time)
        {
            Spherical moon = Astronomy.EclipticGeoMoon(time);
            return Direction * moon.lat;
        }
    }

    internal class SearchContext_Altitude: SearchContext
    {
        private readonly Body body;
        private readonly int direction;
        private readonly Observer observer;
        private readonly double bodyRadiusAu;
        private readonly double targetAltitude;

        public SearchContext_Altitude(Body body, Direction direction, Observer observer, double bodyRadiusAu, double targetAltitude)
        {
            this.body = body;
            this.direction = (int)direction;
            this.observer = observer;
            this.bodyRadiusAu = bodyRadiusAu;
            this.targetAltitude = targetAltitude;
        }

        public override double Eval(AstroTime time)
        {
            Equatorial ofdate = Astronomy.Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
            Topocentric hor = Astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None);
            double altitude = hor.altitude + Astronomy.RAD2DEG*Math.Asin(bodyRadiusAu / ofdate.dist);
            return direction*(altitude - targetAltitude);
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
            // sine, of phi in revolutions, not radians
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

        void Term(int p, int q, int r, int s, out double x, out double y)
        {
            x = 1.0;
            y = 0.0;
            if (p != 0) AddThe(x, y, CO[p, 1], SI[p, 1], out x, out y);
            if (q != 0) AddThe(x, y, CO[q, 2], SI[q, 2], out x, out y);
            if (r != 0) AddThe(x, y, CO[r, 3], SI[r, 3], out x, out y);
            if (s != 0) AddThe(x, y, CO[s, 4], SI[s, 4], out x, out y);
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

        void ADDN(double coeffn, int p, int q, int r, int s)
        {
            double x, y;
            Term(p, q, r, s, out x, out y);
            N += coeffn * y;
        }

        void SolarN()
        {
            N = 0.0;
            ADDN(-526.069,  0, 0, 1, -2);
            ADDN(  -3.352,  0, 0, 1, -4);
            ADDN( +44.297, +1, 0, 1, -2);
            ADDN(  -6.000, +1, 0, 1, -4);
            ADDN( +20.599, -1, 0, 1,  0);
            ADDN( -30.598, -1, 0, 1, -2);
            ADDN( -24.649, -2, 0, 1,  0);
            ADDN(  -2.000, -2, 0, 1, -2);
            ADDN( -22.571,  0,+1, 1, -2);
            ADDN( +10.985,  0,-1, 1, -2);
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
            SINPI = 3422.7;
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
    /// A simulation of zero or more small bodies moving through the Solar System.
    /// </summary>
    /// <remarks>
    /// This class calculates the movement of arbitrary small bodies,
    /// such as asteroids or comets, that move through the Solar System.
    /// It does so by calculating the gravitational forces on the small bodies
    /// from the Sun and planets. The user of this class supplies an enumeration
    /// of initial positions and velocities for the small bodies.
    /// Then the class can update the positions and velocities over small time steps.
    /// </remarks>
    public class GravitySimulator
    {
        /// <summary>
        /// The origin of the reference frame. See constructor for more info.
        /// </summary>
        public readonly Body OriginBody;

        private GravSimEndpoint prev;
        private GravSimEndpoint curr;
        private const int GravitatorArraySize = 1 + (int)Body.Sun;
        private static readonly int[] PlanetIndexes = new int[] {
            (int)Body.Mercury,
            (int)Body.Venus,
            (int)Body.Earth,
            (int)Body.Mars,
            (int)Body.Jupiter,
            (int)Body.Saturn,
            (int)Body.Uranus,
            (int)Body.Neptune
        };

        /// <summary>Creates a gravity simulation object.</summary>
        /// <param name="originBody">
        /// Specifies the origin of the reference frame.
        /// All position vectors and velocity vectors will use `originBody`
        /// as the origin of the coordinate system.
        /// This origin applies to all the input vectors provided in the
        /// `bodyStates` parameter of this function, along with all
        /// output vectors returned by #GravitySimulator.Update.
        /// Most callers will want to provide one of the following:
        /// `Body.Sun` for heliocentric coordinates,
        /// `Body.SSB` for solar system barycentric coordinates,
        /// or `Body.Earth` for geocentric coordinates. Note that the
        /// gravity simulator does not correct for light travel time;
        /// all state vectors are tied to a Newtonian "instantaneous" time.
        /// </param>
        /// <param name="time">
        /// The initial time at which to start the simulation.
        /// </param>
        /// <param name="bodyStates">
        /// An enumeration of zero or more initial state vectors (positions and velocities)
        /// of the small bodies to be simulated.
        /// The caller must know the positions and velocities of the small bodies at an initial moment in time.
        /// Their positions and velocities are expressed with respect to `originBody`, using equatorial
        /// J2000 orientation (EQJ).
        /// Positions are expressed in astronomical units (AU).
        /// Velocities are expressed in AU/day.
        /// All the times embedded within the state vectors must be exactly equal to `time`,
        /// or this constructor will throw an exception.
        /// If `bodyStates` is null, the gravity simulator will contain zero small bodies.
        /// </param>
        public GravitySimulator(
            Body originBody,
            AstroTime time,
            IEnumerable<StateVector> bodyStates)
        {
            OriginBody = originBody;

            // Verify that all the state vectors have matching times.
            StateVector[] bodyStateArray = (bodyStates == null) ? new StateVector[0] : bodyStates.ToArray();
            foreach (StateVector b in bodyStateArray)
                if (b.t.tt != time.tt)
                    throw new ArgumentException("Inconsistent time(s) in bodyStates");

            prev = new GravSimEndpoint
            {
                time = time,
                gravitators = new body_state_t[GravitatorArraySize],
                bodies = new body_grav_calc_t[bodyStateArray.Length],
            };

            curr = new GravSimEndpoint
            {
                time = time,
                gravitators = new body_state_t[GravitatorArraySize],
                bodies = bodyStateArray.Select(b =>
                    new body_grav_calc_t(
                        time.tt,
                        new TerseVector(b.x, b.y, b.z),
                        new TerseVector(b.vx, b.vy, b.vz),
                        TerseVector.Zero
                    )
                ).ToArray(),
            };

            // Calculate the states of the Sun and planets.
            CalcSolarSystem();

            // We need to do all the physics calculations in barycentric coordinates.
            // But the caller provides the input vectors with respect to `originBody`.
            // Correct the input body state vectors for the specified origin.
            if (originBody != Body.SSB)
            {
                // Determine the barycentric state of the origin body.
                body_state_t ostate = InternalBodyState(originBody);

                // Add barycentric origin to origin-centric bodies to obtain barycentric bodies.
                for (int i = 0; i < curr.bodies.Length; ++i)
                {
                    curr.bodies[i].r += ostate.r;
                    curr.bodies[i].v += ostate.v;
                }
            }

            // Calculate the net acceleration experienced by the small bodies.
            CalcBodyAccelerations();

            // To prepare for a possible swap operation, duplicate the current state into the previous state.
            Duplicate();
        }

        /// <summary>
        /// The number of small bodies that are included in this gravity simulation.
        /// </summary>
        /// <remarks>
        /// #GravitySimulator.Update requres the caller to pass in an array to
        /// receive updated state vectors for the small bodies. This array must
        /// have the same number of elements as the bodies that are being simulated.
        /// `NumSmallBodies` returns this number as a convenience.
        /// </remarks>
        public int NumSmallBodies => curr.bodies.Length;

        /// <summary>
        /// The time represented by the current step of the gravity simulation.
        /// </summary>
        public AstroTime Time => curr.time;

        /// <summary>
        /// Advances a gravity simulation by a small time step.
        /// </summary>
        /// <remarks>
        /// Updates the simulation of the user-supplied small bodies
        /// to the time indicated by the `time` parameter.
        /// Updates the supplied array `bodyStates` of state vectors for the small bodies.
        /// This array must be the same size as the number of bodies supplied
        /// to the constructor of this object.
        /// The positions and velocities in the returned array are referenced
        /// to the `originBody` that was used to construct this simulator.
        /// </remarks>
        /// <param name="time">
        /// A time that is a small increment away from the current simulation time.
        /// It is up to the developer to figure out an appropriate time increment.
        /// Depending on the trajectories, a smaller or larger increment
        /// may be needed for the desired accuracy. Some experimentation may be needed.
        /// Generally, bodies that stay in the outer Solar System and move slowly can
        /// use larger time steps.  Bodies that pass into the inner Solar System and
        /// move faster will need a smaller time step to maintain accuracy.
        /// The `time` value may be after or before the current simulation time
        /// to move forward or backward in time.
        /// </param>
        /// <param name="bodyStates">
        /// If this array is not null, it must contain exactly the same number
        /// of elements as the number of small bodies that were added when this
        /// simulator was created. The non-null array receives updated state vectors
        /// for the simulated small bodies.
        /// If `bodyStates` is null, the simulation is updated but without returning
        /// the state vectors.
        /// </param>
        public void Update(AstroTime time, StateVector[] bodyStates)
        {
            int nbodies = NumSmallBodies;

            if (bodyStates != null && bodyStates.Length != nbodies)
                throw new ArgumentException($"This simulation contains {nbodies} small bodies, but the {nameof(bodyStates)} array has length {bodyStates.Length}. The array must either be null, or it must have the same number of elements.");

            double dt = time.tt - curr.time.tt;
            if (dt == 0.0)
            {
                // Special case: the time has not changed, so skip the usual physics calculations.
                // This allows another way for the caller to query the current body states.
                // It is also necessary to avoid dividing by `dt` if `dt` is zero.
                // To prepare for a possible swap operation, duplicate the current state into the previous state.
                Duplicate();
            }
            else
            {
                // Exchange the current state with the previous state. Then calculate the new current state.
                Swap();

                // Update the current time.
                curr.time = time;

                // Now that the time is set, it is safe to update the Solar System.
                CalcSolarSystem();

                // Estimate the positions of the small bodies as if their existing
                // accelerations apply across the whole time interval.
                for (int i = 0; i < nbodies; ++i)
                    curr.bodies[i].r = Astronomy.UpdatePosition(dt, prev.bodies[i].r, prev.bodies[i].v, prev.bodies[i].a);

                // Calculate the acceleration experienced by the small bodies at
                // their respective approximate next locations.
                CalcBodyAccelerations();

                for (int i = 0; i < nbodies; ++i)
                {
                    // Calculate the average of the acceleration vectors
                    // experienced by the previous body positions and
                    // their estimated next positions.
                    // These become estimates of the mean effective accelerations
                    // over the whole interval.
                    TerseVector acc = (curr.bodies[i].a + prev.bodies[i].a) / 2.0;

                    // Refine the estimates of position and velocity at the next time step,
                    // using the mean acceleration as a better approximation of the
                    // continuously changing acceleration acting on each body.
                    curr.bodies[i].tt = time.tt;
                    curr.bodies[i].r = Astronomy.UpdatePosition(dt, prev.bodies[i].r, prev.bodies[i].v, acc);
                    curr.bodies[i].v = Astronomy.UpdateVelocity(dt, prev.bodies[i].v, acc);
                }

                // Re-calculate accelerations experienced by each body.
                // These will be needed for the next simulation step (if any).
                // Also, they will be potentially useful if some day we add
                // a function to query the acceleration vectors for the bodies.
                CalcBodyAccelerations();
            }

            if (bodyStates != null)
            {
                // Translate our internal calculations of body positions and velocities
                // into state vectors that the caller can understand.
                // We have to convert the internal type body_grav_calc_t to the public
                // type StateVector.
                // Also convert from barycentric coordinates to coordinates based on the
                // selected origin body.
                body_state_t ostate = InternalBodyState(OriginBody);
                for (int i = 0; i < nbodies; ++i)
                {
                    bodyStates[i] = new StateVector(
                        curr.bodies[i].r.x - ostate.r.x,
                        curr.bodies[i].r.y - ostate.r.y,
                        curr.bodies[i].r.z - ostate.r.z,
                        curr.bodies[i].v.x - ostate.v.x,
                        curr.bodies[i].v.y - ostate.v.y,
                        curr.bodies[i].v.z - ostate.v.z,
                        time
                    );
                }
            }
        }

        /// <summary>
        /// Exchange the current time step with the previous time step.
        /// </summary>
        /// <remarks>
        /// Sometimes it is helpful to "explore" various times near a given
        /// simulation time step, while repeatedly returning to the original
        /// time step. For example, when backdating a position for light travel
        /// time, the caller may wish to repeatedly try different amounts of
        /// backdating. When the backdating solver has converged, the caller
        /// wants to leave the simulation in its original state.
        ///
        /// This function allows a single "undo" of a simulation, and does so
        /// very efficiently.
        ///
        /// Usually this function will be called immediately after a matching
        /// call to #GravitySimulator.Update. It has the effect of rolling
        /// back the most recent update. If called twice in a row, it reverts
        /// the swap and thus has no net effect.
        ///
        /// The constructor initializes the current state and previous
        /// state to be identical. Both states represent the `time` parameter that was
        /// passed into the constructor. Therefore, `Swap` will
        /// have no effect from the caller's point of view when passed a simulator
        /// that has not yet been updated by a call to #GravitySimulator.Update.
        /// </remarks>
        public void Swap()
        {
            var swap = curr;
            curr = prev;
            prev = swap;
        }

        /// <summary>
        /// Get the position and velocity of a Solar System body included in the simulation.
        /// </summary>
        /// <remarks>
        /// In order to simulate the movement of small bodies through the Solar System,
        /// the simulator needs to calculate the state vectors for the Sun and planets.
        ///
        /// If an application wants to know the positions of one or more of the planets
        /// in addition to the small bodies, this function provides a way to obtain
        /// their state vectors. This is provided for the sake of efficiency, to avoid
        /// redundant calculations.
        ///
        /// The state vector is returned relative to the position and velocity
        /// of the `originBody` parameter that was passed to this object's constructor.
        /// </remarks>
        ///
        /// <param name="body">
        /// The Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, or Neptune.
        /// </param>
        public StateVector SolarSystemBodyState(Body body)
        {
            body_state_t bstate = InternalBodyState(body);
            body_state_t ostate = InternalBodyState(OriginBody);
            return Astronomy.ExportState(bstate - ostate, curr.time);
        }

        private void CalcSolarSystem()
        {
            double tt = curr.time.tt;

            // Initialize the Sun's position/velocity as zero vectors, then adjust from pulls of the planets.
            var ssb = new body_state_t(tt, TerseVector.Zero, TerseVector.Zero);

            curr.gravitators[(int)Body.Mercury] = Astronomy.AdjustBarycenterPosVel(ref ssb, tt, Body.Mercury, Astronomy.MERCURY_GM);
            curr.gravitators[(int)Body.Venus  ] = Astronomy.AdjustBarycenterPosVel(ref ssb, tt, Body.Venus,   Astronomy.VENUS_GM);
            curr.gravitators[(int)Body.Earth  ] = Astronomy.AdjustBarycenterPosVel(ref ssb, tt, Body.Earth,   Astronomy.EARTH_GM + Astronomy.MOON_GM);
            curr.gravitators[(int)Body.Mars   ] = Astronomy.AdjustBarycenterPosVel(ref ssb, tt, Body.Mars,    Astronomy.MARS_GM);
            curr.gravitators[(int)Body.Jupiter] = Astronomy.AdjustBarycenterPosVel(ref ssb, tt, Body.Jupiter, Astronomy.JUPITER_GM);
            curr.gravitators[(int)Body.Saturn ] = Astronomy.AdjustBarycenterPosVel(ref ssb, tt, Body.Saturn,  Astronomy.SATURN_GM);
            curr.gravitators[(int)Body.Uranus ] = Astronomy.AdjustBarycenterPosVel(ref ssb, tt, Body.Uranus,  Astronomy.URANUS_GM);
            curr.gravitators[(int)Body.Neptune] = Astronomy.AdjustBarycenterPosVel(ref ssb, tt, Body.Neptune, Astronomy.NEPTUNE_GM);

            // Convert planets states from heliocentric to barycentric.
            foreach (int bindex in PlanetIndexes)
            {
                curr.gravitators[bindex].r -= ssb.r;
                curr.gravitators[bindex].v -= ssb.v;
            }

            // Convert heliocentric SSB to barycentric Sun.
            curr.gravitators[(int)Body.Sun] = new body_state_t(tt, -ssb.r, -ssb.v);
        }

        private body_state_t InternalBodyState(Body body)
        {
            if (body == Body.Sun || (body >= Body.Mercury && body <= Body.Neptune))
                return curr.gravitators[(int)body];

            if (body == Body.SSB)
                return new body_state_t(curr.time.tt, TerseVector.Zero, TerseVector.Zero);

            throw new InvalidBodyException(body);
        }

        private void CalcBodyAccelerations()
        {
            // Calculate the gravitational acceleration experienced by the simulated bodies.
            const double EMB_GM = Astronomy.EARTH_GM + Astronomy.MOON_GM;
            for (int i = 0; i < curr.bodies.Length; ++i)
            {
                TerseVector a = TerseVector.Zero;
                a += Acceleration(curr.bodies[i].r, curr.gravitators[(int)Body.Sun].r,      Astronomy.SUN_GM    );
                a += Acceleration(curr.bodies[i].r, curr.gravitators[(int)Body.Mercury].r,  Astronomy.MERCURY_GM);
                a += Acceleration(curr.bodies[i].r, curr.gravitators[(int)Body.Venus].r,    Astronomy.VENUS_GM  );
                a += Acceleration(curr.bodies[i].r, curr.gravitators[(int)Body.Earth].r,    EMB_GM              );
                a += Acceleration(curr.bodies[i].r, curr.gravitators[(int)Body.Mars].r,     Astronomy.MARS_GM   );
                a += Acceleration(curr.bodies[i].r, curr.gravitators[(int)Body.Jupiter].r,  Astronomy.JUPITER_GM);
                a += Acceleration(curr.bodies[i].r, curr.gravitators[(int)Body.Saturn].r,   Astronomy.SATURN_GM );
                a += Acceleration(curr.bodies[i].r, curr.gravitators[(int)Body.Uranus].r,   Astronomy.URANUS_GM );
                a += Acceleration(curr.bodies[i].r, curr.gravitators[(int)Body.Neptune].r,  Astronomy.NEPTUNE_GM);
                curr.bodies[i].a = a;
            }
        }

        private static TerseVector Acceleration(TerseVector smallPos, TerseVector majorPos, double gm)
        {
            double dx = majorPos.x - smallPos.x;
            double dy = majorPos.y - smallPos.y;
            double dz = majorPos.z - smallPos.z;
            double r2 = dx*dx + dy*dy + dz*dz;
            double pull = gm / (r2 * Math.Sqrt(r2));
            return new TerseVector(dx * pull, dy * pull, dz * pull);
        }

        private void Duplicate()
        {
            // Copy the current state into the previous state, so that both become the same moment in time.
            prev.time = curr.time;

            for (int i = 0; i < curr.gravitators.Length; ++i)
                prev.gravitators[i] = curr.gravitators[i];

            for (int i = 0; i < curr.bodies.Length; ++i)
                prev.bodies[i] = curr.bodies[i];
        }
    }

    internal class GravSimEndpoint
    {
        public AstroTime time;
        public body_state_t[] gravitators;
        public body_grav_calc_t[] bodies;
    }

    internal struct body_state_t
    {
        public double tt;       // Terrestrial Time in J2000 days
        public TerseVector r;   // position [au]
        public TerseVector v;   // velocity [au/day]

        public body_state_t(double tt, TerseVector r, TerseVector v)
        {
            this.tt = tt;
            this.r = r;
            this.v = v;
        }

        public static body_state_t operator -(body_state_t s)
        {
            return new body_state_t(s.tt, -s.r, -s.v);
        }

        public static body_state_t operator -(body_state_t a, body_state_t b)
        {
            return new body_state_t(a.tt, a.r - b.r, a.v - b.v);
        }
    }

    internal struct body_grav_calc_t
    {
        public double tt;       // J2000 terrestrial time [days]
        public TerseVector r;   // position [au]
        public TerseVector v;   // velocity [au/day]
        public TerseVector a;   // acceleration [au/day^2]

        public body_grav_calc_t(double tt, TerseVector r, TerseVector v, TerseVector a)
        {
            this.tt = tt;
            this.r = r;
            this.v = v;
            this.a = a;
        }
    }

    /// <summary>
    /// A function for which to solve a light-travel time problem.
    /// </summary>
    /// <remarks>
    /// The function #Astronomy.CorrectLightTravel solves a generalized
    /// problem of deducing how far in the past light must have left
    /// a target object to be seen by an observer at a specified time.
    /// This interface expresses an arbitrary position vector as
    /// function of time that is passed to #Astronomy.CorrectLightTravel.
    /// </remarks>
    public interface IPositionFunction
    {
        /// <summary>
        /// Returns a relative position vector for a given time.
        /// </summary>
        /// <param name="time">The time at which to evaluate a relative position vector.</param>
        AstroVector Position(AstroTime time);
    }

    /// <summary>
    /// The wrapper class that holds Astronomy Engine functions.
    /// </summary>
    public static class Astronomy
    {
        /// <summary>
        /// The number of kilometers in one astronomical unit (AU).
        /// </summary>
        public const double KM_PER_AU = 1.4959787069098932e+8;

        /// <summary>
        /// The factor to convert radians to degrees = 180/pi.
        /// </summary>
        public const double RAD2DEG = 57.295779513082321;

        /// <summary>
        /// The factor to convert radians to sidereal hours = 12/pi.
        /// </summary>
        public const double RAD2HOUR  = 3.819718634205488;

        /// <summary>
        /// The factor to convert degrees to radians = pi/180.
        /// </summary>
        public const double DEG2RAD = 0.017453292519943296;

        /// <summary>
        /// The factor to convert sidereal hours to radians = pi/12.
        /// </summary>
        public const double HOUR2RAD = 0.2617993877991494365;


        // Jupiter radius data are nominal values obtained from:
        // https://www.iau.org/static/resolutions/IAU2015_English.pdf
        // https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html

        /// <summary>
        /// The equatorial radius of Jupiter, expressed in kilometers.
        /// </summary>
        public const double JUPITER_EQUATORIAL_RADIUS_KM = 71492.0;

        /// <summary>
        /// The polar radius of Jupiter, expressed in kilometers.
        /// </summary>
        public const double JUPITER_POLAR_RADIUS_KM = 66854.0;

        /// <summary>
        /// The volumetric mean radius of Jupiter, expressed in kilometers.
        /// </summary>
        public const double JUPITER_MEAN_RADIUS_KM = 69911.0;

        // The radii of Jupiter's four major moons are obtained from:
        // https://ssd.jpl.nasa.gov/?sat_phys_par

        /// <summary>
        /// The mean radius of Jupiter's moon Io, expressed in kilometers.
        /// </summary>
        public const double IO_RADIUS_KM = 1821.6;

        /// <summary>
        /// The mean radius of Jupiter's moon Europa, expressed in kilometers.
        /// </summary>
        public const double EUROPA_RADIUS_KM = 1560.8;

        /// <summary>
        /// The mean radius of Jupiter's moon Ganymede, expressed in kilometers.
        /// </summary>
        public const double GANYMEDE_RADIUS_KM = 2631.2;

        /// <summary>
        /// The mean radius of Jupiter's moon Callisto, expressed in kilometers.
        /// </summary>
        public const double CALLISTO_RADIUS_KM = 2410.3;

        /// <summary>
        /// The speed of light in AU/day.
        /// </summary>
        public const double C_AUDAY = 173.1446326846693;

        /// <summary>
        /// The number of astronomical units in one light-year.
        /// </summary>
        public const double AU_PER_LY = 63241.07708807546;

        private const double DAYS_PER_TROPICAL_YEAR = 365.24217;
        private const double ASEC360 = 1296000.0;
        private const double ASEC2RAD = 4.848136811095359935899141e-6;
        internal const double PI2 = 2.0 * Math.PI;
        internal const double ARC = 3600.0 * 180.0 / Math.PI;       // arcseconds per radian

        internal const double SUN_RADIUS_KM  = 695700.0;
        internal const double SUN_RADIUS_AU  = SUN_RADIUS_KM / KM_PER_AU;

        internal const double EARTH_FLATTENING = 0.996647180302104;
        internal const double EARTH_EQUATORIAL_RADIUS_KM = 6378.1366;
        internal const double EARTH_EQUATORIAL_RADIUS_AU = EARTH_EQUATORIAL_RADIUS_KM / KM_PER_AU;
        internal const double EARTH_POLAR_RADIUS_KM = EARTH_EQUATORIAL_RADIUS_KM * EARTH_FLATTENING;
        internal const double EARTH_MEAN_RADIUS_KM = 6371.0;    // mean radius of the Earth's geoid, without atmosphere
        internal const double EARTH_ATMOSPHERE_KM = 88.0;       // effective atmosphere thickness for lunar eclipses
        internal const double EARTH_ECLIPSE_RADIUS_KM = EARTH_MEAN_RADIUS_KM + EARTH_ATMOSPHERE_KM;

        internal const double MOON_EQUATORIAL_RADIUS_KM = 1738.1;
        internal const double MOON_MEAN_RADIUS_KM       = 1737.4;
        internal const double MOON_POLAR_RADIUS_KM      = 1736.0;
        internal const double MOON_POLAR_RADIUS_AU      = (MOON_POLAR_RADIUS_KM / KM_PER_AU);
        internal const double MOON_EQUATORIAL_RADIUS_AU = (MOON_EQUATORIAL_RADIUS_KM / KM_PER_AU);

        private const double ANGVEL = 7.2921150e-5;
        private const double SECONDS_PER_DAY = 24.0 * 3600.0;
        private const double SOLAR_DAYS_PER_SIDEREAL_DAY = 0.9972695717592592;
        private const double MEAN_SYNODIC_MONTH = 29.530588;     // average number of days for Moon to return to the same phase
        private const double EARTH_ORBITAL_PERIOD = 365.256;
        private const double NEPTUNE_ORBITAL_PERIOD = 60189.0;
        internal const double REFRACTION_NEAR_HORIZON = 34.0 / 60.0;   //  degrees of refractive "lift" seen for objects near horizon
        private const double ASEC180 = 180.0 * 60.0 * 60.0;         // arcseconds per 180 degrees (or pi radians)
        private const double AU_PER_PARSEC = (ASEC180 / Math.PI);   // exact definition of how many AU = one parsec
        private const double EARTH_MOON_MASS_RATIO = 81.30056;

        //  Masses of the Sun and outer planets, used for:
        //  (1) Calculating the Solar System Barycenter
        //  (2) Integrating the movement of Pluto
        //
        //  https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
        //
        //  Page 10 in the above document describes the constants used in the DE405 ephemeris.
        //  The following are G*M values (gravity constant * mass) in [au^3 / day^2].
        //  This side-steps issues of not knowing the exact values of G and masses M[i];
        //  the products GM[i] are known extremely accurately.
        internal const double SUN_GM     = 0.2959122082855911e-03;
        internal const double MERCURY_GM = 0.4912547451450812e-10;
        internal const double VENUS_GM   = 0.7243452486162703e-09;
        internal const double EARTH_GM   = 0.8887692390113509e-09;
        internal const double MARS_GM    = 0.9549535105779258e-10;
        internal const double JUPITER_GM = 0.2825345909524226e-06;
        internal const double SATURN_GM  = 0.8459715185680659e-07;
        internal const double URANUS_GM  = 0.1292024916781969e-07;
        internal const double NEPTUNE_GM = 0.1524358900784276e-07;
        internal const double PLUTO_GM   = 0.2188699765425970e-11;

        internal const double MOON_GM = EARTH_GM / EARTH_MOON_MASS_RATIO;

        private static bool isfinite(double x)
        {
            return !double.IsNaN(x) && !double.IsInfinity(x);
        }

        internal static double hypot(double x, double y)
        {
            return Math.Sqrt(x*x + y*y);
        }

        internal static double hypot(double x, double y, double z)
        {
            return Math.Sqrt(x*x + y*y + z*z);
        }

        private class StarDef
        {
            public double ra;       // heliocentric right ascension in EQJ
            public double dec;      // heliocentric declination in EQJ
            public double dist;     // heliocentric distance in AU
        };

        private static readonly StarDef[] StarTable = InitStarTable();

        private static StarDef[] InitStarTable()
        {
            var table = new StarDef[8];
            for (int i = 0; i < table.Length; ++i)
                table[i] = new StarDef();
            return table;
        }

        private static StarDef GetStar(Body body) =>
            ((body >= Body.Star1) && (body <= Body.Star8)) ? StarTable[(int)body - (int)Body.Star1] : null;

        private static StarDef UserDefinedStar(Body body)
        {
            if (GetStar(body) is StarDef star)
                if (star.dist > 0.0)        // has the star been defined yet?
                    return star;

            return null;
        }

        /// <summary>
        /// Assign equatorial coordinates to a user-defined star.
        /// </summary>
        /// <remarks>
        /// Some Astronomy Engine functions allow their `body` parameter to
        /// be a user-defined fixed point in the sky, loosely called a "star".
        /// This function assigns a right ascension, declination, and distance
        /// to one of the eight user-defined stars `Body.Star1`..`Body.Star8`.
        ///
        /// Stars are not valid until defined. Once defined, they retain their
        /// definition until re-defined by another call to `DefineStar`.
        /// </remarks>
        /// <param name="body">
        /// One of the eight user-defined star identifiers: `Body.Star1`, `Body.Star2`, ..., `Body.Star8`.
        /// </param>
        /// <param name="ra">
        /// The right ascension to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ).
        /// The value is in units of sidereal hours, and must be within the half-open range [0, 24).
        /// </param>
        /// <param name="dec">
        /// The declination to be assigned to the star, expressed in J2000 equatorial coordinates (EQJ).
        /// The value is in units of degrees north (positive) or south (negative) of the J2000 equator,
        /// and must be within the closed range [-90, +90].
        /// </param>
        /// <param name="distanceLightYears">
        /// The distance between the star and the Sun, expressed in light-years.
        /// This value is used to calculate the tiny parallax shift as seen by an observer on Earth.
        /// If you don't know the distance to the star, using a large value like 1000 will generally work well.
        /// The minimum allowed distance is 1 light-year, which is required to provide certain internal optimizations.
        /// </param>
        public static void DefineStar(Body body, double ra, double dec, double distanceLightYears)
        {
            StarDef star = GetStar(body) ?? throw new InvalidBodyException(body);

            if (!isfinite(ra) || ra < 0.0 || ra >= 24.0)
                throw new ArgumentException($"Invalid right ascension for star: {ra}");

            if (!isfinite(dec) || dec < -90.0 || dec > +90.0)
                throw new ArgumentException($"Invalid declination for star: {dec}");

            if (!isfinite(distanceLightYears) || distanceLightYears < 1.0)
                throw new ArgumentException($"Invalid heliocentric distance for star: {distanceLightYears}");

            star.ra = ra;
            star.dec = dec;
            star.dist = distanceLightYears * AU_PER_LY;
        }

        /// <summary>
        /// Returns the product of mass and universal gravitational constant of a Solar System body.
        /// </summary>
        /// <remarks>
        /// For problems involving the gravitational interactions of Solar System bodies,
        /// it is helpful to know the product GM, where G = the universal gravitational constant
        /// and M = the mass of the body. In practice, GM is known to a higher precision than
        /// either G or M alone, and thus using the product results in the most accurate results.
        /// This function returns the product GM in the units au^3/day^2.
        /// The values come from page 10 of a
        /// [JPL memorandum regarding the DE405/LE405 ephemeris](https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf).
        /// </remarks>
        /// <param name="body">
        /// The body for which to find the GM product.
        /// Allowed to be the Sun, Moon, EMB (Earth/Moon Barycenter), or any planet.
        /// Any other value will cause an exception to be thrown.
        /// </param>
        /// <returns>The mass product of the given body in au^3/day^2.</returns>
        public static double MassProduct(Body body)
        {
            switch (body)
            {
                case Body.Sun:      return SUN_GM;
                case Body.Mercury:  return MERCURY_GM;
                case Body.Venus:    return VENUS_GM;
                case Body.Earth:    return EARTH_GM;
                case Body.Moon:     return MOON_GM;
                case Body.EMB:      return EARTH_GM + MOON_GM;
                case Body.Mars:     return MARS_GM;
                case Body.Jupiter:  return JUPITER_GM;
                case Body.Saturn:   return SATURN_GM;
                case Body.Uranus:   return URANUS_GM;
                case Body.Neptune:  return NEPTUNE_GM;
                case Body.Pluto:    return PLUTO_GM;
                default:
                    throw new InvalidBodyException(body);
            }
        }

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
            public vsop_formula_t lon;
            public vsop_formula_t lat;
            public vsop_formula_t rad;

            public vsop_model_t(vsop_series_t[] lon, vsop_series_t[] lat, vsop_series_t[] rad)
            {
                this.lon = new vsop_formula_t(lon);
                this.lat = new vsop_formula_t(lat);
                this.rad = new vsop_formula_t(rad);
            }
        };

        private static readonly vsop_term_t[] vsop_lon_Mercury_0 = new vsop_term_t[]
        {
            new vsop_term_t(4.40250710144, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.40989414977, 1.48302034195, 26087.90314157420),
            new vsop_term_t(0.05046294200, 4.47785489551, 52175.80628314840),
            new vsop_term_t(0.00855346844, 1.16520322459, 78263.70942472259),
            new vsop_term_t(0.00165590362, 4.11969163423, 104351.61256629678),
            new vsop_term_t(0.00034561897, 0.77930768443, 130439.51570787099),
            new vsop_term_t(0.00007583476, 3.71348404924, 156527.41884944518)
        };

        private static readonly vsop_term_t[] vsop_lon_Mercury_1 = new vsop_term_t[]
        {
            new vsop_term_t(26087.90313685529, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.01131199811, 6.21874197797, 26087.90314157420),
            new vsop_term_t(0.00292242298, 3.04449355541, 52175.80628314840),
            new vsop_term_t(0.00075775081, 6.08568821653, 78263.70942472259),
            new vsop_term_t(0.00019676525, 2.80965111777, 104351.61256629678)
        };

        private static readonly vsop_series_t[] vsop_lon_Mercury = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Mercury_0),
            new vsop_series_t(vsop_lon_Mercury_1)
        };

        private static readonly vsop_term_t[] vsop_lat_Mercury_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.11737528961, 1.98357498767, 26087.90314157420),
            new vsop_term_t(0.02388076996, 5.03738959686, 52175.80628314840),
            new vsop_term_t(0.01222839532, 3.14159265359, 0.00000000000),
            new vsop_term_t(0.00543251810, 1.79644363964, 78263.70942472259),
            new vsop_term_t(0.00129778770, 4.83232503958, 104351.61256629678),
            new vsop_term_t(0.00031866927, 1.58088495658, 130439.51570787099),
            new vsop_term_t(0.00007963301, 4.60972126127, 156527.41884944518)
        };

        private static readonly vsop_term_t[] vsop_lat_Mercury_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00274646065, 3.95008450011, 26087.90314157420),
            new vsop_term_t(0.00099737713, 3.14159265359, 0.00000000000)
        };

        private static readonly vsop_series_t[] vsop_lat_Mercury = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Mercury_0),
            new vsop_series_t(vsop_lat_Mercury_1)
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


        private static readonly vsop_term_t[] vsop_lon_Venus_0 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lon_Venus_1 = new vsop_term_t[]
        {
            new vsop_term_t(10213.28554621638, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00095617813, 2.46406511110, 10213.28554621100),
            new vsop_term_t(0.00007787201, 0.62478482220, 20426.57109242200)
        };

        private static readonly vsop_series_t[] vsop_lon_Venus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Venus_0),
            new vsop_series_t(vsop_lon_Venus_1)
        };

        private static readonly vsop_term_t[] vsop_lat_Venus_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.05923638472, 0.26702775812, 10213.28554621100),
            new vsop_term_t(0.00040107978, 1.14737178112, 20426.57109242200),
            new vsop_term_t(0.00032814918, 3.14159265359, 0.00000000000)
        };

        private static readonly vsop_term_t[] vsop_lat_Venus_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00287821243, 1.88964962838, 10213.28554621100)
        };

        private static readonly vsop_series_t[] vsop_lat_Venus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Venus_0),
            new vsop_series_t(vsop_lat_Venus_1)
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


        private static readonly vsop_term_t[] vsop_lon_Earth_0 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lon_Earth_1 = new vsop_term_t[]
        {
            new vsop_term_t(6283.07584999140, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00206058863, 2.67823455808, 6283.07584999140),
            new vsop_term_t(0.00004303419, 2.63512233481, 12566.15169998280)
        };

        private static readonly vsop_term_t[] vsop_lon_Earth_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00008721859, 1.07253635559, 6283.07584999140)
        };

        private static readonly vsop_series_t[] vsop_lon_Earth = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Earth_0),
            new vsop_series_t(vsop_lon_Earth_1),
            new vsop_series_t(vsop_lon_Earth_2)
        };

        private static readonly vsop_term_t[] vsop_lat_Earth_0 = new vsop_term_t[]
        {
        };

        private static readonly vsop_term_t[] vsop_lat_Earth_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00227777722, 3.41376620530, 6283.07584999140),
            new vsop_term_t(0.00003805678, 3.37063423795, 12566.15169998280)
        };

        private static readonly vsop_series_t[] vsop_lat_Earth = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Earth_0),
            new vsop_series_t(vsop_lat_Earth_1)
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


        private static readonly vsop_term_t[] vsop_lon_Mars_0 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lon_Mars_1 = new vsop_term_t[]
        {
            new vsop_term_t(3340.61242700512, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.01457554523, 3.60433733236, 3340.61242669980),
            new vsop_term_t(0.00168414711, 3.92318567804, 6681.22485339960),
            new vsop_term_t(0.00020622975, 4.26108844583, 10021.83728009940),
            new vsop_term_t(0.00003452392, 4.73210393190, 3.52311834900),
            new vsop_term_t(0.00002586332, 4.60670058555, 13362.44970679920),
            new vsop_term_t(0.00000841535, 4.45864030426, 2281.23049651060)
        };

        private static readonly vsop_term_t[] vsop_lon_Mars_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00058152577, 2.04961712429, 3340.61242669980),
            new vsop_term_t(0.00013459579, 2.45738706163, 6681.22485339960)
        };

        private static readonly vsop_series_t[] vsop_lon_Mars = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Mars_0),
            new vsop_series_t(vsop_lon_Mars_1),
            new vsop_series_t(vsop_lon_Mars_2)
        };

        private static readonly vsop_term_t[] vsop_lat_Mars_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.03197134986, 3.76832042431, 3340.61242669980),
            new vsop_term_t(0.00298033234, 4.10616996305, 6681.22485339960),
            new vsop_term_t(0.00289104742, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00031365539, 4.44651053090, 10021.83728009940),
            new vsop_term_t(0.00003484100, 4.78812549260, 13362.44970679920)
        };

        private static readonly vsop_term_t[] vsop_lat_Mars_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00217310991, 6.04472194776, 3340.61242669980),
            new vsop_term_t(0.00020976948, 3.14159265359, 0.00000000000),
            new vsop_term_t(0.00012834709, 1.60810667915, 6681.22485339960)
        };

        private static readonly vsop_series_t[] vsop_lat_Mars = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Mars_0),
            new vsop_series_t(vsop_lat_Mars_1)
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


        private static readonly vsop_term_t[] vsop_lon_Jupiter_0 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lon_Jupiter_1 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lon_Jupiter_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00047233601, 4.32148536482, 7.11354700080),
            new vsop_term_t(0.00030649436, 2.92977788700, 529.69096509460),
            new vsop_term_t(0.00014837605, 3.14159265359, 0.00000000000)
        };

        private static readonly vsop_series_t[] vsop_lon_Jupiter = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Jupiter_0),
            new vsop_series_t(vsop_lon_Jupiter_1),
            new vsop_series_t(vsop_lon_Jupiter_2)
        };

        private static readonly vsop_term_t[] vsop_lat_Jupiter_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.02268615702, 3.55852606721, 529.69096509460),
            new vsop_term_t(0.00109971634, 3.90809347197, 1059.38193018920),
            new vsop_term_t(0.00110090358, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00008101428, 3.60509572885, 522.57741809380),
            new vsop_term_t(0.00006043996, 4.25883108339, 1589.07289528380),
            new vsop_term_t(0.00006437782, 0.30627119215, 536.80451209540)
        };

        private static readonly vsop_term_t[] vsop_lat_Jupiter_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00078203446, 1.52377859742, 529.69096509460)
        };

        private static readonly vsop_series_t[] vsop_lat_Jupiter = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Jupiter_0),
            new vsop_series_t(vsop_lat_Jupiter_1)
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


        private static readonly vsop_term_t[] vsop_lon_Saturn_0 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lon_Saturn_1 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lon_Saturn_2 = new vsop_term_t[]
        {
            new vsop_term_t(0.00116441330, 1.17988132879, 7.11354700080),
            new vsop_term_t(0.00091841837, 0.07325195840, 213.29909543800),
            new vsop_term_t(0.00036661728, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00015274496, 4.06493179167, 206.18554843720)
        };

        private static readonly vsop_series_t[] vsop_lon_Saturn = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Saturn_0),
            new vsop_series_t(vsop_lon_Saturn_1),
            new vsop_series_t(vsop_lon_Saturn_2)
        };

        private static readonly vsop_term_t[] vsop_lat_Saturn_0 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lat_Saturn_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00198927992, 4.93901017903, 213.29909543800),
            new vsop_term_t(0.00036947916, 3.14159265359, 0.00000000000),
            new vsop_term_t(0.00017966989, 0.51979431110, 426.59819087600)
        };

        private static readonly vsop_series_t[] vsop_lat_Saturn = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Saturn_0),
            new vsop_series_t(vsop_lat_Saturn_1)
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


        private static readonly vsop_term_t[] vsop_lon_Uranus_0 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lon_Uranus_1 = new vsop_term_t[]
        {
            new vsop_term_t(74.78159860910, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00154332863, 5.24158770553, 74.78159856730),
            new vsop_term_t(0.00024456474, 1.71260334156, 1.48447270830),
            new vsop_term_t(0.00009258442, 0.42829732350, 11.04570026390),
            new vsop_term_t(0.00008265977, 1.50218091379, 63.73589830340),
            new vsop_term_t(0.00009150160, 1.41213765216, 149.56319713460)
        };

        private static readonly vsop_series_t[] vsop_lon_Uranus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Uranus_0),
            new vsop_series_t(vsop_lon_Uranus_1)
        };

        private static readonly vsop_term_t[] vsop_lat_Uranus_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.01346277648, 2.61877810547, 74.78159856730),
            new vsop_term_t(0.00062341400, 5.08111189648, 149.56319713460),
            new vsop_term_t(0.00061601196, 3.14159265359, 0.00000000000),
            new vsop_term_t(0.00009963722, 1.61603805646, 76.26607127560),
            new vsop_term_t(0.00009926160, 0.57630380333, 73.29712585900)
        };

        private static readonly vsop_term_t[] vsop_lat_Uranus_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00034101978, 0.01321929936, 74.78159856730)
        };

        private static readonly vsop_series_t[] vsop_lat_Uranus = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Uranus_0),
            new vsop_series_t(vsop_lat_Uranus_1)
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


        private static readonly vsop_term_t[] vsop_lon_Neptune_0 = new vsop_term_t[]
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

        private static readonly vsop_term_t[] vsop_lon_Neptune_1 = new vsop_term_t[]
        {
            new vsop_term_t(38.13303563957, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00016604172, 4.86323329249, 1.48447270830),
            new vsop_term_t(0.00015744045, 2.27887427527, 38.13303563780)
        };

        private static readonly vsop_series_t[] vsop_lon_Neptune = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lon_Neptune_0),
            new vsop_series_t(vsop_lon_Neptune_1)
        };

        private static readonly vsop_term_t[] vsop_lat_Neptune_0 = new vsop_term_t[]
        {
            new vsop_term_t(0.03088622933, 1.44104372644, 38.13303563780),
            new vsop_term_t(0.00027780087, 5.91271884599, 76.26607127560),
            new vsop_term_t(0.00027623609, 0.00000000000, 0.00000000000),
            new vsop_term_t(0.00015355489, 2.52123799551, 36.64856292950),
            new vsop_term_t(0.00015448133, 3.50877079215, 39.61750834610)
        };

        private static readonly vsop_series_t[] vsop_lat_Neptune = new vsop_series_t[]
        {
            new vsop_series_t(vsop_lat_Neptune_0)
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
            new vsop_model_t(vsop_lon_Mercury,  vsop_lat_Mercury,   vsop_rad_Mercury),
            new vsop_model_t(vsop_lon_Venus,    vsop_lat_Venus,     vsop_rad_Venus  ),
            new vsop_model_t(vsop_lon_Earth,    vsop_lat_Earth,     vsop_rad_Earth  ),
            new vsop_model_t(vsop_lon_Mars,     vsop_lat_Mars,      vsop_rad_Mars   ),
            new vsop_model_t(vsop_lon_Jupiter,  vsop_lat_Jupiter,   vsop_rad_Jupiter),
            new vsop_model_t(vsop_lon_Saturn,   vsop_lat_Saturn,    vsop_rad_Saturn ),
            new vsop_model_t(vsop_lon_Uranus,   vsop_lat_Uranus,    vsop_rad_Uranus ),
            new vsop_model_t(vsop_lon_Neptune,  vsop_lat_Neptune,   vsop_rad_Neptune)
        };

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
            // Fred Espenak writes about Delta-T generically here:
            // https://eclipse.gsfc.nasa.gov/SEhelp/deltaT.html
            // https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html
            //
            // He provides polynomial approximations for distant years here:
            // https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
            //
            // They start with a year value 'y' such that y=2000 corresponds
            // to the UTC Date 15-January-2000. Convert difference in days
            // to mean tropical years.
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

            // all years after 2150
            u = (y - 1820) / 100;
            return -20 + (32 * u*u);
        }

        private static DeltaTimeFunc DeltaT = DeltaT_EspenakMeeus;

        internal static double TerrestrialTime(double ut)
        {
            return ut + DeltaT(ut)/86400.0;
        }

        internal static double UniversalTime(double tt)
        {
            // This is the inverse function of TerrestrialTime.
            // This is an iterative numerical solver, but because
            // the relationship between UT and TT is almost perfectly linear,
            // it converges extremely fast (never more than 3 iterations).

            // dt = tt - ut
            double dt = TerrestrialTime(tt) - tt;
            for(;;)
            {
                double ut = tt - dt;
                double tt_check = TerrestrialTime(ut);
                double err = tt_check - tt;
                if (Math.Abs(err) < 1.0e-12)
                    return ut;
                dt += err;
            }
        }

        private static double VsopFormulaCalc(vsop_formula_t formula, double t, bool clamp_angle)
        {
            double coord = 0.0;
            double tpower = 1.0;
            foreach (vsop_series_t series in formula.series)
            {
                double sum = 0.0;
                foreach (vsop_term_t term in series.term)
                    sum += term.amplitude * Math.Cos(term.phase + (t * term.frequency));
                double incr = tpower * sum;
                if (clamp_angle)
                    incr %= PI2;    // improve precision: longitude angles can be hundreds of radians
                coord += incr;
                tpower *= t;
            }
            return coord;
        }

        private static TerseVector VsopRotate(TerseVector eclip)
        {
            return new TerseVector(
                eclip.x + 0.000000440360*eclip.y - 0.000000190919*eclip.z,
                -0.000000479966*eclip.x + 0.917482137087*eclip.y - 0.397776982902*eclip.z,
                0.397776982902*eclip.y + 0.917482137087*eclip.z
            );
        }

        private static TerseVector VsopSphereToRect(double lon, double lat, double radius)
        {
            double r_coslat = radius * Math.Cos(lat);
            return new TerseVector(
                r_coslat * Math.Cos(lon),
                r_coslat * Math.Sin(lon),
                radius * Math.Sin(lat)
            );
        }

        private const double DAYS_PER_MILLENNIUM = 365250.0;

        private static AstroVector CalcVsop(vsop_model_t model, AstroTime time)
        {
            double t = time.tt / DAYS_PER_MILLENNIUM;    // millennia since 2000

            // Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates.
            double lon = VsopFormulaCalc(model.lon, t, true);
            double lat = VsopFormulaCalc(model.lat, t, false);
            double rad = VsopFormulaCalc(model.rad, t, false);

            // Convert ecliptic spherical coordinates to ecliptic Cartesian coordinates.
            TerseVector eclip = VsopSphereToRect(lon, lat, rad);

            // Convert ecliptic Cartesian coordinates to equatorial Cartesian coordinates.
            return VsopRotate(eclip).ToAstroVector(time);
        }

        private static double VsopDerivCalc(vsop_formula_t formula, double t)
        {
            double tpower = 1.0;        // t^s
            double dpower = 0.0;        // t^(s-1)
            double deriv = 0.0;
            for (int s=0; s < formula.series.Length; ++s)
            {
                double sin_sum = 0.0;
                double cos_sum = 0.0;
                vsop_series_t series = formula.series[s];
                foreach (vsop_term_t term in series.term)
                {
                    double angle = term.phase + (t * term.frequency);
                    sin_sum += term.amplitude * term.frequency * Math.Sin(angle);
                    if (s > 0)
                        cos_sum += term.amplitude * Math.Cos(angle);
                }
                deriv += (s * dpower * cos_sum) - (tpower * sin_sum);
                dpower = tpower;
                tpower *= t;
            }
            return deriv;
        }

        private struct major_bodies_t
        {
            public body_state_t Sun;
            public body_state_t Jupiter;
            public body_state_t Saturn;
            public body_state_t Uranus;
            public body_state_t Neptune;

            public TerseVector Acceleration(TerseVector small_pos)
            {
                // Use barycentric coordinates of the Sun and major planets to calculate
                // the gravitational acceleration vector experienced by a small body at location 'small_pos'.
                return
                    AccelerationIncrement(small_pos, SUN_GM,      Sun.r) +
                    AccelerationIncrement(small_pos, JUPITER_GM,  Jupiter.r) +
                    AccelerationIncrement(small_pos, SATURN_GM,   Saturn.r) +
                    AccelerationIncrement(small_pos, URANUS_GM,   Uranus.r) +
                    AccelerationIncrement(small_pos, NEPTUNE_GM,  Neptune.r);
            }

            private static TerseVector AccelerationIncrement(TerseVector small_pos, double gm, TerseVector major_pos)
            {
                TerseVector delta = major_pos - small_pos;
                double r2 = delta.Quadrature();
                return (gm / (r2 * Math.Sqrt(r2))) * delta;
            }
        }

        private static body_state_t CalcVsopPosVel(vsop_model_t model, double tt)
        {
            double t = tt / DAYS_PER_MILLENNIUM;    // millennia since 2000

            // Calculate the VSOP "B" trigonometric series to obtain ecliptic spherical coordinates.
            double lon = VsopFormulaCalc(model.lon, t, true);
            double lat = VsopFormulaCalc(model.lat, t, false);
            double rad = VsopFormulaCalc(model.rad, t, false);

            TerseVector eclip_pos = VsopSphereToRect(lon, lat, rad);

            double dlon_dt = VsopDerivCalc(model.lon, t);
            double dlat_dt = VsopDerivCalc(model.lat, t);
            double drad_dt = VsopDerivCalc(model.rad, t);

            // Use spherical coords and spherical derivatives to calculate
            // the velocity vector in rectangular coordinates.

            double coslon = Math.Cos(lon);
            double sinlon = Math.Sin(lon);
            double coslat = Math.Cos(lat);
            double sinlat = Math.Sin(lat);

            double vx =
                + (drad_dt * coslat * coslon)
                - (rad * sinlat * coslon * dlat_dt)
                - (rad * coslat * sinlon * dlon_dt);

            double vy =
                + (drad_dt * coslat * sinlon)
                - (rad * sinlat * sinlon * dlat_dt)
                + (rad * coslat * coslon * dlon_dt);

            double vz =
                + (drad_dt * sinlat)
                + (rad * coslat * dlat_dt);

            // Convert speed units from [AU/millennium] to [AU/day].
            var eclip_vel = new TerseVector(
                vx / DAYS_PER_MILLENNIUM,
                vy / DAYS_PER_MILLENNIUM,
                vz / DAYS_PER_MILLENNIUM);

            // Rotate the vectors from ecliptic to equatorial coordinates.
            TerseVector equ_pos = VsopRotate(eclip_pos);
            TerseVector equ_vel = VsopRotate(eclip_vel);
            return new body_state_t(tt, equ_pos, equ_vel);
        }

#region Pluto

        private const int PLUTO_NUM_STATES = 51;
        private const int PLUTO_TIME_STEP  = 29200;
        private const int PLUTO_DT         = 146;
        private const int PLUTO_NSTEPS     = 201;

        private static readonly body_state_t[] PlutoStateTable = new body_state_t[]
        {
            new body_state_t( -730000.0, new TerseVector(-26.118207232108, -14.376168177825,   3.384402515299), new TerseVector( 1.6339372163656e-03, -2.7861699588508e-03, -1.3585880229445e-03))
        ,   new body_state_t( -700800.0, new TerseVector( 41.974905202127,  -0.448502952929, -12.770351505989), new TerseVector( 7.3458569351457e-04,  2.2785014891658e-03,  4.8619778602049e-04))
        ,   new body_state_t( -671600.0, new TerseVector( 14.706930780744,  44.269110540027,   9.353698474772), new TerseVector(-2.1000147999800e-03,  2.2295915939915e-04,  7.0143443551414e-04))
        ,   new body_state_t( -642400.0, new TerseVector(-29.441003929957,  -6.430161530570,   6.858481011305), new TerseVector( 8.4495803960544e-04, -3.0783914758711e-03, -1.2106305981192e-03))
        ,   new body_state_t( -613200.0, new TerseVector( 39.444396946234,  -6.557989760571, -13.913760296463), new TerseVector( 1.1480029005873e-03,  2.2400006880665e-03,  3.5168075922288e-04))
        ,   new body_state_t( -584000.0, new TerseVector( 20.230380950700,  43.266966657189,   7.382966091923), new TerseVector(-1.9754081700585e-03,  5.3457141292226e-04,  7.5929169129793e-04))
        ,   new body_state_t( -554800.0, new TerseVector(-30.658325364620,   2.093818874552,   9.880531138071), new TerseVector( 6.1010603013347e-05, -3.1326500935382e-03, -9.9346125151067e-04))
        ,   new body_state_t( -525600.0, new TerseVector( 35.737703251673, -12.587706024764, -14.677847247563), new TerseVector( 1.5802939375649e-03,  2.1347678412429e-03,  1.9074436384343e-04))
        ,   new body_state_t( -496400.0, new TerseVector( 25.466295188546,  41.367478338417,   5.216476873382), new TerseVector(-1.8054401046468e-03,  8.3283083599510e-04,  8.0260156912107e-04))
        ,   new body_state_t( -467200.0, new TerseVector(-29.847174904071,  10.636426313081,  12.297904180106), new TerseVector(-6.3257063052907e-04, -2.9969577578221e-03, -7.4476074151596e-04))
        ,   new body_state_t( -438000.0, new TerseVector( 30.774692107687, -18.236637015304, -14.945535879896), new TerseVector( 2.0113162005465e-03,  1.9353827024189e-03, -2.0937793168297e-06))
        ,   new body_state_t( -408800.0, new TerseVector( 30.243153324028,  38.656267888503,   2.938501750218), new TerseVector(-1.6052508674468e-03,  1.1183495337525e-03,  8.3333973416824e-04))
        ,   new body_state_t( -379600.0, new TerseVector(-27.288984772533,  18.643162147874,  14.023633623329), new TerseVector(-1.1856388898191e-03, -2.7170609282181e-03, -4.9015526126399e-04))
        ,   new body_state_t( -350400.0, new TerseVector( 24.519605196774, -23.245756064727, -14.626862367368), new TerseVector( 2.4322321483154e-03,  1.6062008146048e-03, -2.3369181613312e-04))
        ,   new body_state_t( -321200.0, new TerseVector( 34.505274805875,  35.125338586954,   0.557361475637), new TerseVector(-1.3824391637782e-03,  1.3833397561817e-03,  8.4823598806262e-04))
        ,   new body_state_t( -292000.0, new TerseVector(-23.275363915119,  25.818514298769,  15.055381588598), new TerseVector(-1.6062295460975e-03, -2.3395961498533e-03, -2.4377362639479e-04))
        ,   new body_state_t( -262800.0, new TerseVector( 17.050384798092, -27.180376290126, -13.608963321694), new TerseVector( 2.8175521080578e-03,  1.1358749093955e-03, -4.9548725258825e-04))
        ,   new body_state_t( -233600.0, new TerseVector( 38.093671910285,  30.880588383337,  -1.843688067413), new TerseVector(-1.1317697153459e-03,  1.6128814698472e-03,  8.4177586176055e-04))
        ,   new body_state_t( -204400.0, new TerseVector(-18.197852930878,  31.932869934309,  15.438294826279), new TerseVector(-1.9117272501813e-03, -1.9146495909842e-03, -1.9657304369835e-05))
        ,   new body_state_t( -175200.0, new TerseVector(  8.528924039997, -29.618422200048, -11.805400994258), new TerseVector( 3.1034370787005e-03,  5.1393633292430e-04, -7.7293066202546e-04))
        ,   new body_state_t( -146000.0, new TerseVector( 40.946857258640,  25.904973592021,  -4.256336240499), new TerseVector(-8.3652705194051e-04,  1.8129497136404e-03,  8.1564228273060e-04))
        ,   new body_state_t( -116800.0, new TerseVector(-12.326958895325,  36.881883446292,  15.217158258711), new TerseVector(-2.1166103705038e-03, -1.4814420035990e-03,  1.7401209844705e-04))
        ,   new body_state_t(  -87600.0, new TerseVector( -0.633258375909, -30.018759794709,  -9.171932874950), new TerseVector( 3.2016994581737e-03, -2.5279858672148e-04, -1.0411088271861e-03))
        ,   new body_state_t(  -58400.0, new TerseVector( 42.936048423883,  20.344685584452,  -6.588027007912), new TerseVector(-5.0525450073192e-04,  1.9910074335507e-03,  7.7440196540269e-04))
        ,   new body_state_t(  -29200.0, new TerseVector( -5.975910552974,  40.611809958460,  14.470131723673), new TerseVector(-2.2184202156107e-03, -1.0562361130164e-03,  3.3652250216211e-04))
        ,   new body_state_t(       0.0, new TerseVector( -9.875369580774, -27.978926224737,  -5.753711824704), new TerseVector( 3.0287533248818e-03, -1.1276087003636e-03, -1.2651326732361e-03))
        ,   new body_state_t(   29200.0, new TerseVector( 43.958831986165,  14.214147973292,  -8.808306227163), new TerseVector(-1.4717608981871e-04,  2.1404187242141e-03,  7.1486567806614e-04))
        ,   new body_state_t(   58400.0, new TerseVector(  0.678136763520,  43.094461639362,  13.243238780721), new TerseVector(-2.2358226110718e-03, -6.3233636090933e-04,  4.7664798895648e-04))
        ,   new body_state_t(   87600.0, new TerseVector(-18.282602096834, -23.305039586660,  -1.766620508028), new TerseVector( 2.5567245263557e-03, -1.9902940754171e-03, -1.3943491701082e-03))
        ,   new body_state_t(  116800.0, new TerseVector( 43.873338744526,   7.700705617215, -10.814273666425), new TerseVector( 2.3174803055677e-04,  2.2402163127924e-03,  6.2988756452032e-04))
        ,   new body_state_t(  146000.0, new TerseVector(  7.392949027906,  44.382678951534,  11.629500214854), new TerseVector(-2.1932815453830e-03, -2.1751799585364e-04,  5.9556516201114e-04))
        ,   new body_state_t(  175200.0, new TerseVector(-24.981690229261, -16.204012851426,   2.466457544298), new TerseVector( 1.8193989149580e-03, -2.6765419531201e-03, -1.3848283502247e-03))
        ,   new body_state_t(  204400.0, new TerseVector( 42.530187039511,   0.845935508021, -12.554907527683), new TerseVector( 6.5059779150669e-04,  2.2725657282262e-03,  5.1133743202822e-04))
        ,   new body_state_t(  233600.0, new TerseVector( 13.999526486822,  44.462363044894,   9.669418486465), new TerseVector(-2.1079296569252e-03,  1.7533423831993e-04,  6.9128485798076e-04))
        ,   new body_state_t(  262800.0, new TerseVector(-29.184024803031,  -7.371243995762,   6.493275957928), new TerseVector( 9.3581363109681e-04, -3.0610357109184e-03, -1.2364201089345e-03))
        ,   new body_state_t(  292000.0, new TerseVector( 39.831980671753,  -6.078405766765, -13.909815358656), new TerseVector( 1.1117769689167e-03,  2.2362097830152e-03,  3.6230548231153e-04))
        ,   new body_state_t(  321200.0, new TerseVector( 20.294955108476,  43.417190420251,   7.450091985932), new TerseVector(-1.9742157451535e-03,  5.3102050468554e-04,  7.5938408813008e-04))
        ,   new body_state_t(  350400.0, new TerseVector(-30.669992302160,   2.318743558955,   9.973480913858), new TerseVector( 4.5605107450676e-05, -3.1308219926928e-03, -9.9066533301924e-04))
        ,   new body_state_t(  379600.0, new TerseVector( 35.626122155983, -12.897647509224, -14.777586508444), new TerseVector( 1.6015684949743e-03,  2.1171931182284e-03,  1.8002516202204e-04))
        ,   new body_state_t(  408800.0, new TerseVector( 26.133186148561,  41.232139187599,   5.006401326220), new TerseVector(-1.7857704419579e-03,  8.6046232702817e-04,  8.0614690298954e-04))
        ,   new body_state_t(  438000.0, new TerseVector(-29.576740229230,  11.863535943587,  12.631323039872), new TerseVector(-7.2292830060955e-04, -2.9587820140709e-03, -7.0824296450300e-04))
        ,   new body_state_t(  467200.0, new TerseVector( 29.910805787391, -19.159019294000, -15.013363865194), new TerseVector( 2.0871080437997e-03,  1.8848372554514e-03, -3.8528655083926e-05))
        ,   new body_state_t(  496400.0, new TerseVector( 31.375957451819,  38.050372720763,   2.433138343754), new TerseVector(-1.5546055556611e-03,  1.1699815465629e-03,  8.3565439266001e-04))
        ,   new body_state_t(  525600.0, new TerseVector(-26.360071336928,  20.662505904952,  14.414696258958), new TerseVector(-1.3142373118349e-03, -2.6236647854842e-03, -4.2542017598193e-04))
        ,   new body_state_t(  554800.0, new TerseVector( 22.599441488648, -24.508879898306, -14.484045731468), new TerseVector( 2.5454108304806e-03,  1.4917058755191e-03, -3.0243665086079e-04))
        ,   new body_state_t(  584000.0, new TerseVector( 35.877864013014,  33.894226366071,  -0.224524636277), new TerseVector(-1.2941245730845e-03,  1.4560427668319e-03,  8.4762160640137e-04))
        ,   new body_state_t(  613200.0, new TerseVector(-21.538149762417,  28.204068269761,  15.321973799534), new TerseVector(-1.7312117409010e-03, -2.1939631314577e-03, -1.6316913275180e-04))
        ,   new body_state_t(  642400.0, new TerseVector( 13.971521374415, -28.339941764789, -13.083792871886), new TerseVector( 2.9334630526035e-03,  9.1860931752944e-04, -5.9939422488627e-04))
        ,   new body_state_t(  671600.0, new TerseVector( 39.526942044143,  28.939897360110,  -2.872799527539), new TerseVector(-1.0068481658095e-03,  1.7021132888090e-03,  8.3578230511981e-04))
        ,   new body_state_t(  700800.0, new TerseVector(-15.576200701394,  34.399412961275,  15.466033737854), new TerseVector(-2.0098814612884e-03, -1.7191109825989e-03,  7.0414782780416e-05))
        ,   new body_state_t(  730000.0, new TerseVector(  4.243252837090, -30.118201690825, -10.707441231349), new TerseVector( 3.1725847067411e-03,  1.6098461202270e-04, -9.0672150593868e-04))
        };

        internal static TerseVector UpdatePosition(double dt, TerseVector r, TerseVector v, TerseVector a)
        {
            return new TerseVector(
                r.x + dt*(v.x + dt*a.x/2),
                r.y + dt*(v.y + dt*a.y/2),
                r.z + dt*(v.z + dt*a.z/2)
            );
        }

        internal static TerseVector UpdateVelocity(double dt, TerseVector v, TerseVector a)
        {
            return new TerseVector(
                v.x + dt*a.x,
                v.y + dt*a.y,
                v.z + dt*a.z
            );
        }

        internal static body_state_t AdjustBarycenterPosVel(ref body_state_t ssb, double tt, Body body, double planet_gm)
        {
            double shift = planet_gm / (planet_gm + SUN_GM);
            body_state_t planet = CalcVsopPosVel(vsop[(int)body], tt);
            ssb.r += shift * planet.r;
            ssb.v += shift * planet.v;
            return planet;
        }

        private static major_bodies_t MajorBodyBary(double tt)
        {
            var bary = new major_bodies_t();
            var ssb = new body_state_t(tt, TerseVector.Zero, TerseVector.Zero);
            bary.Jupiter = AdjustBarycenterPosVel(ref ssb, tt, Body.Jupiter, JUPITER_GM);
            bary.Saturn  = AdjustBarycenterPosVel(ref ssb, tt, Body.Saturn,  SATURN_GM);
            bary.Uranus  = AdjustBarycenterPosVel(ref ssb, tt, Body.Uranus,  URANUS_GM);
            bary.Neptune = AdjustBarycenterPosVel(ref ssb, tt, Body.Neptune, NEPTUNE_GM);

            // Convert planets' [pos, vel] vectors from heliocentric to barycentric.
            bary.Jupiter.r -= ssb.r;    bary.Jupiter.v -= ssb.v;
            bary.Saturn.r  -= ssb.r;    bary.Saturn.v  -= ssb.v;
            bary.Uranus.r  -= ssb.r;    bary.Uranus.v  -= ssb.v;
            bary.Neptune.r -= ssb.r;    bary.Neptune.v -= ssb.v;

            // Convert heliocentric SSB to barycentric Sun.
            bary.Sun = -ssb;

            return bary;
        }

        private static body_grav_calc_t GravSim(    // out: [pos, vel, acc] of the simulated body at time tt2
            out major_bodies_t bary2,               // out: major body barycentric positions at tt2
            double tt2,                             // in:  a target time to be calculated (either before or after tt1
            body_grav_calc_t calc1)                 // in:  [pos, vel, acc] of the simulated body at time tt1
        {
            double dt = tt2 - calc1.tt;

            // Calculate where the major bodies (Sun, Jupiter...Neptune) will be at the next time step.
            bary2 = MajorBodyBary(tt2);

            // Estimate position of small body as if current acceleration applies across the whole time interval.
            // approx_pos = pos1 + vel1*dt + (1/2)acc*dt^2
            TerseVector approx_pos = UpdatePosition(dt, calc1.r, calc1.v, calc1.a);

            // Calculate acceleration experienced by small body at approximate next location.
            TerseVector acc = bary2.Acceleration(approx_pos);

            // Calculate the average acceleration of the endpoints.
            // This becomes our estimate of the mean effective acceleration over the whole interval.
            acc = (acc + calc1.a) / 2.0;

            // Refine the estimates of [pos, vel, acc] at tt2 using the mean acceleration.
            TerseVector pos = UpdatePosition(dt, calc1.r, calc1.v, acc);
            TerseVector vel = calc1.v + (dt * acc);
            acc = bary2.Acceleration(pos);
            return new body_grav_calc_t(tt2, pos, vel, acc);
        }

        private static readonly body_grav_calc_t[][] pluto_cache = new body_grav_calc_t[PLUTO_NUM_STATES-1][];

        private static int ClampIndex(double frac, int nsteps)
        {
            int index = (int) Math.Floor(frac);
            if (index < 0)
                return 0;
            if (index >= nsteps)
                return nsteps-1;
            return index;
        }

        private static body_grav_calc_t GravFromState(out major_bodies_t bary, body_state_t state)
        {
            bary = MajorBodyBary(state.tt);
            TerseVector r = state.r + bary.Sun.r;
            TerseVector v = state.v + bary.Sun.v;
            TerseVector a = bary.Acceleration(r);
            return new body_grav_calc_t(state.tt, r, v, a);
        }

        private static body_grav_calc_t[] GetSegment(body_grav_calc_t[][] cache, double tt)
        {
            if (tt < PlutoStateTable[0].tt || tt > PlutoStateTable[PLUTO_NUM_STATES-1].tt)
                return null;  // Don't bother calculating a segment. Let the caller crawl backward/forward to this time.

            int seg_index = ClampIndex((tt - PlutoStateTable[0].tt) / PLUTO_TIME_STEP, PLUTO_NUM_STATES-1);
            lock (cache)
            {
                if (cache[seg_index] == null)
                {
                    var seg = cache[seg_index] = new body_grav_calc_t[PLUTO_NSTEPS];

                    // Each endpoint is exact.
                    major_bodies_t bary;
                    seg[0] = GravFromState(out bary, PlutoStateTable[seg_index]);
                    seg[PLUTO_NSTEPS-1] = GravFromState(out bary, PlutoStateTable[seg_index + 1]);

                    // Simulate forwards from the lower time bound.
                    int i;
                    double step_tt = seg[0].tt;
                    for (i=1; i < PLUTO_NSTEPS-1; ++i)
                        seg[i] = GravSim(out bary, step_tt += PLUTO_DT, seg[i-1]);

                    // Simulate backwards from the upper time bound.
                    step_tt = seg[PLUTO_NSTEPS-1].tt;
                    var reverse = new body_grav_calc_t[PLUTO_NSTEPS];
                    reverse[PLUTO_NSTEPS-1] = seg[PLUTO_NSTEPS-1];
                    for (i=PLUTO_NSTEPS-2; i > 0; --i)
                        reverse[i] = GravSim(out bary, step_tt -= PLUTO_DT, reverse[i+1]);

                    // Fade-mix the two series so that there are no discontinuities.
                    for (i=PLUTO_NSTEPS-2; i > 0; --i)
                    {
                        double ramp = (double)i / (PLUTO_NSTEPS-1);
                        seg[i].r = (1 - ramp)*seg[i].r + ramp*reverse[i].r;
                        seg[i].v = (1 - ramp)*seg[i].v + ramp*reverse[i].v;
                        seg[i].a = (1 - ramp)*seg[i].a + ramp*reverse[i].a;
                    }
                }
                return cache[seg_index];
            }
        }

        private static body_grav_calc_t CalcPlutoOneWay(
            out major_bodies_t bary,
            body_state_t init_state,
            double target_tt,
            double dt)
        {
            body_grav_calc_t calc = GravFromState(out bary, init_state);
            int n = (int) Math.Ceiling((target_tt - calc.tt) / dt);
            for (int i=0; i < n; ++i)
                calc = GravSim(out bary, (i+1 == n) ? target_tt : (calc.tt + dt), calc);
            return calc;
        }

        private static StateVector CalcPluto(AstroTime time, bool helio)
        {
            body_grav_calc_t calc;
            body_grav_calc_t[] seg = GetSegment(pluto_cache, time.tt);
            var bary = new major_bodies_t();
            if (seg == null)
            {
                // The target time is outside the year range 0000..4000.
                // Calculate it by crawling backward from 0000 or forward from 4000.
                // FIXFIXFIX - This is super slow. Could optimize this with extra caching if needed.
                if (time.tt < PlutoStateTable[0].tt)
                    calc = CalcPlutoOneWay(out bary, PlutoStateTable[0], time.tt, -PLUTO_DT);
                else
                    calc = CalcPlutoOneWay(out bary, PlutoStateTable[PLUTO_NUM_STATES-1], time.tt, +PLUTO_DT);
            }
            else
            {
                int left = ClampIndex((time.tt - seg[0].tt) / PLUTO_DT, PLUTO_NSTEPS-1);
                body_grav_calc_t s1 = seg[left];
                body_grav_calc_t s2 = seg[left+1];

                // Find mean acceleration vector over the interval.
                TerseVector acc = (s1.a + s2.a) / 2.0;

                // Use Newtonian mechanics to extrapolate away from t1 in the positive time direction.
                TerseVector ra = UpdatePosition(time.tt - s1.tt, s1.r, s1.v, acc);
                TerseVector va = UpdateVelocity(time.tt - s1.tt, s1.v, acc);

                // Use Newtonian mechanics to extrapolate away from t2 in the negative time direction.
                TerseVector rb = UpdatePosition(time.tt - s2.tt, s2.r, s2.v, acc);
                TerseVector vb = UpdateVelocity(time.tt - s2.tt, s2.v, acc);

                // Use fade in/out idea to blend the two position estimates.
                double ramp = (time.tt - s1.tt)/PLUTO_DT;
                calc.r = (1 - ramp)*ra + ramp*rb;
                calc.v = (1 - ramp)*va + ramp*vb;
                if (helio)
                    bary = MajorBodyBary(time.tt);
            }

            if (helio)
            {
                // Convert barycentric vectors to heliocentric vectors
                calc.r -= bary.Sun.r;
                calc.v -= bary.Sun.v;
            }

            return new StateVector
            {
                t  = time,
                x  = calc.r.x,
                y  = calc.r.y,
                z  = calc.r.z,
                vx = calc.v.x,
                vy = calc.v.y,
                vz = calc.v.z,
            };
        }

#endregion  // Pluto

#region Jupiter's Moons

        private struct jupiter_moon_t
        {
            public double mu;
            public double al0, al1;
            public vsop_term_t[] a;
            public vsop_term_t[] l;
            public vsop_term_t[] z;
            public vsop_term_t[] zeta;
        }

        private static readonly RotationMatrix Rotation_JUP_EQJ = new RotationMatrix(
            new double[3,3]
            {
                {  9.99432765338654e-01, -3.36771074697641e-02,  0.00000000000000e+00 },
                {  3.03959428906285e-02,  9.02057912352809e-01,  4.30543388542295e-01 },
                { -1.44994559663353e-02, -4.30299169409101e-01,  9.02569881273754e-01 }
            }
        );

        private static readonly jupiter_moon_t[] JupiterMoonModel = new jupiter_moon_t[] {
            // [0] Io
            new jupiter_moon_t {
                mu =  2.8248942843381399e-07,
                al0 =  1.4462132960212239e+00,
                al1 =  3.5515522861824000e+00,
                a = new vsop_term_t[] {
                    new vsop_term_t( 0.0028210960212903,  0.0000000000000000e+00,  0.0000000000000000e+00)
                },
                l = new vsop_term_t[] {
                    new vsop_term_t(-0.0001925258348666,  4.9369589722644998e+00,  1.3584836583050000e-02),
                    new vsop_term_t(-0.0000970803596076,  4.3188796477322002e+00,  1.3034138432430000e-02),
                    new vsop_term_t(-0.0000898817416500,  1.9080016428616999e+00,  3.0506486715799999e-03),
                    new vsop_term_t(-0.0000553101050262,  1.4936156681568999e+00,  1.2938928911549999e-02)
                },
                z = new vsop_term_t[] {
                    new vsop_term_t( 0.0041510849668155,  4.0899396355450000e+00, -1.2906864146660001e-02),
                    new vsop_term_t( 0.0006260521444113,  1.4461888986270000e+00,  3.5515522949801999e+00),
                    new vsop_term_t( 0.0000352747346169,  2.1256287034577999e+00,  1.2727416566999999e-04)
                },
                zeta = new vsop_term_t[] {
                    new vsop_term_t( 0.0003142172466014,  2.7964219722923001e+00, -2.3150960980000000e-03),
                    new vsop_term_t( 0.0000904169207946,  1.0477061879627001e+00, -5.6920638196000003e-04)
                }
            },

            // [1] Europa
            new jupiter_moon_t {
                mu =  2.8248327439289299e-07,
                al0 = -3.7352634374713622e-01,
                al1 =  1.7693227111234699e+00,
                a = new vsop_term_t[] {
                    new vsop_term_t( 0.0044871037804314,  0.0000000000000000e+00,  0.0000000000000000e+00),
                    new vsop_term_t( 0.0000004324367498,  1.8196456062910000e+00,  1.7822295777568000e+00)
                },
                l = new vsop_term_t[] {
                    new vsop_term_t( 0.0008576433172936,  4.3188693178264002e+00,  1.3034138308049999e-02),
                    new vsop_term_t( 0.0004549582875086,  1.4936531751079001e+00,  1.2938928819619999e-02),
                    new vsop_term_t( 0.0003248939825174,  1.8196494533458001e+00,  1.7822295777568000e+00),
                    new vsop_term_t(-0.0003074250079334,  4.9377037005910998e+00,  1.3584832867240000e-02),
                    new vsop_term_t( 0.0001982386144784,  1.9079869054759999e+00,  3.0510121286900001e-03),
                    new vsop_term_t( 0.0001834063551804,  2.1402853388529000e+00,  1.4500978933800000e-03),
                    new vsop_term_t(-0.0001434383188452,  5.6222140366630002e+00,  8.9111478887838003e-01),
                    new vsop_term_t(-0.0000771939140944,  4.3002724372349999e+00,  2.6733443704265998e+00)
                },
                z = new vsop_term_t[] {
                    new vsop_term_t(-0.0093589104136341,  4.0899396509038999e+00, -1.2906864146660001e-02),
                    new vsop_term_t( 0.0002988994545555,  5.9097265185595003e+00,  1.7693227079461999e+00),
                    new vsop_term_t( 0.0002139036390350,  2.1256289300016000e+00,  1.2727418406999999e-04),
                    new vsop_term_t( 0.0001980963564781,  2.7435168292649998e+00,  6.7797343008999997e-04),
                    new vsop_term_t( 0.0001210388158965,  5.5839943711203004e+00,  3.2056614899999997e-05),
                    new vsop_term_t( 0.0000837042048393,  1.6094538368039000e+00, -9.0402165808846002e-01),
                    new vsop_term_t( 0.0000823525166369,  1.4461887708689001e+00,  3.5515522949801999e+00)
                },
                zeta = new vsop_term_t[] {
                    new vsop_term_t( 0.0040404917832303,  1.0477063169425000e+00, -5.6920640539999997e-04),
                    new vsop_term_t( 0.0002200421034564,  3.3368857864364001e+00, -1.2491307306999999e-04),
                    new vsop_term_t( 0.0001662544744719,  2.4134862374710999e+00,  0.0000000000000000e+00),
                    new vsop_term_t( 0.0000590282470983,  5.9719930968366004e+00, -3.0561602250000000e-05)
                }
            },

            // [2] Ganymede
            new jupiter_moon_t {
                mu =  2.8249818418472298e-07,
                al0 =  2.8740893911433479e-01,
                al1 =  8.7820792358932798e-01,
                a = new vsop_term_t[] {
                    new vsop_term_t( 0.0071566594572575,  0.0000000000000000e+00,  0.0000000000000000e+00),
                    new vsop_term_t( 0.0000013930299110,  1.1586745884981000e+00,  2.6733443704265998e+00)
                },
                l = new vsop_term_t[] {
                    new vsop_term_t( 0.0002310797886226,  2.1402987195941998e+00,  1.4500978438400001e-03),
                    new vsop_term_t(-0.0001828635964118,  4.3188672736968003e+00,  1.3034138282630000e-02),
                    new vsop_term_t( 0.0001512378778204,  4.9373102372298003e+00,  1.3584834812520000e-02),
                    new vsop_term_t(-0.0001163720969778,  4.3002659861490002e+00,  2.6733443704265998e+00),
                    new vsop_term_t(-0.0000955478069846,  1.4936612842567001e+00,  1.2938928798570001e-02),
                    new vsop_term_t( 0.0000815246854464,  5.6222137132535002e+00,  8.9111478887838003e-01),
                    new vsop_term_t(-0.0000801219679602,  1.2995922951532000e+00,  1.0034433456728999e+00),
                    new vsop_term_t(-0.0000607017260182,  6.4978769669238001e-01,  5.0172167043264004e-01)
                },
                z = new vsop_term_t[] {
                    new vsop_term_t( 0.0014289811307319,  2.1256295942738999e+00,  1.2727413029000001e-04),
                    new vsop_term_t( 0.0007710931226760,  5.5836330003496002e+00,  3.2064341100000001e-05),
                    new vsop_term_t( 0.0005925911780766,  4.0899396636447998e+00, -1.2906864146660001e-02),
                    new vsop_term_t( 0.0002045597496146,  5.2713683670371996e+00, -1.2523544076106000e-01),
                    new vsop_term_t( 0.0001785118648258,  2.8743156721063001e-01,  8.7820792442520001e-01),
                    new vsop_term_t( 0.0001131999784893,  1.4462127277818000e+00,  3.5515522949801999e+00),
                    new vsop_term_t(-0.0000658778169210,  2.2702423990985001e+00, -1.7951364394536999e+00),
                    new vsop_term_t( 0.0000497058888328,  5.9096792204858000e+00,  1.7693227129285001e+00)
                },
                zeta = new vsop_term_t[] {
                    new vsop_term_t( 0.0015932721570848,  3.3368862796665000e+00, -1.2491307058000000e-04),
                    new vsop_term_t( 0.0008533093128905,  2.4133881688166001e+00,  0.0000000000000000e+00),
                    new vsop_term_t( 0.0003513347911037,  5.9720789850126996e+00, -3.0561017709999999e-05),
                    new vsop_term_t(-0.0001441929255483,  1.0477061764435001e+00, -5.6920632124000004e-04)
                }
            },

            // [3] Callisto
            new jupiter_moon_t {
                mu =  2.8249214488990899e-07,
                al0 = -3.6203412913757038e-01,
                al1 =  3.7648623343382798e-01,
                a = new vsop_term_t[] {
                    new vsop_term_t( 0.0125879701715314,  0.0000000000000000e+00,  0.0000000000000000e+00),
                    new vsop_term_t( 0.0000035952049470,  6.4965776007116005e-01,  5.0172168165034003e-01),
                    new vsop_term_t( 0.0000027580210652,  1.8084235781510001e+00,  3.1750660413359002e+00)
                },
                l = new vsop_term_t[] {
                    new vsop_term_t( 0.0005586040123824,  2.1404207189814999e+00,  1.4500979323100001e-03),
                    new vsop_term_t(-0.0003805813868176,  2.7358844897852999e+00,  2.9729650620000000e-05),
                    new vsop_term_t( 0.0002205152863262,  6.4979652596399995e-01,  5.0172167243580001e-01),
                    new vsop_term_t( 0.0001877895151158,  1.8084787604004999e+00,  3.1750660413359002e+00),
                    new vsop_term_t( 0.0000766916975242,  6.2720114319754998e+00,  1.3928364636651001e+00),
                    new vsop_term_t( 0.0000747056855106,  1.2995916202344000e+00,  1.0034433456728999e+00)
                },
                z = new vsop_term_t[] {
                    new vsop_term_t( 0.0073755808467977,  5.5836071576083999e+00,  3.2065099140000001e-05),
                    new vsop_term_t( 0.0002065924169942,  5.9209831565786004e+00,  3.7648624194703001e-01),
                    new vsop_term_t( 0.0001589869764021,  2.8744006242622999e-01,  8.7820792442520001e-01),
                    new vsop_term_t(-0.0001561131605348,  2.1257397865089001e+00,  1.2727441285000001e-04),
                    new vsop_term_t( 0.0001486043380971,  1.4462134301023000e+00,  3.5515522949801999e+00),
                    new vsop_term_t( 0.0000635073108731,  5.9096803285953996e+00,  1.7693227129285001e+00),
                    new vsop_term_t( 0.0000599351698525,  4.1125517584797997e+00, -2.7985797954588998e+00),
                    new vsop_term_t( 0.0000540660842731,  5.5390350845569003e+00,  2.8683408228299999e-03),
                    new vsop_term_t(-0.0000489596900866,  4.6218149483337996e+00, -6.2695712529518999e-01)
                },
                zeta = new vsop_term_t[] {
                    new vsop_term_t( 0.0038422977898495,  2.4133922085556998e+00,  0.0000000000000000e+00),
                    new vsop_term_t( 0.0022453891791894,  5.9721736773277003e+00, -3.0561255249999997e-05),
                    new vsop_term_t(-0.0002604479450559,  3.3368746306408998e+00, -1.2491309972000001e-04),
                    new vsop_term_t( 0.0000332112143230,  5.5604137742336999e+00,  2.9003768850700000e-03)
                }
            }
        };

        private static StateVector JupiterMoon_elem2pv(
            AstroTime time,
            double mu,
            double A, double AL, double K, double H, double Q, double P)
        {
            // Translation of FORTRAN subroutine ELEM2PV from:
            // https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

            double AN = Math.Sqrt(mu / (A*A*A));

            double CE, SE, DE;
            double EE = AL + K*Math.Sin(AL) - H*Math.Cos(AL);
            do
            {
                CE = Math.Cos(EE);
                SE = Math.Sin(EE);
                DE = (AL - EE + K*SE - H*CE) / (1.0 - K*CE - H*SE);
                EE += DE;
            }
            while (Math.Abs(DE) >= 1.0e-12);

            CE = Math.Cos(EE);
            SE = Math.Sin(EE);
            double DLE = H*CE - K*SE;
            double RSAM1 = -K*CE - H*SE;
            double ASR = 1.0/(1.0 + RSAM1);
            double PHI = Math.Sqrt(1.0 - K*K - H*H);
            double PSI = 1.0/(1.0 + PHI);
            double X1 = A*(CE - K - PSI*H*DLE);
            double Y1 = A*(SE - H + PSI*K*DLE);
            double VX1 = AN*ASR*A*(-SE - PSI*H*RSAM1);
            double VY1 = AN*ASR*A*(+CE + PSI*K*RSAM1);
            double F2 = 2.0*Math.Sqrt(1.0 - Q*Q - P*P);
            double P2 = 1.0 - 2.0*P*P;
            double Q2 = 1.0 - 2.0*Q*Q;
            double PQ = 2.0*P*Q;

            return new StateVector(
                X1*P2 + Y1*PQ,
                X1*PQ + Y1*Q2,
                (Q*Y1 - X1*P)*F2,
                VX1*P2 + VY1*PQ,
                VX1*PQ + VY1*Q2,
                (Q*VY1 - VX1*P)*F2,
                time
            );
        }

        private static StateVector CalcJupiterMoon(AstroTime time, jupiter_moon_t m)
        {
            // This is a translation of FORTRAN code by Duriez, Lainey, and Vienne:
            // https://ftp.imcce.fr/pub/ephem/satel/galilean/L1/L1.2/

            double t = time.tt + 18262.5;     // number of days since 1950-01-01T00:00:00Z

            // Calculate 6 orbital elements at the given time t.
            double elem0 = 0.0;
            foreach (vsop_term_t term in m.a)
                elem0 += term.amplitude * Math.Cos(term.phase + (t * term.frequency));

            double elem1 = m.al0 + (t * m.al1);
            foreach (vsop_term_t term in m.l)
                elem1 += term.amplitude * Math.Sin(term.phase + (t * term.frequency));

            elem1 %= PI2;
            if (elem1 < 0)
                elem1 += PI2;

            double elem2 = 0.0;
            double elem3 = 0.0;
            foreach (vsop_term_t term in m.z)
            {
                double arg = term.phase + (t * term.frequency);
                elem2 += term.amplitude * Math.Cos(arg);
                elem3 += term.amplitude * Math.Sin(arg);
            }

            double elem4 = 0.0;
            double elem5 = 0.0;
            foreach (vsop_term_t term in m.zeta)
            {
                double arg = term.phase + (t * term.frequency);
                elem4 += term.amplitude * Math.Cos(arg);
                elem5 += term.amplitude * Math.Sin(arg);
            }

            // Convert the oribital elements into position vectors in the Jupiter equatorial system (JUP).
            StateVector state = JupiterMoon_elem2pv(time, m.mu, elem0, elem1, elem2, elem3, elem4, elem5);

            // Re-orient position and velocity vectors from Jupiter-equatorial (JUP) to Earth-equatorial in J2000 (EQJ).
            return RotateState(Rotation_JUP_EQJ, state);
        }

        /// <summary>
        /// Calculates jovicentric positions and velocities of Jupiter's largest 4 moons.
        /// </summary>
        /// <remarks>
        /// Calculates position and velocity vectors for Jupiter's moons
        /// Io, Europa, Ganymede, and Callisto, at the given date and time.
        /// The vectors are jovicentric (relative to the center of Jupiter).
        /// Their orientation is the Earth's equatorial system at the J2000 epoch (EQJ).
        /// The position components are expressed in astronomical units (AU), and the
        /// velocity components are in AU/day.
        ///
        /// To convert to heliocentric position vectors, call #Astronomy.HelioVector
        /// with `Body.Jupiter` to get Jupiter's heliocentric position, then
        /// add the jovicentric positions. Likewise, you can call #Astronomy.GeoVector
        /// to convert to geocentric positions; however, you will have to manually
        /// correct for light travel time from the Jupiter system to Earth to
        /// figure out what time to pass to `jupiterMoons` to get an accurate picture
        /// of how Jupiter and its moons look from Earth.
        /// </remarks>
        /// <param name="time">The date and time for which to calculate the position vectors.</param>
        /// <returns>Position and velocity vectors of Jupiter's largest 4 moons.</returns>
        public static JupiterMoonsInfo JupiterMoons(AstroTime time) =>
            new JupiterMoonsInfo
            {
                io       = CalcJupiterMoon(time, JupiterMoonModel[0]),
                europa   = CalcJupiterMoon(time, JupiterMoonModel[1]),
                ganymede = CalcJupiterMoon(time, JupiterMoonModel[2]),
                callisto = CalcJupiterMoon(time, JupiterMoonModel[3]),
            };

#endregion  // Jupiter's Moons

        private enum PrecessDirection
        {
            From2000,
            Into2000,
        }

        private static RotationMatrix precession_rot(AstroTime time, PrecessDirection dir)
        {
            double t = time.tt / 36525;
            double eps0 = 84381.406;

            double psia   = (((((-    0.0000000951  * t
                                +    0.000132851 ) * t
                                -    0.00114045  ) * t
                                -    1.0790069   ) * t
                                + 5038.481507    ) * t);

            double omegaa = (((((+    0.0000003337  * t
                                 -    0.000000467 ) * t
                                 -    0.00772503  ) * t
                                 +    0.0512623   ) * t
                                 -    0.025754    ) * t + eps0);

            double chia   = (((((-    0.0000000560  * t
                                 +    0.000170663 ) * t
                                 -    0.00121197  ) * t
                                 -    2.3814292   ) * t
                                 +   10.556403    ) * t);

            eps0   *= ASEC2RAD;
            psia   *= ASEC2RAD;
            omegaa *= ASEC2RAD;
            chia   *= ASEC2RAD;

            double sa = Math.Sin(eps0);
            double ca = Math.Cos(eps0);
            double sb = Math.Sin(-psia);
            double cb = Math.Cos(-psia);
            double sc = Math.Sin(-omegaa);
            double cc = Math.Cos(-omegaa);
            double sd = Math.Sin(chia);
            double cd = Math.Cos(chia);

            double xx =  cd*cb - sb*sd*cc;
            double yx =  cd*sb*ca + sd*cc*cb*ca - sa*sd*sc;
            double zx =  cd*sb*sa + sd*cc*cb*sa + ca*sd*sc;
            double xy = -sd*cb - sb*cd*cc;
            double yy = -sd*sb * ca + cd*cc*cb*ca - sa*cd*sc;
            double zy = -sd*sb * sa + cd*cc*cb*sa + ca*cd*sc;
            double xz =  sb*sc;
            double yz = -sc*cb*ca - sa*cc;
            double zz = -sc*cb*sa + cc*ca;

            var rot = new double[3,3];
            if (dir == PrecessDirection.Into2000)
            {
                // Perform rotation from other epoch to J2000.0.
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
            else if (dir == PrecessDirection.From2000)
            {
                // Perform rotation from J2000.0 to other epoch.
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
                throw new ArgumentException("Unsupported precess direction: " + dir);
            }

            return new RotationMatrix(rot);
        }

        private static AstroVector precession(AstroVector pos, PrecessDirection dir)
        {
            RotationMatrix r = precession_rot(pos.t, dir);
            return RotateVector(r, pos);
        }

        private static StateVector precession_posvel(StateVector state, PrecessDirection dir)
        {
            RotationMatrix rot = precession_rot(state.t, dir);
            return RotateState(rot, state);
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

        private static void iau2000b(AstroTime time)
        {
            // Adapted from the NOVAS C 3.1 function of the same name.

            if (double.IsNaN(time.psi))
            {
                double t = time.tt / 36525.0;
                double elp = ((1287104.79305 + t * 129596581.0481)  % ASEC360) * ASEC2RAD;
                double f   = ((335779.526232 + t * 1739527262.8478) % ASEC360) * ASEC2RAD;
                double d   = ((1072260.70369 + t * 1602961601.2090) % ASEC360) * ASEC2RAD;
                double om  = ((450160.398036 - t * 6962890.5431)    % ASEC360) * ASEC2RAD;

                double sarg = Math.Sin(om);
                double carg = Math.Cos(om);
                double dp = (-172064161.0 - 174666.0*t)*sarg + 33386.0*carg;
                double de = (92052331.0 + 9086.0*t)*carg + 15377.0*sarg;

                double arg = 2.0*(f - d + om);
                sarg = Math.Sin(arg);
                carg = Math.Cos(arg);
                dp += (-13170906.0 - 1675.0*t)*sarg - 13696.0*carg;
                de += (5730336.0 - 3015.0*t)*carg - 4587.0*sarg;

                arg = 2.0*(f + om);
                sarg = Math.Sin(arg);
                carg = Math.Cos(arg);
                dp += (-2276413.0 - 234.0*t)*sarg + 2796.0*carg;
                de += (978459.0 - 485.0*t)*carg + 1374.0*sarg;

                arg = 2.0*om;
                sarg = Math.Sin(arg);
                carg = Math.Cos(arg);
                dp += (2074554.0 + 207.0*t)*sarg - 698.0*carg;
                de += (-897492.0 + 470.0*t)*carg - 291.0*sarg;

                sarg = Math.Sin(elp);
                carg = Math.Cos(elp);
                dp += (1475877.0 - 3633.0*t)*sarg + 11817.0*carg;
                de += (73871.0 - 184.0*t)*carg - 1924.0*sarg;

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

        private static double era(double ut)        // Earth Rotation Angle
        {
            double thet1 = 0.7790572732640 + 0.00273781191135448 * ut;
            double thet3 = ut % 1.0;
            double theta = 360.0 *((thet1 + thet3) % 1.0);
            if (theta < 0.0)
                theta += 360.0;

            return theta;
        }

        /// <summary>
        /// Calculates Greenwich Apparent Sidereal Time (GAST).
        /// </summary>
        /// <remarks>
        /// Given a date and time, this function calculates the rotation of the
        /// Earth, represented by the equatorial angle of the Greenwich prime meridian
        /// with respect to distant stars (not the Sun, which moves relative to background
        /// stars by almost one degree per day).
        /// This angle is called Greenwich Apparent Sidereal Time (GAST).
        /// GAST is measured in sidereal hours in the half-open range [0, 24).
        /// When GAST = 0, it means the prime meridian is aligned with the of-date equinox,
        /// corrected at that time for precession and nutation of the Earth's axis.
        /// In this context, the "equinox" is the direction in space where the Earth's
        /// orbital plane (the ecliptic) intersects with the plane of the Earth's equator,
        /// at the location on the Earth's orbit of the (seasonal) March equinox.
        /// As the Earth rotates, GAST increases from 0 up to 24 sidereal hours,
        /// then starts over at 0.
        /// To convert to degrees, multiply the return value by 15.
        /// </remarks>
        /// <param name="time">
        /// The date and time for which to find GAST.
        /// As an optimization, this function caches the sidereal time value in `time`,
        /// unless it has already been cached, in which case the cached value is reused.
        /// </param>
        /// <returns>GAST in sidereal hours.</returns>
        public static double SiderealTime(AstroTime time)
        {
            if (double.IsNaN(time.st))
            {
                double t = time.tt / 36525.0;
                double eqeq = 15.0 * e_tilt(time).ee;    // Replace with eqeq=0 to get GMST instead of GAST (if we ever need it)
                double theta = era(time.ut);
                double st = (eqeq + 0.014506 +
                    (((( -    0.0000000368  * t
                        -    0.000029956  ) * t
                        -    0.00000044   ) * t
                        +    1.3915817    ) * t
                        + 4612.156534     ) * t);

                double gst = ((st/3600.0 + theta) % 360.0) / 15.0;
                if (gst < 0.0)
                    gst += 24.0;

                time.st = gst;
            }
            return time.st;     // return sidereal hours in the half-open range [0, 24).
        }

        private static Observer inverse_terra(AstroVector ovec)
        {
            double lon_deg, lat_deg, height_km;

            // Convert from AU to kilometers.
            double x = ovec.x * KM_PER_AU;
            double y = ovec.y * KM_PER_AU;
            double z = ovec.z * KM_PER_AU;
            double p = hypot(x, y);
            if (p < 1.0e-6)
            {
                // Special case: within 1 millimeter of a pole!
                // Use arbitrary longitude, and latitude determined by polarity of z.
                lon_deg = 0.0;
                lat_deg = (z > 0.0) ? +90.0 : -90.0;
                // Elevation is calculated directly from z
                height_km = Math.Abs(z) - EARTH_POLAR_RADIUS_KM;
            }
            else
            {
                double stlocl = Math.Atan2(y, x);
                double st = SiderealTime(ovec.t);
                // Calculate exact longitude.
                lon_deg = RAD2DEG*stlocl - (15.0 * st);
                // Normalize longitude to the range (-180, +180].
                while (lon_deg <= -180.0)
                    lon_deg += 360.0;
                while (lon_deg > +180.0)
                    lon_deg -= 360.0;
                // Numerically solve for exact latitude, using Newton's Method.
                double F = EARTH_FLATTENING * EARTH_FLATTENING;
                // Start with initial latitude estimate, based on a spherical Earth.
                double lat = Math.Atan2(z, p);
                double c, s, denom;
                int count = 0;
                for(;;)
                {
                    if (++count > 10)
                        throw new InternalError("inverse_terra solver failed to converge.");
                    // Calculate the error function W(lat).
                    // We try to find the root of W, meaning where the error is 0.
                    c = Math.Cos(lat);
                    s = Math.Sin(lat);
                    double factor = (F-1)*EARTH_EQUATORIAL_RADIUS_KM;
                    double c2 = c*c;
                    double s2 = s*s;
                    double radicand = c2 + F*s2;
                    denom = Math.Sqrt(radicand);
                    double W = (factor*s*c)/denom - z*c + p*s;
                    if (Math.Abs(W) < 1.0e-8)
                        break;  // The error is now negligible.
                    // Error is still too large. Find the next estimate.
                    // Calculate D = the derivative of W with respect to lat.
                    double D = factor*((c2 - s2)/denom - s2*c2*(F-1)/(factor*radicand)) + z*s + p*c;
                    lat -= W/D;
                }
                // We now have a solution for the latitude in radians.
                lat_deg = lat * RAD2DEG;
                // Solve for exact height in kilometers.
                // There are two formulas I can use. Use whichever has the less risky denominator.
                double adjust = EARTH_EQUATORIAL_RADIUS_KM / denom;
                if (Math.Abs(s) > Math.Abs(c))
                    height_km = z/s - F*adjust;
                else
                    height_km = p/c - adjust;
            }

            return new Observer(lat_deg, lon_deg, 1000.0 * height_km);
        }

        private static StateVector terra(Observer observer, AstroTime time)
        {
            double st = SiderealTime(time);
            double phi = observer.latitude * DEG2RAD;
            double sinphi = Math.Sin(phi);
            double cosphi = Math.Cos(phi);
            double c = 1.0 / hypot(cosphi, EARTH_FLATTENING * sinphi);
            double s = (EARTH_FLATTENING * EARTH_FLATTENING) * c;
            double ht_km = observer.height / 1000.0;
            double ach = EARTH_EQUATORIAL_RADIUS_KM*c + ht_km;
            double ash = EARTH_EQUATORIAL_RADIUS_KM*s + ht_km;
            double stlocl = (15.0*st + observer.longitude) * DEG2RAD;
            double sinst = Math.Sin(stlocl);
            double cosst = Math.Cos(stlocl);

            return new StateVector(
                ach * cosphi * cosst / KM_PER_AU,
                ach * cosphi * sinst / KM_PER_AU,
                ash * sinphi / KM_PER_AU,
                -(ANGVEL * 86400.0 / KM_PER_AU) * ach * cosphi * sinst,
                +(ANGVEL * 86400.0 / KM_PER_AU) * ach * cosphi * cosst,
                0.0,
                time
            );
        }

        private static RotationMatrix nutation_rot(AstroTime time, PrecessDirection dir)
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

            if (dir == PrecessDirection.From2000)
            {
                // convert J2000 to of-date
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
            else if (dir == PrecessDirection.Into2000)
            {
                // convert of-date to J2000
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
                throw new ArgumentException("Unsupported nutation direction: " + dir);
            }

            return new RotationMatrix(rot);
        }


        private static AstroVector nutation(AstroVector pos, PrecessDirection dir)
        {
            RotationMatrix rot = nutation_rot(pos.t, dir);
            return RotateVector(rot, pos);
        }

        private static StateVector nutation_posvel(StateVector state, PrecessDirection dir)
        {
            RotationMatrix rot = nutation_rot(state.t, dir);
            return RotateState(rot, state);
        }

        private static AstroVector gyration(AstroVector pos, PrecessDirection dir)
        {
            // Combine nutation and precession into a single operation I call "gyration".
            // The order they are composed depends on the direction,
            // because both directions are mutual inverse functions.
            return (dir == PrecessDirection.Into2000) ?
                precession(nutation(pos, dir), dir) :
                nutation(precession(pos, dir), dir);
        }

        private static StateVector gyration_posvel(StateVector state, PrecessDirection dir)
        {
            // Combine nutation and precession into a single operation I call "gyration".
            // The order they are composed depends on the direction,
            // because both directions are mutual inverse functions.
            return (dir == PrecessDirection.Into2000) ?
                precession_posvel(nutation_posvel(state, dir), dir) :
                nutation_posvel(precession_posvel(state, dir), dir);
        }

        private static AstroVector geo_pos(AstroTime time, Observer observer)
        {
            AstroVector pos = terra(observer, time).Position();
            return gyration(pos, PrecessDirection.Into2000);
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
                pos.t
            );
        }

        private static AstroVector ecl2equ_vec(AstroVector ecl, double obl)
        {
            double cos_obl = Math.Cos(obl);
            double sin_obl = Math.Sin(obl);

            return new AstroVector(
                ecl.x,
                ecl.y*cos_obl - ecl.z*sin_obl,
                ecl.y*sin_obl + ecl.z*cos_obl,
                ecl.t
            );
        }

        private static AstroVector ecl2equ_vec(AstroVector ecl)
        {
            return ecl2equ_vec(ecl, mean_obliq(ecl.t.tt) * DEG2RAD);
        }

        /// <summary>
        /// Calculates equatorial geocentric position of the Moon at a given time.
        /// </summary>
        /// <remarks>
        /// Given a time of observation, calculates the Moon's position vector.
        /// The vector indicates the Moon's center relative to the Earth's center.
        /// The vector components are expressed in AU (astronomical units).
        /// The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
        /// In Astronomy Engine, this orientation is called EQJ.
        /// </remarks>
        /// <param name="time">The date and time for which to calculate the Moon's position.</param>
        /// <returns>The Moon's position vector in J2000 equatorial coordinates (EQJ).</returns>
        public static AstroVector GeoMoon(AstroTime time)
        {
            var context = new MoonContext(time.tt / 36525.0);
            MoonResult moon = context.CalcMoon();

            // Convert geocentric ecliptic spherical coordinates to Cartesian coordinates.
            double dist_cos_lat = moon.distance_au * Math.Cos(moon.geo_eclip_lat);

            var gepos = new AstroVector(
                dist_cos_lat * Math.Cos(moon.geo_eclip_lon),
                dist_cos_lat * Math.Sin(moon.geo_eclip_lon),
                moon.distance_au * Math.Sin(moon.geo_eclip_lat),
                time
            );

            // Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date.
            AstroVector mpos1 = ecl2equ_vec(gepos);

            // Convert from mean equinox of date to J2000.
            AstroVector mpos2 = precession(mpos1, PrecessDirection.Into2000);

            return mpos2;
        }

        /// <summary>
        /// Calculates spherical ecliptic geocentric position of the Moon.
        /// </summary>
        /// <remarks>
        /// Given a time of observation, calculates the Moon's geocentric position
        /// in ecliptic spherical coordinates. Provides the ecliptic latitude and
        /// longitude in degrees, and the geocentric distance in astronomical units (AU).
        ///
        /// The ecliptic angles are measured in "ECT": relative to the true ecliptic plane and
        /// equatorial plane at the specified time. This means the Earth's equator
        /// is corrected for precession and nutation, and the plane of the Earth's
        /// orbit is corrected for gradual obliquity drift.
        ///
        /// This algorithm is based on the Nautical Almanac Office's *Improved Lunar Ephemeris* of 1954,
        /// which in turn derives from E. W. Brown's lunar theories from the early twentieth century.
        /// It is adapted from Turbo Pascal code from the book
        /// [Astronomy on the Personal Computer](https://www.springer.com/us/book/9783540672210)
        /// by Montenbruck and Pfleger.
        ///
        /// To calculate a J2000 mean equator vector instead, use #Astronomy.GeoMoon.
        /// </remarks>
        /// <param name="time">
        /// The date and time for which to calculate the Moon's position.
        /// </param>
        public static Spherical EclipticGeoMoon(AstroTime time)
        {
            // Find ecliptic coordinates of the Moon in mean equinox of date (ECM).
            var context = new MoonContext(time.tt / 36525.0);
            MoonResult moon = context.CalcMoon();

            // Convert spherical coordinates to a vector.
            // The MoonResult angles are already expressed in radians.
            double dist_cos_lat = moon.distance_au * Math.Cos(moon.geo_eclip_lat);
            var ecm = new AstroVector(
                dist_cos_lat * Math.Cos(moon.geo_eclip_lon),
                dist_cos_lat * Math.Sin(moon.geo_eclip_lon),
                moon.distance_au * Math.Sin(moon.geo_eclip_lat),
                time
            );

            // Obtain true and mean obliquity angles for the given time.
            // This serves to pre-calculate the nutation also, and cache it in `time`.
            earth_tilt_t et = e_tilt(time);

            // Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date.
            AstroVector eqm = ecl2equ_vec(ecm, et.mobl * DEG2RAD);

            // Add nutation to convert ECM to true equatorial coordinates of date (EQD).
            AstroVector eqd = nutation(eqm, PrecessDirection.From2000);

            // Convert back to ecliptic, this time in true equinox of date (ECT).
            Ecliptic eclip = RotateEquatorialToEcliptic(eqd, et.tobl * DEG2RAD);

            return new Spherical(eclip.elat, eclip.elon, moon.distance_au);
        }

        /// <summary>
        /// Calculates equatorial geocentric position and velocity of the Moon at a given time.
        /// </summary>
        /// <remarks>
        /// Given a time of observation, calculates the Moon's position and velocity vectors.
        /// The position and velocity are of the Moon's center relative to the Earth's center.
        /// The position (x, y, z) components are expressed in AU (astronomical units).
        /// The velocity (vx, vy, vz) components are expressed in AU/day.
        /// The coordinates are oriented with respect to the Earth's equator at the J2000 epoch.
        /// In Astronomy Engine, this orientation is called EQJ.
        /// If you need the Moon's position only, and not its velocity,
        /// it is much more efficient to use #Astronomy.GeoMoon instead.
        /// </remarks>
        /// <param name="time">The date and time for which to calculate the Moon's position and velocity.</param>
        /// <returns>The Moon's position and velocity vectors in J2000 equatorial coordinates (EQJ).</returns>
        public static StateVector GeoMoonState(AstroTime time)
        {
            // This is a hack, because trying to figure out how to derive a time
            // derivative for CalcMoon() would be extremely painful!
            // Calculate just before and just after the given time.
            // Average to find position, subtract to find velocity.
            const double dt = 1.0e-5;   // 0.864 seconds

            AstroTime t1 = time.AddDays(-dt);
            AstroTime t2 = time.AddDays(+dt);

            AstroVector r1 = GeoMoon(t1);
            AstroVector r2 = GeoMoon(t2);

            // The desired position is the average of the two calculated positions.
            StateVector s;
            s.x = (r1.x + r2.x) / 2;
            s.y = (r1.y + r2.y) / 2;
            s.z = (r1.z + r2.z) / 2;

            // The difference of the position vectors divided by the time span gives the velocity vector.
            s.vx = (r2.x - r1.x) / (2 * dt);
            s.vy = (r2.y - r1.y) / (2 * dt);
            s.vz = (r2.z - r1.z) / (2 * dt);
            s.t = time;

            return s;
        }

        /// <summary>
        /// Calculates the geocentric position and velocity of the Earth/Moon barycenter.
        /// </summary>
        /// <remarks>
        /// Given a time of observation, calculates the geocentric position and velocity vectors
        /// of the Earth/Moon barycenter (EMB).
        /// The position (x, y, z) components are expressed in AU (astronomical units).
        /// The velocity (vx, vy, vz) components are expressed in AU/day.
        /// </remarks>
        /// <param name="time">The date and time for which to calculate the EMB vectors.</param>
        /// <returns>The EMB's position and velocity vectors in geocentric J2000 equatorial coordinates.</returns>
        public static StateVector GeoEmbState(AstroTime time)
        {
            StateVector s = GeoMoonState(time);
            const double d = 1.0 + EARTH_MOON_MASS_RATIO;
            s.x /= d;
            s.y /= d;
            s.z /= d;
            s.vx /= d;
            s.vy /= d;
            s.vz /= d;
            return s;
        }

        /// <summary>
        /// Calculates the Moon's libration angles at a given moment in time.
        /// </summary>
        /// <remarks>
        /// Libration is an observed back-and-forth wobble of the portion of the
        /// Moon visible from the Earth. It is caused by the imperfect tidal locking
        /// of the Moon's fixed rotation rate, compared to its variable angular speed
        /// of orbit around the Earth.
        ///
        /// This function calculates a pair of perpendicular libration angles,
        /// one representing rotation of the Moon in ecliptic longitude `elon`, the other
        /// in ecliptic latitude `elat`, both relative to the Moon's mean Earth-facing position.
        ///
        /// This function also returns the geocentric position of the Moon
        /// expressed in ecliptic longitude `mlon`, ecliptic latitude `mlat`, the
        /// distance `dist_km` between the centers of the Earth and Moon expressed in kilometers,
        /// and the apparent angular diameter of the Moon `diam_deg`.
        /// </remarks>
        /// <param name="time">The date and time for which to calculate lunar libration.</param>
        /// <returns>The Moon's ecliptic position and libration angles as seen from the Earth.</returns>
        public static LibrationInfo Libration(AstroTime time)
        {
            double t = time.tt / 36525.0;
            double t2 = t * t;
            double t3 = t2 * t;
            double t4 = t2 * t2;

            var context = new MoonContext(t);
            MoonResult moon = context.CalcMoon();

            LibrationInfo lib;
            lib.mlon = RAD2DEG * moon.geo_eclip_lon;
            lib.mlat = RAD2DEG * moon.geo_eclip_lat;
            lib.dist_km = moon.distance_au * KM_PER_AU;
            lib.diam_deg = (2.0 * RAD2DEG) * Math.Atan(MOON_MEAN_RADIUS_KM / Math.Sqrt(lib.dist_km*lib.dist_km - MOON_MEAN_RADIUS_KM*MOON_MEAN_RADIUS_KM));

            // Inclination angle
            const double I = DEG2RAD * 1.543;

            // Moon's argument of latitude in radians.
            double f = DEG2RAD * NormalizeLongitude(93.2720950 + 483202.0175233*t - 0.0036539*t2 - t3/3526000 + t4/863310000);

            // Moon's ascending node's mean longitude in radians.
            double omega = DEG2RAD * NormalizeLongitude(125.0445479 - 1934.1362891*t + 0.0020754*t2 + t3/467441 - t4/60616000);

            // Sun's mean anomaly.
            double m = DEG2RAD * NormalizeLongitude(357.5291092 + 35999.0502909*t - 0.0001536*t2 + t3/24490000);

            // Moon's mean anomaly.
            double mdash = DEG2RAD * NormalizeLongitude(134.9633964 + 477198.8675055*t + 0.0087414*t2 + t3/69699 - t4/14712000);

            // Moon's mean elongation.
            double d = DEG2RAD * NormalizeLongitude(297.8501921 + 445267.1114034*t - 0.0018819*t2 + t3/545868 - t4/113065000);

            // Eccentricity of the Earth's orbit.
            double e = 1.0 - 0.002516*t - 0.0000074*t2;

            // Optical librations
            double w = moon.geo_eclip_lon - omega;
            double a = Math.Atan2(Math.Sin(w)*Math.Cos(moon.geo_eclip_lat)*Math.Cos(I) - Math.Sin(moon.geo_eclip_lat)*Math.Sin(I), Math.Cos(w)*Math.Cos(moon.geo_eclip_lat));
            double ldash = LongitudeOffset(RAD2DEG * (a - f));
            double bdash = Math.Asin(-Math.Sin(w)*Math.Cos(moon.geo_eclip_lat)*Math.Sin(I) - Math.Sin(moon.geo_eclip_lat)*Math.Cos(I));

            // Physical librations
            double k1 = DEG2RAD*(119.75 + 131.849*t);
            double k2 = DEG2RAD*(72.56 + 20.186*t);

            double rho = (
                -0.02752*Math.Cos(mdash) +
                -0.02245*Math.Sin(f) +
                +0.00684*Math.Cos(mdash - 2*f) +
                -0.00293*Math.Cos(2*f) +
                -0.00085*Math.Cos(2*f - 2*d) +
                -0.00054*Math.Cos(mdash - 2*d) +
                -0.00020*Math.Sin(mdash + f) +
                -0.00020*Math.Cos(mdash + 2*f) +
                -0.00020*Math.Cos(mdash - f) +
                +0.00014*Math.Cos(mdash + 2*f - 2*d)
            );

            double sigma = (
                -0.02816*Math.Sin(mdash) +
                +0.02244*Math.Cos(f) +
                -0.00682*Math.Sin(mdash - 2*f) +
                -0.00279*Math.Sin(2*f) +
                -0.00083*Math.Sin(2*f - 2*d) +
                +0.00069*Math.Sin(mdash - 2*d) +
                +0.00040*Math.Cos(mdash + f) +
                -0.00025*Math.Sin(2*mdash) +
                -0.00023*Math.Sin(mdash + 2*f) +
                +0.00020*Math.Cos(mdash - f) +
                +0.00019*Math.Sin(mdash - f) +
                +0.00013*Math.Sin(mdash + 2*f - 2*d) +
                -0.00010*Math.Cos(mdash - 3*f)
            );

            double tau = (
                +0.02520*e*Math.Sin(m) +
                +0.00473*Math.Sin(2*mdash - 2*f) +
                -0.00467*Math.Sin(mdash) +
                +0.00396*Math.Sin(k1) +
                +0.00276*Math.Sin(2*mdash - 2*d) +
                +0.00196*Math.Sin(omega) +
                -0.00183*Math.Cos(mdash - f) +
                +0.00115*Math.Sin(mdash - 2*d) +
                -0.00096*Math.Sin(mdash - d) +
                +0.00046*Math.Sin(2*f - 2*d) +
                -0.00039*Math.Sin(mdash - f) +
                -0.00032*Math.Sin(mdash - m - d) +
                +0.00027*Math.Sin(2*mdash - m - 2*d) +
                +0.00023*Math.Sin(k2) +
                -0.00014*Math.Sin(2*d) +
                +0.00014*Math.Cos(2*mdash - 2*f) +
                -0.00012*Math.Sin(mdash - 2*f) +
                -0.00012*Math.Sin(2*mdash) +
                +0.00011*Math.Sin(2*mdash - 2*m - 2*d)
            );

            double ldash2 = -tau + (rho*Math.Cos(a) + sigma*Math.Sin(a))*Math.Tan(bdash);
            bdash *= RAD2DEG;
            double bdash2 = sigma*Math.Cos(a) - rho*Math.Sin(a);

            lib.elon = ldash + ldash2;
            lib.elat = bdash + bdash2;

            return lib;
        }

        private static AstroVector BarycenterContrib(AstroTime time, Body body, double planet_gm)
        {
            AstroVector p = CalcVsop(vsop[(int)body], time);
            return (planet_gm / (planet_gm + SUN_GM)) * p;
        }

        private static AstroVector CalcSolarSystemBarycenter(AstroTime time)
        {
            AstroVector j = BarycenterContrib(time, Body.Jupiter, JUPITER_GM);
            AstroVector s = BarycenterContrib(time, Body.Saturn,  SATURN_GM);
            AstroVector u = BarycenterContrib(time, Body.Uranus,  URANUS_GM);
            AstroVector n = BarycenterContrib(time, Body.Neptune, NEPTUNE_GM);
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
        /// If given an invalid value for `body`, this function will throw an #InvalidBodyException.
        /// </remarks>
        /// <param name="body">
        /// A body for which to calculate a heliocentric position:
        /// the Sun, Moon, EMB, SSB, or any of the planets.
        /// Also allowed to be a user-defined star created by #Astronomy.DefineStar.
        /// </param>
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
                    StateVector planet = CalcPluto(time, true);
                    return new AstroVector(planet.x, planet.y, planet.z, time);

                case Body.Moon:
                    geomoon = GeoMoon(time);
                    earth = CalcEarth(time);
                    return earth + geomoon;

                case Body.EMB:
                    geomoon = GeoMoon(time);
                    earth = CalcEarth(time);
                    return earth + (geomoon / (1.0 + EARTH_MOON_MASS_RATIO));

                case Body.SSB:
                    return CalcSolarSystemBarycenter(time);

                default:
                    if (UserDefinedStar(body) is StarDef star)
                        return VectorFromSphere(new Spherical(star.dec, 15*star.ra, star.dist), time);
                    throw new InvalidBodyException(body);
            }
        }

        /// <summary>
        /// Calculates the distance between a body and the Sun at a given time.
        /// </summary>
        /// <remarks>
        /// Given a date and time, this function calculates the distance between
        /// the center of `body` and the center of the Sun, expressed in AU.
        /// For the planets Mercury through Neptune, this function is significantly
        /// more efficient than calling #Astronomy.HelioVector followed by taking the length
        /// of the resulting vector.
        /// </remarks>
        /// <param name="body">
        /// A body for which to calculate a heliocentric distance:
        /// the Sun, Moon, EMB, SSB, any of the planets, or a user-defined star.
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
                    return VsopFormulaCalc(vsop[(int)body].rad, time.tt / DAYS_PER_MILLENNIUM, false);

                default:
                    if (UserDefinedStar(body) is StarDef star)
                        return star.dist;

                    // Fall back to taking the length of the heliocentric vector.
                    return HelioVector(body, time).Length();
            }
        }

        private static AstroVector CalcEarth(AstroTime time)
        {
            return CalcVsop(vsop[(int)Body.Earth], time);
        }

        /// <summary>
        /// Solve for light travel time of a vector function.
        /// </summary>
        /// <remarks>
        /// When observing a distant object, for example Jupiter as seen from Earth,
        /// the amount of time it takes for light to travel from the object to the
        /// observer can significantly affect the object's apparent position.
        /// This function is a generic solver that figures out how long in the
        /// past light must have left the observed object to reach the observer
        /// at the specified observation time. It uses #IPositionFunction
        /// to express an arbitrary position vector as a function of time.
        ///
        /// This function repeatedly calls `func.Position`, passing a series of time
        /// estimates in the past. Then `func.Position` must return a relative state vector between
        /// the observer and the target. `CorrectLightTravel` keeps calling
        /// `func.Position` with more and more refined estimates of the time light must have
        /// left the target to arrive at the observer.
        ///
        /// For common use cases, it is simpler to use #Astronomy.BackdatePosition
        /// for calculating the light travel time correction of one body observing another body.
        ///
        /// For geocentric calculations, #Astronomy.GeoVector also backdates the returned
        /// position vector for light travel time, only it returns the observation time in
        /// the returned vector's `t` field rather than the backdated time.
        ///
        /// </remarks>
        /// <param name="func">An arbitrary position vector as a function of time.</param>
        /// <param name="time">The observation time for which to solve for light travel delay.</param>
        /// <returns>
        /// The position vector at the solved backdated time.
        /// The `t` field holds the time that light left the observed
        /// body to arrive at the observer at the observation time.
        /// </returns>
        public static AstroVector CorrectLightTravel(IPositionFunction func, AstroTime time)
        {
            AstroTime ltime = time;
            for (int iter = 0; iter < 10; ++iter)
            {
                AstroVector pos = func.Position(ltime);
                // This solver does not support more than one light-day of distance,
                // because that would cause convergence problems and inaccurate
                // values for stellar aberration angles.
                double lt = pos.Length() / C_AUDAY;
                if (lt > 1.0)
                    throw new ArgumentException("Object is too distant for light-travel solver.");
                AstroTime ltime2 = time.AddDays(-lt);
                double dt = Math.Abs(ltime2.tt - ltime.tt);
                if (dt < 1.0e-9)        // 86.4 microseconds
                    return pos;
                ltime = ltime2;
            }
            throw new InternalError("Light travel time correction did not converge.");
        }

        internal struct BodyPosition: IPositionFunction
        {
            private Body observerBody;
            private Body targetBody;
            private Aberration aberration;
            private AstroVector observerPos;    // used only when aberration == Aberration.None

            public BodyPosition(Body observerBody, Body targetBody, Aberration aberration, AstroVector observerPos)
            {
                this.observerBody = observerBody;
                this.targetBody = targetBody;
                this.aberration = aberration;
                this.observerPos = observerPos;
            }

            public AstroVector Position(AstroTime time)
            {
                if (aberration == Aberration.None)
                {
                    // No aberration, so use the pre-calculated initial position of
                    // the observer body that is already stored in `observerPos`.
                    // To avoid an exception in the subtraction below, patch the time.
                    observerPos.t = time;
                }
                else
                {
                    // The following discussion is worded with the observer body being the Earth,
                    // which is often the case. However, the same reasoning applies to any observer body
                    // without loss of generality.
                    //
                    // To include aberration, make a good first-order approximation
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

                    observerPos = Astronomy.HelioVector(observerBody, time);
                }

                // Subtract the bodies' heliocentric positions to obtain a relative position vector.
                return Astronomy.HelioVector(targetBody, time) - observerPos;
            }
        }

        /// <summary>
        /// Solve for light travel time correction of apparent position.
        /// </summary>
        /// <remarks>
        /// When observing a distant object, for example Jupiter as seen from Earth,
        /// the amount of time it takes for light to travel from the object to the
        /// observer can significantly affect the object's apparent position.
        ///
        /// This function solves the light travel time correction for the apparent
        /// relative position vector of a target body as seen by an observer body
        /// at a given observation time.
        ///
        /// For geocentric calculations, #Astronomy.GeoVector also includes light
        /// travel time correction, but the time `t` embedded in its returned vector
        /// refers to the observation time, not the backdated time that light left
        /// the observed body. Thus `BackdatePosition` provides direct
        /// access to the light departure time for callers that need it.
        ///
        /// For a more generalized light travel correction solver, see #Astronomy.CorrectLightTravel.
        /// </remarks>
        /// <param name="time">The time of observation.</param>
        /// <param name="observerBody">The body to be used as the observation location.</param>
        /// <param name="targetBody">The body to be observed.</param>
        /// <param name="aberration">`Aberration.Corrected` to correct for aberration, or `Aberration.None` to leave uncorrected.</param>
        /// <returns>
        /// The position vector at the solved backdated time.
        /// Its `t` field holds the time that light left the observed
        /// body to arrive at the observer at the observation time.
        /// </returns>
        public static AstroVector BackdatePosition(
            AstroTime time,
            Body observerBody,
            Body targetBody,
            Aberration aberration)
        {
            if (null != UserDefinedStar(targetBody))
            {
                // This is a user-defined star, which must be treated as a special case.
                // First, we assume its heliocentric position does not change with time.
                // Second, we assume its heliocentric position has already been corrected
                // for light-travel time, its coordinates given as it appears on Earth at the present.
                // Therefore, no backdating is applied.
                AstroVector tvec = HelioVector(targetBody, time);
                switch (aberration)
                {
                    case Aberration.None:
                        // No correction is needed. Simply return the star's current position as seen from the observer.
                        return tvec - HelioVector(observerBody, time);

                    case Aberration.Corrected:
                        // (Observer velocity) - (light vector) = (Aberration-corrected direction to target body).
                        // Note that this is an approximation, because technically the light vector should
                        // be measured in barycentric coordinates, not heliocentric. The error is very small.
                        StateVector ostate = HelioState(observerBody, time);
                        AstroVector rvec = tvec - ostate.Position();
                        double s = C_AUDAY / rvec.Length();    // conversion factor from relative distance to speed of light
                        return rvec + ostate.Velocity()/s;

                    default:
                        throw new ArgumentException($"Unsupported aberration option: {aberration}");
                }
            }

            AstroVector observerPos;
            switch (aberration)
            {
                case Aberration.None:
                    // Without aberration, we need the observer body position at the observation time only.
                    // For efficiency, calculate it once and hold onto it, so `BodyPosition` can keep using it.
                    observerPos = HelioVector(observerBody, time);
                    break;

                case Aberration.Corrected:
                    // With aberration, `BackdatePosition` will calculate `observerPos` at different times.
                    // Therefore, do not waste time calculating it now.
                    // Provide a placeholder value.
                    observerPos = new AstroVector();
                    break;

                default:
                    throw new ArgumentException($"Unsupported aberration option: {aberration}");
            }
            var func = new BodyPosition(observerBody, targetBody, aberration, observerPos);
            return CorrectLightTravel(func, time);
        }

        /// <summary>
        /// Calculates geocentric Cartesian coordinates of a body in the J2000 equatorial system.
        /// </summary>
        /// <remarks>
        /// This function calculates the position of the given celestial body as a vector,
        /// using the center of the Earth as the origin.  The result is expressed as a Cartesian
        /// vector in the J2000 equatorial system: the coordinates are based on the mean equator
        /// of the Earth at noon UTC on 1 January 2000.
        ///
        /// If given an invalid value for `body`, this function will throw an exception.
        ///
        /// Unlike #Astronomy.HelioVector, this function corrects for light travel time.
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
            switch (body)
            {
            case Body.Earth:
                // The Earth's geocentric coordinates are always (0,0,0).
                return new AstroVector(0.0, 0.0, 0.0, time);

            case Body.Moon:
                // The moon is so close, aberration and light travel time don't matter.
                return GeoMoon(time);

            default:
                // For all other bodies, apply light travel time correction.
                AstroVector vector = BackdatePosition(time, Body.Earth, body, aberration);
                vector.t = time;    // tricky: return the observation time, not the backdated time.
                return vector;
            }
        }

        internal static StateVector ExportState(body_state_t terse, AstroTime time)
        {
            return new StateVector(
                terse.r.x, terse.r.y, terse.r.z,
                terse.v.x, terse.v.y, terse.v.z,
                time
            );
        }

        /// <summary>
        /// Calculates barycentric position and velocity vectors for the given body.
        /// </summary>
        /// <remarks>
        /// Given a body and a time, calculates the barycentric position and velocity
        /// vectors for the center of that body at that time.
        /// The vectors are expressed in J2000 mean equator coordinates (EQJ).
        /// </remarks>
        /// <param name="body">
        /// The celestial body whose barycentric state vector is to be calculated.
        /// Supported values are `Body.Sun`, `Body.Moon`, `Body.EMB`, `Body.SSB`, and all planets:
        /// `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`,
        /// `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
        /// </param>
        /// <param name="time">
        /// The date and time for which to calculate position and velocity.
        /// </param>
        /// <returns>
        /// A structure that contains barycentric position and velocity vectors.
        /// </returns>
        public static StateVector BaryState(Body body, AstroTime time)
        {
            // Trivial case: the solar system barycenter itself.
            if (body == Body.SSB)
                return new StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time);

            if (body == Body.Pluto)
                return CalcPluto(time, false);

            // Find the barycentric positions and velocities for the 5 major bodies.
            major_bodies_t bary = MajorBodyBary(time.tt);

            // If the caller is asking for one of the major bodies, we can immediately return the answer.
            switch (body)
            {
                case Body.Sun:      return ExportState(bary.Sun, time);
                case Body.Jupiter:  return ExportState(bary.Jupiter, time);
                case Body.Saturn:   return ExportState(bary.Saturn, time);
                case Body.Uranus:   return ExportState(bary.Uranus, time);
                case Body.Neptune:  return ExportState(bary.Neptune, time);

                case Body.Moon:
                case Body.EMB:
                    body_state_t earth = CalcVsopPosVel(vsop[(int)Body.Earth], time.tt);
                    StateVector state;
                    if (body == Body.Moon)
                        state = GeoMoonState(time);
                    else
                        state = GeoEmbState(time);

                    return new StateVector(
                        state.x  + bary.Sun.r.x + earth.r.x,
                        state.y  + bary.Sun.r.y + earth.r.y,
                        state.z  + bary.Sun.r.z + earth.r.z,
                        state.vx + bary.Sun.v.x + earth.v.x,
                        state.vy + bary.Sun.v.y + earth.v.y,
                        state.vz + bary.Sun.v.z + earth.v.z,
                        time
                    );
            }

            // Handle the remaining VSOP bodies: Mercury, Venus, Earth, Mars.
            // BarySun + HelioBody = BaryBody
            int bindex = (int)body;
            if (bindex >= 0 && bindex < vsop.Length)
            {
                body_state_t planet = CalcVsopPosVel(vsop[bindex], time.tt);
                return new StateVector(
                    bary.Sun.r.x + planet.r.x,
                    bary.Sun.r.y + planet.r.y,
                    bary.Sun.r.z + planet.r.z,
                    bary.Sun.v.x + planet.v.x,
                    bary.Sun.v.y + planet.v.y,
                    bary.Sun.v.z + planet.v.z,
                    time
                );
            }

            throw new InvalidBodyException(body);
        }

        /// <summary>
        /// Calculates heliocentric position and velocity vectors for the given body.
        /// </summary>
        /// <remarks>
        /// Given a body and a time, calculates the position and velocity
        /// vectors for the center of that body at that time, relative to the center of the Sun.
        /// The vectors are expressed in J2000 mean equator coordinates (EQJ).
        /// If you need the position vector only, it is more efficient to call #Astronomy.HelioVector.
        /// The Sun's center is a non-inertial frame of reference. In other words, the Sun
        /// experiences acceleration due to gravitational forces, mostly from the larger
        /// planets (Jupiter, Saturn, Uranus, and Neptune). If you want to calculate momentum,
        /// kinetic energy, or other quantities that require a non-accelerating frame
        /// of reference, consider using #Astronomy.BaryState instead.
        /// </remarks>
        /// <param name="body">
        /// The celestial body whose heliocentric state vector is to be calculated.
        /// Supported values are `Body.Sun`, `Body.Moon`, `Body.EMB`, `Body.SSB`, and all planets:
        /// `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`, `Body.Jupiter`,
        /// `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
        /// Also allowed to be a user-defined star created by #Astronomy.DefineStar.
        /// </param>
        /// <param name="time">
        /// The date and time for which to calculate position and velocity.
        /// </param>
        /// <returns>
        /// A structure that contains heliocentric position and velocity vectors.
        /// </returns>
        public static StateVector HelioState(Body body, AstroTime time)
        {
            switch (body)
            {
                case Body.Sun:
                    // Trivial case: the Sun is the origin of the heliocentric frame.
                    return new StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time);

                case Body.SSB:
                    // Calculate the barycentric Sun. Then the negative of that is the heliocentric SSB.
                    major_bodies_t bary = MajorBodyBary(time.tt);
                    return new StateVector(
                        -bary.Sun.r.x,
                        -bary.Sun.r.y,
                        -bary.Sun.r.z,
                        -bary.Sun.v.x,
                        -bary.Sun.v.y,
                        -bary.Sun.v.z,
                        time
                    );

                case Body.Mercury:
                case Body.Venus:
                case Body.Earth:
                case Body.Mars:
                case Body.Jupiter:
                case Body.Saturn:
                case Body.Uranus:
                case Body.Neptune:
                    // Planets included in the VSOP87 model. */
                    body_state_t planet = CalcVsopPosVel(vsop[(int)body], time.tt);
                    return ExportState(planet, time);

                case Body.Pluto:
                    return CalcPluto(time, true);

                case Body.Moon:
                case Body.EMB:
                    body_state_t earth = CalcVsopPosVel(vsop[(int)Body.Earth], time.tt);
                    StateVector state = (body == Body.Moon) ? GeoMoonState(time) : GeoEmbState(time);
                    return new StateVector(
                        state.x  + earth.r.x,
                        state.y  + earth.r.y,
                        state.z  + earth.r.z,
                        state.vx + earth.v.x,
                        state.vy + earth.v.y,
                        state.vz + earth.v.z,
                        time
                    );

                default:
                    if (null != UserDefinedStar(body))
                    {
                        AstroVector vec = HelioVector(body, time);
                        return new StateVector(vec.x, vec.y, vec.z, 0.0, 0.0, 0.0, time);
                    }
                    throw new InvalidBodyException(body);
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
        /// <returns>Topocentric equatorial coordinates of the celestial body.</returns>
        public static Equatorial Equator(
            Body body,
            AstroTime time,
            Observer observer,
            EquatorEpoch equdate,
            Aberration aberration)
        {
            AstroVector gc_observer = geo_pos(time, observer);
            AstroVector gc = GeoVector(body, time, aberration);
            AstroVector j2000 = gc - gc_observer;

            switch (equdate)
            {
                case EquatorEpoch.OfDate:
                    AstroVector datevect = gyration(j2000, PrecessDirection.From2000);
                    return EquatorFromVector(datevect);

                case EquatorEpoch.J2000:
                    return EquatorFromVector(j2000);

                default:
                    throw new ArgumentException(string.Format("Unsupported equator epoch {0}", equdate));
            }
        }

        /// <summary>
        /// Calculates geocentric equatorial coordinates of an observer on the surface of the Earth.
        /// </summary>
        ///
        /// <remarks>
        /// This function calculates a vector from the center of the Earth to
        /// a point on or near the surface of the Earth, expressed in equatorial
        /// coordinates. It takes into account the rotation of the Earth at the given
        /// time, along with the given latitude, longitude, and elevation of the observer.
        ///
        /// The caller may pass a value in `equdate` to select either `EquatorEpoch.J2000`
        /// for using J2000 coordinates, or `EquatorEpoch.OfDate` for using coordinates relative
        /// to the Earth's equator at the specified time.
        ///
        /// The returned vector has components expressed in astronomical units (AU).
        /// To convert to kilometers, multiply the `x`, `y`, and `z` values by
        /// the constant value #Astronomy.KM_PER_AU.
        ///
        /// The inverse of this function is also available: #Astronomy.VectorObserver.
        /// </remarks>
        ///
        /// <param name="time">
        /// The date and time for which to calculate the observer's position vector.
        /// </param>
        ///
        /// <param name="observer">
        /// The geographic location of a point on or near the surface of the Earth.
        /// </param>
        ///
        /// <param name="equdate">
        /// Selects the date of the Earth's equator in which to express the equatorial coordinates.
        /// The caller may select `EquatorEpoch.J2000` to use the orientation of the Earth's equator
        /// at noon UTC on January 1, 2000, in which case this function corrects for precession
        /// and nutation of the Earth as it was at the moment specified by the `time` parameter.
        /// Or the caller may select `EquatorEpoch.OfDate` to use the Earth's equator at `time`
        /// as the orientation.
        /// </param>
        ///
        /// <returns>
        /// An equatorial vector from the center of the Earth to the specified location
        /// on (or near) the Earth's surface.
        /// </returns>
        public static AstroVector ObserverVector(
            AstroTime time,
            Observer observer,
            EquatorEpoch equdate)
        {
            return ObserverState(time, observer, equdate).Position();
        }

        /// <summary>
        /// Calculates geocentric equatorial position and velocity of an observer on the surface of the Earth.
        /// </summary>
        ///
        /// <remarks>
        /// This function calculates position and velocity vectors of an observer
        /// on or near the surface of the Earth, expressed in equatorial
        /// coordinates. It takes into account the rotation of the Earth at the given
        /// time, along with the given latitude, longitude, and elevation of the observer.
        ///
        /// The caller may pass a value in `equdate` to select either `EquatorEpoch.J2000`
        /// for using J2000 coordinates, or `EquatorEpoch.OfDate` for using coordinates relative
        /// to the Earth's equator at the specified time.
        ///
        /// The returned position vector has components expressed in astronomical units (AU).
        /// To convert to kilometers, multiply the `x`, `y`, and `z` values by
        /// the constant value #Astronomy.KM_PER_AU.
        ///
        /// The returned velocity vector is measured in AU/day.
        /// </remarks>
        ///
        /// <param name="time">
        /// The date and time for which to calculate the observer's geocentric state vector.
        /// </param>
        ///
        /// <param name="observer">
        /// The geographic location of a point on or near the surface of the Earth.
        /// </param>
        ///
        /// <param name="equdate">
        /// Selects the date of the Earth's equator in which to express the equatorial coordinates.
        /// The caller may select `EquatorEpoch.J2000` to use the orientation of the Earth's equator
        /// at noon UTC on January 1, 2000, in which case this function corrects for precession
        /// and nutation of the Earth as it was at the moment specified by the `time` parameter.
        /// Or the caller may select `EquatorEpoch.OfDate` to use the Earth's equator at `time`
        /// as the orientation.
        /// </param>
        ///
        /// <returns>
        /// The position and velocity of the given geographic location, relative to the center of the Earth.
        /// </returns>
        public static StateVector ObserverState(
            AstroTime time,
            Observer observer,
            EquatorEpoch equdate)
        {
            StateVector state = terra(observer, time);

            if (equdate == EquatorEpoch.OfDate)
                return state;

            if (equdate == EquatorEpoch.J2000)
                return gyration_posvel(state, PrecessDirection.Into2000);

            throw new ArgumentException(string.Format("Unsupported equator epoch {0}", equdate));
        }

        /// <summary>
        /// Calculates the geographic location corresponding to an equatorial vector.
        /// </summary>
        ///
        /// <remarks>
        /// This is the inverse function of #Astronomy.ObserverVector.
        /// Given a geocentric equatorial vector, it returns the geographic
        /// latitude, longitude, and elevation for that vector.
        /// </remarks>
        ///
        /// <param name="vector">
        /// The geocentric equatorial position vector for which to find geographic coordinates.
        /// The components are expressed in Astronomical Units (AU).
        /// You can calculate AU by dividing kilometers by the constant #Astronomy.KM_PER_AU.
        /// The time `vector.t` determines the Earth's rotation.
        /// </param>
        ///
        /// <param name="equdate">
        /// Selects the date of the Earth's equator in which `vector` is expressed.
        /// The caller may select `EquatorEpoch.J2000` to use the orientation of the Earth's equator
        /// at noon UTC on January 1, 2000, in which case this function corrects for precession
        /// and nutation of the Earth as it was at the moment specified by `vector.t`.
        /// Or the caller may select `EquatorEpoch.OfDate` to use the Earth's equator at `vector.t`
        /// as the orientation.
        /// </param>
        ///
        /// <returns>
        /// The geographic latitude, longitude, and elevation above sea level
        /// that corresponds to the given equatorial vector.
        /// </returns>
        public static Observer VectorObserver(
            AstroVector vector,
            EquatorEpoch equdate)
        {
            if (equdate == EquatorEpoch.J2000)
                vector = gyration(vector, PrecessDirection.From2000);
            return inverse_terra(vector);
        }

        /// <summary>
        /// Calculates the gravitational acceleration experienced by an observer on the Earth.
        /// </summary>
        /// <remarks>
        /// This function implements the WGS 84 Ellipsoidal Gravity Formula.
        /// The result is a combination of inward gravitational acceleration
        /// with outward centrifugal acceleration, as experienced by an observer
        /// in the Earth's rotating frame of reference.
        /// The resulting value increases toward the Earth's poles and decreases
        /// toward the equator, consistent with changes of the weight measured
        /// by a spring scale of a fixed mass moved to different latitudes and heights
        /// on the Earth.
        /// </remarks>
        /// <param name="latitude">
        /// The latitude of the observer in degrees north or south of the equator.
        /// By formula symmetry, positive latitudes give the same answer as negative
        /// latitudes, so the sign does not matter.
        /// </param>
        /// <param name="height">
        /// The height above the sea level geoid in meters.
        /// No range checking is done; however, accuracy is only valid in the
        /// range 0 to 100000 meters.
        /// </param>
        /// <returns>
        /// The effective gravitational acceleration expressed in meters per second squared [m/s^2].
        /// </returns>
        public static double ObserverGravity(double latitude, double height)
        {
            double s = Math.Sin(latitude * DEG2RAD);
            double s2 = s*s;
            double g0 = 9.7803253359 * (1.0 + 0.00193185265241*s2) / Math.Sqrt(1.0 - 0.00669437999013*s2);
            return g0 * (1.0 - (3.15704e-07 - 2.10269e-09*s2)*height + 7.37452e-14*height*height);
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
            double sinra = Math.Sin(ra * HOUR2RAD);
            double cosra = Math.Cos(ra * HOUR2RAD);

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
            var uze = new AstroVector(coslat * coslon, coslat * sinlon, sinlat, time);
            var une = new AstroVector(-sinlat * coslon, -sinlat * sinlon, coslat, time);
            var uwe = new AstroVector(sinlon, -coslon, 0.0, time);

            // Correct the vectors uze, une, uwe for the Earth's rotation by calculating
            // sidereal time. Call spin() for each uncorrected vector to rotate about
            // the Earth's axis to yield corrected unit vectors uz, un, uw.
            // Multiply sidereal hours by -15 to convert to degrees and flip eastward
            // rotation of the Earth to westward apparent movement of objects with time.
            double angle = -15.0 * SiderealTime(time);
            AstroVector uz = spin(angle, uze);
            AstroVector un = spin(angle, une);
            AstroVector uw = spin(angle, uwe);

            // Convert angular equatorial coordinates (RA, DEC) to
            // cartesian equatorial coordinates in 'p', using the
            // same orientation system as uze, une, uwe.
            var p = new AstroVector(cosdc * cosra, cosdc * sinra, sindc, time);

            // Use dot products of p with the zenith, north, and west
            // vectors to obtain the cartesian coordinates of the body in
            // the observer's horizontal orientation system.
            // pz = zenith component [-1, +1]
            // pn = north  component [-1, +1]
            // pw = west   component [-1, +1]
            double pz = p.x*uz.x + p.y*uz.y + p.z*uz.z;
            double pn = p.x*un.x + p.y*un.y + p.z*un.z;
            double pw = p.x*uw.x + p.y*uw.y + p.z*uw.z;

            // proj is the "shadow" of the body vector along the observer's flat ground.
            double proj = hypot(pn, pw);

            // Calculate az = azimuth (compass direction clockwise from East.)
            double az;
            if (proj > 0.0)
            {
                // If the body is not exactly straight up/down, it has an azimuth.
                // Invert the angle to produce degrees eastward from north.
                az = -Math.Atan2(pw, pn) * RAD2DEG;
                if (az < 0.0)
                    az += 360.0;
            }
            else
            {
                // The body is straight up/down, so it does not have an azimuth.
                // Report an arbitrary but reasonable value.
                az = 0.0;
            }

            // zd = the angle of the body away from the observer's zenith, in degrees.
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

                    proj = hypot(prx, pry);
                    if (proj > 0.0)
                    {
                        hor_ra = RAD2HOUR * Math.Atan2(pry, prx);
                        if (hor_ra < 0.0)
                            hor_ra += 24.0;
                    }
                    else
                    {
                        hor_ra = 0.0;
                    }
                    hor_dec = RAD2DEG * Math.Atan2(prz, proj);
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
            // Correct for light travel time from the Sun.
            // Otherwise season calculations (equinox, solstice) will all be early by about 8 minutes!
            AstroTime adjusted_time = time.AddDays(-1.0 / C_AUDAY);

            AstroVector earth2000 = CalcEarth(adjusted_time);

            // Convert heliocentric location of Earth to geocentric location of Sun.
            AstroVector sun2000 = new AstroVector(-earth2000.x, -earth2000.y, -earth2000.z, adjusted_time);

            // Convert to equatorial Cartesian coordinates of date.
            AstroVector sun_ofdate = gyration(sun2000, PrecessDirection.From2000);

            // Convert equatorial coordinates to ecliptic coordinates.
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

            double xyproj = hypot(ex, ey);
            double elon = 0.0;
            if (xyproj > 0.0)
            {
                elon = RAD2DEG * Math.Atan2(ey, ex);
                if (elon < 0.0)
                    elon += 360.0;
            }

            double elat = RAD2DEG * Math.Atan2(ez, xyproj);

            var vec = new AstroVector(ex, ey, ez, pos.t);
            return new Ecliptic(vec, elat, elon);
        }

        /// <summary>
        /// Converts a J2000 mean equator (EQJ) vector to a true ecliptic of date (ETC) vector and angles.
        /// </summary>
        /// <remarks>
        /// Given coordinates relative to the Earth's equator at J2000 (the instant of noon UTC
        /// on 1 January 2000), this function converts those coordinates to true ecliptic coordinates of date,
        /// which are relative to the plane of the Earth's orbit around the Sun.
        /// </remarks>
        /// <param name="eqj">
        /// Equatorial coordinates in the J2000 frame of reference.
        /// You can call #Astronomy.GeoVector to obtain suitable equatorial coordinates.
        /// </param>
        /// <returns>Spherical and vector coordinates expressed in true ecliptic coordinates of date (ECT).</returns>
        public static Ecliptic EquatorialToEcliptic(AstroVector eqj)
        {
            // Calculate nutation and obliquity for this time.
            // As an optimization, the nutation angles are cached in `eqj.t`,
            // and reused below when the `nutation` function is called.
            earth_tilt_t et = e_tilt(eqj.t);

            // Convert J2000 mean equator (EQJ) to true equator of date (EQD).
            AstroVector mean_pos = precession(eqj, PrecessDirection.From2000);
            AstroVector eqd = nutation(mean_pos, PrecessDirection.From2000);

            // Rotate from EQD to true ecliptic of date (ECT).
            return RotateEquatorialToEcliptic(eqd, et.tobl * DEG2RAD);
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
            // https://github.com/cosinekitty/astronomy/issues/187
            // Solstices and equinoxes drift over long spans of time,
            // due to precession of the Earth's axis.
            // Therefore, we have to search a wider range of time than
            // one might expect. It turns out this has very little
            // effect on efficiency, thanks to the quick convergence
            // of quadratic interpolation inside the `Search` function.

            return new SeasonsInfo(
                FindSeasonChange(  0, year,  3, 10),
                FindSeasonChange( 90, year,  6, 10),
                FindSeasonChange(180, year,  9, 10),
                FindSeasonChange(270, year, 12, 10)
            );
        }

        private static AstroTime FindSeasonChange(double targetLon, int year, int month, int day)
        {
            var startTime = new AstroTime(year, month, day, 0, 0, 0);
            return SearchSunLongitude(targetLon, startTime, 20.0) ??
                throw new InternalError($"Cannot find solution for Sun longitude {targetLon} in year {year}");
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
        /// The search will return `null` if the Sun never reaches the longitude `targetLon` or
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
            return Search(sun_offset, startTime, t2, 0.01);
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
                    throw new InternalError(string.Format("Search did not converge within {0} iterations.", iter_limit));

                double dt = (t2.tt - t1.tt) / 2.0;
                AstroTime tmid = t1.AddDays(dt);
                if (Math.Abs(dt) < dt_days)
                {
                    // We are close enough to the event to stop the search.
                    return tmid;
                }

                if (calc_fmid)
                    fmid = func.Eval(tmid);
                else
                    calc_fmid = true;   // we already have the correct value of fmid from the previous loop

                // Quadratic interpolation:
                // Try to find a parabola that passes through the 3 points we have sampled:
                // (t1,f1), (tmid,fmid), (t2,f2)

                double q_ut, q_df_dt;
                if (QuadInterp(tmid.ut, t2.ut - tmid.ut, f1, fmid, f2, out q_ut, out q_df_dt))
                {
                    var tq = new AstroTime(q_ut);
                    double fq = func.Eval(tq);
                    if (q_df_dt != 0.0)
                    {
                        double dt_guess = Math.Abs(fq / q_df_dt);
                        if (dt_guess < dt_days)
                        {
                            // The estimated time error is small enough that we can quit now.
                            return tq;
                        }

                        // Try guessing a tighter boundary with the interpolated root at the center.
                        dt_guess *= 1.2;
                        if (dt_guess < dt/10.0)
                        {
                            AstroTime tleft = tq.AddDays(-dt_guess);
                            AstroTime tright = tq.AddDays(+dt_guess);
                            if ((tleft.ut - t1.ut)*(tleft.ut - t2.ut) < 0.0)
                            {
                                if ((tright.ut - t1.ut)*(tright.ut - t2.ut) < 0.0)
                                {
                                    double fleft = func.Eval(tleft);
                                    double fright = func.Eval(tright);
                                    if (fleft<0.0 && fright>=0.0)
                                    {
                                        f1 = fleft;
                                        f2 = fright;
                                        t1 = tleft;
                                        t2 = tright;
                                        fmid = fq;
                                        calc_fmid = false;  // save a little work -- no need to re-calculate fmid next time around the loop
                                        continue;
                                    }
                                }
                            }
                        }
                    }
                }

                // After quadratic interpolation attempt.
                // Now just divide the region in two parts and pick whichever one appears to contain a root.
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

                // Either there is no ascending zero-crossing in this range
                // or the search window is too wide (more than one zero-crossing).
                return null;
            }
        }

        private static bool QuadInterp(
            double tm, double dt, double fa, double fm, double fb,
            out double out_t, out double out_df_dt)
        {
            double Q, R, S;
            double u, ru;
            double x, x1, x2;

            out_t = out_df_dt = 0.0;

            Q = (fb + fa)/2.0 - fm;
            R = (fb - fa)/2.0;
            S = fm;

            if (Q == 0.0)
            {
                // This is a line, not a parabola.
                if (R == 0.0)
                    return false;       // This is a HORIZONTAL line... can't make progress!
                x = -S / R;
                if (x < -1.0 || x > +1.0)
                    return false;   // out of bounds
            }
            else
            {
                // This really is a parabola. Find roots x1, x2.
                u = R*R - 4*Q*S;
                if (u <= 0.0)
                    return false;   // can't solve if imaginary, or if vertex of parabola is tangent.

                ru = Math.Sqrt(u);
                x1 = (-R + ru) / (2.0 * Q);
                x2 = (-R - ru) / (2.0 * Q);
                if (-1.0 <= x1 && x1 <= +1.0)
                {
                    if (-1.0 <= x2 && x2 <= +1.0)
                        return false;   // two roots are within bounds; we require a unique zero-crossing.
                    x = x1;
                }
                else if (-1.0 <= x2 && x2 <= +1.0)
                    x = x2;
                else
                    return false;   // neither root is within bounds
            }

            out_t = tm + x*dt;
            out_df_dt = (2*Q*x + R) / dt;
            return true;   // success
        }

        ///
        /// <summary>
        /// Returns one body's ecliptic longitude with respect to another, as seen from the Earth.
        /// </summary>
        /// <remarks>
        /// This function determines where one body appears around the ecliptic plane
        /// (the plane of the Earth's orbit around the Sun) as seen from the Earth,
        /// relative to the another body's apparent position.
        /// The function returns an angle in the half-open range [0, 360) degrees.
        /// The value is the ecliptic longitude of `body1` relative to the ecliptic
        /// longitude of `body2`.
        ///
        /// The angle is 0 when the two bodies are at the same ecliptic longitude
        /// as seen from the Earth. The angle increases in the prograde direction
        /// (the direction that the planets orbit the Sun and the Moon orbits the Earth).
        ///
        /// When the angle is 180 degrees, it means the two bodies appear on opposite sides
        /// of the sky for an Earthly observer.
        ///
        /// Neither `body1` nor `body2` is allowed to be `Body.Earth`.
        /// If this happens, the function throws an exception.
        /// </remarks>
        /// <param name="body1">The first body, whose longitude is to be found relative to the second body.</param>
        /// <param name="body2">The second body, relative to which the longitude of the first body is to be found.</param>
        /// <param name="time">The date and time of the observation.</param>
        /// <returns>
        /// An angle in the range [0, 360), expressed in degrees.
        /// </returns>
        public static double PairLongitude(Body body1, Body body2, AstroTime time)
        {
            if (body1 == Body.Earth || body2 == Body.Earth)
                throw new EarthNotAllowedException();

            AstroVector vector1 = GeoVector(body1, time, Aberration.None);
            Ecliptic eclip1 = EquatorialToEcliptic(vector1);

            AstroVector vector2 = GeoVector(body2, time, Aberration.None);
            Ecliptic eclip2 = EquatorialToEcliptic(vector2);

            return NormalizeLongitude(eclip1.elon - eclip2.elon);
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
            return PairLongitude(Body.Moon, Body.Sun, time);
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
        ///
        /// See #Astronomy.MoonQuartersAfter for a convenient enumerator.
        /// </remarks>
        /// <param name="startTime">The date and time at which to start the search.</param>
        /// <returns>
        /// A #MoonQuarterInfo structure reporting the next quarter phase and the time it will occur.
        /// </returns>
        public static MoonQuarterInfo SearchMoonQuarter(AstroTime startTime)
        {
            double currentPhaseAngle = MoonPhase(startTime);
            int quarter = (1 + (int)Math.Floor(currentPhaseAngle / 90.0)) % 4;
            AstroTime qtime = SearchMoonPhase(90.0 * quarter, startTime, 10.0) ??
                throw new InternalError($"Unable to find moon quarter {quarter} for startTime={startTime}.");
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
        ///
        /// See #Astronomy.MoonQuartersAfter for a convenient enumerator.
        /// </remarks>
        /// <param name="mq">The previous moon quarter found by a call to #Astronomy.SearchMoonQuarter or `Astronomy.NextMoonQuarter`.</param>
        /// <returns>The moon quarter that occurs next in time after the one passed in `mq`.</returns>
        public static MoonQuarterInfo NextMoonQuarter(MoonQuarterInfo mq)
        {
            // Skip 6 days past the previous found moon quarter to find the next one.
            // This is less than the minimum possible increment.
            // So far I have seen the interval well contained by the range (6.5, 8.3) days.

            AstroTime time = mq.time.AddDays(6.0);
            MoonQuarterInfo next_mq = SearchMoonQuarter(time);
            // Verify that we found the expected moon quarter.
            if (next_mq.quarter != (1 + mq.quarter) % 4)
                throw new InternalError("found the wrong moon quarter.");
            return next_mq;
        }


        /// <summary>Enumerates a series of lunar quarter phases that occur after a specified time.</summary>
        /// <remarks>
        /// This is a convenience wrapper around the functions
        /// #Astronomy.SearchMoonQuarter and #Astronomy.NextMoonQuarter.
        /// </remarks>
        /// <param name="startTime">
        /// Specifies the time to begin searching for consecutive lunar quarter phases.
        /// </param>
        public static IEnumerable<MoonQuarterInfo> MoonQuartersAfter(AstroTime startTime)
        {
            MoonQuarterInfo mq = SearchMoonQuarter(startTime);
            yield return mq;
            while (true)
            {
                mq = NextMoonQuarter(mq);
                yield return mq;
            }
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
        /// The number of days away from `startTime` that limits the time window for the search.
        /// If the value is negative, the search is performed into the past from `startTime`.
        /// Otherwise, the search is performed into the future from `startTime`.
        /// </param>
        /// <returns>
        /// If successful, returns the date and time the moon reaches the phase specified by
        /// `targetlon`. This function will return `null` if the phase does not
        /// occur within `limitDays` of `startTime`; that is, if the search window is too small.
        /// </returns>
        public static AstroTime SearchMoonPhase(double targetLon, AstroTime startTime, double limitDays)
        {
            // To avoid discontinuities in the moon_offset function causing problems,
            // we need to approximate when that function will next return 0.
            // We probe it with the start time and take advantage of the fact
            // that every lunar phase repeats roughly every 29.5 days.
            // There is a surprising uncertainty in the quarter timing,
            // due to the eccentricity of the moon's orbit.
            // I have seen more than 0.9 days away from the simple prediction.
            // To be safe, we take the predicted time of the event and search
            // +/-1.5 days around it (a 3-day wide window).

            const double uncertainty = 1.5;
            var moon_offset = new SearchContext_MoonOffset(targetLon);

            double ya = moon_offset.Eval(startTime);

            double est_dt, dt1, dt2;
            if (limitDays < 0.0)
            {
                // Search backward in time.
                if (ya < 0.0) ya += 360.0;
                est_dt = -(MEAN_SYNODIC_MONTH * ya) / 360.0;
                dt1 = est_dt - uncertainty;
                dt2 = est_dt + uncertainty;
                if (dt2 < limitDays)
                    return null;    // not possible for moon phase to occur within specified window.
                if (dt1 < limitDays)
                    dt1 = limitDays;
            }
            else
            {
                // Search forward in time.
                if (ya > 0.0) ya -= 360.0;
                est_dt = -(MEAN_SYNODIC_MONTH * ya) / 360.0;
                dt1 = est_dt - uncertainty;
                dt2 = est_dt + uncertainty;
                if (dt1 > limitDays)
                    return null;    // not possible for moon phase to occur within specified window.
                if (dt2 > limitDays)
                    dt2 = limitDays;
            }
            AstroTime t1 = startTime.AddDays(dt1);
            AstroTime t2 = startTime.AddDays(dt2);
            return Search(moon_offset, t1, t2, 0.1);
        }


        /// <summary>
        /// Calculates U.S. Standard Atmosphere (1976) variables as a function of elevation.
        /// </summary>
        /// <remarks>
        /// This function calculates idealized values of pressure, temperature, and density
        /// using the U.S. Standard Atmosphere (1976) model.
        /// 1. COESA, U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, DC, 1976.
        /// 2. Jursa, A. S., Ed., Handbook of Geophysics and the Space Environment, Air Force Geophysics Laboratory, 1985.
        /// See:
        /// https://hbcp.chemnetbase.com/faces/documents/14_12/14_12_0001.xhtml
        /// https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf
        /// https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf
        /// </remarks>
        /// <param name="elevationMeters">
        /// The elevation above sea level at which to calculate atmospheric variables.
        /// Must be in the range -500 to +100000, or an exception will occur.
        /// </param>
        public static AtmosphereInfo Atmosphere(double elevationMeters)
        {
            const double P0 = 101325.0;     // pressure at sea level [pascals]
            const double T0 = 288.15;       // temperature at sea level [kelvins]
            const double T1 = 216.65;       // temperature between 20 km and 32 km [kelvins]

            if (!isfinite(elevationMeters) || elevationMeters < -500.0 || elevationMeters > 100000.0)
                throw new ArgumentOutOfRangeException(nameof(elevationMeters));

            double temperature;
            double pressure;
            if (elevationMeters <= 11000.0)
            {
                temperature = T0 - 0.0065*elevationMeters;
                pressure = P0 * Math.Pow(T0 / temperature, -5.25577);
            }
            else if (elevationMeters <= 20000.0)
            {
                temperature = T1;
                pressure = 22632.0 * Math.Exp(-0.00015768832 * (elevationMeters - 11000.0));
            }
            else
            {
                temperature = T1 + 0.001*(elevationMeters - 20000.0);
                pressure = 5474.87 * Math.Pow(T1 / temperature, 34.16319);
            }
            // The density is calculated relative to the sea level value.
            // Using the ideal gas law PV=nRT, we deduce that density is proportional to P/T.
            double density = (pressure / temperature) / (P0 / T0);
            return new AtmosphereInfo { pressure = pressure, temperature = temperature, density = density };
        }


        private static double HorizonDipAngle(Observer observer, double metersAboveGround)
        {
            // Calculate the effective radius of the Earth at ground level below the observer.
            // Correct for the Earth's oblateness.
            double phi = observer.latitude * DEG2RAD;
            double sinphi = Math.Sin(phi);
            double cosphi = Math.Cos(phi);
            double c = 1.0 / hypot(cosphi, sinphi*EARTH_FLATTENING);
            double s = c * (EARTH_FLATTENING * EARTH_FLATTENING);
            double ht_km = (observer.height - metersAboveGround) / 1000.0;     // height of ground above sea level
            double ach = EARTH_EQUATORIAL_RADIUS_KM*c + ht_km;
            double ash = EARTH_EQUATORIAL_RADIUS_KM*s + ht_km;
            double radius_m = 1000.0 * hypot(ach*cosphi, ash*sinphi);

            // Correct refraction of a ray of light traveling tangent to the Earth's surface.
            // Based on: https://www.largeformatphotography.info/sunmooncalc/SMCalc.js
            // which in turn derives from:
            // Sweer, John. 1938.  The Path of a Ray of Light Tangent to the Surface of the Earth.
            // Journal of the Optical Society of America 28 (September):327-329.

            // k = refraction index
            double k = 0.175 * Math.Pow(1.0 - (6.5e-3/283.15)*(observer.height - (2.0/3.0)*metersAboveGround), 3.256);

            // Calculate how far below the observer's horizontal plane the observed horizon dips.
            return RAD2DEG * -(Math.Sqrt(2*(1 - k)*metersAboveGround / radius_m) / (1 - k));
        }


        private struct AscentInfo
        {
            public bool valid;
            public AstroTime tx;
            public AstroTime ty;
            public double ax;
            public double ay;

            public static readonly AscentInfo Fail = new AscentInfo { valid = false };

            public override string ToString() =>
                valid ? $"AscentInfo(tx={tx}, ty={ty}, ax={ax}, ay={ay})" : "AscentInfo.Fail";
        }


        private static AscentInfo FindAscent(
            int depth,
            SearchContext_Altitude context,
            double max_deriv_alt,
            AstroTime t1,
            AstroTime t2,
            double a1,
            double a2)
        {
            // See if we can find any time interval where the altitude-diff function
            // rises from non-positive to positive.

            if (a1 < 0.0 && a2 >= 0.0)
            {
                // Trivial success case: the endpoints already rise through zero.
                return new AscentInfo { valid = true, tx = t1, ty = t2, ax = a1, ay = a2 };
            }

            if (a1 >= 0.0 && a2 < 0.0)
            {
                // Trivial failure case: Assume Nyquist condition prevents an ascent.
                return AscentInfo.Fail;
            }

            if (depth > 17)
            {
                // Safety valve: do not allow unlimited recursion.
                // This should never happen if the rest of the logic is working correctly,
                // so fail the whole search if it does happen. It's a bug!
                throw new InternalError("Excessive recursion in rise/set ascent search.");
            }

            // Both altitudes are on the same side of zero: both are negative, or both are non-negative.
            // There could be a convex "hill" or a concave "valley" that passes through zero.
            // In polar regions sometimes there is a rise/set or set/rise pair within minutes of each other.
            // For example, the Moon can be below the horizon, then the very top of it becomes
            // visible (moonrise) for a few minutes, then it moves sideways and down below
            // the horizon again (moonset). We want to catch these cases.
            // However, for efficiency and practicality concerns, because the rise/set search itself
            // has a 0.1 second threshold, we do not worry about rise/set pairs that are less than
            // one second apart. These are marginal cases that are rendered highly uncertain
            // anyway, due to unpredictable atmospheric refraction conditions (air temperature and pressure).

            double dt = t2.ut - t1.ut;
            if (dt * SECONDS_PER_DAY < 1.0)
                return new AscentInfo { valid = false };

            // Is it possible to reach zero from the altitude that is closer to zero?
            double da = Math.Min(Math.Abs(a1), Math.Abs(a2));

            // Without loss of generality, assume |a1| <= |a2|.
            // (Reverse the argument in the case |a2| < |a1|.)
            // Imagine you have to "drive" from a1 to 0, then back to a2.
            // You can't go faster than max_deriv_alt. If you can't reach 0 in half the time,
            // you certainly don't have time to reach 0, turn around, and still make your way
            // back up to a2 (which is at least as far from 0 than a1 is) in the time interval dt.
            // Therefore, the time threshold is half the time interval, or dt/2.
            if (da > max_deriv_alt*(dt / 2))
            {
                // Prune: the altitude cannot change fast enough to reach zero.
                return AscentInfo.Fail;
            }

            // Bisect the time interval and evaluate the altitude at the midpoint.
            var tmid = new AstroTime((t1.ut + t2.ut)/2);
            double amid = context.Eval(tmid);

            // Recurse to the left interval.
            AscentInfo ascent = FindAscent(1+depth, context, max_deriv_alt, t1, tmid, a1, amid);
            if (!ascent.valid)
            {
                // Recurse to the right interval.
                ascent = FindAscent(1+depth, context, max_deriv_alt, tmid, t2, amid, a2);
            }

            return ascent;
        }


        private static double MaxAltitudeSlope(Body body, double latitude)
        {
            // Calculate the maximum possible rate that this body's altitude
            // could change [degrees/day] as seen by this observer.
            // First use experimentally determined extreme bounds for this body
            // of how much topocentric RA and DEC can ever change per rate of time.
            // We need minimum possible d(RA)/dt, and maximum possible magnitude of d(DEC)/dt.
            // Conservatively, we round d(RA)/dt down, d(DEC)/dt up.
            // Then calculate the resulting maximum possible altitude change rate.

            double deriv_ra;
            double deriv_dec;

            switch (body)
            {
            case Body.Moon:
                deriv_ra  = +4.5;
                deriv_dec = +8.2;
                break;

            case Body.Sun:
                deriv_ra  = +0.8;
                deriv_dec = +0.5;
                break;

            case Body.Mercury:
                deriv_ra  = -1.6;
                deriv_dec = +1.0;
                break;

            case Body.Venus:
                deriv_ra  = -0.8;
                deriv_dec = +0.6;
                break;

            case Body.Mars:
                deriv_ra  = -0.5;
                deriv_dec = +0.4;
                break;

            case Body.Jupiter:
            case Body.Saturn:
            case Body.Uranus:
            case Body.Neptune:
            case Body.Pluto:
                deriv_ra  = -0.2;
                deriv_dec = +0.2;
                break;

            case Body.Star1:
            case Body.Star2:
            case Body.Star3:
            case Body.Star4:
            case Body.Star5:
            case Body.Star6:
            case Body.Star7:
            case Body.Star8:
                // The minimum allowed heliocentric distance of a user-defined star
                // is one light-year. This can cause a tiny amount of parallax (about 0.001 degrees).
                // Also, including stellar aberration (22 arcsec = 0.006 degrees), we provide a
                // generous safety buffer of 0.008 degrees.
                deriv_ra  = -0.008;
                deriv_dec = +0.008;
                break;

            case Body.Earth:
                throw new EarthNotAllowedException();

            default:
                throw new InvalidBodyException(body);
            }

            double latrad = DEG2RAD * latitude;
            return Math.Abs(((360.0 / SOLAR_DAYS_PER_SIDEREAL_DAY) - deriv_ra)*Math.Cos(latrad)) + Math.Abs(deriv_dec*Math.Sin(latrad));
        }

        private const double RISE_SET_DT = 0.42;    // 10.08 hours: Nyquist-safe for 22-hour period.

        private static AstroTime InternalSearchAltitude(
            Body body,
            Observer observer,
            Direction direction,
            AstroTime startTime,
            double limitDays,
            double bodyRadiusAu,
            double targetAltitude)
        {
            double max_deriv_alt = MaxAltitudeSlope(body, observer.latitude);
            var context = new SearchContext_Altitude(body, direction, observer, bodyRadiusAu, targetAltitude);

            // We allow searching forward or backward in time.
            // But we want to keep t1 < t2, so we need a few if/else statements.
            AstroTime t1 = startTime;
            AstroTime t2 = t1;
            double a1 = context.Eval(t1);
            double a2 = a1;

            for(;;)
            {
                if (limitDays < 0.0)
                {
                    t1 = t2.AddDays(-RISE_SET_DT);
                    a1 = context.Eval(t1);
                }
                else
                {
                    t2 = t1.AddDays(+RISE_SET_DT);
                    a2 = context.Eval(t2);
                }

                AscentInfo ascent = FindAscent(0, context, max_deriv_alt, t1, t2, a1, a2);
                if (ascent.valid)
                {
                    // We found a time interval [t1, t2] that contains an alt-diff
                    // rising from negative a1 to non-negative a2.
                    // Search for the time where the root occurs.
                    AstroTime time = Search(context, ascent.tx, ascent.ty, 0.1);
                    if (time != null)
                    {
                        // Now that we have a solution, we have to check whether it goes outside the time bounds.
                        if (limitDays < 0.0)
                        {
                            if (time.ut < startTime.ut + limitDays)
                                return null;
                        }
                        else
                        {
                            if (time.ut > startTime.ut + limitDays)
                                return null;
                        }
                        return time;    // success!
                    }

                    // The search should have succeeded. Something is wrong with the ascent finder!
                    throw new InternalError($"Rise/set search failed after finding {ascent}");
                }

                // There is no ascent in this interval, so keep searching.
                if (limitDays < 0.0)
                {
                    if (t1.ut < startTime.ut + limitDays)
                        return null;
                    t2 = t1;
                    a2 = a1;
                }
                else
                {
                    if (t2.ut > startTime.ut + limitDays)
                        return null;
                    t1 = t2;
                    a1 = a2;
                }
            }
        }


        /// <summary>
        /// Searches for the next time a celestial body rises or sets as seen by an observer on the Earth.
        /// </summary>
        ///
        /// <remarks>
        /// This function finds the next rise or set time of the Sun, Moon, or planet other than the Earth.
        /// Rise time is when the body first starts to be visible above the horizon.
        /// For example, sunrise is the moment that the top of the Sun first appears to peek above the horizon.
        /// Set time is the moment when the body appears to vanish below the horizon.
        /// Therefore, this function adjusts for the apparent angular radius of the observed body
        /// (significant only for the Sun and Moon).
        ///
        /// This function corrects for a typical value of atmospheric refraction, which causes celestial
        /// bodies to appear higher above the horizon than they would if the Earth had no atmosphere.
        /// Astronomy Engine uses a correction of 34 arcminutes. Real-world refraction varies based
        /// on air temperature, pressure, and humidity; such weather-based conditions are outside
        /// the scope of Astronomy Engine.
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
        /// <param name="body">
        /// The Sun, Moon, any planet other than the Earth,
        /// or a user-defined star that was created by a call to #Astronomy.DefineStar.
        /// </param>
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
        /// Limits how many days to search for a rise or set time, and defines
        /// the direction in time to search. When `limitDays` is positive, the
        /// search is performed into the future, after `startTime`.
        /// When negative, the search is performed into the past, before `startTime`.
        /// To limit a rise or set time to the same day, you can use a value of 1 day.
        /// In cases where you want to find the next rise or set time no matter how far
        /// in the future (for example, for an observer near the south pole), you can
        /// pass in a larger value like 365.
        /// </param>
        ///
        /// <param name="metersAboveGround">
        /// Usually the observer is located at ground level. Then this parameter
        /// should be zero. But if the observer is significantly higher than ground
        /// level, for example in an airplane, this parameter should be a positive
        /// number indicating how far above the ground the observer is.
        /// An exception occurs if `metersAboveGround` is negative.
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
            double limitDays,
            double metersAboveGround = 0.0)
        {
            if (!isfinite(metersAboveGround) || metersAboveGround < 0.0)
                throw new ArgumentOutOfRangeException(nameof(metersAboveGround));

            double bodyRadiusAu;
            switch (body)
            {
                case Body.Sun:
                    bodyRadiusAu = SUN_RADIUS_AU;
                    break;

                case Body.Moon:
                    bodyRadiusAu = MOON_EQUATORIAL_RADIUS_AU;
                    break;

                default:
                    bodyRadiusAu = 0.0;
                    break;
            }

            // Calculate atmospheric density at ground level.
            AtmosphereInfo atmos = Astronomy.Atmosphere(observer.height - metersAboveGround);

            // Calculate the apparent angular dip of the horizon.
            double dip = HorizonDipAngle(observer, metersAboveGround);

            // Correct refraction for objects near the horizon, using atmospheric density at the ground.
            double altitude = dip - (REFRACTION_NEAR_HORIZON * atmos.density);

            // Search for the top of the body crossing the corrected altitude angle.
            return InternalSearchAltitude(body, observer, direction, startTime, limitDays, bodyRadiusAu, altitude);
        }

        /// <summary>
        /// Finds the next time the center of a body passes through a given altitude.
        /// </summary>
        /// <remarks>
        /// Finds when the center of the given body ascends or descends through a given
        /// altitude angle, as seen by an observer at the specified location on the Earth.
        /// By using the appropriate combination of `direction` and `altitude` parameters,
        /// this function can be used to find when civil, nautical, or astronomical twilight
        /// begins (dawn) or ends (dusk).
        ///
        /// Civil dawn begins before sunrise when the Sun ascends through 6 degrees below
        /// the horizon. To find civil dawn, pass `Direction.Rise` for `direction` and -6 for `altitude`.
        ///
        /// Civil dusk ends after sunset when the Sun descends through 6 degrees below the horizon.
        /// To find civil dusk, pass `Direction.Set` for `direction` and -6 for `altitude`.
        ///
        /// Nautical twilight is similar to civil twilight, only the `altitude` value should be -12 degrees.
        ///
        /// Astronomical twilight uses -18 degrees as the `altitude` value.
        ///
        /// By convention for twilight time calculations, the altitude is not corrected for
        /// atmospheric refraction. This is because the target altitudes are below the horizon,
        /// and refraction is not directly observable.
        ///
        /// `SearchAltitude` is not intended to find rise/set times of a body for two reasons:
        /// (1) Rise/set times of the Sun or Moon are defined by their topmost visible portion, not their centers.
        /// (2) Rise/set times are affected significantly by atmospheric refraction.
        /// Therefore, it is better to use #Astronomy.SearchRiseSet to find rise/set times, which
        /// corrects for both of these considerations.
        ///
        /// `SearchAltitude` will not work reliably for altitudes at or near the body's
        /// maximum or minimum altitudes. To find the time a body reaches minimum or maximum altitude
        /// angles, use #Astronomy.SearchHourAngle.
        /// </remarks>
        ///
        /// <param name="body">
        /// The Sun, Moon, any planet other than the Earth,
        /// or a user-defined star that was created by a call to #Astronomy.DefineStar.
        /// </param>
        ///
        /// <param name="observer">The location where observation takes place.</param>
        ///
        /// <param name="direction">
        /// Either `Direction.Rise` to find an ascending altitude event
        /// or `Direction.Set` to find a descending altitude event.
        /// </param>
        ///
        /// <param name="startTime">The date and time at which to start the search.</param>
        ///
        /// <param name="limitDays">
        /// Limits how many days to search for the body reaching the altitude angle,
        /// and defines the direction in time to search. When `limitDays` is positive, the
        /// search is performed into the future, after `startTime`.
        /// When negative, the search is performed into the past, before `startTime`.
        /// To limit the search to the same day, you can use a value of 1 day.
        /// In cases where you want to find the altitude event no matter how far
        /// in the future (for example, for an observer near the south pole), you can
        /// pass in a larger value like 365.
        /// </param>
        ///
        /// <param name="altitude">
        /// The desired altitude angle of the body's center above (positive)
        /// or below (negative) the observer's local horizon, expressed in degrees.
        /// Must be in the range [-90, +90].
        /// </param>
        ///
        /// <returns>
        /// The date and time of the altitude event, or `null` if no such event
        /// occurs within the specified time window.
        /// </returns>
        public static AstroTime SearchAltitude(
            Body body,
            Observer observer,
            Direction direction,
            AstroTime startTime,
            double limitDays,
            double altitude)
        {
            if (altitude < -90.0 || altitude > +90.0)
                throw new ArgumentOutOfRangeException(nameof(altitude));
            return InternalSearchAltitude(body, observer, direction, startTime, limitDays, 0.0, altitude);
        }

        /// <summary>
        /// Finds the hour angle of a body for a given observer and time.
        /// </summary>
        /// <remarks>
        /// The *hour angle* of a celestial body indicates its position in the sky with respect
        /// to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
        /// The hour angle is 0 when the body's center reaches its highest angle above the horizon in a given day.
        /// The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
        /// to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
        /// the number of hours that have passed since the most recent time that the body has culminated,
        /// or reached its highest point.
        /// </remarks>
        /// <param name="body">The body whose observed hour angle is to be found.</param>
        /// <param name="time">The time of the observation.</param>
        /// <param name="observer">The geographic location where the observation takes place.</param>
        /// <returns>The real-valued hour angle of the body in the half-open range [0, 24).</returns>
        public static double HourAngle(
            Body body,
            AstroTime time,
            Observer observer)
        {
            double gast = SiderealTime(time);
            Equatorial ofdate = Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
            double hourAngle = (observer.longitude/15 + gast - ofdate.ra) % 24.0;
            if (hourAngle < 0.0)
                hourAngle += 24.0;
            return hourAngle;
        }

        /// <summary>
        /// Searches for the time when the center of a body reaches a specified hour angle as seen by an observer on the Earth.
        /// </summary>
        ///
        /// <remarks>
        /// The *hour angle* of a celestial body indicates its position in the sky with respect
        /// to the Earth's rotation. The hour angle depends on the location of the observer on the Earth.
        /// The hour angle is 0 when the body's center reaches its highest angle above the horizon in a given day.
        /// The hour angle increases by 1 unit for every sidereal hour that passes after that point, up
        /// to 24 sidereal hours when it reaches the highest point again. So the hour angle indicates
        /// the number of hours that have passed since the most recent time that the body has culminated,
        /// or reached its highest point.
        ///
        /// This function searches for the next or previous time a celestial body reaches the given hour angle
        /// relative to the date and time specified by `startTime`.
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
        /// The Sun, Moon, any planet other than the Earth,
        /// or a user-defined star that was created by a call to #Astronomy.DefineStar.
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
        /// <param name="direction">
        /// The direction in time to perform the search: a positive value
        /// searches forward in time, a negative value searches backward in time.
        /// The function throws an exception if `direction` is zero.
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
            AstroTime startTime,
            int direction = +1)
        {
            int iter = 0;

            if (body == Body.Earth)
                throw new EarthNotAllowedException();

            if (hourAngle < 0.0 || hourAngle >= 24.0)
                throw new ArgumentException("hourAngle is out of the allowed range [0, 24).");

            if (direction == 0)
                throw new ArgumentException("direction must be positive or negative.");

            AstroTime time = startTime;
            for(;;)
            {
                ++iter;

                // Calculate Greenwich Apparent Sidereal Time (GAST) at the given time.
                double gast = SiderealTime(time);

                // Obtain equatorial coordinates of date for the body.
                Equatorial ofdate = Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);

                // Calculate the adjustment needed in sidereal time
                // to bring the hour angle to the desired value.

                double delta_sidereal_hours = ((hourAngle + ofdate.ra - observer.longitude/15.0) - gast) % 24.0;
                if (iter == 1)
                {
                    // On the first iteration, always search in the requested time direction.
                    if (direction > 0)
                    {
                        // Search forward in time.
                        if (delta_sidereal_hours < 0.0)
                            delta_sidereal_hours += 24.0;
                    }
                    else
                    {
                        // Search backward in time.
                        if (delta_sidereal_hours > 0.0)
                            delta_sidereal_hours -= 24.0;
                    }
                }
                else
                {
                    // On subsequent iterations, we make the smallest possible adjustment,
                    // either forward or backward in time.
                    if (delta_sidereal_hours < -12.0)
                        delta_sidereal_hours += 24.0;
                    else if (delta_sidereal_hours > +12.0)
                        delta_sidereal_hours -= 24.0;
                }

                // If the error is tolerable (less than 0.1 seconds), the search has succeeded.
                if (Math.Abs(delta_sidereal_hours) * 3600.0 < 0.1)
                {
                    Topocentric hor = Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.Normal);
                    return new HourAngleInfo(time, hor);
                }

                // We need to loop another time to get more accuracy.
                // Update the terrestrial time (in solar days) adjusting by sidereal time (sidereal hours).
                time = time.AddDays((delta_sidereal_hours / 24.0) * SOLAR_DAYS_PER_SIDEREAL_DAY);
            }
        }

        /// <summary>
        /// Searches for the time when the Earth and another planet are separated by a specified angle
        /// in ecliptic longitude, as seen from the Sun.
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
        /// <returns>The date and time of the relative longitude event.</returns>
        public static AstroTime SearchRelativeLongitude(Body body, double targetRelLon, AstroTime startTime)
        {
            if (body == Body.Earth || body == Body.Sun || body == Body.Moon)
                throw new InvalidBodyException(body);

            double syn = SynodicPeriod(body);
            int direction = IsSuperiorPlanet(body) ? +1 : -1;

            // Iterate until we converge on the desired event.
            // Calculate the error angle, which will be a negative number of degrees,
            // meaning we are "behind" the target relative longitude.

            double error_angle = RelativeLongitudeOffset(body, startTime, direction, targetRelLon);
            if (error_angle > 0.0)
                error_angle -= 360.0;    // force searching forward in time

            AstroTime time = startTime;
            for (int iter = 0; iter < 100; ++iter)
            {
                // Estimate how many days in the future (positive) or past (negative)
                // we have to go to get closer to the target relative longitude.
                double day_adjust = (-error_angle/360.0) * syn;
                time = time.AddDays(day_adjust);
                if (Math.Abs(day_adjust) * SECONDS_PER_DAY < 1.0)
                    return time;

                double prev_angle = error_angle;
                error_angle = RelativeLongitudeOffset(body, time, direction, targetRelLon);
                if (Math.Abs(prev_angle) < 30.0 && (prev_angle != error_angle))
                {
                    // Improve convergence for Mercury/Mars (eccentric orbits)
                    // by adjusting the synodic period to more closely match the
                    // variable speed of both planets in this part of their respective orbits.
                    double ratio = prev_angle / (prev_angle - error_angle);
                    if (ratio > 0.5 && ratio < 2.0)
                        syn *= ratio;
                }
            }

            throw new InternalError("Relative longitude search failed to converge.");
        }

        private static double RelativeLongitudeOffset(Body body, AstroTime time, int direction, double targetRelLon)
        {
            double plon = EclipticLongitude(body, time);
            double elon = EclipticLongitude(Body.Earth, time);
            double diff = direction * (elon - plon);
            return LongitudeOffset(diff - targetRelLon);
        }

        private static double SynodicPeriod(Body body)
        {
            // The Earth does not have a synodic period as seen from itself.
            if (body == Body.Earth)
                throw new EarthNotAllowedException();

            if (body == Body.Moon)
                return MEAN_SYNODIC_MONTH;

            double Tp = PlanetOrbitalPeriod(body);
            return Math.Abs(EARTH_ORBITAL_PERIOD / (EARTH_ORBITAL_PERIOD/Tp - 1.0));
        }

        /// <summary>Calculates heliocentric ecliptic longitude of a body.</summary>
        /// <remarks>
        /// This function calculates the angle around the plane of the Earth's orbit
        /// of a celestial body, as seen from the center of the Sun.
        /// The angle is measured prograde (in the direction of the Earth's orbit around the Sun)
        /// in degrees from the true equinox of date. The ecliptic longitude is always in the range [0, 360).
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

        /// <summary>Returns the average number of days it takes for a planet to orbit the Sun.</summary>
        /// <param name="body">One of the planets: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, or Pluto.</param>
        /// <returns>The mean orbital period of the body in days.</returns>
        public static double PlanetOrbitalPeriod(Body body)
        {
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
            double ecliptic_separation = PairLongitude(body, Body.Sun, time);
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
                double rlon = LongitudeOffset(plon - elon);     // clamp to (-180, +180]

                // The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees
                // because there is a cusp there that causes a discontinuity in the derivative.
                // So we need to guard against searching near such times.
                double adjust_days, rlon_lo, rlon_hi;
                if (rlon >= -s1 && rlon < +s1)
                {
                    // Seek to the window [+s1, +s2].
                    adjust_days = 0.0;
                    // Search forward for the time t1 when rel lon = +s1.
                    rlon_lo = +s1;
                    // Search forward for the time t2 when rel lon = +s2.
                    rlon_hi = +s2;
                }
                else if (rlon > +s2 || rlon < -s2)
                {
                    // Seek to the next search window at [-s2, -s1].
                    adjust_days = 0.0;
                    // Search forward for the time t1 when rel lon = -s2.
                    rlon_lo = -s2;
                    // Search forward for the time t2 when rel lon = -s1.
                    rlon_hi = -s1;
                }
                else if (rlon >= 0.0)
                {
                    // rlon must be in the middle of the window [+s1, +s2].
                    // Search BACKWARD for the time t1 when rel lon = +s1.
                    adjust_days = -syn / 4.0;
                    rlon_lo = +s1;
                    rlon_hi = +s2;
                    // Search forward from t1 to find t2 such that rel lon = +s2.
                }
                else
                {
                    // rlon must be in the middle of the window [-s2, -s1].
                    // Search BACKWARD for the time t1 when rel lon = -s2.
                    adjust_days = -syn / 4.0;
                    rlon_lo = -s2;
                    // Search forward from t1 to find t2 such that rel lon = -s1.
                    rlon_hi = -s1;
                }

                AstroTime t_start = startTime.AddDays(adjust_days);

                AstroTime t1 = SearchRelativeLongitude(body, rlon_lo, t_start);
                AstroTime t2 = SearchRelativeLongitude(body, rlon_hi, t1);

                // Now we have a time range [t1,t2] that brackets a maximum elongation event.
                // Confirm the bracketing.
                double m1 = neg_elong_slope.Eval(t1);
                if (m1 >= 0.0)
                    throw new InternalError("There is a bug in the bracketing algorithm! m1 = " + m1);

                double m2 = neg_elong_slope.Eval(t2);
                if (m2 <= 0.0)
                    throw new InternalError("There is a bug in the bracketing algorithm! m2 = " + m2);

                // Use the generic search algorithm to home in on where the slope crosses from negative to positive.
                AstroTime searchx = Search(neg_elong_slope, t1, t2, 10.0) ??
                    throw new InternalError("Maximum elongation search failed.");

                if (searchx.tt >= startTime.tt)
                    return Elongation(body, searchx);

                // This event is in the past (earlier than startTime).
                // We need to search forward from t2 to find the next possible window.
                // We never need to search more than twice.
                startTime = t2.AddDays(1.0);
            }

            throw new InternalError("Maximum elongation search iterated too many times.");
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

        /// <summary>
        /// Calculates the angle in degrees between two vectors.
        /// </summary>
        /// <remarks>
        /// Given a pair of vectors, this function returns the angle in degrees
        /// between the two vectors in 3D space.
        /// The angle is measured in the plane that contains both vectors.
        /// </remarks>
        /// <param name="a">The first of a pair of vectors between which to measure an angle.</param>
        /// <param name="b">The second of a pair of vectors between which to measure an angle.</param>
        /// <returns>
        /// The angle between the two vectors expressed in degrees.
        /// The value is in the range [0, 180].
        /// </returns>
        public static double AngleBetween(AstroVector a, AstroVector b)
        {
            double r = a.Length() * b.Length();
            if (r < 1.0e-8)
                throw new ArgumentException("Cannot find angle between vectors because they are too short.");

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
        ///
        /// See #Astronomy.LunarApsidesAfter for a convenient enumerator.
        /// </remarks>
        /// <param name="startTime">
        ///      The date and time at which to start searching for the next perigee or apogee.
        /// </param>
        /// <returns>
        /// Returns an #ApsisInfo structure containing information about the next lunar apsis.
        /// </returns>
        public static ApsisInfo SearchLunarApsis(AstroTime startTime)
        {
            const double increment = 5.0;   // number of days to skip in each iteration
            var positive_slope = new SearchContext_MoonDistanceSlope(+1);
            var negative_slope = new SearchContext_MoonDistanceSlope(-1);

            // Check the rate of change of the distance dr/dt at the start time.
            // If it is positive, the Moon is currently getting farther away,
            // so start looking for apogee.
            // Conversely, if dr/dt < 0, start looking for perigee.
            // Either way, the polarity of the slope will change, so the product will be negative.
            // Handle the crazy corner case of exactly touching zero by checking for m1*m2 <= 0.
            AstroTime t1 = startTime;
            double m1 = positive_slope.Eval(t1);
            for (int iter=0; iter * increment < 2.0 * Astronomy.MEAN_SYNODIC_MONTH; ++iter)
            {
                AstroTime t2 = t1.AddDays(increment);
                double m2 = positive_slope.Eval(t2);
                if (m1 * m2 <= 0.0)
                {
                    // There is a change of slope polarity within the time range [t1, t2].
                    // Therefore this time range contains an apsis.
                    // Figure out whether it is perigee or apogee.

                    AstroTime search;
                    ApsisKind kind;
                    if (m1 < 0.0 || m2 > 0.0)
                    {
                        // We found a minimum-distance event: perigee.
                        // Search the time range for the time when the slope goes from negative to positive.
                        search = Search(positive_slope, t1, t2, 1.0);
                        kind = ApsisKind.Pericenter;
                    }
                    else if (m1 > 0.0 || m2 < 0.0)
                    {
                        // We found a maximum-distance event: apogee.
                        // Search the time range for the time when the slope goes from positive to negative.
                        search = Search(negative_slope, t1, t2, 1.0);
                        kind = ApsisKind.Apocenter;
                    }
                    else
                    {
                        // This should never happen. It should not be possible for both slopes to be zero.
                        throw new InternalError("both slopes are zero in SearchLunarApsis.");
                    }

                    if (search == null)
                        throw new InternalError("Failed to find slope transition in lunar apsis search.");

                    double dist_au = SearchContext_MoonDistanceSlope.MoonDistance(search);
                    return new ApsisInfo(search, kind, dist_au);
                }
                // We have not yet found a slope polarity change. Keep searching.
                t1 = t2;
                m1 = m2;
            }

            // It should not be possible to fail to find an apsis within 2 synodic months.
            throw new InternalError("should have found lunar apsis within 2 synodic months.");
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
        /// See #Astronomy.LunarApsidesAfter for a convenient enumerator.
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
                throw new InternalError($"Internal error: previous apsis was {apsis.kind}, but found {next.kind} for next apsis.");
            return next;
        }


        /// <summary>Enumerates a series of apogees/perigees that occur after a specified time.</summary>
        /// <remarks>
        /// This is a convenience wrapper around the functions
        /// #Astronomy.SearchLunarApsis and #Astronomy.NextLunarApsis.
        /// </remarks>
        /// <param name="startTime">
        /// Specifies the time to begin searching for consecutive lunar apsides.
        /// </param>
        public static IEnumerable<ApsisInfo> LunarApsidesAfter(AstroTime startTime)
        {
            ApsisInfo apsis = SearchLunarApsis(startTime);
            yield return apsis;
            while (true)
            {
                apsis = NextLunarApsis(apsis);
                yield return apsis;
            }
        }


        private static ApsisInfo PlanetExtreme(Body body, ApsisKind kind, AstroTime start_time, double dayspan)
        {
            double direction = (kind == ApsisKind.Apocenter) ? +1.0 : -1.0;
            const int npoints = 10;

            for(;;)
            {
                double interval = dayspan / (npoints - 1);

                if (interval < 1.0 / 1440.0)    // iterate until uncertainty is less than one minute
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

                // Narrow in on the extreme point.
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

            // Neptune is a special case for two reasons:
            // 1. Its orbit is nearly circular (low orbital eccentricity).
            // 2. It is so distant from the Sun that the orbital period is very long.
            // Put together, this causes wobbling of the Sun around the Solar System Barycenter (SSB)
            // to be so significant that there are 3 local minima in the distance-vs-time curve
            // near each apsis. Therefore, unlike for other planets, we can't use an optimized
            // algorithm for finding dr/dt = 0.
            // Instead, we use a dumb, brute-force algorithm of sampling and finding min/max
            // heliocentric distance.
            //
            // There is a similar problem in the TOP2013 model for Pluto:
            // Its position vector has high-frequency oscillations that confuse the
            // slope-based determination of apsides.
            //
            // Rewind approximately 30 degrees in the orbit,
            // then search forward for 270 degrees.
            // This is a very cautious way to prevent missing an apsis.
            // Typically we will find two apsides, and we pick whichever
            // apsis is ealier, but after startTime.
            // Sample points around this orbital arc and find when the distance
            // is greatest and smallest.

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
                    // Perihelion and aphelion are both valid. Pick the one that comes first.
                    if (aphelion.time.tt < perihelion.time.tt)
                        return aphelion;
                }
                return perihelion;
            }

            if (aphelion.time.tt >= startTime.tt)
                return aphelion;

            throw new InternalError("failed to find planet apsis.");
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
        ///
        /// See #Astronomy.PlanetApsidesAfter for a convenient enumerator.
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
                    // There is a change of slope polarity within the time range [t1, t2].
                    // Therefore this time range contains an apsis.
                    // Figure out whether it is perihelion or aphelion.

                    SearchContext_PlanetDistanceSlope slope_func;
                    ApsisKind kind;
                    if (m1 < 0.0 || m2 > 0.0)
                    {
                        // We found a minimum-distance event: perihelion.
                        // Search the time range for the time when the slope goes from negative to positive.
                        slope_func = positive_slope;
                        kind = ApsisKind.Pericenter;
                    }
                    else if (m1 > 0.0 || m2 < 0.0)
                    {
                        // We found a maximum-distance event: aphelion.
                        // Search the time range for the time when the slope goes from positive to negative.
                        slope_func = negative_slope;
                        kind = ApsisKind.Apocenter;
                    }
                    else
                    {
                        // This should never happen. It should not be possible for both slopes to be zero.
                        throw new InternalError("Both slopes were zero in SearchPlanetApsis");
                    }

                    AstroTime search = Search(slope_func, t1, t2, 1.0) ??
                        throw new InternalError("Failed to find slope transition in planetary apsis search.");

                    double dist = HelioDistance(body, search);
                    return new ApsisInfo(search, kind, dist);
                }
                // We have not yet found a slope polarity change. Keep searching.
                t1 = t2;
                m1 = m2;
            }
            // It should not be possible to fail to find an apsis within 2 planet orbits.
            throw new InternalError("should have found planetary apsis within 2 orbital periods.");
        }

        /// <summary>
        /// Finds the next planetary perihelion or aphelion event in a series.
        /// </summary>
        /// <remarks>
        /// This function requires an #ApsisInfo value obtained from a call
        /// to #Astronomy.SearchPlanetApsis or `Astronomy.NextPlanetApsis`.
        /// Given an aphelion event, this function finds the next perihelion event, and vice versa.
        ///
        /// See #Astronomy.SearchPlanetApsis for more details.
        /// See #Astronomy.PlanetApsidesAfter for a convenient enumerator.
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

            // skip 1/4 of an orbit before starting search again
            double skip = 0.25 * PlanetOrbitalPeriod(body);
            if (skip <= 0.0)
                throw new InvalidBodyException(body);

            AstroTime time = apsis.time.AddDays(skip);
            ApsisInfo next = SearchPlanetApsis(body, time);

            // Verify that we found the opposite apsis from the previous one.
            if ((int)next.kind + (int)apsis.kind != 1)
                throw new InternalError($"previous apsis was {apsis.kind}, but found {next.kind} for next apsis.");

            return next;
        }


        /// <summary>Enumerates a series of planet aphelia/perihelia that occur after a specified time.</summary>
        /// <remarks>
        /// This is a convenience wrapper around the functions
        /// #Astronomy.SearchPlanetApsis and #Astronomy.NextPlanetApsis.
        /// </remarks>
        /// <param name="body">
        /// The planet for which to find a series of consecutive aphelia/perihelia.
        /// Not allowed to be `Body.Sun` or `Body.Moon`.
        /// </param>
        /// <param name="startTime">
        /// Specifies the time to begin searching for consecutive planetary apsides.
        /// </param>
        public static IEnumerable<ApsisInfo> PlanetApsidesAfter(Body body, AstroTime startTime)
        {
            ApsisInfo apsis = SearchPlanetApsis(body, startTime);
            yield return apsis;
            while (true)
            {
                apsis = NextPlanetApsis(body, apsis);
                yield return apsis;
            }
        }


        // We can get away with creating a single EarthShadowSlope context
        // because it contains no state and it has no side-effects.
        // This reduces memory allocation overhead.
        private static readonly SearchContext_EarthShadowSlope earthShadowSlopeContext = new SearchContext_EarthShadowSlope();

        private static ShadowInfo PeakEarthShadow(AstroTime search_center_time)
        {
            const double window = 0.03;        // initial search window, in days, before/after given time
            AstroTime t1 = search_center_time.AddDays(-window);
            AstroTime t2 = search_center_time.AddDays(+window);
            AstroTime tx = Search(earthShadowSlopeContext, t1, t2, 1.0) ??
                throw new InternalError("Failed to find Earth peak shadow event.");
            return EarthShadow(tx);
        }


        private static double Obscuration(  // returns area of intersection of two discs, divided by area of first disc
            double a,       // radius of first disc
            double b,       // radius of second disc
            double c)       // distance between the centers of the discs
        {
            if (a <= 0.0)
                throw new ArgumentException("Radius of first disc must be positive.");

            if (b <= 0.0)
                throw new ArgumentException("Radius of second disc must be positive.");

            if (c < 0.0)
                throw new ArgumentException("Distance between discs is not allowed to be negative.");

            if (c >= a + b)
            {
                // The discs are too far apart to have any overlapping area.
                return 0.0;
            }

            if (c == 0.0)
            {
                // The discs have a common center. Therefore, one disc is inside the other.
                return (a <= b) ? 1.0 : (b*b)/(a*a);
            }

            double x = (a*a - b*b + c*c) / (2*c);
            double radicand = a*a - x*x;
            if (radicand <= 0.0)
            {
                // The circumferences do not intersect, or are tangent.
                // We already ruled out the case of non-overlapping discs.
                // Therefore, one disc is inside the other.
                return (a <= b) ? 1.0 : (b*b)/(a*a);
            }

            // The discs overlap fractionally in a pair of lens-shaped areas.

            double y = Math.Sqrt(radicand);

            // Return the overlapping fractional area.
            // There are two lens-shaped areas, one to the left of x, the other to the right of x.
            // Each part is calculated by subtracting a triangular area from a sector's area.
            double lens1 = a*a*Math.Acos(x/a) - x*y;
            double lens2 = b*b*Math.Acos((c-x)/b) - (c-x)*y;

            // Find the fractional area with respect to the first disc.
            return (lens1 + lens2) / (Math.PI*a*a);
        }


        private static double SolarEclipseObscuration(
            AstroVector hm,     // heliocentric Moon
            AstroVector lo)     // lunacentric observer
        {
            // Find heliocentric observer.
            AstroVector ho = hm + lo;

            // Calculate the apparent angular radius of the Sun for the observer.
            double sun_radius = Math.Asin(SUN_RADIUS_AU / ho.Length());

            // Calculate the apparent angular radius of the Moon for the observer.
            double moon_radius = Math.Asin(MOON_POLAR_RADIUS_AU / lo.Length());

            // Calculate the apparent angular separation between the Sun's center and the Moon's center.
            double sun_moon_separation = AngleBetween(lo, ho);

            // Find the fraction of the Sun's apparent disc area that is covered by the Moon.
            double obscuration = Obscuration(sun_radius, moon_radius, sun_moon_separation * DEG2RAD);

            // HACK: In marginal cases, we need to clamp obscuration to less than 1.0.
            // This function is never called for total eclipses, so it should never return 1.0.
            return Math.Min(0.9999, obscuration);
        }


        /// <summary>Searches for a lunar eclipse.</summary>
        /// <remarks>
        /// This function finds the first lunar eclipse that occurs after `startTime`.
        /// A lunar eclipse may be penumbral, partial, or total.
        /// See #LunarEclipseInfo for more information.
        /// To find a series of lunar eclipses, call this function once,
        /// then keep calling #Astronomy.NextLunarEclipse as many times as desired,
        /// passing in the `center` value returned from the previous call.
        ///
        /// See #Astronomy.LunarEclipsesAfter for a convenient enumerator.
        /// </remarks>
        /// <param name="startTime">
        /// The date and time for starting the search for a lunar eclipse.
        /// </param>
        /// <returns>
        /// A #LunarEclipseInfo structure containing information about the lunar eclipse.
        /// </returns>
        public static LunarEclipseInfo SearchLunarEclipse(AstroTime startTime)
        {
            const double PruneLatitude = 1.8;   // full Moon's ecliptic latitude above which eclipse is impossible
            // Iterate through consecutive full moons until we find any kind of lunar eclipse.
            AstroTime fmtime = startTime;
            for (int fmcount=0; fmcount < 12; ++fmcount)
            {
                // Search for the next full moon. Any eclipse will be near it.
                AstroTime fullmoon = SearchMoonPhase(180.0, fmtime, 40.0) ??
                    throw new InternalError("Failed to find next full moon.");

                // Pruning: if the full Moon's ecliptic latitude is too large,
                // a lunar eclipse is not possible. Avoid needless work searching for
                // the minimum moon distance.
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
                        double obscuration = 0.0;
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
                                obscuration = 1.0;
                                sd_total = ShadowSemiDurationMinutes(shadow.time, shadow.k - MOON_MEAN_RADIUS_KM, sd_partial);
                            }
                            else
                            {
                                obscuration = Obscuration(MOON_MEAN_RADIUS_KM, shadow.k, shadow.r);
                            }
                        }
                        return new LunarEclipseInfo(kind, obscuration, shadow.time, sd_penum, sd_partial, sd_total);
                    }
                }

                // We didn't find an eclipse on this full moon, so search for the next one.
                fmtime = fullmoon.AddDays(10.0);
            }

            // This should never happen, because there should be at least 2 lunar eclipses per year.
            throw new InternalError("failed to find lunar eclipse within 12 full moons.");
        }


        /// <summary>Searches for the next lunar eclipse in a series.</summary>
        /// <remarks>
        /// After using #Astronomy.SearchLunarEclipse to find the first lunar eclipse
        /// in a series, you can call this function to find the next consecutive lunar eclipse.
        /// Pass in the `center` value from the #LunarEclipseInfo returned by the
        /// previous call to `Astronomy.SearchLunarEclipse` or `Astronomy.NextLunarEclipse`
        /// to find the next lunar eclipse.
        ///
        /// See #Astronomy.LunarEclipsesAfter for a convenient enumerator.
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


        /// <summary>Enumerates a series of lunar eclipses that occur after a specified time.</summary>
        /// <remarks>
        /// This is a convenience wrapper around the functions
        /// #Astronomy.SearchLunarEclipse and #Astronomy.NextLunarEclipse.
        /// </remarks>
        /// <param name="startTime">
        /// Specifies the time to begin searching for consecutive lunar eclipses.
        /// </param>
        public static IEnumerable<LunarEclipseInfo> LunarEclipsesAfter(AstroTime startTime)
        {
            LunarEclipseInfo eclipse = SearchLunarEclipse(startTime);
            yield return eclipse;
            while (true)
            {
                eclipse = NextLunarEclipse(eclipse.peak);
                yield return eclipse;
            }
        }


        private static double ShadowSemiDurationMinutes(AstroTime center_time, double radius_limit, double window_minutes)
        {
            // Search backwards and forwards from the center time until shadow axis distance crosses radius limit.
            double window = window_minutes / (24.0 * 60.0);
            AstroTime before = center_time.AddDays(-window);
            AstroTime after  = center_time.AddDays(+window);
            AstroTime t1 = Search(new SearchContext_EarthShadow(radius_limit, -1.0), before, center_time, 1.0) ??
                throw new InternalError("Failed to find start of shadow event.");
            AstroTime t2 = Search(new SearchContext_EarthShadow(radius_limit, +1.0), center_time, after, 1.0) ??
                throw new InternalError("Failed to find end of shadow event.");
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
        ///
        /// See #Astronomy.GlobalSolarEclipsesAfter for a convenient enumerator.
        /// </remarks>
        /// <param name="startTime">The date and time for starting the search for a solar eclipse.</param>
        public static GlobalSolarEclipseInfo SearchGlobalSolarEclipse(AstroTime startTime)
        {
            const double PruneLatitude = 1.8;   // Moon's ecliptic latitude beyond which eclipse is impossible

            // Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth.
            AstroTime nmtime = startTime;
            for (int nmcount=0; nmcount < 12; ++nmcount)
            {
                // Search for the next new moon. Any eclipse will be near it.
                AstroTime newmoon = SearchMoonPhase(0.0, nmtime, 40.0) ??
                    throw new InternalError("Failed to find next new moon.");

                // Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible.
                double eclip_lat = MoonEclipticLatitudeDegrees(newmoon);
                if (Math.Abs(eclip_lat) < PruneLatitude)
                {
                    // Search near the new moon for the time when the center of the Earth
                    // is closest to the line passing through the centers of the Sun and Moon.
                    ShadowInfo shadow = PeakMoonShadow(newmoon);
                    if (shadow.r < shadow.p + EARTH_MEAN_RADIUS_KM)
                    {
                        // This is at least a partial solar eclipse visible somewhere on Earth.
                        // Try to find an intersection between the shadow axis and the Earth's oblate geoid.
                        return GeoidIntersect(shadow);
                    }
                }

                // We didn't find an eclipse on this new moon, so search for the next one.
                nmtime = newmoon.AddDays(10.0);
            }

            // Safety valve to prevent infinite loop.
            // This should never happen, because at least 2 solar eclipses happen per year.
            throw new InternalError("Failure to find global solar eclipse.");
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
        ///
        /// See #Astronomy.GlobalSolarEclipsesAfter for a convenient enumerator.
        /// </remarks>
        /// <param name="prevEclipseTime">
        /// A date and time near a new moon. Solar eclipse search will start at the next new moon.
        /// </param>
        public static GlobalSolarEclipseInfo NextGlobalSolarEclipse(AstroTime prevEclipseTime)
        {
            AstroTime startTime = prevEclipseTime.AddDays(10.0);
            return SearchGlobalSolarEclipse(startTime);
        }


        /// <summary>Enumerates a series of global solar eclipses that occur after a specified time.</summary>
        /// <remarks>
        /// This is a convenience wrapper around the functions
        /// #Astronomy.SearchGlobalSolarEclipse and #Astronomy.NextGlobalSolarEclipse.
        /// </remarks>
        /// <param name="startTime">
        /// Specifies the time to begin searching for consecutive solar eclipses.
        /// </param>
        public static IEnumerable<GlobalSolarEclipseInfo> GlobalSolarEclipsesAfter(AstroTime startTime)
        {
            GlobalSolarEclipseInfo eclipse = SearchGlobalSolarEclipse(startTime);
            yield return eclipse;
            while (true)
            {
                eclipse = NextGlobalSolarEclipse(eclipse.peak);
                yield return eclipse;
            }
        }


        private static GlobalSolarEclipseInfo GeoidIntersect(ShadowInfo shadow)
        {
            var eclipse = new GlobalSolarEclipseInfo();
            eclipse.kind = EclipseKind.Partial;
            eclipse.peak = shadow.time;
            eclipse.distance = shadow.r;
            eclipse.latitude = eclipse.longitude = double.NaN;

            // We want to calculate the intersection of the shadow axis with the Earth's geoid.
            // First we must convert EQJ (equator of J2000) coordinates to EQD (equator of date)
            // coordinates that are perfectly aligned with the Earth's equator at this
            // moment in time.
            RotationMatrix rot = Rotation_EQJ_EQD(shadow.time);

            AstroVector v = RotateVector(rot, shadow.dir);        // shadow-axis vector in equator-of-date coordinates
            AstroVector e = RotateVector(rot, shadow.target);     // lunacentric Earth in equator-of-date coordinates

            // Convert all distances from AU to km.
            // But dilate the z-coordinates so that the Earth becomes a perfect sphere.
            // Then find the intersection of the vector with the sphere.
            // See p 184 in Montenbruck & Pfleger's "Astronomy on the Personal Computer", second edition.
            v.x *= KM_PER_AU;
            v.y *= KM_PER_AU;
            v.z *= KM_PER_AU / EARTH_FLATTENING;

            e.x *= KM_PER_AU;
            e.y *= KM_PER_AU;
            e.z *= KM_PER_AU / EARTH_FLATTENING;

            // Solve the quadratic equation that finds whether and where
            // the shadow axis intersects with the Earth in the dilated coordinate system.
            double R = EARTH_EQUATORIAL_RADIUS_KM;
            double A = v.x*v.x + v.y*v.y + v.z*v.z;
            double B = -2.0 * (v.x*e.x + v.y*e.y + v.z*e.z);
            double C = (e.x*e.x + e.y*e.y + e.z*e.z) - R*R;
            double radic = B*B - 4*A*C;

            if (radic > 0.0)
            {
                // Calculate the closer of the two intersection points.
                // This will be on the day side of the Earth.
                double u = (-B - Math.Sqrt(radic)) / (2 * A);

                // Convert lunacentric dilated coordinates to geocentric coordinates.
                double px = u*v.x - e.x;
                double py = u*v.y - e.y;
                double pz = (u*v.z - e.z) * EARTH_FLATTENING;

                // Convert cartesian coordinates into geodetic latitude/longitude.
                double proj = hypot(px, py) * (EARTH_FLATTENING * EARTH_FLATTENING);
                if (proj == 0.0)
                    eclipse.latitude = (pz > 0.0) ? +90.0 : -90.0;
                else
                    eclipse.latitude = RAD2DEG * Math.Atan(pz / proj);

                // Adjust longitude for Earth's rotation at the given UT.
                double gast = SiderealTime(eclipse.peak);
                eclipse.longitude = ((RAD2DEG*Math.Atan2(py, px)) - (15*gast)) % 360.0;
                if (eclipse.longitude <= -180.0)
                    eclipse.longitude += 360.0;
                else if (eclipse.longitude > +180.0)
                    eclipse.longitude -= 360.0;

                // We want to determine whether the observer sees a total eclipse or an annular eclipse.
                // We need to perform a series of vector calculations...
                // Calculate the inverse rotation matrix, so we can convert EQD to EQJ.
                RotationMatrix inv = InverseRotation(rot);

                // Put the EQD geocentric coordinates of the observer into the vector 'o'.
                // Also convert back from kilometers to astronomical units.
                var o = new AstroVector(px / KM_PER_AU, py / KM_PER_AU, pz / KM_PER_AU, shadow.time);

                // Rotate the observer's geocentric EQD back to the EQJ system.
                o = RotateVector(inv, o);

                // Convert geocentric vector to lunacentric vector.
                o.x += shadow.target.x;
                o.y += shadow.target.y;
                o.z += shadow.target.z;

                // Recalculate the shadow using a vector from the Moon's center toward the observer.
                ShadowInfo surface = CalcShadow(MOON_POLAR_RADIUS_KM, shadow.time, o, shadow.dir);

                // If we did everything right, the shadow distance should be very close to zero.
                // That's because we already determined the observer 'o' is on the shadow axis!
                if (surface.r > 1.0e-9 || surface.r < 0.0)
                    throw new InternalError("Invalid surface distance from intersection.");

                eclipse.kind = EclipseKindFromUmbra(surface.k);
                eclipse.obscuration = (eclipse.kind == EclipseKind.Total) ? 1.0 : SolarEclipseObscuration(shadow.dir, o);
            }
            else
            {
                // This is a partial solar eclipse. It does not make practical sense to calculate obscuration.
                // Anyone who wants obscuration should use Astronomy.SearchLocalSolarEclipse for a specific location on the Earth.
                eclipse.obscuration = double.NaN;
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
            // Search for when the Moon's shadow axis is closest to the center of the Earth.

            const double window = 0.03;     // days before/after new moon to search for minimum shadow distance
            AstroTime t1 = search_center_time.AddDays(-window);
            AstroTime t2 = search_center_time.AddDays(+window);
            AstroTime time = Search(moonShadowSlopeContext, t1, t2, 1.0) ??
                throw new InternalError("Failed to find Moon shadow event.");
            return MoonShadow(time);
        }

        private static ShadowInfo PeakLocalMoonShadow(AstroTime search_center_time, Observer observer)
        {
            // Search for the time near search_center_time that the Moon's shadow comes
            // closest to the given observer.
            const double window = 0.2;
            AstroTime t1 = search_center_time.AddDays(-window);
            AstroTime t2 = search_center_time.AddDays(+window);
            var context = new SearchContext_LocalMoonShadowSlope(observer);
            AstroTime time = Search(context, t1, t2, 1.0) ??
                throw new InternalError("Failed to find local Moon peak shadow event.");
            return LocalMoonShadow(time, observer);
        }

        private static ShadowInfo PeakPlanetShadow(Body body, double planet_radius_km, AstroTime search_center_time)
        {
            // Search for when the body's shadow is closest to the center of the Earth.
            const double window = 1.0;     // days before/after inferior conjunction to search for minimum shadow distance
            AstroTime t1 = search_center_time.AddDays(-window);
            AstroTime t2 = search_center_time.AddDays(+window);
            var context = new SearchContext_PlanetShadowSlope(body, planet_radius_km);
            AstroTime time = Search(context, t1, t2, 1.0) ??
                throw new InternalError("Failed to find peak planet shadow event.");
            return PlanetShadow(body, planet_radius_km, time);
        }

        private static ShadowInfo CalcShadow(
            double body_radius_km,
            AstroTime time,
            AstroVector target,
            AstroVector dir)
        {
            double u = (dir * target) / (dir * dir);
            double dx = (u * dir.x) - target.x;
            double dy = (u * dir.y) - target.y;
            double dz = (u * dir.z) - target.z;
            double r = KM_PER_AU * hypot(dx, dy, dz);
            double k = +SUN_RADIUS_KM - (1.0 + u)*(SUN_RADIUS_KM - body_radius_km);
            double p = -SUN_RADIUS_KM + (1.0 + u)*(SUN_RADIUS_KM + body_radius_km);
            return new ShadowInfo(time, u, r, k, p, target, dir);
        }


        internal static ShadowInfo EarthShadow(AstroTime time)
        {
            // This function helps find when the Earth's shadow falls upon the Moon.

            // Light-travel and aberration corrected vector from the Earth to the Sun.
            // The negative vector -s is thus the path of sunlight through the center of the Earth.
            AstroVector s = GeoVector(Body.Sun, time, Aberration.Corrected);

            // Geocentric Moon.
            AstroVector m = GeoMoon(time);

            return CalcShadow(EARTH_ECLIPSE_RADIUS_KM, time, m, -s);
        }


        internal static ShadowInfo MoonShadow(AstroTime time)
        {
            // This function helps find when the Moon's shadow falls upon the Earth.

            AstroVector s = GeoVector(Body.Sun, time, Aberration.Corrected);
            AstroVector m = GeoMoon(time);      // geocentric Moon

            // -m  = lunacentric Earth
            // m-s = heliocentric Moon
            return CalcShadow(MOON_MEAN_RADIUS_KM, time, -m, m-s);
        }


        internal static ShadowInfo LocalMoonShadow(AstroTime time, Observer observer)
        {
            // Calculate observer's geocentric position.
            AstroVector o = geo_pos(time, observer);

            // Calculate light-travel and aberration corrected Sun.
            AstroVector s = GeoVector(Body.Sun, time, Aberration.Corrected);

            // Calculate geocentric Moon.
            AstroVector m = GeoMoon(time);

            // o-m = lunacentric observer
            // m-s = heliocentric Moon
            return CalcShadow(MOON_MEAN_RADIUS_KM, time, o-m, m-s);
        }


        internal static ShadowInfo PlanetShadow(Body body, double planet_radius_km, AstroTime time)
        {
            // Calculate light-travel-corrected vector from Earth to planet.
            AstroVector p = GeoVector(body, time, Aberration.Corrected);

            // Calculate light-travel-corrected vector from Earth to Sun.
            AstroVector s = GeoVector(Body.Sun, time, Aberration.Corrected);

            // -p  = planetcentric Earth
            // p-s = heliocentric planet
            return CalcShadow(planet_radius_km, time, -p, p-s);
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
        ///
        /// See #LocalSolarEclipseInfo for more information about this topic.
        /// See #Astronomy.LocalSolarEclipsesAfter for a convenient enumerator.
        /// </remarks>
        ///
        /// <param name="startTime">The date and time for starting the search for a solar eclipse.</param>
        /// <param name="observer">The geographic location of the observer.</param>
        public static LocalSolarEclipseInfo SearchLocalSolarEclipse(AstroTime startTime, Observer observer)
        {
            const double PruneLatitude = 1.8;   // Moon's ecliptic latitude beyond which eclipse is impossible

            // Iterate through consecutive new moons until we find a solar eclipse visible somewhere on Earth.
            AstroTime nmtime = startTime;
            for(;;)
            {
                // Search for the next new moon. Any eclipse will be near it.
                AstroTime newmoon = SearchMoonPhase(0.0, nmtime, 40.0) ??
                    throw new InternalError("Failed to find next new moon.");

                // Pruning: if the new moon's ecliptic latitude is too large, a solar eclipse is not possible.
                double eclip_lat = MoonEclipticLatitudeDegrees(newmoon);
                if (Math.Abs(eclip_lat) < PruneLatitude)
                {
                    // Search near the new moon for the time when the observer
                    // is closest to the line passing through the centers of the Sun and Moon.
                    ShadowInfo shadow = PeakLocalMoonShadow(newmoon, observer);
                    if (shadow.r < shadow.p)
                    {
                        // This is at least a partial solar eclipse for the observer.
                        LocalSolarEclipseInfo eclipse = LocalEclipse(shadow, observer);

                        // Ignore any eclipse that happens completely at night.
                        // More precisely, the center of the Sun must be above the horizon
                        // at the beginning or the end of the eclipse, or we skip the event.
                        if (eclipse.partial_begin.altitude > 0.0 || eclipse.partial_end.altitude > 0.0)
                            return eclipse;
                    }
                }

                // We didn't find an eclipse on this new moon, so search for the next one.
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
        ///
        /// See #Astronomy.LocalSolarEclipsesAfter for a convenient enumerator.
        /// </remarks>
        ///
        /// <param name="prevEclipseTime">
        /// A date and time near a new moon. Solar eclipse search will start at the next new moon.
        /// </param>
        ///
        /// <param name="observer">
        /// The geographic location of the observer.
        /// </param>
        public static LocalSolarEclipseInfo NextLocalSolarEclipse(AstroTime prevEclipseTime, Observer observer)
        {
            AstroTime startTime = prevEclipseTime.AddDays(10.0);
            return SearchLocalSolarEclipse(startTime, observer);
        }


        /// <summary>Enumerates a series of local solar eclipses that occur after a specified time.</summary>
        /// <remarks>
        /// This is a convenience wrapper around the functions
        /// #Astronomy.SearchLocalSolarEclipse and #Astronomy.NextLocalSolarEclipse.
        /// </remarks>
        /// <param name="startTime">
        /// Specifies the time to begin searching for consecutive solar eclipses.
        /// </param>
        /// <param name="observer">
        /// The geographic location of the observer.
        /// </param>
        public static IEnumerable<LocalSolarEclipseInfo> LocalSolarEclipsesAfter(
            AstroTime startTime,
            Observer observer)
        {
            LocalSolarEclipseInfo eclipse = SearchLocalSolarEclipse(startTime, observer);
            yield return eclipse;
            while (true)
            {
                eclipse = NextLocalSolarEclipse(eclipse.peak.time, observer);
                yield return eclipse;
            }
        }


        private static double local_partial_distance(ShadowInfo shadow)
        {
            return shadow.p - shadow.r;
        }

        private static double local_total_distance(ShadowInfo shadow)
        {
            // Must take the absolute value of the umbra radius 'k'
            // because it can be negative for an annular eclipse.
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

            if (shadow.r < Math.Abs(shadow.k))      // take absolute value of 'k' to handle annular eclipses too.
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

            eclipse.obscuration = (eclipse.kind == EclipseKind.Total) ? 1.0 : SolarEclipseObscuration(shadow.dir, shadow.target);
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
            AstroTime search = Search(context, t1, t2, 1.0) ??
                throw new InternalError("Local eclipse transition search failed.");
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
            // Search for the time the planet's penumbra begins/ends making contact with the center of the Earth.
            var context = new SearchContext_PlanetShadowBoundary(body, planet_radius_km, direction);
            AstroTime time = Search(context, t1, t2, 1.0) ??
                throw new InternalError("Planet transit boundary search failed");
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
        ///
        /// See #Astronomy.TransitsAfter for a convenient enumerator.
        /// </remarks>
        /// <param name="body">
        /// The planet whose transit is to be found. Must be `Body.Mercury` or `Body.Venus`.
        /// </param>
        /// <param name="startTime">
        /// The date and time for starting the search for a transit.
        /// </param>
        public static TransitInfo SearchTransit(Body body, AstroTime startTime)
        {
            const double threshold_angle = 0.4;     // maximum angular separation to attempt transit calculation
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
                // Search for the next inferior conjunction of the given planet.
                // This is the next time the Earth and the other planet have the same
                // ecliptic longitude as seen from the Sun.
                AstroTime conj = SearchRelativeLongitude(body, 0.0, search_time);

                // Calculate the angular separation between the body and the Sun at this time.
                double separation = AngleFromSun(body, conj);

                if (separation < threshold_angle)
                {
                    // The planet's angular separation from the Sun is small enough
                    // to consider it a transit candidate.
                    // Search for the moment when the line passing through the Sun
                    // and planet are closest to the Earth's center.
                    ShadowInfo shadow = PeakPlanetShadow(body, planet_radius_km, conj);

                    if (shadow.r < shadow.p)        // does the planet's penumbra touch the Earth's center?
                    {
                        var transit = new TransitInfo();

                        // Find the beginning and end of the penumbral contact.
                        AstroTime tx = shadow.time.AddDays(-dt_days);
                        transit.start = PlanetTransitBoundary(body, planet_radius_km, tx, shadow.time, -1.0);

                        tx = shadow.time.AddDays(+dt_days);
                        transit.finish = PlanetTransitBoundary(body, planet_radius_km, shadow.time, tx, +1.0);

                        transit.peak = shadow.time;
                        transit.separation = 60.0 * AngleFromSun(body, shadow.time);
                        return transit;
                    }
                }

                // This inferior conjunction was not a transit. Try the next inferior conjunction.
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
        ///
        /// See #Astronomy.TransitsAfter for a convenient enumerator.
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

        /// <summary>Enumerates a series of transits of Mercury or Venus.</summary>
        /// <remarks>
        /// This is a convenience wrapper around the functions
        /// #Astronomy.SearchTransit and #Astronomy.NextTransit.
        /// </remarks>
        /// <param name="body">
        /// The planet whose transits are to be enumerated. Must be `Body.Mercury` or `Body.Venus`.
        /// </param>
        /// <param name="startTime">
        /// Specifies the time to begin searching for consecutive transits.
        /// </param>
        public static IEnumerable<TransitInfo> TransitsAfter(Body body, AstroTime startTime)
        {
            TransitInfo transit = SearchTransit(body, startTime);
            yield return transit;
            while (true)
            {
                transit = NextTransit(body, transit.peak);
                yield return transit;
            }
        }


        private const double MoonNodeStepDays = +10.0; // a safe number of days to step without missing a Moon node

        /// <summary>
        /// Searches for a time when the Moon's center crosses through the ecliptic plane.
        /// </summary>
        /// <remarks>
        /// Searches for the first ascending or descending node of the Moon after `startTime`.
        /// An ascending node is when the Moon's center passes through the ecliptic plane
        /// (the plane of the Earth's orbit around the Sun) from south to north.
        /// A descending node is when the Moon's center passes through the ecliptic plane
        /// from north to south. Nodes indicate possible times of solar or lunar eclipses,
        /// if the Moon also happens to be in the correct phase (new or full, respectively).
        /// Call `Astronomy.SearchMoonNode` to find the first of a series of nodes.
        /// Then call #Astronomy.NextMoonNode to find as many more consecutive nodes as desired.
        ///
        /// See #Astronomy.MoonNodesAfter for a convenient enumerator.
        /// </remarks>
        /// <param name="startTime">
        /// The date and time for starting the search for an ascending or descending node of the Moon.
        /// </param>
        public static NodeEventInfo SearchMoonNode(AstroTime startTime)
        {
            // Start at the given moment in time and sample the Moon's ecliptic latitude.
            // Step 10 days at a time, searching for an interval where that latitude crosses zero.
            AstroTime time1 = startTime;
            Spherical eclip1 = EclipticGeoMoon(time1);
            var context = new SearchContext_MoonNode();

            for(;;)
            {
                AstroTime time2 = time1.AddDays(MoonNodeStepDays);
                Spherical eclip2 = EclipticGeoMoon(time2);
                if (eclip1.lat * eclip2.lat <= 0.0)
                {
                    // There is a node somewhere inside this closed time interval.
                    // Figure out whether it is an ascending node or a descending node.
                    NodeEventKind kind;
                    if (eclip2.lat > eclip1.lat)
                    {
                        context.Direction = +1;
                        kind = NodeEventKind.Ascending;
                    }
                    else
                    {
                        context.Direction = -1;
                        kind = NodeEventKind.Descending;
                    }
                    AstroTime result = Search(context, time1, time2, 1.0) ??
                        throw new InternalError("Could not find Moon node.");

                    return new NodeEventInfo { time = result, kind = kind };
                }
                time1 = time2;
                eclip1 = eclip2;
            }
        }

        /// <summary>
        /// Searches for the next time when the Moon's center crosses through the ecliptic plane.
        /// </summary>
        /// <remarks>
        /// Call #Astronomy.SearchMoonNode to find the first of a series of nodes.
        /// Then call `Astronomy.NextMoonNode` to find as many more consecutive nodes as desired.
        ///
        /// See #Astronomy.MoonNodesAfter for a convenient enumerator.
        /// </remarks>
        /// <param name="prevNode">
        /// The previous node found from calling #Astronomy.SearchMoonNode or `Astronomy.NextMoonNode`.
        /// </param>
        public static NodeEventInfo NextMoonNode(NodeEventInfo prevNode)
        {
            AstroTime time = prevNode.time.AddDays(MoonNodeStepDays);
            NodeEventInfo node = SearchMoonNode(time);
            switch (prevNode.kind)
            {
                case NodeEventKind.Ascending:
                    if (node.kind != NodeEventKind.Descending)
                        throw new InternalError("previous node was ascending, but this node was: " + node.kind);
                    break;

                case NodeEventKind.Descending:
                    if (node.kind != NodeEventKind.Ascending)
                        throw new InternalError("previous node was descending, but this node was: " + node.kind);
                    break;

                default:
                    throw new ArgumentException("Previous node has an invalid node kind.");
            }
            return node;
        }


        /// <summary>Enumerates a series of ascending/descending nodes of the Moon that occur after a specified time.</summary>
        /// <remarks>
        /// This is a convenience wrapper around the functions
        /// #Astronomy.SearchMoonNode and #Astronomy.NextMoonNode.
        /// </remarks>
        /// <param name="startTime">
        /// Specifies the time to begin searching for consecutive lunar apsides.
        /// </param>
        public static IEnumerable<NodeEventInfo> MoonNodesAfter(AstroTime startTime)
        {
            NodeEventInfo node = SearchMoonNode(startTime);
            yield return node;
            while (true)
            {
                node = NextMoonNode(node);
                yield return node;
            }
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
                gc = -earth;
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
                    hc = earth + gc;
                }
                else
                {
                    // For planets, the heliocentric vector is more direct to calculate.
                    hc = HelioVector(body, time);
                    gc = hc - earth;
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
            // https://astronomy.stackexchange.com/questions/10246/is-there-a-simple-analytical-formula-for-the-lunar-phase-brightness-curve
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
            // For Mercury and Venus, see:  https://iopscience.iop.org/article/10.1086/430212
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
            // Based on formulas by Paul Schlyter found here:
            // http://www.stjarnhimlen.se/comp/ppcomp.html#15

            // We must handle Saturn's rings as a major component of its visual magnitude.
            // Find geocentric ecliptic coordinates of Saturn.
            Ecliptic eclip = EquatorialToEcliptic(gc);

            double ir = DEG2RAD * 28.06;   // tilt of Saturn's rings to the ecliptic, in radians
            double Nr = DEG2RAD * (169.51 + (3.82e-5 * time.tt));    // ascending node of Saturn's rings, in radians

            // Find tilt of Saturn's rings, as seen from Earth.
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
            // s1 and s2 are relative longitudes within which peak magnitude of Venus can occur.
            const double s1 = 10.0;
            const double s2 = 30.0;

            if (body != Body.Venus)
                throw new InvalidBodyException(body);

            var mag_slope = new SearchContext_MagnitudeSlope(body);

            int iter = 0;
            while (++iter <= 2)
            {
                // Find current heliocentric relative longitude between the
                // inferior planet and the Earth.
                double plon = EclipticLongitude(body, startTime);
                double elon = EclipticLongitude(Body.Earth, startTime);
                double rlon = LongitudeOffset(plon - elon);     // clamp to (-180, +180].

                // The slope function is not well-behaved when rlon is near 0 degrees or 180 degrees
                // because there is a cusp there that causes a discontinuity in the derivative.
                // So we need to guard against searching near such times.

                double rlon_lo, rlon_hi, adjust_days, syn;
                if (rlon >= -s1 && rlon < +s1)
                {
                    // Seek to the window [+s1, +s2].
                    adjust_days = 0.0;
                    // Search forward for the time t1 when rel lon = +s1.
                    rlon_lo = +s1;
                    // Search forward for the time t2 when rel lon = +s2.
                    rlon_hi = +s2;
                }
                else if (rlon >= +s2 || rlon < -s2)
                {
                    // Seek to the next search window at [-s2, -s1].
                    adjust_days = 0.0;
                    // Search forward for the time t1 when rel lon = -s2.
                    rlon_lo = -s2;
                    // Search forward for the time t2 when rel lon = -s1.
                    rlon_hi = -s1;
                }
                else if (rlon >= 0)
                {
                    // rlon must be in the middle of the window [+s1, +s2].
                    // Search BACKWARD for the time t1 when rel lon = +s1.
                    syn = SynodicPeriod(body);
                    adjust_days = -syn / 4;
                    rlon_lo = +s1;
                    // Search forward from t1 to find t2 such that rel lon = +s2.
                    rlon_hi = +s2;
                }
                else
                {
                    // rlon must be in the middle of the window [-s2, -s1].
                    // Search BACKWARD for the time t1 when rel lon = -s2.
                    syn = SynodicPeriod(body);
                    adjust_days = -syn / 4;
                    rlon_lo = -s2;
                    // Search forward from t1 to find t2 such that rel lon = -s1.
                    rlon_hi = -s1;
                }
                AstroTime t_start = startTime.AddDays(adjust_days);
                AstroTime t1 = SearchRelativeLongitude(body, rlon_lo, t_start);
                AstroTime t2 = SearchRelativeLongitude(body, rlon_hi, t1);

                // Now we have a time range [t1,t2] that brackets a maximum magnitude event.
                // Confirm the bracketing.
                double m1 = mag_slope.Eval(t1);
                if (m1 >= 0.0)
                    throw new InternalError("m1 >= 0");    // should never happen!

                double m2 = mag_slope.Eval(t2);
                if (m2 <= 0.0)
                    throw new InternalError("m2 <= 0");    // should never happen!

                // Use the generic search algorithm to home in on where the slope crosses from negative to positive.
                AstroTime tx = Search(mag_slope, t1, t2, 10.0) ??
                    throw new InternalError("Failed to find magnitude slope transition.");

                if (tx.tt >= startTime.tt)
                    return Illumination(body, tx);

                // This event is in the past (earlier than startTime).
                // We need to search forward from t2 to find the next possible window.
                // We never need to search more than twice.
                startTime = t2.AddDays(1.0);
            }
            // This should never happen. If it does, please report as a bug in Astronomy Engine.
            throw new InternalError("Peak magnitude search failed.");
        }


        private static double CubeRoot(double x)
        {
            // Astronomy Engine is targeted at .NET Standard 2.0.
            // That means it supports the older Framework 4+ platform
            // as well as .NET Core 5+.
            // .NET Core has a Math.Cbrt function, but .NET Framework doesn't.
            // Therefore, I have to implement my own cube root function.

            if (x < 0.0)
                return -CubeRoot(-x);

            if (x == 0.0)
                return 0.0;

            return Math.Pow(x, (1.0 / 3.0));
        }

        /// <summary>
        /// Calculates one of the 5 Lagrange points for a pair of co-orbiting bodies.
        /// </summary>
        /// <remarks>
        /// Given a more massive "major" body and a much less massive "minor" body,
        /// calculates one of the five Lagrange points in relation to the minor body's
        /// orbit around the major body. The parameter `point` is an integer that
        /// selects the Lagrange point as follows:
        ///
        /// 1 = the Lagrange point between the major body and minor body.
        /// 2 = the Lagrange point on the far side of the minor body.
        /// 3 = the Lagrange point on the far side of the major body.
        /// 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
        /// 5 = the Lagrange point 60 degrees behind the minor body's orbital position.
        ///
        /// The function returns the state vector for the selected Lagrange point
        /// in J2000 mean equator coordinates (EQJ), with respect to the center of the
        /// major body.
        ///
        /// To calculate Sun/Earth Lagrange points, pass in `Body.Sun` for `major_body`
        /// and `Body.EMB` (Earth/Moon barycenter) for `minor_body`.
        /// For Lagrange points of the Sun and any other planet, pass in just that planet
        /// (e.g. `Body.Jupiter`) for `minor_body`.
        /// To calculate Earth/Moon Lagrange points, pass in `Body.Earth` and `Body.Moon`
        /// for the major and minor bodies respectively.
        ///
        /// In some cases, it may be more efficient to call #Astronomy.LagrangePointFast,
        /// especially when the state vectors have already been calculated, or are needed
        /// for some other purpose.
        /// </remarks>
        /// <param name="point">An integer 1..5 that selects which of the Lagrange points to calculate.</param>
        /// <param name="time">The time for which the Lagrange point is to be calculated.</param>
        /// <param name="major_body">The more massive of the co-orbiting bodies: `Body.Sun` or `Body.Earth`.</param>
        /// <param name="minor_body">The less massive of the co-orbiting bodies. See main remarks.</param>
        /// <returns>The position and velocity of the selected Lagrange point with respect to the major body's center.</returns>
        public static StateVector LagrangePoint(
            int point,
            AstroTime time,
            Body major_body,
            Body minor_body)
        {
            double major_mass = MassProduct(major_body);
            double minor_mass = MassProduct(minor_body);

            StateVector major_state;
            StateVector minor_state;

            // Calculate the state vectors for the major and minor bodies.
            if (major_body == Body.Earth && minor_body == Body.Moon)
            {
                // Use geocentric calculations for more precision.

                // The Earth's geocentric state is trivial.
                major_state.t = time;
                major_state.x = major_state.y = major_state.z = 0.0;
                major_state.vx = major_state.vy = major_state.vz = 0.0;

                minor_state = GeoMoonState(time);
            }
            else
            {
                major_state = HelioState(major_body, time);
                minor_state = HelioState(minor_body, time);
            }

            return LagrangePointFast(
                point,
                major_state,
                major_mass,
                minor_state,
                minor_mass
            );
        }

        /// <summary>
        /// Calculates one of the 5 Lagrange points from body masses and state vectors.
        /// </summary>
        /// <remarks>
        /// Given a more massive "major" body and a much less massive "minor" body,
        /// calculates one of the five Lagrange points in relation to the minor body's
        /// orbit around the major body. The parameter `point` is an integer that
        /// selects the Lagrange point as follows:
        ///
        /// 1 = the Lagrange point between the major body and minor body.
        /// 2 = the Lagrange point on the far side of the minor body.
        /// 3 = the Lagrange point on the far side of the major body.
        /// 4 = the Lagrange point 60 degrees ahead of the minor body's orbital position.
        /// 5 = the Lagrange point 60 degrees behind the minor body's orbital position.
        ///
        /// The caller passes in the state vector and mass for both bodies.
        /// The state vectors can be in any orientation and frame of reference.
        /// The body masses are expressed as GM products, where G = the universal
        /// gravitation constant and M = the body's mass. Thus the units for
        /// `major_mass` and `minor_mass` must be au^3/day^2.
        /// Use #Astronomy.MassProduct to obtain GM values for various solar system bodies.
        ///
        /// The function returns the state vector for the selected Lagrange point
        /// using the same orientation as the state vector parameters `major_state` and `minor_state`,
        /// and the position and velocity components are with respect to the major body's center.
        ///
        /// Consider calling #Astronomy.LagrangePoint, instead of this function, for simpler usage in most cases.
        /// </remarks>
        /// <param name="point">An integer 1..5 that selects which of the Lagrange points to calculate.</param>
        /// <param name="major_state">The state vector of the major (more massive) of the pair of bodies.</param>
        /// <param name="major_mass">The mass product GM of the major body.</param>
        /// <param name="minor_state">The state vector of the minor (less massive) of the pair of bodies.</param>
        /// <param name="minor_mass">The mass product GM of the minor body.</param>
        /// <returns>The position and velocity of the selected Lagrange point with respect to the major body's center.</returns>
        public static StateVector LagrangePointFast(
            int point,
            StateVector major_state,
            double major_mass,
            StateVector minor_state,
            double minor_mass)
        {
            const double cos_60 = 0.5;
            const double sin_60 = 0.8660254037844386;   // sqrt(3) / 2

            if (point < 1 || point > 5)
                throw new ArgumentException($"Invalid lagrange point {point}");

            if (!isfinite(major_mass) || major_mass <= 0.0)
                throw new ArgumentException("Major mass must be a positive number.");

            if (!isfinite(minor_mass) || minor_mass <= 0.0)
                throw new ArgumentException("Minor mass must be a positive number.");

            // Find the relative position vector <dx, dy, dz>.
            double dx = minor_state.x - major_state.x;
            double dy = minor_state.y - major_state.y;
            double dz = minor_state.z - major_state.z;
            double R2 = (dx*dx + dy*dy + dz*dz);

            // R = Total distance between the bodies.
            double R = Math.Sqrt(R2);

            // Find the velocity vector <vx, vy, vz>.
            double vx = minor_state.vx - major_state.vx;
            double vy = minor_state.vy - major_state.vy;
            double vz = minor_state.vz - major_state.vz;

            StateVector p;
            if (point == 4 || point == 5)
            {
                // For L4 and L5, we need to find points 60 degrees away from the
                // line connecting the two bodies and in the instantaneous orbital plane.
                // Define the instantaneous orbital plane as the unique plane that contains
                // both the relative position vector and the relative velocity vector.

                // Take the cross product of position and velocity to find a normal vector <nx, ny, nz>.
                double nx = dy*vz - dz*vy;
                double ny = dz*vx - dx*vz;
                double nz = dx*vy - dy*vx;

                // Take the cross product normal*position to get a tangential vector <ux, uy, uz>.
                double ux = ny*dz - nz*dy;
                double uy = nz*dx - nx*dz;
                double uz = nx*dy - ny*dx;

                // Convert the tangential direction vector to a unit vector.
                double U = hypot(ux, uy, uz);
                ux /= U;
                uy /= U;
                uz /= U;

                // Convert the relative position vector into a unit vector.
                dx /= R;
                dy /= R;
                dz /= R;

                // Now we have two perpendicular unit vectors in the orbital plane: 'd' and 'u'.

                // Create new unit vectors rotated (+/-)60 degrees from the radius/tangent directions.
                double vert = (point == 4) ? +sin_60 : -sin_60;

                // Rotated radial vector
                double Dx = cos_60*dx + vert*ux;
                double Dy = cos_60*dy + vert*uy;
                double Dz = cos_60*dz + vert*uz;

                // Rotated tangent vector
                double Ux = cos_60*ux - vert*dx;
                double Uy = cos_60*uy - vert*dy;
                double Uz = cos_60*uz - vert*dz;

                // Calculate L4/L5 positions relative to the major body.
                p.x = R * Dx;
                p.y = R * Dy;
                p.z = R * Dz;

                // Use dot products to find radial and tangential components of the relative velocity.
                double vrad = vx*dx + vy*dy + vz*dz;
                double vtan = vx*ux + vy*uy + vz*uz;

                // Calculate L4/L5 velocities.
                p.vx = vrad*Dx + vtan*Ux;
                p.vy = vrad*Dy + vtan*Uy;
                p.vz = vrad*Dz + vtan*Uz;
            }
            else
            {
                // Calculate the distances of each body from their mutual barycenter.
                // r1 = negative distance of major mass from barycenter (e.g. Sun to the left of barycenter)
                // r2 = positive distance of minor mass from barycenter (e.g. Earth to the right of barycenter)
                double r1 = -R * (minor_mass / (major_mass + minor_mass));
                double r2 = +R * (major_mass / (major_mass + minor_mass));

                // Calculate the square of the angular orbital speed in [rad^2 / day^2].
                double omega2 = (major_mass + minor_mass) / (R2*R);

                // Use Newton's Method to numerically solve for the location where
                // outward centrifugal acceleration in the rotating frame of reference
                // is equal to net inward gravitational acceleration.
                // First derive a good initial guess based on approximate analysis.
                double scale, numer1, numer2;
                if (point == 1 || point == 2)
                {
                    scale = (major_mass / (major_mass + minor_mass)) * CubeRoot(minor_mass / (3.0 * major_mass));
                    numer1 = -major_mass;    // The major mass is to the left of L1 and L2
                    if (point == 1)
                    {
                        scale = 1.0 - scale;
                        numer2 = +minor_mass;    // The minor mass is to the right of L1.
                    }
                    else
                    {
                        scale = 1.0 + scale;
                        numer2 = -minor_mass;    // The minor mass is to the left of L2.
                    }
                }
                else // point == 3
                {
                    scale = ((7.0/12.0)*minor_mass - major_mass) / (minor_mass + major_mass);
                    numer1 = +major_mass;    // major mass is to the right of L3.
                    numer2 = +minor_mass;    // minor mass is to the right of L3.
                }

                // Iterate Newton's Method until it converges.
                double x = R*scale - r1;
                double deltax;
                do
                {
                    double dr1 = x - r1;
                    double dr2 = x - r2;
                    double accel = omega2*x + numer1/(dr1*dr1) + numer2/(dr2*dr2);
                    double deriv = omega2 - 2*numer1/(dr1*dr1*dr1) - 2*numer2/(dr2*dr2*dr2);
                    deltax = accel/deriv;
                    x -= deltax;
                }
                while (Math.Abs(deltax/R) > 1.0e-14);
                scale = (x - r1) / R;

                p.x  = scale * dx;
                p.y  = scale * dy;
                p.z  = scale * dz;
                p.vx = scale * vx;
                p.vy = scale * vy;
                p.vz = scale * vz;
            }
            p.t = major_state.t;
            return p;
        }

        /// <summary>Calculates the inverse of a rotation matrix.</summary>
        /// <remarks>
        /// Given a rotation matrix that performs some coordinate transform,
        /// this function returns the matrix that reverses that transform.
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

        /// <summary>Creates an identity rotation matrix.</summary>
        /// <remarks>
        /// Returns a rotation matrix that has no effect on orientation.
        /// This matrix can be the starting point for other operations,
        /// such as using a series of calls to #Astronomy.Pivot to
        /// create a custom rotation matrix.
        /// </remarks>
        /// <returns>The identity matrix.</returns>
        public static RotationMatrix IdentityMatrix()
        {
            var rot = new double[3, 3]
            {
                { 1, 0, 0 },
                { 0, 1, 0 },
                { 0, 0, 1 }
            };

            return new RotationMatrix(rot);
        }

        /// <summary>Re-orients a rotation matrix by pivoting it by an angle around one of its axes.</summary>
        /// <remarks>
        /// Given a rotation matrix, a selected coordinate axis, and an angle in degrees,
        /// this function pivots the rotation matrix by that angle around that coordinate axis.
        ///
        /// For example, if you have rotation matrix that converts ecliptic coordinates (ECL)
        /// to horizontal coordinates (HOR), but you really want to convert ECL to the orientation
        /// of a telescope camera pointed at a given body, you can use `Astronomy.Pivot` twice:
        /// (1) pivot around the zenith axis by the body's azimuth, then (2) pivot around the
        /// western axis by the body's altitude angle. The resulting rotation matrix will then
        /// reorient ECL coordinates to the orientation of your telescope camera.
        /// </remarks>
        ///
        /// <param name="rotation">The input rotation matrix.</param>
        ///
        /// <param name="axis">
        ///      An integer that selects which coordinate axis to rotate around:
        ///      0 = x, 1 = y, 2 = z. Any other value will cause an ArgumentException to be thrown.
        /// </param>
        ///
        /// <param name="angle">
        ///      An angle in degrees indicating the amount of rotation around the specified axis.
        ///      Positive angles indicate rotation counterclockwise as seen from the positive
        ///      direction along that axis, looking towards the origin point of the orientation system.
        ///      Any finite number of degrees is allowed, but best precision will result from keeping
        ///      `angle` in the range [-360, +360].
        /// </param>
        ///
        /// <returns>A pivoted matrix object.</returns>
        public static RotationMatrix Pivot(RotationMatrix rotation, int axis, double angle)
        {
            // Check for an invalid coordinate axis.
            if (axis < 0 || axis > 2)
                throw new ArgumentException($"Invalid coordinate axis = {axis}. Must be 0..2.");

            // Check for an invalid angle value.
            if (!isfinite(angle))
                throw new ArgumentException("Angle is not a finite number.");

            double radians = angle * DEG2RAD;
            double c = Math.Cos(radians);
            double s = Math.Sin(radians);

            // We need to maintain the "right-hand" rule, no matter which
            // axis was selected. That means we pick (i, j, k) axis order
            // such that the following vector cross product is satisfied:
            // i x j = k
            int i = (axis + 1) % 3;
            int j = (axis + 2) % 3;
            int k = axis;

            var rot = new double[3, 3];

            rot[i, i] = c*rotation.rot[i, i] - s*rotation.rot[i, j];
            rot[i, j] = s*rotation.rot[i, i] + c*rotation.rot[i, j];
            rot[i, k] = rotation.rot[i, k];

            rot[j, i] = c*rotation.rot[j, i] - s*rotation.rot[j, j];
            rot[j, j] = s*rotation.rot[j, i] + c*rotation.rot[j, j];
            rot[j, k] = rotation.rot[j, k];

            rot[k, i] = c*rotation.rot[k, i] - s*rotation.rot[k, j];
            rot[k, j] = s*rotation.rot[k, i] + c*rotation.rot[k, j];
            rot[k, k] = rotation.rot[k, k];

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

        /// <summary>Applies a rotation to a state vector, yielding a rotated state vector.</summary>
        /// <remarks>
        /// This function transforms a state vector in one orientation to a state vector in another orientation.
        /// </remarks>
        /// <param name="rotation">A rotation matrix that specifies how the orientation of the state vector is to be changed.</param>
        /// <param name="state">The state vector whose orientation is to be changed.</param>
        /// <returns>A state vector in the orientation specified by `rotation`.</returns>
        public static StateVector RotateState(RotationMatrix rotation, StateVector state)
        {
            return new StateVector(
                rotation.rot[0, 0]*state.x + rotation.rot[1, 0]*state.y + rotation.rot[2, 0]*state.z,
                rotation.rot[0, 1]*state.x + rotation.rot[1, 1]*state.y + rotation.rot[2, 1]*state.z,
                rotation.rot[0, 2]*state.x + rotation.rot[1, 2]*state.y + rotation.rot[2, 2]*state.z,
                rotation.rot[0, 0]*state.vx + rotation.rot[1, 0]*state.vy + rotation.rot[2, 0]*state.vz,
                rotation.rot[0, 1]*state.vx + rotation.rot[1, 1]*state.vy + rotation.rot[2, 1]*state.vz,
                rotation.rot[0, 2]*state.vx + rotation.rot[1, 2]*state.vy + rotation.rot[2, 2]*state.vz,
                state.t
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
                    // Indeterminate coordinates; pos vector has zero length.
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


        /// <summary>Given an equatorial vector, calculates equatorial angular coordinates.</summary>
        /// <param name="vector">A vector in an equatorial coordinate system.</param>
        /// <returns>Angular coordinates expressed in the same equatorial system as `vector`.</returns>
        public static Equatorial EquatorFromVector(AstroVector vector)
        {
            Spherical sphere = SphereFromVector(vector);
            return new Equatorial(sphere.lon / 15.0, sphere.lat, sphere.dist, vector);
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
        /// `lat` holds the refracted altitude angle,
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
        /// This function works best near sea level.
        /// To correct for higher elevations, call #Astronomy.Atmosphere for that
        /// elevation and multiply the refraction angle by the resulting relative density.
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
                return 0.0;     // no attempt to correct an invalid altitude

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
                // No refraction, or the refraction option is invalid.
                refr = 0.0;
            }

            return refr;
        }

        /// <summary>
        /// Calculates the inverse of an atmospheric refraction angle.
        /// </summary>
        /// <remarks>
        /// Given an observed altitude angle that includes atmospheric refraction,
        /// calculates the negative angular correction to obtain the unrefracted
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
                return 0.0;     // no attempt to correct an invalid altitude

            // Find the pre-adjusted altitude whose refraction correction leads to 'altitude'.
            double altitude = bent_altitude - RefractionAngle(refraction, bent_altitude);
            for(;;)
            {
                // See how close we got.
                double diff = (altitude + RefractionAngle(refraction, altitude)) - bent_altitude;
                if (Math.Abs(diff) < 1.0e-14)
                    return altitude - bent_altitude;

                altitude -= diff;
            }
        }


        private static AxisInfo EarthRotationAxis(AstroTime time)
        {
            AxisInfo axis;

            // Unlike the other planets, we have a model of precession and nutation
            // for the Earth's axis that provides a north pole vector.
            // So calculate the vector first, then derive the (RA,DEC) angles from the vector.

            // Start with a north pole vector in equator-of-date coordinates: (0,0,1).
            var pos1 = new AstroVector(0.0, 0.0, 1.0, time);

            // Convert the vector into J2000 coordinates.
            AstroVector pos2 = nutation(pos1, PrecessDirection.Into2000);
            axis.north = precession(pos2, PrecessDirection.Into2000);

            // Derive angular values: right ascension and declination.
            Equatorial equ = Astronomy.EquatorFromVector(axis.north);
            axis.ra = equ.ra;
            axis.dec = equ.dec;

            // Use a modified version of the era() function that does not trim to 0..360 degrees.
            // This expression is also corrected to give the correct angle at the J2000 epoch.
            axis.spin = 190.41375788700253 + (360.9856122880876 * time.ut);

            return axis;
        }


        /// <summary>
        /// Calculates information about a body's rotation axis at a given time.
        /// </summary>
        /// <remarks>
        /// Calculates the orientation of a body's rotation axis, along with
        /// the rotation angle of its prime meridian, at a given moment in time.
        ///
        /// This function uses formulas standardized by the IAU Working Group
        /// on Cartographics and Rotational Elements 2015 report, as described
        /// in the following document:
        ///
        /// https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf
        ///
        /// See #AxisInfo for more detailed information.
        /// </remarks>
        /// <param name="body">
        /// One of the following values:
        /// `Body.Sun`, `Body.Moon`, `Body.Mercury`, `Body.Venus`, `Body.Earth`, `Body.Mars`,
        /// `Body.Jupiter`, `Body.Saturn`, `Body.Uranus`, `Body.Neptune`, `Body.Pluto`.
        /// </param>
        /// <param name="time">The time at which to calculate the body's rotation axis.</param>
        /// <returns>North pole orientation and body spin angle.</returns>
        public static AxisInfo RotationAxis(Body body, AstroTime time)
        {
            double d = time.tt;
            double T = d / 36525.0;
            double ra, dec, w;

            switch (body)
            {
            case Body.Sun:
                ra = 286.13;
                dec = 63.87;
                w = 84.176 + (14.1844 * d);
                break;

            case Body.Mercury:
                ra = 281.0103 - (0.0328 * T);
                dec = 61.4155 - (0.0049 * T);
                w = (
                    329.5988
                    + (6.1385108 * d)
                    + (0.01067257 * Math.Sin(DEG2RAD*(174.7910857 + 4.092335*d)))
                    - (0.00112309 * Math.Sin(DEG2RAD*(349.5821714 + 8.184670*d)))
                    - (0.00011040 * Math.Sin(DEG2RAD*(164.3732571 + 12.277005*d)))
                    - (0.00002539 * Math.Sin(DEG2RAD*(339.1643429 + 16.369340*d)))
                    - (0.00000571 * Math.Sin(DEG2RAD*(153.9554286 + 20.461675*d)))
                );
                break;

            case Body.Venus:
                ra = 272.76;
                dec = 67.16;
                w = 160.20 - (1.4813688 * d);
                break;

            case Body.Earth:
                return EarthRotationAxis(time);

            case Body.Moon:
                // See page 8, Table 2 in:
                // https://astropedia.astrogeology.usgs.gov/alfresco/d/d/workspace/SpacesStore/28fd9e81-1964-44d6-a58b-fbbf61e64e15/WGCCRE2009reprint.pdf
                double E1  = DEG2RAD * (125.045 -  0.0529921*d);
                double E2  = DEG2RAD * (250.089 -  0.1059842*d);
                double E3  = DEG2RAD * (260.008 + 13.0120009*d);
                double E4  = DEG2RAD * (176.625 + 13.3407154*d);
                double E5  = DEG2RAD * (357.529 +  0.9856003*d);
                double E6  = DEG2RAD * (311.589 + 26.4057084*d);
                double E7  = DEG2RAD * (134.963 + 13.0649930*d);
                double E8  = DEG2RAD * (276.617 +  0.3287146*d);
                double E9  = DEG2RAD * (34.226  +  1.7484877*d);
                double E10 = DEG2RAD * (15.134  -  0.1589763*d);
                double E11 = DEG2RAD * (119.743 +  0.0036096*d);
                double E12 = DEG2RAD * (239.961 +  0.1643573*d);
                double E13 = DEG2RAD * (25.053  + 12.9590088*d);

                ra = (
                    269.9949 + 0.0031*T
                    - 3.8787*Math.Sin(E1)
                    - 0.1204*Math.Sin(E2)
                    + 0.0700*Math.Sin(E3)
                    - 0.0172*Math.Sin(E4)
                    + 0.0072*Math.Sin(E6)
                    - 0.0052*Math.Sin(E10)
                    + 0.0043*Math.Sin(E13)
                );

                dec = (
                    66.5392 + 0.0130*T
                    + 1.5419*Math.Cos(E1)
                    + 0.0239*Math.Cos(E2)
                    - 0.0278*Math.Cos(E3)
                    + 0.0068*Math.Cos(E4)
                    - 0.0029*Math.Cos(E6)
                    + 0.0009*Math.Cos(E7)
                    + 0.0008*Math.Cos(E10)
                    - 0.0009*Math.Cos(E13)
                );

                w = (
                    38.3213 + (13.17635815 - 1.4e-12*d)*d
                    + 3.5610*Math.Sin(E1)
                    + 0.1208*Math.Sin(E2)
                    - 0.0642*Math.Sin(E3)
                    + 0.0158*Math.Sin(E4)
                    + 0.0252*Math.Sin(E5)
                    - 0.0066*Math.Sin(E6)
                    - 0.0047*Math.Sin(E7)
                    - 0.0046*Math.Sin(E8)
                    + 0.0028*Math.Sin(E9)
                    + 0.0052*Math.Sin(E10)
                    + 0.0040*Math.Sin(E11)
                    + 0.0019*Math.Sin(E12)
                    - 0.0044*Math.Sin(E13)
                );
                break;

            case Body.Mars:
                ra = (
                    317.269202 - 0.10927547*T
                    + 0.000068 * Math.Sin(DEG2RAD*(198.991226 + 19139.4819985*T))
                    + 0.000238 * Math.Sin(DEG2RAD*(226.292679 + 38280.8511281*T))
                    + 0.000052 * Math.Sin(DEG2RAD*(249.663391 + 57420.7251593*T))
                    + 0.000009 * Math.Sin(DEG2RAD*(266.183510 + 76560.6367950*T))
                    + 0.419057 * Math.Sin(DEG2RAD*(79.398797 + 0.5042615*T))
                );

                dec = (
                    54.432516 - 0.05827105*T
                    + 0.000051*Math.Cos(DEG2RAD*(122.433576 + 19139.9407476*T))
                    + 0.000141*Math.Cos(DEG2RAD*(43.058401 + 38280.8753272*T))
                    + 0.000031*Math.Cos(DEG2RAD*(57.663379 + 57420.7517205*T))
                    + 0.000005*Math.Cos(DEG2RAD*(79.476401 + 76560.6495004*T))
                    + 1.591274*Math.Cos(DEG2RAD*(166.325722 + 0.5042615*T))
                );

                w = (
                    176.049863 + 350.891982443297*d
                    + 0.000145*Math.Sin(DEG2RAD*(129.071773 + 19140.0328244*T))
                    + 0.000157*Math.Sin(DEG2RAD*(36.352167 + 38281.0473591*T))
                    + 0.000040*Math.Sin(DEG2RAD*(56.668646 + 57420.9295360*T))
                    + 0.000001*Math.Sin(DEG2RAD*(67.364003 + 76560.2552215*T))
                    + 0.000001*Math.Sin(DEG2RAD*(104.792680 + 95700.4387578*T))
                    + 0.584542*Math.Sin(DEG2RAD*(95.391654 + 0.5042615*T))
                );
                break;

            case Body.Jupiter:
                double Ja = DEG2RAD*(99.360714 + 4850.4046*T);
                double Jb = DEG2RAD*(175.895369 + 1191.9605*T);
                double Jc = DEG2RAD*(300.323162 + 262.5475*T);
                double Jd = DEG2RAD*(114.012305 + 6070.2476*T);
                double Je = DEG2RAD*(49.511251 + 64.3000*T);

                ra = (
                    268.056595 - 0.006499*T
                    + 0.000117*Math.Sin(Ja)
                    + 0.000938*Math.Sin(Jb)
                    + 0.001432*Math.Sin(Jc)
                    + 0.000030*Math.Sin(Jd)
                    + 0.002150*Math.Sin(Je)
                );

                dec = (
                    64.495303 + 0.002413*T
                    + 0.000050*Math.Cos(Ja)
                    + 0.000404*Math.Cos(Jb)
                    + 0.000617*Math.Cos(Jc)
                    - 0.000013*Math.Cos(Jd)
                    + 0.000926*Math.Cos(Je)
                );

                w = 284.95 + 870.536*d;
                break;

            case Body.Saturn:
                ra = 40.589 - 0.036*T;
                dec = 83.537 - 0.004*T;
                w = 38.90 + 810.7939024*d;
                break;

            case Body.Uranus:
                ra = 257.311;
                dec = -15.175;
                w = 203.81 - 501.1600928*d;
                break;

            case Body.Neptune:
                double N = DEG2RAD*(357.85 + 52.316*T);
                ra = 299.36 + 0.70*Math.Sin(N);
                dec = 43.46 - 0.51*Math.Cos(N);
                w = 249.978 + 541.1397757*d - 0.48*Math.Sin(N);
                break;

            case Body.Pluto:
                ra = 132.993;
                dec = -6.163;
                w = 302.695 + 56.3625225*d;
                break;

            default:
                throw new InvalidBodyException(body);
            }

            AxisInfo axis;
            axis.ra = ra / 15.0;      // convert degrees to sidereal hours
            axis.dec = dec;
            axis.spin = w;

            // Calculate the north pole vector using the given angles.
            double radlat = dec * DEG2RAD;
            double radlon = ra * DEG2RAD;
            double rcoslat = Math.Cos(radlat);
            axis.north = new AstroVector(
                rcoslat * Math.Cos(radlon),
                rcoslat * Math.Sin(radlon),
                Math.Sin(radlat),
                time
            );

            return axis;
        }

        /// <summary>Calculates a rotation matrix from J2000 mean equator (EQJ) to J2000 mean ecliptic (ECL).</summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQJ = equatorial system, using equator at J2000 epoch.
        /// Target: ECL = ecliptic system, using equator at J2000 epoch.
        /// </remarks>
        /// <returns>A rotation matrix that converts EQJ to ECL.</returns>
        public static RotationMatrix Rotation_EQJ_ECL()
        {
            // ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
            const double c = 0.9174821430670688;    // cos(ob)
            const double s = 0.3977769691083922;    // sin(ob)
            var r = new RotationMatrix(new double[3,3]);

            r.rot[0, 0] = 1.0;  r.rot[1, 0] = 0.0;  r.rot[2, 0] = 0.0;
            r.rot[0, 1] = 0.0;  r.rot[1, 1] = +c;   r.rot[2, 1] = +s;
            r.rot[0, 2] = 0.0;  r.rot[1, 2] = -s;   r.rot[2, 2] = +c;

            return r;
        }


        /// <summary>Calculates a rotation matrix from J2000 mean ecliptic (ECL) to J2000 mean equator (EQJ).</summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: ECL = ecliptic system, using equator at J2000 epoch.
        /// Target: EQJ = equatorial system, using equator at J2000 epoch.
        /// </remarks>
        /// <returns>A rotation matrix that converts ECL to EQJ.</returns>
        public static RotationMatrix Rotation_ECL_EQJ()
        {
            // ob = mean obliquity of the J2000 ecliptic = 0.40909260059599012 radians.
            const double c = 0.9174821430670688;    // cos(ob)
            const double s = 0.3977769691083922;    // sin(ob)
            var r = new RotationMatrix(new double[3,3]);

            r.rot[0, 0] = 1.0;  r.rot[1, 0] = 0.0;  r.rot[2, 0] = 0.0;
            r.rot[0, 1] = 0.0;  r.rot[1, 1] = +c;   r.rot[2, 1] = -s;
            r.rot[0, 2] = 0.0;  r.rot[1, 2] = +s;   r.rot[2, 2] = +c;

            return r;
        }


        /// <summary>
        /// Calculates a rotation matrix from J2000 mean equator (EQJ) to equatorial of-date (EQD).
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
            RotationMatrix prec = precession_rot(time, PrecessDirection.From2000);
            RotationMatrix nut = nutation_rot(time, PrecessDirection.From2000);
            return CombineRotation(prec, nut);
        }


        /// <summary>
        /// Calculates a rotation matrix from J2000 mean equator (EQJ) to true ecliptic of date (ECT).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQJ = equatorial system, using equator at J2000 epoch.
        /// Target: ECT = ecliptic system, using true equinox of the specified date/time.
        /// </remarks>
        /// <param name="time">
        /// The date and time at which the Earth's equator defines the target orientation.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts EQJ to ECT at `time`.
        /// </returns>
        public static RotationMatrix Rotation_EQJ_ECT(AstroTime time)
        {
            RotationMatrix rot  = Rotation_EQJ_EQD(time);
            RotationMatrix step = Rotation_EQD_ECT(time);
            return CombineRotation(rot, step);
        }

        /// <summary>
        /// Calculates a rotation matrix from true ecliptic of date (ECT) to J2000 mean equator (EQJ).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: ECT = ecliptic system, using true equinox of the specified date/time.
        /// Target: EQJ = equatorial system, using equator at J2000 epoch.
        /// </remarks>
        /// <param name="time">
        /// The date and time at which the Earth's equator defines the target orientation.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts ECT to EQJ at `time`.
        /// </returns>
        public static RotationMatrix Rotation_ECT_EQJ(AstroTime time)
        {
            RotationMatrix rot  = Rotation_ECT_EQD(time);
            RotationMatrix step = Rotation_EQD_EQJ(time);
            return CombineRotation(rot, step);
        }

        /// <summary>
        /// Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean equator (EQJ).
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
            RotationMatrix nut = nutation_rot(time, PrecessDirection.Into2000);
            RotationMatrix prec = precession_rot(time, PrecessDirection.Into2000);
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

            var uze = new AstroVector(coslat * coslon, coslat * sinlon, sinlat, time);
            var une = new AstroVector(-sinlat * coslon, -sinlat * sinlon, coslat, time);
            var uwe = new AstroVector(sinlon, -coslon, 0.0, time);

            // Multiply sidereal hours by -15 to convert to degrees and flip eastward
            // rotation of the Earth to westward apparent movement of objects with time.
            double angle = -15.0 * SiderealTime(time);
            AstroVector uz = spin(angle, uze);
            AstroVector un = spin(angle, une);
            AstroVector uw = spin(angle, uwe);

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
        /// A rotation matrix that converts HOR to EQJ at `time` and for `observer`.
        /// </returns>
        public static RotationMatrix Rotation_HOR_EQJ(AstroTime time, Observer observer)
        {
            RotationMatrix hor_eqd = Rotation_HOR_EQD(time, observer);
            RotationMatrix eqd_eqj = Rotation_EQD_EQJ(time);
            return CombineRotation(hor_eqd, eqd_eqj);
        }


        /// <summary>
        /// Calculates a rotation matrix from J2000 mean equator (EQJ) to horizontal (HOR).
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
        /// Calculates a rotation matrix from equatorial of-date (EQD) to J2000 mean ecliptic (ECL).
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
        /// Calculates a rotation matrix from J2000 mean ecliptic (ECL) to equatorial of-date (EQD).
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
        /// Calculates a rotation matrix from J2000 mean ecliptic (ECL) to horizontal (HOR).
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
        /// Calculates a rotation matrix from horizontal (HOR) to J2000 mean ecliptic (ECL).
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

        /// <summary>
        /// Calculates a rotation matrix from J2000 mean equator (EQJ) to galactic (GAL).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQJ = equatorial system, using the equator at the J2000 epoch.
        /// Target: GAL = galactic system (IAU 1958 definition).
        /// </remarks>
        /// <returns>
        /// A rotation matrix that converts EQJ to GAL.
        /// </returns>
        public static RotationMatrix Rotation_EQJ_GAL()
        {
            var rot = new double[3, 3];

            // This rotation matrix was calculated by the following script
            // in this same source code repository:
            // demo/python/galeqj_matrix.py

            rot[0, 0] = -0.0548624779711344;
            rot[0, 1] = +0.4941095946388765;
            rot[0, 2] = -0.8676668813529025;

            rot[1, 0] = -0.8734572784246782;
            rot[1, 1] = -0.4447938112296831;
            rot[1, 2] = -0.1980677870294097;

            rot[2, 0] = -0.4838000529948520;
            rot[2, 1] = +0.7470034631630423;
            rot[2, 2] = +0.4559861124470794;

            return new RotationMatrix(rot);
        }

        /// <summary>
        /// Calculates a rotation matrix from galactic (GAL) to J2000 mean equator (EQJ).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: GAL = galactic system (IAU 1958 definition).
        /// Target: EQJ = equatorial system, using the equator at the J2000 epoch.
        /// </remarks>
        /// <returns>
        /// A rotation matrix that converts GAL to EQJ.
        /// </returns>
        public static RotationMatrix Rotation_GAL_EQJ()
        {
            var rot = new double[3, 3];

            // This rotation matrix was calculated by the following script
            // in this same source code repository:
            // demo/python/galeqj_matrix.py

            rot[0, 0] = -0.0548624779711344;
            rot[0, 1] = -0.8734572784246782;
            rot[0, 2] = -0.4838000529948520;

            rot[1, 0] = +0.4941095946388765;
            rot[1, 1] = -0.4447938112296831;
            rot[1, 2] = +0.7470034631630423;

            rot[2, 0] = -0.8676668813529025;
            rot[2, 1] = -0.1980677870294097;
            rot[2, 2] = +0.4559861124470794;

            return new RotationMatrix(rot);
        }

        /// <summary>
        /// Calculates a rotation matrix from true ecliptic of date (ECT) to equator of date (EQD).
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: ECT = true ecliptic of date.
        /// Target: EQD = equator of date.
        /// </remarks>
        /// <param name="time">
        /// The date and time of the ecliptic/equator conversion.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts ECT to EQD.
        /// </returns>
        public static RotationMatrix Rotation_ECT_EQD(AstroTime time)
        {
            var rot = new double[3, 3];

            // Find true ecliptic obliquity for this time.
            earth_tilt_t et = e_tilt(time);
            double tobl = et.tobl * DEG2RAD;
            double cos_tobl = Math.Cos(tobl);
            double sin_tobl = Math.Sin(tobl);

            // EQD.x = ECT.x
            rot[0, 0] = 1.0;
            rot[1, 0] = 0.0;
            rot[2, 0] = 0.0;

            // EQD.y = +cos*ECT.y - sin*ECT.z
            rot[0, 1] = 0.0;
            rot[1, 1] = +cos_tobl;
            rot[2, 1] = -sin_tobl;

            // EQD.z = +sin*ECT.y + cos*ECT.z
            rot[0, 2] = 0.0;
            rot[1, 2] = +sin_tobl;
            rot[2, 2] = +cos_tobl;

            return new RotationMatrix(rot);
        }

        /// <summary>
        /// Calculates a rotation matrix from equator of date (EQD) to true ecliptic of date (ECT) .
        /// </summary>
        /// <remarks>
        /// This is one of the family of functions that returns a rotation matrix
        /// for converting from one orientation to another.
        /// Source: EQD = equator of date.
        /// Target: ECT = true ecliptic of date.
        /// </remarks>
        /// <param name="time">
        /// The date and time of the equator/ecliptic conversion.
        /// </param>
        /// <returns>
        /// A rotation matrix that converts EQD to ECT.
        /// </returns>
        public static RotationMatrix Rotation_EQD_ECT(AstroTime time)
        {
            var rot = new double[3, 3];

            // Find true ecliptic obliquity for this time.
            earth_tilt_t et = e_tilt(time);
            double tobl = et.tobl * DEG2RAD;
            double cos_tobl = Math.Cos(tobl);
            double sin_tobl = Math.Sin(tobl);

            // ECT.x = EQD.x
            rot[0, 0] = 1.0;
            rot[1, 0] = 0.0;
            rot[2, 0] = 0.0;

            // ECT.y = +cos*EQJ.y + sin*EQJ.z
            rot[0, 1] = 0.0;
            rot[1, 1] = +cos_tobl;
            rot[2, 1] = +sin_tobl;

            // ECT.z = -sin*EQJ.y + cos*EQJ.z
            rot[0, 2] = 0.0;
            rot[1, 2] = -sin_tobl;
            rot[2, 2] = +cos_tobl;

            return new RotationMatrix(rot);
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
        /// The right ascension (RA) of a point in the sky, using the J2000 equatorial system (EQJ).
        /// </param>
        /// <param name="dec">
        /// The declination (DEC) of a point in the sky, using the J2000 equatorial system (EQJ).
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
            var sph2000 = new Spherical(dec, 15.0 * ra, 1.0);
            AstroVector vec2000 = VectorFromSphere(sph2000, Epoch2000);
            AstroVector vec1875 = RotateVector(ConstelRot, vec2000);
            Equatorial equ1875 = EquatorFromVector(vec1875);

            // Convert DEC from degrees and RA from hours, into compact angle units used in the _ConstelBounds table.
            double x_dec = 24.0 * equ1875.dec;
            double x_ra = (24.0 * 15.0) * equ1875.ra;

            // Search for the constellation using the B1875 coordinates.
            foreach (constel_boundary_t b in ConstelBounds)
                if ((b.dec_lo <= x_dec) && (b.ra_hi > x_ra) && (b.ra_lo <= x_ra))
                    return new ConstellationInfo(ConstelNames[b.index].symbol, ConstelNames[b.index].name, equ1875.ra, equ1875.dec);

            // This should never happen!
            throw new InternalError($"Unable to find constellation for coordinates: RA={ra}, DEC={dec}");
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
            new constel_boundary_t(83,      0,   8640,   2112)    // UMi
        ,   new constel_boundary_t(83,   2880,   5220,   2076)    // UMi
        ,   new constel_boundary_t(83,   7560,   8280,   2068)    // UMi
        ,   new constel_boundary_t(83,   6480,   7560,   2064)    // UMi
        ,   new constel_boundary_t(15,      0,   2880,   2040)    // Cep
        ,   new constel_boundary_t(10,   3300,   3840,   1968)    // Cam
        ,   new constel_boundary_t(15,      0,   1800,   1920)    // Cep
        ,   new constel_boundary_t(10,   3840,   5220,   1920)    // Cam
        ,   new constel_boundary_t(83,   6300,   6480,   1920)    // UMi
        ,   new constel_boundary_t(33,   7260,   7560,   1920)    // Dra
        ,   new constel_boundary_t(15,      0,   1263,   1848)    // Cep
        ,   new constel_boundary_t(10,   4140,   4890,   1848)    // Cam
        ,   new constel_boundary_t(83,   5952,   6300,   1800)    // UMi
        ,   new constel_boundary_t(15,   7260,   7440,   1800)    // Cep
        ,   new constel_boundary_t(10,   2868,   3300,   1764)    // Cam
        ,   new constel_boundary_t(33,   3300,   4080,   1764)    // Dra
        ,   new constel_boundary_t(83,   4680,   5952,   1680)    // UMi
        ,   new constel_boundary_t(13,   1116,   1230,   1632)    // Cas
        ,   new constel_boundary_t(33,   7350,   7440,   1608)    // Dra
        ,   new constel_boundary_t(33,   4080,   4320,   1596)    // Dra
        ,   new constel_boundary_t(15,      0,    120,   1584)    // Cep
        ,   new constel_boundary_t(83,   5040,   5640,   1584)    // UMi
        ,   new constel_boundary_t(15,   8490,   8640,   1584)    // Cep
        ,   new constel_boundary_t(33,   4320,   4860,   1536)    // Dra
        ,   new constel_boundary_t(33,   4860,   5190,   1512)    // Dra
        ,   new constel_boundary_t(15,   8340,   8490,   1512)    // Cep
        ,   new constel_boundary_t(10,   2196,   2520,   1488)    // Cam
        ,   new constel_boundary_t(33,   7200,   7350,   1476)    // Dra
        ,   new constel_boundary_t(15, 7393.2,   7416,   1462)    // Cep
        ,   new constel_boundary_t(10,   2520,   2868,   1440)    // Cam
        ,   new constel_boundary_t(82,   2868,   3030,   1440)    // UMa
        ,   new constel_boundary_t(33,   7116,   7200,   1428)    // Dra
        ,   new constel_boundary_t(15,   7200, 7393.2,   1428)    // Cep
        ,   new constel_boundary_t(15,   8232,   8340,   1418)    // Cep
        ,   new constel_boundary_t(13,      0,    876,   1404)    // Cas
        ,   new constel_boundary_t(33,   6990,   7116,   1392)    // Dra
        ,   new constel_boundary_t(13,    612,    687,   1380)    // Cas
        ,   new constel_boundary_t(13,    876,   1116,   1368)    // Cas
        ,   new constel_boundary_t(10,   1116,   1140,   1368)    // Cam
        ,   new constel_boundary_t(15,   8034,   8232,   1350)    // Cep
        ,   new constel_boundary_t(10,   1800,   2196,   1344)    // Cam
        ,   new constel_boundary_t(82,   5052,   5190,   1332)    // UMa
        ,   new constel_boundary_t(33,   5190,   6990,   1332)    // Dra
        ,   new constel_boundary_t(10,   1140,   1200,   1320)    // Cam
        ,   new constel_boundary_t(15,   7968,   8034,   1320)    // Cep
        ,   new constel_boundary_t(15,   7416,   7908,   1316)    // Cep
        ,   new constel_boundary_t(13,      0,    612,   1296)    // Cas
        ,   new constel_boundary_t(50,   2196,   2340,   1296)    // Lyn
        ,   new constel_boundary_t(82,   4350,   4860,   1272)    // UMa
        ,   new constel_boundary_t(33,   5490,   5670,   1272)    // Dra
        ,   new constel_boundary_t(15,   7908,   7968,   1266)    // Cep
        ,   new constel_boundary_t(10,   1200,   1800,   1260)    // Cam
        ,   new constel_boundary_t(13,   8232,   8400,   1260)    // Cas
        ,   new constel_boundary_t(33,   5670,   6120,   1236)    // Dra
        ,   new constel_boundary_t(62,    735,    906,   1212)    // Per
        ,   new constel_boundary_t(33,   6120,   6564,   1212)    // Dra
        ,   new constel_boundary_t(13,      0,    492,   1200)    // Cas
        ,   new constel_boundary_t(62,    492,    600,   1200)    // Per
        ,   new constel_boundary_t(50,   2340,   2448,   1200)    // Lyn
        ,   new constel_boundary_t(13,   8400,   8640,   1200)    // Cas
        ,   new constel_boundary_t(82,   4860,   5052,   1164)    // UMa
        ,   new constel_boundary_t(13,      0,    402,   1152)    // Cas
        ,   new constel_boundary_t(13,   8490,   8640,   1152)    // Cas
        ,   new constel_boundary_t(39,   6543,   6564,   1140)    // Her
        ,   new constel_boundary_t(33,   6564,   6870,   1140)    // Dra
        ,   new constel_boundary_t(30,   6870,   6900,   1140)    // Cyg
        ,   new constel_boundary_t(62,    600,    735,   1128)    // Per
        ,   new constel_boundary_t(82,   3030,   3300,   1128)    // UMa
        ,   new constel_boundary_t(13,     60,    312,   1104)    // Cas
        ,   new constel_boundary_t(82,   4320,   4350,   1080)    // UMa
        ,   new constel_boundary_t(50,   2448,   2652,   1068)    // Lyn
        ,   new constel_boundary_t(30,   7887,   7908,   1056)    // Cyg
        ,   new constel_boundary_t(30,   7875,   7887,   1050)    // Cyg
        ,   new constel_boundary_t(30,   6900,   6984,   1044)    // Cyg
        ,   new constel_boundary_t(82,   3300,   3660,   1008)    // UMa
        ,   new constel_boundary_t(82,   3660,   3882,    960)    // UMa
        ,   new constel_boundary_t( 8,   5556,   5670,    960)    // Boo
        ,   new constel_boundary_t(39,   5670,   5880,    960)    // Her
        ,   new constel_boundary_t(50,   3330,   3450,    954)    // Lyn
        ,   new constel_boundary_t( 0,      0,    906,    882)    // And
        ,   new constel_boundary_t(62,    906,    924,    882)    // Per
        ,   new constel_boundary_t(51,   6969,   6984,    876)    // Lyr
        ,   new constel_boundary_t(62,   1620,   1689,    864)    // Per
        ,   new constel_boundary_t(30,   7824,   7875,    864)    // Cyg
        ,   new constel_boundary_t(44,   7875,   7920,    864)    // Lac
        ,   new constel_boundary_t( 7,   2352,   2652,    852)    // Aur
        ,   new constel_boundary_t(50,   2652,   2790,    852)    // Lyn
        ,   new constel_boundary_t( 0,      0,    720,    840)    // And
        ,   new constel_boundary_t(44,   7920,   8214,    840)    // Lac
        ,   new constel_boundary_t(44,   8214,   8232,    828)    // Lac
        ,   new constel_boundary_t( 0,   8232,   8460,    828)    // And
        ,   new constel_boundary_t(62,    924,    978,    816)    // Per
        ,   new constel_boundary_t(82,   3882,   3960,    816)    // UMa
        ,   new constel_boundary_t(29,   4320,   4440,    816)    // CVn
        ,   new constel_boundary_t(50,   2790,   3330,    804)    // Lyn
        ,   new constel_boundary_t(48,   3330,   3558,    804)    // LMi
        ,   new constel_boundary_t( 0,    258,    507,    792)    // And
        ,   new constel_boundary_t( 8,   5466,   5556,    792)    // Boo
        ,   new constel_boundary_t( 0,   8460,   8550,    770)    // And
        ,   new constel_boundary_t(29,   4440,   4770,    768)    // CVn
        ,   new constel_boundary_t( 0,   8550,   8640,    752)    // And
        ,   new constel_boundary_t(29,   5025,   5052,    738)    // CVn
        ,   new constel_boundary_t(80,    870,    978,    736)    // Tri
        ,   new constel_boundary_t(62,    978,   1620,    736)    // Per
        ,   new constel_boundary_t( 7,   1620,   1710,    720)    // Aur
        ,   new constel_boundary_t(51,   6543,   6969,    720)    // Lyr
        ,   new constel_boundary_t(82,   3960,   4320,    696)    // UMa
        ,   new constel_boundary_t(30,   7080,   7530,    696)    // Cyg
        ,   new constel_boundary_t( 7,   1710,   2118,    684)    // Aur
        ,   new constel_boundary_t(48,   3558,   3780,    684)    // LMi
        ,   new constel_boundary_t(29,   4770,   5025,    684)    // CVn
        ,   new constel_boundary_t( 0,      0,     24,    672)    // And
        ,   new constel_boundary_t(80,    507,    600,    672)    // Tri
        ,   new constel_boundary_t( 7,   2118,   2352,    672)    // Aur
        ,   new constel_boundary_t(37,   2838,   2880,    672)    // Gem
        ,   new constel_boundary_t(30,   7530,   7824,    672)    // Cyg
        ,   new constel_boundary_t(30,   6933,   7080,    660)    // Cyg
        ,   new constel_boundary_t(80,    690,    870,    654)    // Tri
        ,   new constel_boundary_t(25,   5820,   5880,    648)    // CrB
        ,   new constel_boundary_t( 8,   5430,   5466,    624)    // Boo
        ,   new constel_boundary_t(25,   5466,   5820,    624)    // CrB
        ,   new constel_boundary_t(51,   6612,   6792,    624)    // Lyr
        ,   new constel_boundary_t(48,   3870,   3960,    612)    // LMi
        ,   new constel_boundary_t(51,   6792,   6933,    612)    // Lyr
        ,   new constel_boundary_t(80,    600,    690,    600)    // Tri
        ,   new constel_boundary_t(66,    258,    306,    570)    // Psc
        ,   new constel_boundary_t(48,   3780,   3870,    564)    // LMi
        ,   new constel_boundary_t(87,   7650,   7710,    564)    // Vul
        ,   new constel_boundary_t(77,   2052,   2118,    548)    // Tau
        ,   new constel_boundary_t( 0,     24,     51,    528)    // And
        ,   new constel_boundary_t(73,   5730,   5772,    528)    // Ser
        ,   new constel_boundary_t(37,   2118,   2238,    516)    // Gem
        ,   new constel_boundary_t(87,   7140,   7290,    510)    // Vul
        ,   new constel_boundary_t(87,   6792,   6930,    506)    // Vul
        ,   new constel_boundary_t( 0,     51,    306,    504)    // And
        ,   new constel_boundary_t(87,   7290,   7404,    492)    // Vul
        ,   new constel_boundary_t(37,   2811,   2838,    480)    // Gem
        ,   new constel_boundary_t(87,   7404,   7650,    468)    // Vul
        ,   new constel_boundary_t(87,   6930,   7140,    460)    // Vul
        ,   new constel_boundary_t( 6,   1182,   1212,    456)    // Ari
        ,   new constel_boundary_t(75,   6792,   6840,    444)    // Sge
        ,   new constel_boundary_t(59,   2052,   2076,    432)    // Ori
        ,   new constel_boundary_t(37,   2238,   2271,    420)    // Gem
        ,   new constel_boundary_t(75,   6840,   7140,    388)    // Sge
        ,   new constel_boundary_t(77,   1788,   1920,    384)    // Tau
        ,   new constel_boundary_t(39,   5730,   5790,    384)    // Her
        ,   new constel_boundary_t(75,   7140,   7290,    378)    // Sge
        ,   new constel_boundary_t(77,   1662,   1788,    372)    // Tau
        ,   new constel_boundary_t(77,   1920,   2016,    372)    // Tau
        ,   new constel_boundary_t(23,   4620,   4860,    360)    // Com
        ,   new constel_boundary_t(39,   6210,   6570,    344)    // Her
        ,   new constel_boundary_t(23,   4272,   4620,    336)    // Com
        ,   new constel_boundary_t(37,   2700,   2811,    324)    // Gem
        ,   new constel_boundary_t(39,   6030,   6210,    308)    // Her
        ,   new constel_boundary_t(61,      0,     51,    300)    // Peg
        ,   new constel_boundary_t(77,   2016,   2076,    300)    // Tau
        ,   new constel_boundary_t(37,   2520,   2700,    300)    // Gem
        ,   new constel_boundary_t(61,   7602,   7680,    300)    // Peg
        ,   new constel_boundary_t(37,   2271,   2496,    288)    // Gem
        ,   new constel_boundary_t(39,   6570,   6792,    288)    // Her
        ,   new constel_boundary_t(31,   7515,   7578,    284)    // Del
        ,   new constel_boundary_t(61,   7578,   7602,    284)    // Peg
        ,   new constel_boundary_t(45,   4146,   4272,    264)    // Leo
        ,   new constel_boundary_t(59,   2247,   2271,    240)    // Ori
        ,   new constel_boundary_t(37,   2496,   2520,    240)    // Gem
        ,   new constel_boundary_t(21,   2811,   2853,    240)    // Cnc
        ,   new constel_boundary_t(61,   8580,   8640,    240)    // Peg
        ,   new constel_boundary_t( 6,    600,   1182,    238)    // Ari
        ,   new constel_boundary_t(31,   7251,   7308,    204)    // Del
        ,   new constel_boundary_t( 8,   4860,   5430,    192)    // Boo
        ,   new constel_boundary_t(61,   8190,   8580,    180)    // Peg
        ,   new constel_boundary_t(21,   2853,   3330,    168)    // Cnc
        ,   new constel_boundary_t(45,   3330,   3870,    168)    // Leo
        ,   new constel_boundary_t(58,   6570, 6718.4,    150)    // Oph
        ,   new constel_boundary_t( 3, 6718.4,   6792,    150)    // Aql
        ,   new constel_boundary_t(31,   7500,   7515,    144)    // Del
        ,   new constel_boundary_t(20,   2520,   2526,    132)    // CMi
        ,   new constel_boundary_t(73,   6570,   6633,    108)    // Ser
        ,   new constel_boundary_t(39,   5790,   6030,     96)    // Her
        ,   new constel_boundary_t(58,   6570,   6633,     72)    // Oph
        ,   new constel_boundary_t(61,   7728,   7800,     66)    // Peg
        ,   new constel_boundary_t(66,      0,    720,     48)    // Psc
        ,   new constel_boundary_t(73,   6690,   6792,     48)    // Ser
        ,   new constel_boundary_t(31,   7308,   7500,     48)    // Del
        ,   new constel_boundary_t(34,   7500,   7680,     48)    // Equ
        ,   new constel_boundary_t(61,   7680,   7728,     48)    // Peg
        ,   new constel_boundary_t(61,   7920,   8190,     48)    // Peg
        ,   new constel_boundary_t(61,   7800,   7920,     42)    // Peg
        ,   new constel_boundary_t(20,   2526,   2592,     36)    // CMi
        ,   new constel_boundary_t(77,   1290,   1662,      0)    // Tau
        ,   new constel_boundary_t(59,   1662,   1680,      0)    // Ori
        ,   new constel_boundary_t(20,   2592,   2910,      0)    // CMi
        ,   new constel_boundary_t(85,   5280,   5430,      0)    // Vir
        ,   new constel_boundary_t(58,   6420,   6570,      0)    // Oph
        ,   new constel_boundary_t(16,    954,   1182,    -42)    // Cet
        ,   new constel_boundary_t(77,   1182,   1290,    -42)    // Tau
        ,   new constel_boundary_t(73,   5430,   5856,    -78)    // Ser
        ,   new constel_boundary_t(59,   1680,   1830,    -96)    // Ori
        ,   new constel_boundary_t(59,   2100,   2247,    -96)    // Ori
        ,   new constel_boundary_t(73,   6420,   6468,    -96)    // Ser
        ,   new constel_boundary_t(73,   6570,   6690,    -96)    // Ser
        ,   new constel_boundary_t( 3,   6690,   6792,    -96)    // Aql
        ,   new constel_boundary_t(66,   8190,   8580,    -96)    // Psc
        ,   new constel_boundary_t(45,   3870,   4146,   -144)    // Leo
        ,   new constel_boundary_t(85,   4146,   4260,   -144)    // Vir
        ,   new constel_boundary_t(66,      0,    120,   -168)    // Psc
        ,   new constel_boundary_t(66,   8580,   8640,   -168)    // Psc
        ,   new constel_boundary_t(85,   5130,   5280,   -192)    // Vir
        ,   new constel_boundary_t(58,   5730,   5856,   -192)    // Oph
        ,   new constel_boundary_t( 3,   7200,   7392,   -216)    // Aql
        ,   new constel_boundary_t( 4,   7680,   7872,   -216)    // Aqr
        ,   new constel_boundary_t(58,   6180,   6468,   -240)    // Oph
        ,   new constel_boundary_t(54,   2100,   2910,   -264)    // Mon
        ,   new constel_boundary_t(35,   1770,   1830,   -264)    // Eri
        ,   new constel_boundary_t(59,   1830,   2100,   -264)    // Ori
        ,   new constel_boundary_t(41,   2910,   3012,   -264)    // Hya
        ,   new constel_boundary_t(74,   3450,   3870,   -264)    // Sex
        ,   new constel_boundary_t(85,   4260,   4620,   -264)    // Vir
        ,   new constel_boundary_t(58,   6330,   6360,   -280)    // Oph
        ,   new constel_boundary_t( 3,   6792,   7200, -288.8)    // Aql
        ,   new constel_boundary_t(35,   1740,   1770,   -348)    // Eri
        ,   new constel_boundary_t( 4,   7392,   7680,   -360)    // Aqr
        ,   new constel_boundary_t(73,   6180,   6570,   -384)    // Ser
        ,   new constel_boundary_t(72,   6570,   6792,   -384)    // Sct
        ,   new constel_boundary_t(41,   3012,   3090,   -408)    // Hya
        ,   new constel_boundary_t(58,   5856,   5895,   -438)    // Oph
        ,   new constel_boundary_t(41,   3090,   3270,   -456)    // Hya
        ,   new constel_boundary_t(26,   3870,   3900,   -456)    // Crt
        ,   new constel_boundary_t(71,   5856,   5895,   -462)    // Sco
        ,   new constel_boundary_t(47,   5640,   5730,   -480)    // Lib
        ,   new constel_boundary_t(28,   4530,   4620,   -528)    // Crv
        ,   new constel_boundary_t(85,   4620,   5130,   -528)    // Vir
        ,   new constel_boundary_t(41,   3270,   3510,   -576)    // Hya
        ,   new constel_boundary_t(16,    600,    954, -585.2)    // Cet
        ,   new constel_boundary_t(35,    954,   1350, -585.2)    // Eri
        ,   new constel_boundary_t(26,   3900,   4260,   -588)    // Crt
        ,   new constel_boundary_t(28,   4260,   4530,   -588)    // Crv
        ,   new constel_boundary_t(47,   5130,   5370,   -588)    // Lib
        ,   new constel_boundary_t(58,   5856,   6030,   -590)    // Oph
        ,   new constel_boundary_t(16,      0,    600,   -612)    // Cet
        ,   new constel_boundary_t(11,   7680,   7872,   -612)    // Cap
        ,   new constel_boundary_t( 4,   7872,   8580,   -612)    // Aqr
        ,   new constel_boundary_t(16,   8580,   8640,   -612)    // Cet
        ,   new constel_boundary_t(41,   3510,   3690,   -636)    // Hya
        ,   new constel_boundary_t(35,   1692,   1740,   -654)    // Eri
        ,   new constel_boundary_t(46,   1740,   2202,   -654)    // Lep
        ,   new constel_boundary_t(11,   7200,   7680,   -672)    // Cap
        ,   new constel_boundary_t(41,   3690,   3810,   -700)    // Hya
        ,   new constel_boundary_t(41,   4530,   5370,   -708)    // Hya
        ,   new constel_boundary_t(47,   5370,   5640,   -708)    // Lib
        ,   new constel_boundary_t(71,   5640,   5760,   -708)    // Sco
        ,   new constel_boundary_t(35,   1650,   1692,   -720)    // Eri
        ,   new constel_boundary_t(58,   6030,   6336,   -720)    // Oph
        ,   new constel_boundary_t(76,   6336,   6420,   -720)    // Sgr
        ,   new constel_boundary_t(41,   3810,   3900,   -748)    // Hya
        ,   new constel_boundary_t(19,   2202,   2652,   -792)    // CMa
        ,   new constel_boundary_t(41,   4410,   4530,   -792)    // Hya
        ,   new constel_boundary_t(41,   3900,   4410,   -840)    // Hya
        ,   new constel_boundary_t(36,   1260,   1350,   -864)    // For
        ,   new constel_boundary_t(68,   3012,   3372,   -882)    // Pyx
        ,   new constel_boundary_t(35,   1536,   1650,   -888)    // Eri
        ,   new constel_boundary_t(76,   6420,   6900,   -888)    // Sgr
        ,   new constel_boundary_t(65,   7680,   8280,   -888)    // PsA
        ,   new constel_boundary_t(70,   8280,   8400,   -888)    // Scl
        ,   new constel_boundary_t(36,   1080,   1260,   -950)    // For
        ,   new constel_boundary_t( 1,   3372,   3960,   -954)    // Ant
        ,   new constel_boundary_t(70,      0,    600,   -960)    // Scl
        ,   new constel_boundary_t(36,    600,   1080,   -960)    // For
        ,   new constel_boundary_t(35,   1392,   1536,   -960)    // Eri
        ,   new constel_boundary_t(70,   8400,   8640,   -960)    // Scl
        ,   new constel_boundary_t(14,   5100,   5370,  -1008)    // Cen
        ,   new constel_boundary_t(49,   5640,   5760,  -1008)    // Lup
        ,   new constel_boundary_t(71,   5760, 5911.5,  -1008)    // Sco
        ,   new constel_boundary_t( 9,   1740,   1800,  -1032)    // Cae
        ,   new constel_boundary_t(22,   1800,   2370,  -1032)    // Col
        ,   new constel_boundary_t(67,   2880,   3012,  -1032)    // Pup
        ,   new constel_boundary_t(35,   1230,   1392,  -1056)    // Eri
        ,   new constel_boundary_t(71, 5911.5,   6420,  -1092)    // Sco
        ,   new constel_boundary_t(24,   6420,   6900,  -1092)    // CrA
        ,   new constel_boundary_t(76,   6900,   7320,  -1092)    // Sgr
        ,   new constel_boundary_t(53,   7320,   7680,  -1092)    // Mic
        ,   new constel_boundary_t(35,   1080,   1230,  -1104)    // Eri
        ,   new constel_boundary_t( 9,   1620,   1740,  -1116)    // Cae
        ,   new constel_boundary_t(49,   5520,   5640,  -1152)    // Lup
        ,   new constel_boundary_t(63,      0,    840,  -1156)    // Phe
        ,   new constel_boundary_t(35,    960,   1080,  -1176)    // Eri
        ,   new constel_boundary_t(40,   1470,   1536,  -1176)    // Hor
        ,   new constel_boundary_t( 9,   1536,   1620,  -1176)    // Cae
        ,   new constel_boundary_t(38,   7680,   7920,  -1200)    // Gru
        ,   new constel_boundary_t(67,   2160,   2880,  -1218)    // Pup
        ,   new constel_boundary_t(84,   2880,   2940,  -1218)    // Vel
        ,   new constel_boundary_t(35,    870,    960,  -1224)    // Eri
        ,   new constel_boundary_t(40,   1380,   1470,  -1224)    // Hor
        ,   new constel_boundary_t(63,      0,    660,  -1236)    // Phe
        ,   new constel_boundary_t(12,   2160,   2220,  -1260)    // Car
        ,   new constel_boundary_t(84,   2940,   3042,  -1272)    // Vel
        ,   new constel_boundary_t(40,   1260,   1380,  -1276)    // Hor
        ,   new constel_boundary_t(32,   1380,   1440,  -1276)    // Dor
        ,   new constel_boundary_t(63,      0,    570,  -1284)    // Phe
        ,   new constel_boundary_t(35,    780,    870,  -1296)    // Eri
        ,   new constel_boundary_t(64,   1620,   1800,  -1296)    // Pic
        ,   new constel_boundary_t(49,   5418,   5520,  -1296)    // Lup
        ,   new constel_boundary_t(84,   3042,   3180,  -1308)    // Vel
        ,   new constel_boundary_t(12,   2220,   2340,  -1320)    // Car
        ,   new constel_boundary_t(14,   4260,   4620,  -1320)    // Cen
        ,   new constel_boundary_t(49,   5100,   5418,  -1320)    // Lup
        ,   new constel_boundary_t(56,   5418,   5520,  -1320)    // Nor
        ,   new constel_boundary_t(32,   1440,   1560,  -1356)    // Dor
        ,   new constel_boundary_t(84,   3180,   3960,  -1356)    // Vel
        ,   new constel_boundary_t(14,   3960,   4050,  -1356)    // Cen
        ,   new constel_boundary_t( 5,   6300,   6480,  -1368)    // Ara
        ,   new constel_boundary_t(78,   6480,   7320,  -1368)    // Tel
        ,   new constel_boundary_t(38,   7920,   8400,  -1368)    // Gru
        ,   new constel_boundary_t(40,   1152,   1260,  -1380)    // Hor
        ,   new constel_boundary_t(64,   1800,   1980,  -1380)    // Pic
        ,   new constel_boundary_t(12,   2340,   2460,  -1392)    // Car
        ,   new constel_boundary_t(63,      0,    480,  -1404)    // Phe
        ,   new constel_boundary_t(35,    480,    780,  -1404)    // Eri
        ,   new constel_boundary_t(63,   8400,   8640,  -1404)    // Phe
        ,   new constel_boundary_t(32,   1560,   1650,  -1416)    // Dor
        ,   new constel_boundary_t(56,   5520, 5911.5,  -1440)    // Nor
        ,   new constel_boundary_t(43,   7320,   7680,  -1440)    // Ind
        ,   new constel_boundary_t(64,   1980,   2160,  -1464)    // Pic
        ,   new constel_boundary_t(18,   5460,   5520,  -1464)    // Cir
        ,   new constel_boundary_t( 5, 5911.5,   5970,  -1464)    // Ara
        ,   new constel_boundary_t(18,   5370,   5460,  -1526)    // Cir
        ,   new constel_boundary_t( 5,   5970,   6030,  -1526)    // Ara
        ,   new constel_boundary_t(64,   2160,   2460,  -1536)    // Pic
        ,   new constel_boundary_t(12,   2460,   3252,  -1536)    // Car
        ,   new constel_boundary_t(14,   4050,   4260,  -1536)    // Cen
        ,   new constel_boundary_t(27,   4260,   4620,  -1536)    // Cru
        ,   new constel_boundary_t(14,   4620,   5232,  -1536)    // Cen
        ,   new constel_boundary_t(18,   4860,   4920,  -1560)    // Cir
        ,   new constel_boundary_t( 5,   6030,   6060,  -1560)    // Ara
        ,   new constel_boundary_t(40,    780,   1152,  -1620)    // Hor
        ,   new constel_boundary_t(69,   1152,   1650,  -1620)    // Ret
        ,   new constel_boundary_t(18,   5310,   5370,  -1620)    // Cir
        ,   new constel_boundary_t( 5,   6060,   6300,  -1620)    // Ara
        ,   new constel_boundary_t(60,   6300,   6480,  -1620)    // Pav
        ,   new constel_boundary_t(81,   7920,   8400,  -1620)    // Tuc
        ,   new constel_boundary_t(32,   1650,   2370,  -1680)    // Dor
        ,   new constel_boundary_t(18,   4920,   5310,  -1680)    // Cir
        ,   new constel_boundary_t(79,   5310,   6120,  -1680)    // TrA
        ,   new constel_boundary_t(81,      0,    480,  -1800)    // Tuc
        ,   new constel_boundary_t(42,   1260,   1650,  -1800)    // Hyi
        ,   new constel_boundary_t(86,   2370,   3252,  -1800)    // Vol
        ,   new constel_boundary_t(12,   3252,   4050,  -1800)    // Car
        ,   new constel_boundary_t(55,   4050,   4920,  -1800)    // Mus
        ,   new constel_boundary_t(60,   6480,   7680,  -1800)    // Pav
        ,   new constel_boundary_t(43,   7680,   8400,  -1800)    // Ind
        ,   new constel_boundary_t(81,   8400,   8640,  -1800)    // Tuc
        ,   new constel_boundary_t(81,    270,    480,  -1824)    // Tuc
        ,   new constel_boundary_t(42,      0,   1260,  -1980)    // Hyi
        ,   new constel_boundary_t(17,   2760,   4920,  -1980)    // Cha
        ,   new constel_boundary_t( 2,   4920,   6480,  -1980)    // Aps
        ,   new constel_boundary_t(52,   1260,   2760,  -2040)    // Men
        ,   new constel_boundary_t(57,      0,   8640,  -2160)    // Oct
        };



    }
}
