/*
    Astronomy Engine for C# / .NET.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019 Don Cross <cosinekitty@gmail.com>

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
            : this((d - Origin).TotalDays)
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
            return ToUtcDateTime().ToString("yyyy-MM-ddThh:mm:ss.fffZ");
        }
    }

    /// <summary>
    /// A 3D Cartesian vector whose components are expressed in Astronomical Units (AU).
    /// </summary>
    public class AstroVector
    {
        /// <summary>
        /// The Cartesian x-coordinate of the vector in AU.
        /// </summary>
        public readonly double x;

        /// <summary>
        /// The Cartesian y-coordinate of the vector in AU.
        /// </summary>
        public readonly double y;

        /// <summary>
        /// The Cartesian z-coordinate of the vector in AU.
        /// </summary>
        public readonly double z;

        /// <summary>
        /// The date and time at which this vector is valid.
        /// </summary>
        public readonly AstroTime t;

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
    }

    /// <summary>
    /// The location of an observer on (or near) the surface of the Earth.
    /// </summary>
    /// <remarks>
    /// This structure is passed to functions that calculate phenomena as observed
    /// from a particular place on the Earth.
    /// </remarks>
    public class Observer
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
    /// Equatorial angular coordinates.
    /// </summary>
    /// <remarks>
    /// Coordinates of a celestial body as seen from the Earth
    /// (geocentric or topocentric, depending on context),
    /// oriented with respect to the projection of the Earth's equator onto the sky.
    /// </remarks>
    public class Equatorial
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
    public class Ecliptic
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
    }

    /// <summary>
    /// The wrapper class that holds Astronomy Engine functions.
    /// </summary>
    public static class Astronomy
    {
        private struct deltat_entry_t
        {
            public double mjd;
            public double dt;
        }

        private static readonly deltat_entry_t[] DT = new deltat_entry_t[] {
new deltat_entry_t { mjd=-72638.0, dt=38 },
new deltat_entry_t { mjd=-65333.0, dt=26 },
new deltat_entry_t { mjd=-58028.0, dt=21 },
new deltat_entry_t { mjd=-50724.0, dt=21.1 },
new deltat_entry_t { mjd=-43419.0, dt=13.5 },
new deltat_entry_t { mjd=-39766.0, dt=13.7 },
new deltat_entry_t { mjd=-36114.0, dt=14.8 },
new deltat_entry_t { mjd=-32461.0, dt=15.7 },
new deltat_entry_t { mjd=-28809.0, dt=15.6 },
new deltat_entry_t { mjd=-25156.0, dt=13.3 },
new deltat_entry_t { mjd=-21504.0, dt=12.6 },
new deltat_entry_t { mjd=-17852.0, dt=11.2 },
new deltat_entry_t { mjd=-14200.0, dt=11.13 },
new deltat_entry_t { mjd=-10547.0, dt=7.95 },
new deltat_entry_t { mjd=-6895.0, dt=6.22 },
new deltat_entry_t { mjd=-3242.0, dt=6.55 },
new deltat_entry_t { mjd=-1416.0, dt=7.26 },
new deltat_entry_t { mjd=410.0, dt=7.35 },
new deltat_entry_t { mjd=2237.0, dt=5.92 },
new deltat_entry_t { mjd=4063.0, dt=1.04 },
new deltat_entry_t { mjd=5889.0, dt=-3.19 },
new deltat_entry_t { mjd=7715.0, dt=-5.36 },
new deltat_entry_t { mjd=9542.0, dt=-5.74 },
new deltat_entry_t { mjd=11368.0, dt=-5.86 },
new deltat_entry_t { mjd=13194.0, dt=-6.41 },
new deltat_entry_t { mjd=15020.0, dt=-2.70 },
new deltat_entry_t { mjd=16846.0, dt=3.92 },
new deltat_entry_t { mjd=18672.0, dt=10.38 },
new deltat_entry_t { mjd=20498.0, dt=17.19 },
new deltat_entry_t { mjd=22324.0, dt=21.41 },
new deltat_entry_t { mjd=24151.0, dt=23.63 },
new deltat_entry_t { mjd=25977.0, dt=24.02 },
new deltat_entry_t { mjd=27803.0, dt=23.91 },
new deltat_entry_t { mjd=29629.0, dt=24.35 },
new deltat_entry_t { mjd=31456.0, dt=26.76 },
new deltat_entry_t { mjd=33282.0, dt=29.15 },
new deltat_entry_t { mjd=35108.0, dt=31.07 },
new deltat_entry_t { mjd=36934.0, dt=33.150 },
new deltat_entry_t { mjd=38761.0, dt=35.738 },
new deltat_entry_t { mjd=40587.0, dt=40.182 },
new deltat_entry_t { mjd=42413.0, dt=45.477 },
new deltat_entry_t { mjd=44239.0, dt=50.540 },
new deltat_entry_t { mjd=44605.0, dt=51.3808 },
new deltat_entry_t { mjd=44970.0, dt=52.1668 },
new deltat_entry_t { mjd=45335.0, dt=52.9565 },
new deltat_entry_t { mjd=45700.0, dt=53.7882 },
new deltat_entry_t { mjd=46066.0, dt=54.3427 },
new deltat_entry_t { mjd=46431.0, dt=54.8712 },
new deltat_entry_t { mjd=46796.0, dt=55.3222 },
new deltat_entry_t { mjd=47161.0, dt=55.8197 },
new deltat_entry_t { mjd=47527.0, dt=56.3000 },
new deltat_entry_t { mjd=47892.0, dt=56.8553 },
new deltat_entry_t { mjd=48257.0, dt=57.5653 },
new deltat_entry_t { mjd=48622.0, dt=58.3092 },
new deltat_entry_t { mjd=48988.0, dt=59.1218 },
new deltat_entry_t { mjd=49353.0, dt=59.9845 },
new deltat_entry_t { mjd=49718.0, dt=60.7853 },
new deltat_entry_t { mjd=50083.0, dt=61.6287 },
new deltat_entry_t { mjd=50449.0, dt=62.2950 },
new deltat_entry_t { mjd=50814.0, dt=62.9659 },
new deltat_entry_t { mjd=51179.0, dt=63.4673 },
new deltat_entry_t { mjd=51544.0, dt=63.8285 },
new deltat_entry_t { mjd=51910.0, dt=64.0908 },
new deltat_entry_t { mjd=52275.0, dt=64.2998 },
new deltat_entry_t { mjd=52640.0, dt=64.4734 },
new deltat_entry_t { mjd=53005.0, dt=64.5736 },
new deltat_entry_t { mjd=53371.0, dt=64.6876 },
new deltat_entry_t { mjd=53736.0, dt=64.8452 },
new deltat_entry_t { mjd=54101.0, dt=65.1464 },
new deltat_entry_t { mjd=54466.0, dt=65.4573 },
new deltat_entry_t { mjd=54832.0, dt=65.7768 },
new deltat_entry_t { mjd=55197.0, dt=66.0699 },
new deltat_entry_t { mjd=55562.0, dt=66.3246 },
new deltat_entry_t { mjd=55927.0, dt=66.6030 },
new deltat_entry_t { mjd=56293.0, dt=66.9069 },
new deltat_entry_t { mjd=56658.0, dt=67.2810 },
new deltat_entry_t { mjd=57023.0, dt=67.6439 },
new deltat_entry_t { mjd=57388.0, dt=68.1024 },
new deltat_entry_t { mjd=57754.0, dt=68.5927 },
new deltat_entry_t { mjd=58119.0, dt=68.9676 },
new deltat_entry_t { mjd=58484.0, dt=69.2201 },
new deltat_entry_t { mjd=58849.0, dt=69.87 },
new deltat_entry_t { mjd=59214.0, dt=70.39 },
new deltat_entry_t { mjd=59580.0, dt=70.91 },
new deltat_entry_t { mjd=59945.0, dt=71.40 },
new deltat_entry_t { mjd=60310.0, dt=71.88 },
new deltat_entry_t { mjd=60675.0, dt=72.36 },
new deltat_entry_t { mjd=61041.0, dt=72.83 },
new deltat_entry_t { mjd=61406.0, dt=73.32 },
new deltat_entry_t { mjd=61680.0, dt=73.66 }
};
        private const double T0 = 2451545.0;
        private const double MJD_BASIS = 2400000.5;
        private const double Y2000_IN_MJD  =  T0 - MJD_BASIS;

        /// <summary>
        /// The minimum year value supported by Astronomy Engine.
        /// </summary>
        public const int MinYear = 1700;

        /// <summary>
        /// The maximum year value supported by Astronomy Engine.
        /// </summary>
        public const int MaxYear = 2200;

        private static double DeltaT(double mjd)
        {
            int lo, hi, c;
            double frac;

            if (mjd <= DT[0].mjd)
                return DT[0].dt;

            if (mjd >= DT[DT.Length-1].mjd)
                return DT[DT.Length-1].dt;

            // Do a binary search to find the pair of indexes this mjd lies between.

            lo = 0;
            hi = DT.Length-2;   // make sure there is always an array element after the one we are looking at.
            for(;;)
            {
                if (lo > hi)
                {
                    // This should never happen unless there is a bug in the binary search.
                    throw new Exception("Could not find delta-t value");
                }

                c = (lo + hi) / 2;
                if (mjd < DT[c].mjd)
                    hi = c-1;
                else if (mjd > DT[c+1].mjd)
                    lo = c+1;
                else
                {
                    frac = (mjd - DT[c].mjd) / (DT[c+1].mjd - DT[c].mjd);
                    return DT[c].dt + frac*(DT[c+1].dt - DT[c].dt);
                }
            }
        }

        internal static double TerrestrialTime(double ut)
        {
            return ut + DeltaT(ut + Y2000_IN_MJD)/86400.0;
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
        /// This is different from the behavior of #GeoVector.
        ///
        /// If given an invalid value for `body`, or the body is `Body.Pluto` and the `time` is outside
        /// the year range 1700..2200, this function will throw an `ArgumentException`.
        /// </remarks>
        /// <param name="body">A body for which to calculate a heliocentric position: the Sun, Moon, or any of the planets.</param>
        /// <param name="time">The date and time for which to calculate the position.</param>
        /// <returns>A heliocentric position vector of the center of the given body.</returns>
        public static AstroVector HelioVector(Body body, AstroTime time)
        {
            switch (body)
            {
                case Body.Sun:
                    return new AstroVector(0.0, 0.0, 0.0, time);

                default:
                    throw new ArgumentException(string.Format("Invalid body: {0}", body));
            }
        }
    }
}
