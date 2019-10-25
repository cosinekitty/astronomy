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

        private static readonly deltat_entry_t[] DT = new deltat_entry_t[]
        {
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
            new vsop_term_t(0.00001658058, 4.90206728031, 20426.57109242200)
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

        private static readonly vsop_term_t[] vsop_lon_Earth_1 = new vsop_term_t[]
        {
            new vsop_term_t(0.00227777722, 3.41376620530, 6283.07584999140),
            new vsop_term_t(0.00003805678, 3.37063423795, 12566.15169998280)
        };

        private static readonly vsop_series_t[] vsop_lon_Earth = new vsop_series_t[]
        {
            new vsop_series_t(null),
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
            new vsop_term_t(0.00000472110, 3.66100022149, 5884.92684658320)
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
            new vsop_term_t(0.00012749023, 2.71550286592, 1052.26838318840)
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
            new vsop_term_t(0.00020936596, 0.46349251129, 735.87651353180)
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
            new vsop_term_t(0.00029156413, 3.18056336700, 77.75054398390)
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
            new vsop_term_t(0.00274571975, 1.84552258866, 175.16605980020)
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

                case Body.Mercury:
                case Body.Venus:
                case Body.Earth:
                case Body.Mars:
                case Body.Jupiter:
                case Body.Saturn:
                case Body.Uranus:
                case Body.Neptune:
                    return CalcVsop(vsop[(int)body], time);

                default:
                    throw new ArgumentException(string.Format("Invalid body: {0}", body));
            }
        }
    }
}
