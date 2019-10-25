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
    public enum Body
    {
        Invalid = -1,
        Mercury,
        Venus,
        Earth,
        Mars,
        Jupiter,
        Saturn,
        Uranus,
        Neptune,
        Pluto,
        Sun,
        Moon,
    }

    public class AstroTime
    {
        private static readonly DateTime Origin = new DateTime(2000, 1, 1, 12, 0, 0, DateTimeKind.Utc);

        public readonly double ut;
        public readonly double tt;
        internal double psi;
        internal double eps;

        public AstroTime(double ut)
        {
            this.ut = ut;
            this.tt = Astronomy.TerrestrialTime(ut);
            this.psi = this.eps = double.NaN;
        }

        public AstroTime(DateTime d)
            : this((d - Origin).TotalDays)
        {
        }

        public DateTime ToUtcDateTime()
        {
            return Origin.AddDays(ut).ToUniversalTime();
        }

        public override string ToString()
        {
            return ToUtcDateTime().ToString("yyyy-MM-ddThh:mm:ss.fffZ");
        }
    }

    public class AstroVector
    {
        public readonly double x;
        public readonly double y;
        public readonly double z;
        public readonly AstroTime t;

        public AstroVector(double x, double y, double z, AstroTime t)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.t = t;
        }
    }

    public class Observer
    {
        public readonly double Latitude;
        public readonly double Longitude;
        public readonly double Height;

        public Observer(double latitude, double longitude, double height)
        {
            Latitude = latitude;
            Longitude = longitude;
            Height = height;
        }
    }

    public static class Astronomy
    {
        private static readonly deltat_entry_t[] DT = $ASTRO_DELTA_T();
        private const double T0 = 2451545.0;
        private const double MJD_BASIS = 2400000.5;
        private const double Y2000_IN_MJD  =  T0 - MJD_BASIS;

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

        /// <summary>Calculates heliocentric Cartesian coordinates of a body in the J2000 equatorial system.</summary>
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

    internal struct deltat_entry_t
    {
        public double mjd;
        public double dt;
    }
}
