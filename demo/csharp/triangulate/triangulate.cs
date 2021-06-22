using System;
using CosineKitty;

//
//    triangulate.cs  -  by Don Cross - 2021-06-22
//
//    Example C# program for Astronomy Engine:
//    https://github.com/cosinekitty/astronomy
//


namespace triangulate
{
    class Program
    {
        const string UsageText = @"
USAGE:  triangulate  lat1 lon1 elv1 az1 alt1  lat2 lon2 elv2 az2 alt2

Calculate the best-fit location of a point as observed
from two different locations on or near the Earth's surface.

lat1, lat2 = Geographic latitudes in degrees north of the equator.
lon1, lon2 = Geographic longitudes in degrees east of the prime meridian.
elv1, elv2 = Elevations above sea level in meters.
az1,  az2  = Azimuths toward observed object in degrees clockwise from north.
alt1, alt2 = Altitude angles toward observed object in degrees above horizon.

This program extrapolates lines in the given directions from the two
geographic locations and finds the location in space where they
come closest to intersecting. It then prints out the coordinates
of that triangulation point, along with the error radius in meters.
";

        static int Main(string[] args)
        {
            if (args.Length != 10)
            {
                Console.WriteLine(UsageText);
                return 1;
            }

            // Validate and parse command line arguments.
            double lat1 = ParseNumber("lat1", args[0]);
            double lon1 = ParseNumber("lon1", args[1]);
            double elv1 = ParseNumber("elv1", args[2]);
            double  az1 = ParseNumber("az1",  args[3]);
            double alt1 = ParseNumber("alt1", args[4]);
            double lat2 = ParseNumber("lat2", args[5]);
            double lon2 = ParseNumber("lon2", args[6]);
            double elv2 = ParseNumber("elv2", args[7]);
            double  az2 = ParseNumber("az2",  args[8]);
            double alt2 = ParseNumber("alt2", args[9]);

            var obs1 = new Observer(lat1, lon1, elv1);
            var obs2 = new Observer(lat2, lon2, elv2);

            // Use an arbitrary but consistent time for the Earth's rotation.
            AstroTime time = new AstroTime(0.0);

            // Convert geographic coordinates of the observers to vectors.
            AstroVector pos1 = Astronomy.ObserverVector(time, obs1, EquatorEpoch.OfDate);
            AstroVector pos2 = Astronomy.ObserverVector(time, obs2, EquatorEpoch.OfDate);

            // Convert horizontal coordinates into unit direction vectors.
            AstroVector dir1 = DirectionVector(time, obs1, alt1, az1);
            AstroVector dir2 = DirectionVector(time, obs2, alt2, az2);

            // Find the closest point between the skew lines.
            return Intersect(pos1, dir1, pos2, dir2);
        }

        static double ParseNumber(string name, string text)
        {
            if (double.TryParse(text, out double value) && double.IsFinite(value))
                return value;

            throw new ArgumentException($"Invalid value for {name}: {text}");
        }

        static AstroVector DirectionVector(AstroTime time, Observer observer, double altitude, double azimuth)
        {
            // Convert horizontal angles to a horizontal unit vector.
            var hor = new Spherical(altitude, azimuth, 1.0);
            AstroVector hvec = Astronomy.VectorFromHorizon(hor, time, Refraction.None);

            // Find the rotation matrix that converts horizontal vectors to equatorial vectors.
            RotationMatrix rot = Astronomy.Rotation_HOR_EQD(time, observer);

            // Rotate the horizontal (HOR) vector to an equator-of-date (EQD) vector.
            AstroVector evec = Astronomy.RotateVector(rot, hvec);

            return evec;
        }

        static int Intersect(AstroVector pos1, AstroVector dir1, AstroVector pos2, AstroVector dir2)
        {
            double F = dir1 * dir2;
            AstroVector amb = pos1 - pos2;
            double E = dir1 * amb;
            double G = dir2 * amb;
            double denom = 1.0 - F*F;
            if (denom == 0.0)
            {
                Console.WriteLine("ERROR: Cannot solve because directions are parallel.");
                return 1;
            }

            double u = (F*G - E) / denom;
            double v = G + F*u;
            if (u < 0.0 || v < 0.0)
            {
                Console.WriteLine("ERROR: Lines of sight do not converge.");
                return 1;
            }

            AstroVector a = pos1 + u*dir1;
            AstroVector b = pos2 + v*dir2;
            AstroVector c = (a + b) / 2.0;
            AstroVector miss = a - b;

            double dist = (Astronomy.KM_PER_AU * 1000 / 2) * miss.Length();   // error radius in meters
            Observer obs = Astronomy.VectorObserver(c, EquatorEpoch.OfDate);

            Console.WriteLine($"Solution: lat = {obs.latitude:F6}, lon = {obs.longitude:F6}, elv = {obs.height:F3} meters; error = {dist:F3} meters.");
            return 0;
        }
    }
}
