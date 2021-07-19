using System;
using System.Globalization;
using System.Threading;
using CosineKitty;

namespace gravity
{
    class Program
    {
        const string UsageText = @"
    USAGE:

    gravity latitude height

    Calculates the gravitational acceleration experienced
    by an observer on the surface of the Earth at the specified
    latitude (degrees north of the equator) and height
    (meters above sea level).
    The output is the gravitational acceleration in m/s^2.
";

        static int Main(string[] args)
        {
            if (args.Length != 2)
            {
                Console.WriteLine("{0}", UsageText);
                return 1;
            }

            const double MAX_HEIGHT_METERS = 100000.0;

            // Force use of "." for the decimal mark, regardless of local culture settings.
            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

            double latitude = double.Parse(args[0]);
            if (!double.IsFinite(latitude) || latitude < -90.0 || latitude > +90.0)
            {
                Console.WriteLine($"ERROR: Invalid latitude '{args[0]}'. Must be a number between -90 and +90.");
                return 1;
            }

            double height = double.Parse(args[1]);
            if (!double.IsFinite(height) || height < 0.0 || height > MAX_HEIGHT_METERS)
            {
                Console.WriteLine($"ERROR: Invalid height '{args[1]}'. Must be a number between 0 and {MAX_HEIGHT_METERS:F0}.");
                return 1;
            }

            double gravity = Astronomy.ObserverGravity(latitude, height);
            Console.WriteLine($"latitude = {latitude,8:F4},  height = {height,6:F0},  gravity = {gravity,8:F6}");
            return 0;
        }
    }
}
