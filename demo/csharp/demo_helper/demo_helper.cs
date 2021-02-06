using System;
using System.Globalization;
using System.Threading;
using CosineKitty;

namespace demo_helper
{
    public class DemoHelper
    {
        static DemoHelper()
        {
            // Force use of "." for the decimal mark, regardless of local culture settings.
            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
        }

        public static void ParseArgs(string program, string[] args, out Observer observer, out AstroTime time)
        {
            if (args.Length == 2 || args.Length == 3)
            {
                double latitude;
                if (!double.TryParse(args[0], out latitude) || (latitude < -90.0) || (latitude > +90.0))
                    throw new ArgumentException(string.Format("ERROR({0}): Invalid latitude '{1}' on command line.", program, args[0]));

                double longitude;
                if (!double.TryParse(args[1], out longitude) || (longitude < -180.0) || (longitude > +180.0))
                    throw new ArgumentException(string.Format("ERROR({0}): Invalid longitude '{1}' on command line.", program, args[1]));

                observer = new Observer(latitude, longitude, 0.0);

                if (args.Length == 3)
                {
                    // Time is present on the command line, so use it.
                    time = ParseTime(program, args[2]);
                }
                else
                {
                    // Time is absent on the command line, so use the current time.
                    time = new AstroTime(DateTime.UtcNow);
                }
            }
            else
            {
                throw new ArgumentException(string.Format("USAGE: {0} latitude longitude [yyyy-mm-ddThh:mm:ssZ]", program));
            }
        }

        public static AstroTime ParseTime(string program, string text)
        {
            DateTime dt;
            if (!DateTime.TryParse(text, out dt))
                throw new ArgumentException(string.Format("ERROR({0}): Cannot parse date/time string from '{1}'", program, text));
            return new AstroTime(dt);
        }
    }
}
