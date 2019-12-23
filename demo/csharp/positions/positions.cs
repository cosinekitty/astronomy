using System;
using CosineKitty;
using demo_helper;

namespace positions
{
    class Program
    {
        static int Main(string[] args)
        {
            var bodies = new Body[]
            {
                Body.Sun, Body.Moon, Body.Mercury, Body.Venus, Body.Mars,
                Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto
            };

            Observer observer;
            AstroTime time;
            DemoHelper.ParseArgs("positions", args, out observer, out time);
            Console.WriteLine("UTC date = {0}", time);
            Console.WriteLine();
            Console.WriteLine("BODY           RA      DEC       AZ      ALT");
            foreach (Body body in bodies)
            {
                Equatorial equ_2000 = Astronomy.Equator(body, time, observer, EquatorEpoch.J2000, Aberration.Corrected);
                Equatorial equ_ofdate = Astronomy.Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
                Topocentric hor = Astronomy.Horizon(time, observer, equ_ofdate.ra, equ_ofdate.dec, Refraction.Normal);
                Console.WriteLine("{0,-8} {1,8:0.00} {2,8:0.00} {3,8:0.00} {4,8:0.00}", body, equ_2000.ra, equ_2000.dec, hor.azimuth, hor.altitude);
            }

            return 0;
        }
    }
}
