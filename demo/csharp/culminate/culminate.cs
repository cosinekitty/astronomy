using System;
using demo_helper;
using CosineKitty;

namespace culminate
{
    class Program
    {
        static int Main(string[] args)
        {
            Observer observer;
            AstroTime time;
            DemoHelper.ParseArgs("culminate", args, out observer, out time);

            var bodies = new Body[]
            {
                Body.Sun, Body.Moon, Body.Mercury, Body.Venus, Body.Mars,
                Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto
            };

            Console.WriteLine("search   : {0}", time);
            foreach (Body body in bodies)
            {
                HourAngleInfo evt = Astronomy.SearchHourAngle(body, observer, 0.0, time);
                Console.WriteLine("{0,-8} : {1}  altitude={2:##0.00}  azimuth={3:##0.00}",
                    body, evt.time, evt.hor.altitude, evt.hor.azimuth);
            }

            return 0;
        }
    }
}
