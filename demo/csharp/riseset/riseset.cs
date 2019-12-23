using System;
using CosineKitty;
using demo_helper;

namespace riseset
{
    class Program
    {
        static int PrintEvent(string name, AstroTime time)
        {
            if (time == null)
            {
                Console.WriteLine("ERROR: Search failed for {0}", name);
                return 1;
            }
            Console.WriteLine("{0,-8} : {1}", name, time);
            return 0;
        }

        static int Main(string[] args)
        {
            Observer observer;
            AstroTime time;
            DemoHelper.ParseArgs("riseset", args, out observer, out time);
            Console.WriteLine("search   : {0}", time);
            if (0 != PrintEvent("sunrise",  Astronomy.SearchRiseSet(Body.Sun,  observer, Direction.Rise, time, 300.0))) return 1;
            if (0 != PrintEvent("sunset",   Astronomy.SearchRiseSet(Body.Sun,  observer, Direction.Set,  time, 300.0))) return 1;
            if (0 != PrintEvent("moonrise", Astronomy.SearchRiseSet(Body.Moon, observer, Direction.Rise, time, 300.0))) return 1;
            if (0 != PrintEvent("moonset",  Astronomy.SearchRiseSet(Body.Moon, observer, Direction.Set,  time, 300.0))) return 1;
            return 0;
        }
    }
}
