using CosineKitty;
using System;

namespace AstronomyTest
{
    class Program
    {
        static int Main(string[] args)
        {
            var now = new AstroTime(DateTime.Now);
            AstroTime full = Astronomy.SearchMoonPhase(180.0, now, 40.0);
            Console.WriteLine("The next full moon will be: {0}", full);
            return 0;
        }
    }
}
