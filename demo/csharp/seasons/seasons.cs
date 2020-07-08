using System;
using CosineKitty;
using demo_helper;

namespace seasons
{
    class Program
    {
        static int Main(string[] args)
        {
            int year;
            if (args.Length != 1 || !int.TryParse(args[0], out year))
            {
                Console.WriteLine("ERROR: Must provide a year value on the command line.");
                return 1;
            }
            SeasonsInfo seasons = Astronomy.Seasons(year);
            Console.WriteLine("March equinox     : {0}", seasons.mar_equinox);
            Console.WriteLine("June solstice     : {0}", seasons.jun_solstice);
            Console.WriteLine("September equinox : {0}", seasons.sep_equinox);
            Console.WriteLine("December solstice : {0}", seasons.dec_solstice);
            return 0;
        }
    }
}
