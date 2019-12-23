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
            if (args.Length != 1 || !int.TryParse(args[0], out year) || year < Astronomy.MinYear || year > Astronomy.MaxYear)
            {
                Console.WriteLine("ERROR: Must provide year {0}..{1} on command line.", Astronomy.MinYear, Astronomy.MaxYear);
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
