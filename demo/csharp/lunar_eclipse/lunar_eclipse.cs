using System;
using demo_helper;
using CosineKitty;

namespace lunar_eclipse
{
    class Program
    {
        static int Main(string[] args)
        {
            AstroTime time;
            switch (args.Length)
            {
                case 0:
                    time = new AstroTime(DateTime.Now);
                    break;

                case 1:
                    time = DemoHelper.ParseTime("lunar_eclipse", args[0]);
                    break;

                default:
                    Console.WriteLine("USAGE: lunar_eclipse [date]");
                    return 1;
            }

            int count = 0;
            LunarEclipseInfo eclipse = Astronomy.SearchLunarEclipse(time);
            for(;;)
            {
                if (eclipse.kind != EclipseKind.Penumbral)
                {
                    PrintEclipse(eclipse);
                    if (++count == 10)
                        break;
                }
                eclipse = Astronomy.NextLunarEclipse(eclipse.center);
            }

            return 0;
        }

        static void PrintEclipse(LunarEclipseInfo eclipse)
        {
            // Calculate beginning/ending of different phases
            // of an eclipse by subtracting/adding the center time
            // with the number of minutes indicated by the "semi-duration"
            // fields sd_partial and sd_total.
            const double MINUTES_PER_DAY = 24 * 60;

            AstroTime p1 = eclipse.center.AddDays(-eclipse.sd_partial / MINUTES_PER_DAY);
            Console.WriteLine("{0}  Partial eclipse begins.", p1);

            if (eclipse.sd_total > 0.0)
            {
                AstroTime t1 = eclipse.center.AddDays(-eclipse.sd_total / MINUTES_PER_DAY);
                Console.WriteLine("{0}  Total eclipse begins.", t1);
            }

            Console.WriteLine("{0}  Peak of {1} eclipse.", eclipse.center, eclipse.kind.ToString().ToLowerInvariant());

            if (eclipse.sd_total > 0.0)
            {
                AstroTime t2 = eclipse.center.AddDays(+eclipse.sd_total / MINUTES_PER_DAY);
                Console.WriteLine("{0}  Total eclipse ends.", t2);
            }

            AstroTime p2 = eclipse.center.AddDays(+eclipse.sd_partial / MINUTES_PER_DAY);
            Console.WriteLine("{0}  Partial eclipse ends.", p2);
            Console.WriteLine();
        }
    }
}
