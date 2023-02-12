using System;
using System.Linq;
using demo_helper;
using CosineKitty;

namespace solar_time
{
    class Program
    {
        static int Main(string[] args)
        {
            DemoHelper.ParseArgs("solar_time", args, out Observer observer, out AstroTime time);
            double hourAngle = Astronomy.HourAngle(Body.Sun, time, observer);
            double solarTime = (hourAngle + 12.0) % 24.0;
            int millis = (int)Math.Round(solarTime * 3.6e+6);
            int seconds = millis / 1000;
            millis %= 1000;
            int minutes = seconds / 60;
            seconds %= 60;
            int hours = minutes / 60;
            minutes %= 60;
            hours %= 24;
            Console.WriteLine($"True solar time = {solarTime:F4} hours ({hours:00}:{minutes:00}:{seconds:00}.{millis:000})");
            return 0;
        }
    }
}
