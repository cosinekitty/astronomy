using System;
using CosineKitty;

namespace csharp_test
{
    class Program
    {
        static int Main(string[] args)
        {
            AstroTime time = new AstroTime(0.0);
            Console.WriteLine("ut={0:0.000000}, tt={1:0.000000}", time.ut, time.tt);
            return 0;
        }
    }
}
