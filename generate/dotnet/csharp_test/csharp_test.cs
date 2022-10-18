using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading;

using CosineKitty;

namespace csharp_test
{
    class Program
    {
        static bool Verbose;
        const double SECONDS_PER_DAY = 86400.0;

        static void Debug(string format, params object[] args)
        {
            if (Verbose)
                Console.WriteLine(format, args);
        }

        struct Test
        {
            public string Name;
            public Func<int> TestFunc;
            public Test(string name, Func<int> testFunc)
            {
                this.Name = name;
                this.TestFunc = testFunc;
            }
        }

        static Test[] UnitTests = new Test[]
        {
            new Test("time", TestTime),
            new Test("moon", MoonTest),
            new Test("geoid", GeoidTest),
            new Test("constellation", ConstellationTest),
            new Test("dates250", DatesIssue250),
            new Test("elongation", ElongationTest),
            new Test("global_solar_eclipse", GlobalSolarEclipseTest),
            new Test("gravsim", GravitySimulatorTest),
            new Test("jupiter_moons", JupiterMoonsTest),
            new Test("libration", LibrationTest),
            new Test("lagrange", LagrangeTest),
            new Test("local_solar_eclipse", LocalSolarEclipseTest),
            new Test("lunar_apsis", LunarApsisTest),
            new Test("lunar_eclipse", LunarEclipseTest),
            new Test("lunar_eclipse_78", LunarEclipseIssue78),
            new Test("lunar_fraction", LunarFractionTest),
            new Test("magnitude", MagnitudeTest),
            new Test("moonphase", MoonPhaseTest),
            new Test("moon_nodes", MoonNodesTest),
            new Test("moon_reverse", MoonReverseTest),
            new Test("planet_apsis", PlanetApsisTest),
            new Test("pluto", PlutoCheck),
            new Test("refraction", RefractionTest),
            new Test("riseset", RiseSetTest),
            new Test("riseset_reverse", RiseSetReverseTest),
            new Test("rotation", RotationTest),
            new Test("seasons", SeasonsTest),
            new Test("seasons187", SeasonsIssue187),
            new Test("sidereal", SiderealTimeTest),
            new Test("transit", TransitTest),
            new Test("astro_check", AstroCheck),
            new Test("barystate", BaryStateTest),
            new Test("heliostate", HelioStateTest),
            new Test("topostate", TopoStateTest),
            new Test("aberration", AberrationTest),
            new Test("twilight", TwilightTest),
            new Test("axis", AxisTest),
        };

        static int Main(string[] args)
        {
            try
            {
                // Force use of "." for the decimal mark, regardless of local culture settings.
                Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

                if (args.Length > 0 && args[0] == "-v")
                {
                    Verbose = true;
                    args = args.Skip(1).ToArray();
                }

                if (args.Length == 1)
                {
                    string name = args[0];
                    if (name == "all")
                    {
                        Console.WriteLine("csharp_test: starting");
                        foreach (Test t in UnitTests)
                            if (0 != t.TestFunc())
                                return 1;
                        Console.WriteLine("csharp_test: PASS");
                        return 0;
                    }

                    foreach (Test t in UnitTests)
                        if (t.Name == name)
                            return t.TestFunc();
                }

                Console.WriteLine("csharp_test: Invalid command line parameters.");
                return 1;
            }
            catch (Exception ex)
            {
                Console.WriteLine("charp_test: EXCEPTION: {0}", ex);
                return 1;
            }
        }

        static double v(double x)
        {
            if (!double.IsFinite(x))
                throw new ArgumentException("Non-finite result");
            return x;
        }

        static double abs(double x)
        {
            return Math.Abs(v(x));
        }

        static double max(double a, double b)
        {
            return Math.Max(v(a), v(b));
        }

        static double min(double a, double b)
        {
            return Math.Min(v(a), v(b));
        }

        static double sqrt(double x)
        {
            return v(Math.Sqrt(v(x)));
        }

        static double sin(double x)
        {
            return Math.Sin(v(x));
        }

        static double cos(double x)
        {
            return Math.Cos(v(x));
        }

        static readonly char[] TokenSeparators = new char[] { ' ', '\t', '\r', '\n' };

        static string[] Tokenize(string line)
        {
            return line.Split(TokenSeparators, StringSplitOptions.RemoveEmptyEntries);
        }

        static int TestTime()
        {
            const int year = 2018;
            const int month = 12;
            const int day = 2;
            const int hour = 18;
            const int minute = 30;
            const int second = 12;
            const int milli = 543;

            DateTime d = new DateTime(year, month, day, hour, minute, second, milli, DateTimeKind.Utc);
            AstroTime time = new AstroTime(d);
            Console.WriteLine("C# TestTime: text={0}, ut={1}, tt={2}", time.ToString(), time.ut.ToString("F6"), time.tt.ToString("F6"));

            const double expected_ut = 6910.270978506945;
            double diff = time.ut - expected_ut;
            if (abs(diff) > 1.0e-12)
            {
                Console.WriteLine("C# TestTime: ERROR - excessive UT error {0}", diff);
                return 1;
            }

            const double expected_tt = 6910.271800214368;
            diff = time.tt - expected_tt;
            if (abs(diff) > 1.0e-12)
            {
                Console.WriteLine("C# TestTime: ERROR - excessive TT error {0}", diff);
                return 1;
            }

            DateTime utc = time.ToUtcDateTime();
            if (utc.Year != year || utc.Month != month || utc.Day != day || utc.Hour != hour || utc.Minute != minute || utc.Second != second || utc.Millisecond != milli)
            {
                Console.WriteLine("C# TestTime: ERROR - Expected {0:o}, found {1:o}", d, utc);
                return 1;
            }

            return 0;
        }

        static int MoonTest()
        {
            var time = new AstroTime(2019, 6, 24, 15, 45, 37);
            AstroVector vec = Astronomy.GeoVector(Body.Moon, time, Aberration.None);
            Console.WriteLine("C# MoonTest: {0} {1} {2}", vec.x.ToString("f17"), vec.y.ToString("f17"), vec.z.ToString("f17"));

            double dx = vec.x - (+0.002674037026701135);
            double dy = vec.y - (-0.0001531610316600666);
            double dz = vec.z - (-0.0003150159927069429);
            double diff = sqrt(dx*dx + dy*dy + dz*dz);
            Console.WriteLine("C# MoonTest: diff = {0}", diff.ToString("g5"));
            if (diff > 4.34e-19)
            {
                Console.WriteLine("C# MoonTest: EXCESSIVE ERROR");
                return 1;
            }

            return 0;
        }

        static int AstroCheck()
        {
            const string filename = "csharp_check.txt";
            using (StreamWriter outfile = File.CreateText(filename))
            {
                var bodylist = new Body[]
                {
                    Body.Sun, Body.Mercury, Body.Venus, Body.Earth, Body.Mars,
                    Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto,
                    Body.SSB, Body.EMB
                };

                var observer = new Observer(29.0, -81.0, 10.0);
                var time = new AstroTime(new DateTime(1700, 1, 1, 0, 0, 0, DateTimeKind.Utc));
                var stop = new AstroTime(new DateTime(2200, 1, 1, 0, 0, 0, DateTimeKind.Utc));

                AstroVector pos;
                Equatorial j2000, ofdate;
                Topocentric hor;

                outfile.WriteLine("o {0} {1} {2}", observer.latitude, observer.longitude, observer.height);
                while (time.tt < stop.tt)
                {
                    foreach (Body body in bodylist)
                    {
                        pos = Astronomy.HelioVector(body, time);
                        outfile.WriteLine("v {0} {1} {2} {3} {4}", body, pos.t.tt.ToString("G18"), pos.x.ToString("G18"), pos.y.ToString("G18"), pos.z.ToString("G18"));
                        if (body != Body.Earth && body != Body.SSB && body != Body.EMB)
                        {
                            j2000 = Astronomy.Equator(body, time, observer, EquatorEpoch.J2000, Aberration.None);
                            ofdate = Astronomy.Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
                            hor = Astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None);
                            outfile.WriteLine("s {0} {1} {2} {3} {4} {5} {6} {7}",
                                body,
                                time.tt.ToString("G18"), time.ut.ToString("G18"),
                                j2000.ra.ToString("G18"), j2000.dec.ToString("G18"), j2000.dist.ToString("G18"),
                                hor.azimuth.ToString("G18"), hor.altitude.ToString("G18"));
                        }
                    }

                    pos = Astronomy.GeoVector(Body.Moon, time, Aberration.None);
                    outfile.WriteLine("v GM {0} {1} {2} {3}", pos.t.tt.ToString("G18"), pos.x.ToString("G18"), pos.y.ToString("G18"), pos.z.ToString("G18"));
                    j2000 = Astronomy.Equator(Body.Moon, time, observer, EquatorEpoch.J2000, Aberration.None);
                    ofdate = Astronomy.Equator(Body.Moon, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
                    hor = Astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None);
                    outfile.WriteLine("s GM {0} {1} {2} {3} {4} {5} {6}",
                        time.tt.ToString("G18"), time.ut.ToString("G18"),
                        j2000.ra.ToString("G18"), j2000.dec.ToString("G18"), j2000.dist.ToString("G18"),
                        hor.azimuth.ToString("G18"), hor.altitude.ToString("G18"));

                    JupiterMoonsInfo jm = Astronomy.JupiterMoons(time);
                    for (int mindex = 0; mindex < NUM_JUPITER_MOONS; ++mindex)
                    {
                        StateVector moon = SelectJupiterMoon(jm, mindex);
                        outfile.WriteLine($"j {mindex} {time.tt:G18} {time.ut:G18} {moon.x:G18} {moon.y:G18} {moon.z:G18} {moon.vx:G18} {moon.vy:G18} {moon.vz:G18}");
                    }

                    time = time.AddDays(10.0 + Math.PI/100.0);
                }
            }
            Console.WriteLine("C# AstroCheck: finished");
            return 0;
        }

        static int SeasonsTest()
        {
            const string filename = "../../seasons/seasons.txt";
            var re = new Regex(@"^(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z\s+([A-Za-z]+)\s*$");
            using (StreamReader infile = File.OpenText(filename))
            {
                string line;
                int lnum = 0;
                int current_year = 0;
                int mar_count=0, jun_count=0, sep_count=0, dec_count=0;
                double max_minutes = 0.0;
                SeasonsInfo seasons = new SeasonsInfo();
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    // 2019-01-03T05:20Z Perihelion
                    // 2019-03-20T21:58Z Equinox
                    // 2019-06-21T15:54Z Solstice
                    // 2019-07-04T22:11Z Aphelion
                    // 2019-09-23T07:50Z Equinox
                    // 2019-12-22T04:19Z Solstice
                    Match m = re.Match(line);
                    if (!m.Success)
                    {
                        Console.WriteLine("C# SeasonsTest: ERROR {0} line {1}: cannot parse", filename, lnum);
                        return 1;
                    }

                    int year = int.Parse(m.Groups[1].Value);
                    int month = int.Parse(m.Groups[2].Value);
                    int day = int.Parse(m.Groups[3].Value);
                    int hour = int.Parse(m.Groups[4].Value);
                    int minute = int.Parse(m.Groups[5].Value);
                    string name = m.Groups[6].Value;
                    var correct_time = new AstroTime(year, month, day, hour, minute, 0);

                    if (year != current_year)
                    {
                        current_year = year;
                        seasons = Astronomy.Seasons(year);
                    }

                    AstroTime calc_time = null;
                    if (name == "Equinox")
                    {
                        switch (month)
                        {
                            case 3:
                                calc_time = seasons.mar_equinox;
                                ++mar_count;
                                break;

                            case 9:
                                calc_time = seasons.sep_equinox;
                                ++sep_count;
                                break;

                            default:
                                Console.WriteLine("C# SeasonsTest: {0} line {1}: Invalid equinox date in test data.", filename, lnum);
                                return 1;
                        }
                    }
                    else if (name == "Solstice")
                    {
                        switch (month)
                        {
                            case 6:
                                calc_time = seasons.jun_solstice;
                                ++jun_count;
                                break;

                            case 12:
                                calc_time = seasons.dec_solstice;
                                ++dec_count;
                                break;

                            default:
                                Console.WriteLine("C# SeasonsTest: {0} line {1}: Invalid solstice date in test data.", filename, lnum);
                                return 1;
                        }
                    }
                    else if (name == "Aphelion")
                    {
                        // not yet calculated
                        continue;
                    }
                    else if (name == "Perihelion")
                    {
                        // not yet calculated
                        continue;
                    }
                    else
                    {
                        Console.WriteLine("C# SeasonsTest: {0} line {1}: unknown event type {2}", filename, lnum, name);
                        return 1;
                    }

                    // Verify that the calculated time matches the correct time for this event.
                    double diff_minutes = (24.0 * 60.0) * abs(calc_time.tt - correct_time.tt);
                    if (diff_minutes > max_minutes)
                        max_minutes = diff_minutes;

                    if (diff_minutes > 2.37)
                    {
                        Console.WriteLine("C# SeasonsTest: {0} line {1}: excessive error ({2}): {3} minutes.", filename, lnum, name, diff_minutes);
                        return 1;
                    }
                }
                Console.WriteLine("C# SeasonsTest: verified {0} lines from file {1} : max error minutes = {2:0.000}", lnum, filename, max_minutes);
                Console.WriteLine("C# SeasonsTest: Event counts: mar={0}, jun={1}, sep={2}, dec={3}", mar_count, jun_count, sep_count, dec_count);
                return 0;
            }
        }

        static int SeasonsIssue187()
        {
            // This is a regression test for:
            // https://github.com/cosinekitty/astronomy/issues/187
            // For years far from the present, the seasons search was sometimes failing.
            for (int year = 1; year <= 9999; ++year)
                Astronomy.Seasons(year);

            return 0;
        }

        static int CheckDecemberSolstice(int year, string expected)
        {
            SeasonsInfo si = Astronomy.Seasons(year);
            string actual = si.dec_solstice.ToString();
            if (actual != expected)
            {
                Console.WriteLine($"C# DatesIssue250: FAIL: year {year}, expected [{expected}], actual [{actual}]");
                return 1;
            }
            return 0;
        }

        static int DatesIssue250()
        {
            // Make sure we can handle dates outside the range supported by System.DateTime.
            // https://github.com/cosinekitty/astronomy/issues/250
            if (0 != CheckDecemberSolstice( 2022, "2022-12-21T21:47:58.189Z")) return 1;
            if (0 != CheckDecemberSolstice(-2300, "-002300-12-19T16:22:26.325Z")) return 1;
            if (0 != CheckDecemberSolstice(12345, "+012345-12-11T13:30:10.041Z")) return 1;
            Console.WriteLine("C# DatesIssue250: PASS");
            return 0;
        }

        static int MoonPhaseTest()
        {
            const string filename = "../../moonphase/moonphases.txt";
            using (StreamReader infile = File.OpenText(filename))
            {
                const double threshold_seconds = 90.0;
                int lnum = 0;
                string line;
                double max_arcmin = 0.0;
                int prev_year = 0;
                int expected_quarter = 0;
                int quarter_count = 0;
                double maxdiff = 0.0;
                MoonQuarterInfo mq = new MoonQuarterInfo();
                var re = new Regex(@"^([0-3])\s+(\d+)-(\d+)-(\d+)T(\d+):(\d+):(\d+)\.000Z$");
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    // 0 1800-01-25T03:21:00.000Z
                    // 1 1800-02-01T20:40:00.000Z
                    // 2 1800-02-09T17:26:00.000Z
                    // 3 1800-02-16T15:49:00.000Z
                    Match m = re.Match(line);
                    if (!m.Success)
                    {
                        Console.WriteLine("C# MoonPhaseTest: ERROR {0} line {1}: cannot parse", filename, lnum);
                        return 1;
                    }
                    int quarter = int.Parse(m.Groups[1].Value);
                    int year = int.Parse(m.Groups[2].Value);
                    int month = int.Parse(m.Groups[3].Value);
                    int day = int.Parse(m.Groups[4].Value);
                    int hour = int.Parse(m.Groups[5].Value);
                    int minute = int.Parse(m.Groups[6].Value);
                    int second = int.Parse(m.Groups[7].Value);

                    double expected_elong = 90.0 * quarter;
                    AstroTime expected_time = new AstroTime(year, month, day, hour, minute, second);
                    double calc_elong = Astronomy.MoonPhase(expected_time);
                    double degree_error = abs(calc_elong - expected_elong);
                    if (degree_error > 180.0)
                        degree_error = 360.0 - degree_error;
                    double arcmin = 60.0 * degree_error;
                    if (arcmin > 1.0)
                    {
                        Console.WriteLine("C# MoonPhaseTest({0} line {1}): EXCESSIVE ANGULAR ERROR: {2} arcmin", filename, lnum, arcmin);
                        return 1;
                    }
                    if (arcmin > max_arcmin)
                        max_arcmin = arcmin;

                    if (year != prev_year)
                    {
                        prev_year = year;
                        // The test data contains a single year's worth of data for every 10 years.
                        // Every time we see the year value change, it breaks continuity of the phases.
                        // Start the search over again.
                        AstroTime start_time = new AstroTime(year, 1, 1, 0, 0, 0);
                        mq = Astronomy.SearchMoonQuarter(start_time);
                        expected_quarter = -1;  // we have no idea what the quarter should be
                    }
                    else
                    {
                        // Yet another lunar quarter in the same year.
                        expected_quarter = (1 + mq.quarter) % 4;
                        mq = Astronomy.NextMoonQuarter(mq);

                        // Make sure we find the next expected quarter.
                        if (expected_quarter != mq.quarter)
                        {
                            Console.WriteLine("C# MoonPhaseTest({0} line {1}): SearchMoonQuarter returned quarter {2}, but expected {3}", filename, lnum, mq.quarter, expected_quarter);
                            return 1;
                        }
                    }
                    ++quarter_count;
                    // Make sure the time matches what we expect.
                    double diff_seconds = abs(mq.time.tt - expected_time.tt) * SECONDS_PER_DAY;
                    if (diff_seconds > threshold_seconds)
                    {
                        Console.WriteLine("C# MoonPhaseTest({0} line {1}): excessive time error {2:0.000} seconds", filename, lnum, diff_seconds);
                        return 1;
                    }

                    if (diff_seconds > maxdiff)
                        maxdiff = diff_seconds;
                }

                Console.WriteLine("C# MoonPhaseTest: passed {0} lines for file {1} : max_arcmin = {2:0.000000}, maxdiff = {3:0.000} seconds, {4} quarters",
                    lnum, filename, max_arcmin, maxdiff, quarter_count);

                return 0;
            }
        }

        static int MoonNodesTest()
        {
            const string filename = "../../moon_nodes/moon_nodes.txt";
            using (StreamReader infile = File.OpenText(filename))
            {
                var node = new NodeEventInfo();
                double max_angle = 0.0;
                double max_minutes = 0.0;
                int lnum = 0;
                string line;
                string prev_kind = "?";
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;

                    // Parse the line from the test data file.
                    // A 2001-01-09T13:53Z    7.1233   22.5350
                    // D 2001-01-22T22:22Z   19.1250  -21.4683
                    string[] token = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    if (token.Length != 4)
                    {
                        Console.WriteLine($"C# MoonNodesTest({filename} line {lnum}): syntax error");
                        return 1;
                    }

                    string kind = token[0];
                    if (kind != "A" && kind != "D")
                    {
                        Console.WriteLine($"C# MoonNodesTest({filename} line {lnum}): incorrect node kind [{kind}].");
                        return 1;
                    }

                    if (kind == prev_kind)
                    {
                        Console.WriteLine($"C# MoonNodesTest({filename} line {lnum}): duplicate ascending/descending node.");
                        return 1;
                    }

                    AstroTime time = ParseDate(token[1]);

                    double ra;
                    if (!double.TryParse(token[2], out ra) || !double.IsFinite(ra) || (ra < 0.0) || (ra > 24.0))
                    {
                        Console.WriteLine($"C# MoonNodesTest({filename} line {lnum}): invalid RA.");
                        return 1;
                    }

                    double dec;
                    if (!double.TryParse(token[3], out dec) || !double.IsFinite(dec) || (dec < -90.0) || (dec > +90.0))
                    {
                        Console.WriteLine($"C# MoonNodesTest({filename} line {lnum}): invalid DEC.");
                        return 1;
                    }

                    var sphere = new Spherical(dec, 15.0 * ra, 1.0);
                    AstroVector vec_test = Astronomy.VectorFromSphere(sphere, time);

                    // Calculate EQD coordinates of the Moon. Verify against input file.
                    AstroVector vec_eqj = Astronomy.GeoMoon(time);
                    RotationMatrix rot = Astronomy.Rotation_EQJ_EQD(time);
                    AstroVector vec_eqd = Astronomy.RotateVector(rot, vec_eqj);
                    double angle = Astronomy.AngleBetween(vec_test, vec_eqd);
                    double diff_angle = 60.0 * abs(angle);
                    if (diff_angle > max_angle)
                        max_angle = diff_angle;
                    if (diff_angle > 1.54)
                    {
                        Console.WriteLine($"C# MoonNodesTest({filename} line {lnum}): EXCESSIVE equatorial error = {diff_angle:F3} arcmin.");
                        return 1;
                    }

                    // Test the Astronomy Engine moon node searcher.
                    if (lnum == 1)
                    {
                        // The very first time, so search for the first node in the series.
                        // Back up a few days to make sure we really are finding it ourselves.
                        AstroTime earlier = time.AddDays(-6.5472);    // back up by a weird amount of time
                        node = Astronomy.SearchMoonNode(earlier);
                    }
                    else
                    {
                        // Use the previous node to find the next node.
                        node = Astronomy.NextMoonNode(node);
                    }

                    // Verify the ecliptic latitude is very close to zero at the alleged node.
                    Spherical ecl = Astronomy.EclipticGeoMoon(node.time);
                    double diff_lat = 60.0 * abs(ecl.lat);
                    if (diff_lat > 8.1e-4)
                    {
                        Console.WriteLine($"C# MoonNodesTest({filename} line {lnum}): found node has excessive latitude = {diff_lat:F4} arcmin");
                        return 1;
                    }

                    // Verify the time agrees with Espenak's time to within a few minutes.
                    double diff_minutes = (24.0 * 60.0) * abs(node.time.tt - time.tt);
                    if (diff_minutes > max_minutes)
                        max_minutes = diff_minutes;

                    // Verify the kind of node matches what Espenak says (ascending or descending).
                    if (kind == "A" && node.kind != NodeEventKind.Ascending)
                    {
                        Console.WriteLine($"C# MoonNodesTest({filename} line {lnum}): did not find ascending node as expected.");
                        return 1;
                    }

                    if (kind == "D" && node.kind != NodeEventKind.Descending)
                    {
                        Console.WriteLine($"C# MoonNodesTest({filename} line {lnum}): did not find descending node as expected.");
                        return 1;
                    }

                    // Prepare for the next iteration.
                    prev_kind = kind;
                }

                if (max_minutes > 3.681)
                {
                    Console.WriteLine($"C# MoonNodesTest: EXCESSIVE time prediction error = {max_minutes:F3} minutes.");
                    return 1;
                }

                Console.WriteLine($"C# MoonNodesTest: PASS ({lnum} nodes, max equ error = {max_angle:F3}, max time error = {max_minutes:F3} minutes).");
                return 0;
            }
        }

        static int MoonReverseTest()
        {
            return (
                MoonReverse(  0.0) == 0 &&
                MoonReverse( 90.0) == 0 &&
                MoonReverse(180.0) == 0 &&
                MoonReverse(270.0) == 0
            ) ? 0 : 1;
        }

        static int MoonReverse(double longitude)
        {
            // Verify that SearchMoonPhase works both forward and backward in time.

            const int nphases = 5000;
            var utList = new double[nphases];
            double dtMin = +1000.0;
            double dtMax = -1000.0;
            double diff;

            // Search forward in time from 1800 to find consecutive phase events.
            var time = new AstroTime(1800, 1, 1, 0, 0, 0);
            for (int i = 0; i < nphases; ++i)
            {
                AstroTime result = Astronomy.SearchMoonPhase(longitude, time, +40.0);
                if (result == null)
                {
                    Console.WriteLine($"C# MoonReverse(i={i}): failed to find phase {longitude} after {time}");
                    return 1;
                }
                utList[i] = result.ut;
                if (i > 0)
                {
                    // Verify that consecutive events are reasonably close to the synodic period (29.5 days) apart.
                    double dt = v(utList[i] - utList[i-1]);
                    if (dt < dtMin) dtMin = dt;
                    if (dt > dtMax) dtMax = dt;
                }
                time = result.AddDays(+0.1);
            }

            Debug($"C# MoonReverse({longitude}): dtMin={dtMin:F6} days, dtMax={dtMax:F6} days.");
            if (dtMin < 29.175 || dtMax > 29.926)
            {
                Console.WriteLine($"C# MoonReverse({longitude}): Time between consecutive phases is suspicious.");
                return 1;
            }

            // Do a reverse chronological search and make sure the results are consistent with the forward search.
            time = time.AddDays(20.0);
            double maxDiff = 0.0;
            for (int i = nphases-1; i >= 0; --i)
            {
                AstroTime result = Astronomy.SearchMoonPhase(longitude, time, -40.0);
                if (result == null)
                {
                    Console.WriteLine($"C# MoonReverse(i={i}): failed to find phase {longitude} before {time}");
                    return 1;
                }
                diff = SECONDS_PER_DAY * abs(result.ut - utList[i]);
                if (diff > maxDiff) maxDiff = diff;
                time = result.AddDays(-0.1);
            }

            Debug($"C# MoonReverse({longitude}): Maximum discrepancy in reverse search = {maxDiff:F6} seconds.");
            if (maxDiff > 0.164)
            {
                Console.WriteLine($"C# MoonReverse({longitude}): EXCESSIVE DISCREPANCY in reverse search.");
                return 1;
            }

            // Pick a pair of consecutive events from the middle of the list.
            // Verify forward and backward searches work correctly from many intermediate times.
            const int nslots = 100;
            int k = nphases / 2;
            double ut1 = utList[k];
            double ut2 = utList[k+1];
            for (int i = 1; i < nslots; ++i)
            {
                double ut = ut1 + ((double)i/nslots)*(ut2 - ut1);
                time = new AstroTime(ut);

                AstroTime before = Astronomy.SearchMoonPhase(longitude, time, -40.0);
                if (before == null)
                {
                    Console.WriteLine($"C# MoonReverse({longitude}): backward search from {time} failed.");
                    return 1;
                }
                diff = SECONDS_PER_DAY * abs(before.ut - ut1);
                if (diff > 0.07)
                {
                    Console.WriteLine($"C# MoonReverse({longitude}): backward search error = {diff:E4} seconds from {time}.");
                    return 1;
                }

                AstroTime after = Astronomy.SearchMoonPhase(longitude, time, +40.0);
                if (after == null)
                {
                    Console.WriteLine($"C# MoonReverse({longitude}): forward search from {time} failed.");
                    return 1;
                }
                diff = SECONDS_PER_DAY * abs(after.ut - ut2);
                if (diff > 0.07)
                {
                    Console.WriteLine($"C# MoonReverse({longitude}): forward search error = {diff:E4} seconds from {time}.");
                    return 1;
                }
            }

            Console.WriteLine($"C# MoonReverse({longitude}): PASS");
            return 0;
        }

        static int RiseSetTest()
        {
            const string filename = "../../riseset/riseset.txt";
            using (StreamReader infile = File.OpenText(filename))
            {
                int lnum = 0;
                string line;
                var re = new Regex(@"^([A-Za-z]+)\s+([\-\+]?\d+\.?\d*)\s+([\-\+]?\d+\.?\d*)\s+(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z\s+([rs])\s*$");
                Body current_body = Body.Invalid;
                Observer observer = new Observer();
                bool foundObserver = false;
                AstroTime r_search_date = null, s_search_date = null;
                AstroTime r_evt = null, s_evt = null;     // rise event, set event: search results
                AstroTime a_evt = null, b_evt = null;     // chronologically first and second events
                Direction a_dir = Direction.Rise, b_dir = Direction.Rise;
                const double nudge_days = 0.01;
                double sum_minutes = 0.0;
                double max_minutes = 0.0;

                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;

                    // Moon  103 -61 1944-01-02T17:08Z s
                    // Moon  103 -61 1944-01-03T05:47Z r
                    Match m = re.Match(line);
                    if (!m.Success)
                    {
                        Console.WriteLine("C# RiseSetTest({0} line {1}): invalid input format", filename, lnum);
                        return 1;
                    }
                    Body body = Enum.Parse<Body>(m.Groups[1].Value);
                    double longitude = double.Parse(m.Groups[2].Value);
                    double latitude = double.Parse(m.Groups[3].Value);
                    int year = int.Parse(m.Groups[4].Value);
                    int month = int.Parse(m.Groups[5].Value);
                    int day = int.Parse(m.Groups[6].Value);
                    int hour = int.Parse(m.Groups[7].Value);
                    int minute = int.Parse(m.Groups[8].Value);
                    Direction direction = (m.Groups[9].Value == "r") ? Direction.Rise : Direction.Set;
                    var correct_date = new AstroTime(year, month, day, hour, minute, 0);

                    // Every time we see a new geographic location or body, start a new iteration
                    // of finding all rise/set times for that UTC calendar year.
                    if (!foundObserver || observer.latitude != latitude || observer.longitude != longitude || current_body != body)
                    {
                        current_body = body;
                        observer = new Observer(latitude, longitude, 0.0);
                        foundObserver = true;
                        r_search_date = s_search_date = new AstroTime(year, 1, 1, 0, 0, 0);
                        b_evt = null;
                        Debug("C# RiseSetTest: {0} lat={1} lon={2}", body, latitude, longitude);
                    }

                    if (b_evt != null)
                    {
                        // The previous iteration found two events.
                        // We already processed the earlier event (a_evt).
                        // Now it is time to process the later event (b_evt).
                        a_evt = b_evt;
                        a_dir = b_dir;
                        b_evt = null;
                    }
                    else
                    {
                        r_evt = Astronomy.SearchRiseSet(body, observer, Direction.Rise, r_search_date, 366.0);
                        if (r_evt == null)
                        {
                            Console.WriteLine("C# RiseSetTest({0} line {1}): Did not find {2} rise event.", filename, lnum, body);
                            return 1;
                        }

                        s_evt = Astronomy.SearchRiseSet(body, observer, Direction.Set, s_search_date, 366.0);
                        if (s_evt == null)
                        {
                            Console.WriteLine("C# RiseSetTest({0} line {1}): Did not find {2} set event.", filename, lnum, body);
                            return 1;
                        }

                        // Sort the two events chronologically.
                        // We will check the earlier event in this iteration,
                        // and check the later event in the next iteration.
                        if (r_evt.tt < s_evt.tt)
                        {
                            a_evt = r_evt;
                            b_evt = s_evt;
                            a_dir = Direction.Rise;
                            b_dir = Direction.Set;
                        }
                        else
                        {
                            a_evt = s_evt;
                            b_evt = r_evt;
                            a_dir = Direction.Set;
                            b_dir = Direction.Rise;
                        }

                        // Nudge the event times forward a tiny amount.
                        // This prevents us from getting stuck in a loop, finding the same event repeatedly.
                        r_search_date = r_evt.AddDays(nudge_days);
                        s_search_date = s_evt.AddDays(nudge_days);
                    }

                    // Expect the current search result to match the earlier of the found dates.

                    if (a_dir != direction)
                    {
                        Console.WriteLine("C# RiseSetTest({0} line {1}): expected dir={2} but found {3}", filename, lnum, a_dir, direction);
                        return 1;
                    }

                    double error_minutes = (24.0 * 60.0) * abs(a_evt.tt - correct_date.tt);
                    sum_minutes += error_minutes * error_minutes;
                    if (error_minutes > max_minutes)
                        max_minutes = error_minutes;

                    if (error_minutes > 0.57)
                    {
                        Console.WriteLine("C# RiseSetTest({0} line {1}): excessive prediction time error = {2} minutes.", filename, lnum, error_minutes);
                        return 1;
                    }
                }

                double rms_minutes = sqrt(sum_minutes / lnum);
                Console.WriteLine("C# RiseSetTest: passed {0} lines: time errors in minutes: rms={1}, max={2}", lnum, rms_minutes, max_minutes);
                return 0;
            }
        }

        static Direction Toggle(Direction dir)
        {
            switch (dir)
            {
            case Direction.Rise: return Direction.Set;
            case Direction.Set:  return Direction.Rise;
            default: throw new ArgumentException($"Invalid direction: {dir}");
            }
        }

        static int RiseSetSlot(double ut1, double ut2, Direction dir, Observer observer)
        {
            const int nslots = 100;
            double maxDiff = 0.0;
            AstroTime time, result;
            for (int i = 1; i < nslots; ++i)
            {
                double ut = ut1 + ((double)i / nslots)*(ut2 - ut1);
                time = new AstroTime(ut);

                result = Astronomy.SearchRiseSet(Body.Sun, observer, dir, time, -1.0);
                if (result == null)
                {
                    Console.WriteLine($"C# RiseSetSlot({dir}): backward slot search failed before {time}");
                    return 1;
                }
                double diff = SECONDS_PER_DAY * abs(result.ut - ut1);
                if (diff > maxDiff) maxDiff = diff;

                result = Astronomy.SearchRiseSet(Body.Sun, observer, dir, time, +1.0);
                if (result == null)
                {
                    Console.WriteLine($"C# RiseSetSlot({dir}): forward slot search failed after {time}");
                    return 1;
                }
                diff = SECONDS_PER_DAY * abs(result.ut - ut2);
                if (diff > maxDiff) maxDiff = diff;
            }
            if (maxDiff > 0.9)
            {
                Console.WriteLine($"C# RiseSetSlot({dir}): EXCESSIVE slot-test discrepancy = {maxDiff:F6} seconds.");
                return 1;
            }
            Debug($"C# RiseSetSlot({dir}): slot-test discrepancy = {maxDiff:F6} seconds.");
            return 0;
        }

        static int RiseSetReverseTest()
        {
            // Verify that the rise/set search works equally well forwards and backwards in time.
            const int nsamples = 5000;
            const double nudge = 0.1;
            var utList = new double[nsamples];
            var observer = new Observer(30.5, -90.7, 0.0);
            double dtMin = +1000.0;
            double dtMax = -1000.0;
            double maxDiff = 0.0;
            AstroTime result;

            // Find alternating sunrise/sunset events in forward chronological order.
            Direction dir = Direction.Rise;
            var time = new AstroTime(2022, 1, 1, 0, 0, 0);
            for (int i = 0; i < nsamples; ++i)
            {
                result = Astronomy.SearchRiseSet(Body.Sun, observer, dir, time, +1.0);
                if (result == null)
                {
                    Console.WriteLine($"C# RiseSetReverseTest: cannot find {dir} event after {time}.");
                    return 1;
                }
                utList[i] = result.ut;
                if (i > 0)
                {
                    // Check the time between consecutive sunrise/sunset events.
                    // These will vary considerably with the seasons, so just make sure we don't miss any entirely.
                    double dt = v(utList[i] - utList[i-1]);
                    if (dt < dtMin) dtMin = dt;
                    if (dt > dtMax) dtMax = dt;
                }
                dir = Toggle(dir);
                time = result.AddDays(+nudge);
            }

            Debug($"C# RiseSetReverse: dtMin={dtMin:F6} days, dtMax={dtMax:F6} days.");
            if (dtMin < 0.411 || dtMax > 0.589)
            {
                Console.WriteLine($"C# RiseSetReverse: Invalid intervals between sunrise/sunset.");
                return 1;
            }

            // Perform the same search in reverse. Verify we get consistent rise/set times.
            for (int i = nsamples-1; i >= 0; --i)
            {
                dir = Toggle(dir);
                result = Astronomy.SearchRiseSet(Body.Sun, observer, dir, time, -1.0);
                if (result == null)
                {
                    Console.WriteLine($"C# RiseSetReverseTest: cannot find {dir} event before {time}.");
                    return 1;
                }
                double diff = SECONDS_PER_DAY * abs(utList[i] - result.ut);
                if (diff > maxDiff) maxDiff = diff;
                time = result.AddDays(-nudge);
            }

            if (maxDiff > 0.982)
            {
                Console.WriteLine($"C# RiseSetReverse: EXCESSIVE forward/backward discrepancy = {maxDiff:F6} seconds.");
                return 1;
            }
            Debug($"C# RiseSetReverse: forward/backward discrepancy = {maxDiff:F6} seconds.");

            // All even indexes in utList hold sunrise times.
            // All odd indexes in utList hold sunset times.
            // Verify that forward/backward searches for consecutive sunrises/sunsets
            // resolve correctly for 100 time slots between them.
            int k = (nsamples / 2) & ~1;
            if (0 != RiseSetSlot(utList[k+0], utList[k+2], Direction.Rise, observer)) return 1;
            if (0 != RiseSetSlot(utList[k+1], utList[k+3], Direction.Set,  observer)) return 1;

            Console.WriteLine("C# RiseSetReverse: PASS");
            return 0;
        }

        static int TestElongFile(string filename, double targetRelLon)
        {
            using (StreamReader infile = File.OpenText(filename))
            {
                int lnum = 0;
                string line;
                var re = new Regex(@"^(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z\s+([A-Z][a-z]+)\s*$");
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    // 2018-05-09T00:28Z Jupiter
                    Match m = re.Match(line);
                    if (!m.Success)
                    {
                        Console.WriteLine("C# TestElongFile({0} line {1}): invalid data format.", filename, lnum);
                        return 1;
                    }
                    int year = int.Parse(m.Groups[1].Value);
                    int month = int.Parse(m.Groups[2].Value);
                    int day = int.Parse(m.Groups[3].Value);
                    int hour = int.Parse(m.Groups[4].Value);
                    int minute = int.Parse(m.Groups[5].Value);
                    Body body = Enum.Parse<Body>(m.Groups[6].Value);
                    var search_date = new AstroTime(year, 1, 1, 0, 0, 0);
                    var expected_time = new AstroTime(year, month, day, hour, minute, 0);
                    AstroTime search_result = Astronomy.SearchRelativeLongitude(body, targetRelLon, search_date);
                    if (search_result == null)
                    {
                        // This function should NEVER return null.
                        Console.WriteLine("C# TestElongFile({0} line {1}): SearchRelativeLongitude returned null.", filename, lnum);
                        return 1;
                    }
                    double diff_minutes = (24.0 * 60.0) * (search_result.tt - expected_time.tt);
                    Console.WriteLine("{0} error = {1} minutes.", body, diff_minutes.ToString("f3"));
                    if (abs(diff_minutes) > 6.8)
                    {
                        Console.WriteLine("C# TestElongFile({0} line {1}): EXCESSIVE ERROR.", filename, lnum);
                        return 1;
                    }
                }
                Console.WriteLine("C# TestElongFile: passed {0} rows of data.", lnum);
                return 0;
            }
        }

        static int TestPlanetLongitudes(Body body, string outFileName, string zeroLonEventName)
        {
            const int startYear = 1700;
            const int stopYear = 2200;
            int count = 0;
            double rlon = 0.0;
            double min_diff = 1.0e+99;
            double max_diff = 1.0e+99;

            using (StreamWriter outfile = File.CreateText(outFileName))
            {
                var time = new AstroTime(startYear, 1, 1, 0, 0, 0);
                var stopTime = new AstroTime(stopYear, 1, 1, 0, 0, 0);
                while (time.tt < stopTime.tt)
                {
                    ++count;
                    string event_name = (rlon == 0.0) ? zeroLonEventName : "sup";
                    AstroTime search_result = Astronomy.SearchRelativeLongitude(body, rlon, time);
                    if (search_result == null)
                    {
                        Console.WriteLine("C# TestPlanetLongitudes({0}): SearchRelativeLongitude returned null.", body);
                        return 1;
                    }

                    if (count >= 2)
                    {
                        // Check for consistent intervals.
                        // Mainly I don't want to skip over an event!
                        double day_diff = search_result.tt - time.tt;
                        if (count == 2)
                        {
                            min_diff = max_diff = day_diff;
                        }
                        else
                        {
                            if (day_diff < min_diff)
                                min_diff = day_diff;

                            if (day_diff > max_diff)
                                max_diff = day_diff;
                        }
                    }

                    AstroVector geo = Astronomy.GeoVector(body, search_result, Aberration.Corrected);
                    double dist = geo.Length();
                    outfile.WriteLine("e {0} {1} {2} {3}", body, event_name, search_result.tt.ToString("G18"), dist.ToString("G18"));

                    // Search for the opposite longitude event next time.
                    time = search_result;
                    rlon = 180.0 - rlon;
                }
            }

            double thresh;
            switch (body)
            {
                case Body.Mercury:  thresh = 1.65;  break;
                case Body.Mars:     thresh = 1.30;  break;
                default:            thresh = 1.07;  break;
            }

            double ratio = max_diff / min_diff;
            Debug("C# TestPlanetLongitudes({0,7}): {1,5} events, ratio={2,5}, file: {3}", body, count, ratio.ToString("f3"), outFileName);

            if (ratio > thresh)
            {
                Console.WriteLine("C# TestPlanetLongitudes({0}): excessive event interval ratio.", body);
                return 1;
            }
            return 0;
        }

        static int ElongationTest()
        {
            if (0 != TestElongFile("../../longitude/opposition_2018.txt", 0.0)) return 1;
            if (0 != TestPlanetLongitudes(Body.Mercury, "csharp_longitude_Mercury.txt", "inf")) return 1;
            if (0 != TestPlanetLongitudes(Body.Venus,   "csharp_longitude_Venus.txt",   "inf")) return 1;
            if (0 != TestPlanetLongitudes(Body.Mars,    "csharp_longitude_Mars.txt",    "opp")) return 1;
            if (0 != TestPlanetLongitudes(Body.Jupiter, "csharp_longitude_Jupiter.txt", "opp")) return 1;
            if (0 != TestPlanetLongitudes(Body.Saturn,  "csharp_longitude_Saturn.txt",  "opp")) return 1;
            if (0 != TestPlanetLongitudes(Body.Uranus,  "csharp_longitude_Uranus.txt",  "opp")) return 1;
            if (0 != TestPlanetLongitudes(Body.Neptune, "csharp_longitude_Neptune.txt", "opp")) return 1;
            if (0 != TestPlanetLongitudes(Body.Pluto,   "csharp_longitude_Pluto.txt",   "opp")) return 1;

            foreach (elong_test_t et in ElongTestData)
                if (0 != TestMaxElong(et))
                    return 1;

            Console.WriteLine("C# ElongationTest: PASS");
            return 0;
        }

        static readonly Regex regexDate = new Regex(@"^(\d+)-(\d+)-(\d+)T(\d+):(\d+)(:(\d+))?Z$");

        static AstroTime ParseDate(string text)
        {
            Match m = regexDate.Match(text);
            if (!m.Success)
                throw new Exception(string.Format("ParseDate failed for string: '{0}'", text));
            int year = int.Parse(m.Groups[1].Value);
            int month = int.Parse(m.Groups[2].Value);
            int day = int.Parse(m.Groups[3].Value);
            int hour = int.Parse(m.Groups[4].Value);
            int minute = int.Parse(m.Groups[5].Value);
            int second = 0;
            if (!string.IsNullOrEmpty(m.Groups[7].Value))
                second = int.Parse(m.Groups[7].Value);
            return new AstroTime(year, month, day, hour, minute, second);
        }

        static AstroTime OptionalParseDate(string text)
        {
            if (text == "-")
                return null;

            return ParseDate(text);
        }

        static int TestMaxElong(elong_test_t test)
        {
            AstroTime searchTime = ParseDate(test.searchDate);
            AstroTime eventTime = ParseDate(test.eventDate);
            ElongationInfo evt = Astronomy.SearchMaxElongation(test.body, searchTime);
            double hour_diff = 24.0 * abs(evt.time.tt - eventTime.tt);
            double arcmin_diff = 60.0 * abs(evt.elongation - test.angle);
            Debug("C# TestMaxElong: {0,7} {1,7} elong={2,5} ({3} arcmin, {4} hours)", test.body, test.visibility, evt.elongation, arcmin_diff, hour_diff);
            if (hour_diff > 0.6)
            {
                Console.WriteLine("C# TestMaxElong({0} {1}): excessive hour error.", test.body, test.searchDate);
                return 1;
            }

            if (arcmin_diff > 3.4)
            {
                Console.WriteLine("C# TestMaxElong({0} {1}): excessive arcmin error.", test.body, test.searchDate);
                return 1;
            }

            return 0;
        }

        struct elong_test_t
        {
            public Body body;
            public string searchDate;
            public string eventDate;
            public double angle;
            public Visibility visibility;

            public elong_test_t(Body body, string searchDate, string eventDate, double angle, Visibility visibility)
            {
                this.body = body;
                this.searchDate = searchDate;
                this.eventDate = eventDate;
                this.angle = angle;
                this.visibility = visibility;
            }
        }

        static readonly elong_test_t[] ElongTestData = new elong_test_t[]
        {
            // Max elongation data obtained from:
            // http://www.skycaramba.com/greatest_elongations.shtml
            new elong_test_t( Body.Mercury, "2010-01-17T05:22Z", "2010-01-27T05:22Z", 24.80, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2010-05-16T02:15Z", "2010-05-26T02:15Z", 25.10, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2010-09-09T17:24Z", "2010-09-19T17:24Z", 17.90, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2010-12-30T14:33Z", "2011-01-09T14:33Z", 23.30, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2011-04-27T19:03Z", "2011-05-07T19:03Z", 26.60, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2011-08-24T05:52Z", "2011-09-03T05:52Z", 18.10, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2011-12-13T02:56Z", "2011-12-23T02:56Z", 21.80, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2012-04-08T17:22Z", "2012-04-18T17:22Z", 27.50, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2012-08-06T12:04Z", "2012-08-16T12:04Z", 18.70, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2012-11-24T22:55Z", "2012-12-04T22:55Z", 20.60, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2013-03-21T22:02Z", "2013-03-31T22:02Z", 27.80, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2013-07-20T08:51Z", "2013-07-30T08:51Z", 19.60, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2013-11-08T02:28Z", "2013-11-18T02:28Z", 19.50, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2014-03-04T06:38Z", "2014-03-14T06:38Z", 27.60, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2014-07-02T18:22Z", "2014-07-12T18:22Z", 20.90, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2014-10-22T12:36Z", "2014-11-01T12:36Z", 18.70, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2015-02-14T16:20Z", "2015-02-24T16:20Z", 26.70, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2015-06-14T17:10Z", "2015-06-24T17:10Z", 22.50, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2015-10-06T03:20Z", "2015-10-16T03:20Z", 18.10, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2016-01-28T01:22Z", "2016-02-07T01:22Z", 25.60, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2016-05-26T08:45Z", "2016-06-05T08:45Z", 24.20, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2016-09-18T19:27Z", "2016-09-28T19:27Z", 17.90, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2017-01-09T09:42Z", "2017-01-19T09:42Z", 24.10, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2017-05-07T23:19Z", "2017-05-17T23:19Z", 25.80, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2017-09-02T10:14Z", "2017-09-12T10:14Z", 17.90, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2017-12-22T19:48Z", "2018-01-01T19:48Z", 22.70, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2018-04-19T18:17Z", "2018-04-29T18:17Z", 27.00, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2018-08-16T20:35Z", "2018-08-26T20:35Z", 18.30, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2018-12-05T11:34Z", "2018-12-15T11:34Z", 21.30, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2019-04-01T19:40Z", "2019-04-11T19:40Z", 27.70, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2019-07-30T23:08Z", "2019-08-09T23:08Z", 19.00, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2019-11-18T10:31Z", "2019-11-28T10:31Z", 20.10, Visibility.Morning ),
            new elong_test_t( Body.Mercury, "2010-03-29T23:32Z", "2010-04-08T23:32Z", 19.40, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2010-07-28T01:03Z", "2010-08-07T01:03Z", 27.40, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2010-11-21T15:42Z", "2010-12-01T15:42Z", 21.50, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2011-03-13T01:07Z", "2011-03-23T01:07Z", 18.60, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2011-07-10T04:56Z", "2011-07-20T04:56Z", 26.80, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2011-11-04T08:40Z", "2011-11-14T08:40Z", 22.70, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2012-02-24T09:39Z", "2012-03-05T09:39Z", 18.20, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2012-06-21T02:00Z", "2012-07-01T02:00Z", 25.70, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2012-10-16T21:59Z", "2012-10-26T21:59Z", 24.10, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2013-02-06T21:24Z", "2013-02-16T21:24Z", 18.10, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2013-06-02T16:45Z", "2013-06-12T16:45Z", 24.30, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2013-09-29T09:59Z", "2013-10-09T09:59Z", 25.30, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2014-01-21T10:00Z", "2014-01-31T10:00Z", 18.40, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2014-05-15T07:06Z", "2014-05-25T07:06Z", 22.70, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2014-09-11T22:20Z", "2014-09-21T22:20Z", 26.40, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2015-01-04T20:26Z", "2015-01-14T20:26Z", 18.90, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2015-04-27T04:46Z", "2015-05-07T04:46Z", 21.20, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2015-08-25T10:20Z", "2015-09-04T10:20Z", 27.10, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2015-12-19T03:11Z", "2015-12-29T03:11Z", 19.70, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2016-04-08T14:00Z", "2016-04-18T14:00Z", 19.90, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2016-08-06T21:24Z", "2016-08-16T21:24Z", 27.40, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2016-12-01T04:36Z", "2016-12-11T04:36Z", 20.80, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2017-03-22T10:24Z", "2017-04-01T10:24Z", 19.00, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2017-07-20T04:34Z", "2017-07-30T04:34Z", 27.20, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2017-11-14T00:32Z", "2017-11-24T00:32Z", 22.00, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2018-03-05T15:07Z", "2018-03-15T15:07Z", 18.40, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2018-07-02T05:24Z", "2018-07-12T05:24Z", 26.40, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2018-10-27T15:25Z", "2018-11-06T15:25Z", 23.30, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2019-02-17T01:23Z", "2019-02-27T01:23Z", 18.10, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2019-06-13T23:14Z", "2019-06-23T23:14Z", 25.20, Visibility.Evening ),
            new elong_test_t( Body.Mercury, "2019-10-10T04:00Z", "2019-10-20T04:00Z", 24.60, Visibility.Evening ),
            new elong_test_t( Body.Venus,   "2010-12-29T15:57Z", "2011-01-08T15:57Z", 47.00, Visibility.Morning ),
            new elong_test_t( Body.Venus,   "2012-08-05T08:59Z", "2012-08-15T08:59Z", 45.80, Visibility.Morning ),
            new elong_test_t( Body.Venus,   "2014-03-12T19:25Z", "2014-03-22T19:25Z", 46.60, Visibility.Morning ),
            new elong_test_t( Body.Venus,   "2015-10-16T06:57Z", "2015-10-26T06:57Z", 46.40, Visibility.Morning ),
            new elong_test_t( Body.Venus,   "2017-05-24T13:09Z", "2017-06-03T13:09Z", 45.90, Visibility.Morning ),
            new elong_test_t( Body.Venus,   "2018-12-27T04:24Z", "2019-01-06T04:24Z", 47.00, Visibility.Morning ),
            new elong_test_t( Body.Venus,   "2010-08-10T03:19Z", "2010-08-20T03:19Z", 46.00, Visibility.Evening ),
            new elong_test_t( Body.Venus,   "2012-03-17T08:03Z", "2012-03-27T08:03Z", 46.00, Visibility.Evening ),
            new elong_test_t( Body.Venus,   "2013-10-22T08:00Z", "2013-11-01T08:00Z", 47.10, Visibility.Evening ),
            new elong_test_t( Body.Venus,   "2015-05-27T18:46Z", "2015-06-06T18:46Z", 45.40, Visibility.Evening ),
            new elong_test_t( Body.Venus,   "2017-01-02T13:19Z", "2017-01-12T13:19Z", 47.10, Visibility.Evening ),
            new elong_test_t( Body.Venus,   "2018-08-07T17:02Z", "2018-08-17T17:02Z", 45.90, Visibility.Evening )
        };

        static double PlanetOrbitalPeriod(Body body)
        {
            switch (body)
            {
            case Body.Mercury:  return     87.969;
            case Body.Venus:    return    224.701;
            case Body.Earth:    return    365.256;
            case Body.Mars:     return    686.980;
            case Body.Jupiter:  return   4332.589;
            case Body.Saturn:   return  10759.22;
            case Body.Uranus:   return  30685.4;
            case Body.Neptune:  return  60189.0;
            case Body.Pluto:    return  90560.0;
            default:
                throw new ArgumentException(string.Format("Invalid body {0}", body));
            }
        }

        static int PlanetApsisTest()
        {
            const double degree_threshold = 0.1;
            const string testDataPath = "../../apsides";
            var start_time = new AstroTime(1700, 1, 1, 0, 0, 0);
            for (Body body = Body.Mercury; body <= Body.Pluto; ++body)
            {
                double period = PlanetOrbitalPeriod(body);
                double max_dist_ratio = 0.0;
                double max_diff_days = 0.0;
                double min_interval = -1.0;
                double max_interval = -1.0;
                ApsisInfo apsis = Astronomy.SearchPlanetApsis(body, start_time);
                int count = 1;

                string filename = Path.Combine(testDataPath, string.Format("apsis_{0}.txt", (int)body));
                using (StreamReader infile = File.OpenText(filename))
                {
                    string line;
                    while (null != (line = infile.ReadLine()))
                    {
                        // Parse the line of test data.
                        string[] token = Tokenize(line);
                        if (token.Length != 3)
                        {
                            Console.WriteLine("C# PlanetApsisTest({0} line {1}): Invalid data format: {2} tokens", filename, count, token.Length);
                            return 1;
                        }
                        int expected_kind = int.Parse(token[0]);
                        AstroTime expected_time = ParseDate(token[1]);
                        double expected_distance = double.Parse(token[2]);

                        // Compare computed values against expected values.
                        if ((int)apsis.kind != expected_kind)
                        {
                            Console.WriteLine("C# PlanetApsisTest({0} line {1}): WRONG APSIS KIND", filename, count);
                            return 1;
                        }

                        double diff_days = abs(expected_time.tt - apsis.time.tt);
                        max_diff_days = max(max_diff_days, diff_days);
                        double diff_degrees = (diff_days / period) * 360.0;
                        if (diff_degrees > degree_threshold)
                        {
                            Console.WriteLine("C# PlanetApsis: FAIL - {0} exceeded angular threshold ({1} vs {2} degrees)", body, diff_degrees, degree_threshold);
                            return 1;
                        }

                        double diff_dist_ratio = abs(expected_distance - apsis.dist_au) / expected_distance;
                        max_dist_ratio = max(max_dist_ratio, diff_dist_ratio);
                        if (diff_dist_ratio > 1.05e-4)
                        {
                            Console.WriteLine("C# PlanetApsisTest({0} line {1}): distance ratio {2} is too large.", filename, count, diff_dist_ratio);
                            return 1;
                        }

                        // Calculate the next apsis.
                        AstroTime prev_time = apsis.time;
                        apsis = Astronomy.NextPlanetApsis(body, apsis);

                        // Update statistics.
                        ++count;
                        double interval = apsis.time.tt - prev_time.tt;
                        if (min_interval < 0.0)
                        {
                            min_interval = max_interval = interval;
                        }
                        else
                        {
                            min_interval = min(min_interval, interval);
                            max_interval = max(max_interval, interval);
                        }
                    }
                }

                if (count < 2)
                {
                    Console.WriteLine("C# PlanetApsis: FAILED to find apsides for {0}", body);
                    return 1;
                }

                Debug("C# PlanetApsis: {0} apsides for {1,-9} -- intervals: min={2:0.00}, max={3:0.00}, ratio={4:0.000000}; max day={5}, degrees={6:0.000}, dist ratio={7}",
                    count, body,
                    min_interval, max_interval, max_interval / min_interval,
                    max_diff_days,
                    (max_diff_days / period) * 360.0,
                    max_dist_ratio);
            }

            Console.WriteLine("C# PlanetApsis: PASS");
            return 0;
        }

        static int LunarApsisTest()
        {
            const string inFileName = "../../apsides/moon.txt";
            using (StreamReader infile = File.OpenText(inFileName))
            {
                int lnum = 0;
                string line;
                var start_time = new AstroTime(2001,1, 1, 0, 0, 0);
                ApsisInfo apsis = new ApsisInfo();
                double max_minutes = 0.0;
                double max_km = 0.0;
                // 0 2001-01-10T08:59Z 357132
                // 1 2001-01-24T19:02Z 406565
                var regex = new Regex(@"^\s*([01])\s+(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z\s+(\d+)\s*$");
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    Match m = regex.Match(line);
                    if (!m.Success)
                    {
                        Console.WriteLine("C# LunarApsisTest({0} line {1}): invalid data format.", inFileName, lnum);
                        return 1;
                    }
                    ApsisKind kind = (m.Groups[1].Value == "0") ? ApsisKind.Pericenter : ApsisKind.Apocenter;

                    int year = int.Parse(m.Groups[2].Value);
                    int month = int.Parse(m.Groups[3].Value);
                    int day = int.Parse(m.Groups[4].Value);
                    int hour = int.Parse(m.Groups[5].Value);
                    int minute = int.Parse(m.Groups[6].Value);
                    double dist_km = double.Parse(m.Groups[7].Value);

                    var correct_time = new AstroTime(year, month, day, hour, minute, 0);

                    if (lnum == 1)
                        apsis = Astronomy.SearchLunarApsis(start_time);
                    else
                        apsis = Astronomy.NextLunarApsis(apsis);

                    if (kind != apsis.kind)
                    {
                        Console.WriteLine("C# LunarApsisTest({0} line {1}): expected apsis kind {2} but found {3}", inFileName, lnum, kind, apsis.kind);
                        return 1;
                    }
                    double diff_minutes = (24.0 * 60.0) * abs(apsis.time.ut - correct_time.ut);
                    if (diff_minutes > 35.0)
                    {
                        Console.WriteLine("C# LunarApsisTest({0} line {1}): excessive time error: {2} minutes", inFileName, lnum, diff_minutes);
                        return 1;
                    }
                    double diff_km =  abs(apsis.dist_km - dist_km);
                    if (diff_km > 25.0)
                    {
                        Console.WriteLine("C# LunarApsisTest({0} line {1}): excessive distance error: {2} km", inFileName, lnum, diff_km);
                        return 1;
                    }

                    if (diff_minutes > max_minutes)
                        max_minutes = diff_minutes;

                    if (diff_km > max_km)
                        max_km = diff_km;
                }
                Console.WriteLine("C# LunarApsisTest: Found {0} events, max time error = {1} minutes, max distance error = {2} km.", lnum, max_minutes, max_km);
                return 0;
            }
        }

        class JplDateTime
        {
            public string Rest;
            public AstroTime Time;
        }

        static readonly Regex JplRegex = new Regex(@"^\s*(\d{4})-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-(\d{2})\s+(\d{2}):(\d{2})\s+(.*)");

        static int MonthNumber(string mtext)
        {
            switch (mtext)
            {
                case "Jan": return  1;
                case "Feb": return  2;
                case "Mar": return  3;
                case "Apr": return  4;
                case "May": return  5;
                case "Jun": return  6;
                case "Jul": return  7;
                case "Aug": return  8;
                case "Sep": return  9;
                case "Oct": return 10;
                case "Nov": return 11;
                case "Dec": return 12;
                default:
                    throw new Exception(string.Format("Internal error: unexpected month name '{0}'", mtext));
            }
        }

        static JplDateTime ParseJplHorizonsDateTime(string line)
        {
            Match m = JplRegex.Match(line);
            if (!m.Success)
                return null;
            int year = int.Parse(m.Groups[1].Value);
            int month = MonthNumber(m.Groups[2].Value);
            int day = int.Parse(m.Groups[3].Value);
            int hour = int.Parse(m.Groups[4].Value);
            int minute = int.Parse(m.Groups[5].Value);
            string rest = m.Groups[6].Value;
            AstroTime time = new AstroTime(year, month, day, hour, minute, 0);
            return new JplDateTime { Rest=rest, Time=time };
        }

        static int CheckMagnitudeData(Body body, string filename)
        {
            using (StreamReader infile = File.OpenText(filename))
            {
                const double limit = 0.012;
                double diff_lo = 0.0;
                double diff_hi = 0.0;
                double sum_squared_diff = 0.0;
                int lnum = 0;
                int count = 0;
                string line;
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    JplDateTime jpl = ParseJplHorizonsDateTime(line);
                    if (jpl == null)
                        continue;

                    string[] token = Tokenize(jpl.Rest);
                    if (token.Length > 0 && token[0] == "n.a.")
                        continue;

                    if (token.Length != 7)
                    {
                        Console.WriteLine("C# CheckMagnitudeData({0} line {1}): invalid data format", lnum, filename);
                        return 1;
                    }
                    double mag;
                    if (!double.TryParse(token[0], out mag))
                    {
                        Console.WriteLine("C# CheckMagnitudeData({0} line {1}): cannot parse number from '{2}'", filename, lnum, token[0]);
                        return 1;
                    }
                    var illum = Astronomy.Illumination(body, jpl.Time);
                    double diff = illum.mag - mag;
                    if (abs(diff) > limit)
                    {
                        Console.WriteLine("C# CheckMagnitudeData({0} line {1}): EXCESSIVE ERROR: correct mag={0}, calc mag={1}, diff={2}", mag, illum.mag, diff);
                        return 1;
                    }
                    sum_squared_diff += diff * diff;
                    if (count == 0)
                    {
                        diff_lo = diff_hi = diff;
                    }
                    else
                    {
                        if (diff < diff_lo)
                            diff_lo = diff;

                        if (diff > diff_hi)
                            diff_hi = diff;
                    }
                    ++count;
                }
                if (count == 0)
                {
                    Console.WriteLine("C# CheckMagnitudeData: Did not find any data in file: {0}", filename);
                    return 1;
                }
                double rms = sqrt(sum_squared_diff / count);
                Debug("C# CheckMagnitudeData: {0} {1} rows diff_lo={2} diff_hi={3} rms={4}", filename, count, diff_lo, diff_hi, rms);
                return 0;
            }
        }

        struct saturn_test_case
        {
            public readonly string date;
            public readonly double mag;
            public readonly double tilt;

            public saturn_test_case(string date, double mag, double tilt)
            {
                this.date = date;
                this.mag = mag;
                this.tilt = tilt;
            }
        }

        // JPL Horizons does not include Saturn's rings in its magnitude models.
        // I still don't have authoritative test data for Saturn's magnitude.
        // For now, I just test for consistency with Paul Schlyter's formulas at:
        // http://www.stjarnhimlen.se/comp/ppcomp.html#15
        static saturn_test_case[] saturn_data = new saturn_test_case[]
        {
            new saturn_test_case("1972-01-01T00:00Z", -0.31904865,  +24.50061220),
            new saturn_test_case("1980-01-01T00:00Z", +0.85213663,   -1.85761461),
            new saturn_test_case("2009-09-04T00:00Z", +1.01626809,   +0.08380716),
            new saturn_test_case("2017-06-15T00:00Z", -0.12318790,  -26.60871409),
            new saturn_test_case("2019-05-01T00:00Z", +0.32954097,  -23.53880802),
            new saturn_test_case("2025-09-25T00:00Z", +0.51286575,   +1.52327932),
            new saturn_test_case("2032-05-15T00:00Z", -0.04652109,  +26.95717765)
        };

        static int CheckSaturn()
        {
            int error = 0;

            foreach (saturn_test_case data in saturn_data)
            {
                AstroTime time = ParseDate(data.date);

                IllumInfo illum = Astronomy.Illumination(Body.Saturn, time);
                Debug("C# Saturn: date={0}  calc mag={1}  ring_tilt={2}", data.date, illum.mag, illum.ring_tilt);

                double mag_diff = abs(illum.mag - data.mag);
                if (mag_diff > 1.0e-4)
                {
                    Console.WriteLine("C# CheckSaturn ERROR: Excessive magnitude error {0}", mag_diff);
                    error = 1;      // keep going -- print all errors before exiting
                }

                double tilt_diff = abs(illum.ring_tilt - data.tilt);
                if (tilt_diff > 3.0e-5)
                {
                    Console.WriteLine("C# CheckSaturn ERROR: Excessive ring tilt error {0}", tilt_diff);
                    error = 1;      // keep going -- print all errors before exiting
                }
            }

            return error;
        }

        static int TestMaxMag(Body body, string filename)
        {
            // Example of input data:
            //
            // 2001-02-21T08:00Z 2001-02-27T08:00Z 23.17 19.53 -4.84
            //
            // JPL Horizons test data has limited floating point precision in the magnitude values.
            // There is a pair of dates for the beginning and end of the max magnitude period,
            // given the limited precision.
            // We pick the point halfway between as the supposed max magnitude time.

            using (StreamReader infile = File.OpenText(filename))
            {
                int lnum = 0;
                string line;
                var search_time = new AstroTime(2001, 1, 1, 0, 0, 0);
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    string[] token = Tokenize(line);
                    if (token.Length != 5)
                    {
                        Console.WriteLine("C# TestMaxMag({0} line {1}): invalid data format", filename, lnum);
                        return 1;
                    }
                    AstroTime time1 = ParseDate(token[0]);
                    AstroTime time2 = ParseDate(token[1]);
                    double correct_angle1 = double.Parse(token[2]);
                    double correct_angle2 = double.Parse(token[3]);
                    double correct_mag = double.Parse(token[4]);
                    AstroTime center_time = time1.AddDays(0.5*(time2.ut - time1.ut));
                    IllumInfo illum = Astronomy.SearchPeakMagnitude(body, search_time);
                    double mag_diff = abs(illum.mag - correct_mag);
                    double hours_diff = 24.0 * abs(illum.time.ut - center_time.ut);
                    Debug("C# TestMaxMag: mag_diff={0}, hours_diff={1}", mag_diff, hours_diff);
                    if (hours_diff > 7.1)
                    {
                        Console.WriteLine("C# TestMaxMag({0} line {1}): EXCESSIVE TIME DIFFERENCE.", filename, lnum);
                        return 1;
                    }
                    if (mag_diff > 0.005)
                    {
                        Console.WriteLine("C# TestMaxMag({0} line {1}): EXCESSIVE MAGNITUDE DIFFERENCE.", filename, lnum);
                        return 1;
                    }
                    search_time = time2;
                }
                Debug("C# TestMaxMag: Processed {0} lines from file {1}", lnum, filename);
                return 0;
            }
        }

        static int MagnitudeTest()
        {
            int nfailed = 0;
            nfailed += CheckMagnitudeData(Body.Sun, "../../magnitude/Sun.txt");
            nfailed += CheckMagnitudeData(Body.Moon, "../../magnitude/Moon.txt");
            nfailed += CheckMagnitudeData(Body.Mercury, "../../magnitude/Mercury.txt");
            nfailed += CheckMagnitudeData(Body.Venus, "../../magnitude/Venus.txt");
            nfailed += CheckMagnitudeData(Body.Mars, "../../magnitude/Mars.txt");
            nfailed += CheckMagnitudeData(Body.Jupiter, "../../magnitude/Jupiter.txt");
            nfailed += CheckSaturn();
            nfailed += CheckMagnitudeData(Body.Uranus, "../../magnitude/Uranus.txt");
            nfailed += CheckMagnitudeData(Body.Neptune, "../../magnitude/Neptune.txt");
            nfailed += CheckMagnitudeData(Body.Pluto, "../../magnitude/Pluto.txt");
            nfailed += TestMaxMag(Body.Venus, "../../magnitude/maxmag_Venus.txt");
            if (nfailed == 0)
                Console.WriteLine("C# MagnitudeTest: PASS");
            else
                Console.WriteLine("C# MagnitudeTest: FAILED {0} test(s).", nfailed);
            return nfailed;
        }

        static double VectorDiff(AstroVector a, AstroVector b)
        {
            double dx = a.x - b.x;
            double dy = a.y - b.y;
            double dz = a.z - b.z;
            return sqrt(dx*dx + dy*dy + dz*dz);
        }

        static int CompareVectors(string caller, AstroVector a, AstroVector b, double tolerance)
        {
            double diff;

            diff = abs(a.x - b.x);
            if (diff > tolerance)
            {
                Console.WriteLine("C# CompareVectors ERROR({0}): x={1}, expected {2}, diff {3}", caller, a.x, b.x, diff);
                return 1;
            }

            diff = abs(a.y - b.y);
            if (diff > tolerance)
            {
                Console.WriteLine("C# CompareVectors ERROR({0}): y={1}, expected {2}, diff {3}", caller, a.y, b.y, diff);
                return 1;
            }

            diff = abs(a.z - b.z);
            if (diff > tolerance)
            {
                Console.WriteLine("C# CompareVectors ERROR({0}): z={1}, expected {2}, diff {3}", caller, a.z, b.z, diff);
                return 1;
            }

            return 0;
        }

        static int CompareMatrices(string caller, RotationMatrix a, RotationMatrix b, double tolerance)
        {
            for (int i=0; i<3; ++i)
            {
                for (int j=0; j<3; ++j)
                {
                    double diff = abs(a.rot[i,j] - b.rot[i,j]);
                    if (diff > tolerance)
                    {
                        Console.WriteLine("C# CompareMatrices ERROR({0}): matrix[{1},{2}]={3}, expected {4}, diff {5}", caller, i, j, a.rot[i,j], b.rot[i,j], diff);
                        return 1;
                    }
                }
            }
            return 0;
        }

        static int Rotation_MatrixInverse()
        {
            var a = new RotationMatrix(new double[3,3]);
            a.rot[0, 0] = 1.0; a.rot[1, 0] = 2.0; a.rot[2, 0] = 3.0;
            a.rot[0, 1] = 4.0; a.rot[1, 1] = 5.0; a.rot[2, 1] = 6.0;
            a.rot[0, 2] = 7.0; a.rot[1, 2] = 8.0; a.rot[2, 2] = 9.0;

            var v = new RotationMatrix(new double[3,3]);
            v.rot[0, 0] = 1.0; v.rot[1, 0] = 4.0; v.rot[2, 0] = 7.0;
            v.rot[0, 1] = 2.0; v.rot[1, 1] = 5.0; v.rot[2, 1] = 8.0;
            v.rot[0, 2] = 3.0; v.rot[1, 2] = 6.0; v.rot[2, 2] = 9.0;

            RotationMatrix b = Astronomy.InverseRotation(a);
            if (0 != CompareMatrices("Rotation_MatrixInverse", b, v, 0.0)) return 1;
            Console.WriteLine("C# Rotation_MatrixInverse: PASS");
            return 0;
        }

        static int Rotation_MatrixMultiply()
        {
            var a = new RotationMatrix(new double[3,3]);
            a.rot[0, 0] = 1.0; a.rot[1, 0] = 2.0; a.rot[2, 0] = 3.0;
            a.rot[0, 1] = 4.0; a.rot[1, 1] = 5.0; a.rot[2, 1] = 6.0;
            a.rot[0, 2] = 7.0; a.rot[1, 2] = 8.0; a.rot[2, 2] = 9.0;

            var b = new RotationMatrix(new double[3,3]);
            b.rot[0, 0] = 10.0; b.rot[1, 0] = 11.0; b.rot[2, 0] = 12.0;
            b.rot[0, 1] = 13.0; b.rot[1, 1] = 14.0; b.rot[2, 1] = 15.0;
            b.rot[0, 2] = 16.0; b.rot[1, 2] = 17.0; b.rot[2, 2] = 18.0;

            var v = new RotationMatrix(new double[3,3]);
            v.rot[0, 0] =  84.0; v.rot[1, 0] =  90.0; v.rot[2, 0] =  96.0;
            v.rot[0, 1] = 201.0; v.rot[1, 1] = 216.0; v.rot[2, 1] = 231.0;
            v.rot[0, 2] = 318.0; v.rot[1, 2] = 342.0; v.rot[2, 2] = 366.0;

            RotationMatrix c = Astronomy.CombineRotation(b, a);
            if (0 != CompareMatrices("Rotation_MatrixMultiply", c, v, 0.0)) return 1;
            Console.WriteLine("C# Rotation_MatrixMultiply: PASS");
            return 0;
        }

        static int Test_EQJ_ECL()
        {
            RotationMatrix r = Astronomy.Rotation_EQJ_ECL();

            // Calculate heliocentric Earth position at a test time.
            var time = new AstroTime(2019, 12, 8, 19, 39, 15);
            var ev = Astronomy.HelioVector(Body.Earth, time);

            // Use the older function to calculate ecliptic vector and angles.
            Ecliptic ecl = Astronomy.EquatorialToEcliptic(ev);
            Debug("C# Test_EQJ_ECL ecl = ({0}, {1}, {2})", ecl.vec.x, ecl.vec.y, ecl.vec.z);

            // Now compute the same vector via rotation matrix.
            AstroVector ee = Astronomy.RotateVector(r, ev);
            double dx = ee.x - ecl.vec.x;
            double dy = ee.y - ecl.vec.y;
            double dz = ee.z - ecl.vec.z;
            double diff = sqrt(dx*dx + dy*dy + dz*dz);
            Debug("C# Test_EQJ_ECL ee = ({0}, {1}, {2}); diff={3}", ee.x, ee.y, ee.z, diff);
            if (diff > 1.0e-16)
            {
                Console.WriteLine("C# Test_EQJ_ECL: EXCESSIVE VECTOR ERROR");
                return 1;
            }

            // Reverse the test: go from ecliptic back to equatorial.
            r = Astronomy.Rotation_ECL_EQJ();
            AstroVector et = Astronomy.RotateVector(r, ee);
            diff = VectorDiff(et, ev);
            Debug("C# Test_EQJ_ECL  ev diff={0}", diff);
            if (diff > 2.3e-16)
            {
                Console.WriteLine("C# Test_EQJ_ECL: EXCESSIVE REVERSE ROTATION ERROR");
                return 1;
            }

            Debug("C# Test_EQJ_ECL: PASS");
            return 0;
        }

        static int Test_EQJ_EQD(Body body)
        {
            // Verify convresion of equatorial J2000 to equatorial of-date, and back.
            var time = new AstroTime(2019, 12, 8, 20, 50, 0);
            var observer = new Observer(35, -85, 0);
            Equatorial eq2000 = Astronomy.Equator(body, time, observer, EquatorEpoch.J2000, Aberration.Corrected);
            Equatorial eqdate = Astronomy.Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
            AstroVector v2000 = eq2000.vec;
            RotationMatrix r = Astronomy.Rotation_EQJ_EQD(time);
            AstroVector vdate = Astronomy.RotateVector(r, v2000);
            Equatorial eqcheck = Astronomy.EquatorFromVector(vdate);

            double ra_diff = abs(eqcheck.ra - eqdate.ra);
            double dec_diff = abs(eqcheck.dec - eqdate.dec);
            double dist_diff = abs(eqcheck.dist - eqdate.dist);
            Debug("C# Test_EQJ_EQD: {0} ra={1}, dec={2}, dist={3}, ra_diff={4}, dec_diff={5}, dist_diff={6}",
                body, eqdate.ra, eqdate.dec, eqdate.dist, ra_diff, dec_diff, dist_diff);

            if (ra_diff > 1.0e-14 || dec_diff > 1.0e-14 || dist_diff > 4.0e-15)
            {
                Console.WriteLine("C# Test_EQJ_EQD: EXCESSIVE ERROR");
                return 1;
            }

            r = Astronomy.Rotation_EQD_EQJ(time);
            AstroVector t2000 = Astronomy.RotateVector(r, vdate);
            double diff = VectorDiff(t2000, v2000);
            Debug("C# Test_EQJ_EQD: {0} inverse diff = {1}", body, diff);
            if (diff > 5.0e-15)
            {
                Console.WriteLine("C# Test_EQJ_EQD: EXCESSIVE INVERSE ERROR");
                return 1;
            }

            Debug("C# Test_EQJ_EQD: PASS");
            return 0;
        }

        static int Test_EQD_HOR(Body body)
        {
            var time = new AstroTime(1970, 12, 13, 5, 15, 0);
            var observer = new Observer(-37.0, +45.0, 0.0);
            Equatorial eqd = Astronomy.Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
            Topocentric hor = Astronomy.Horizon(time, observer, eqd.ra, eqd.dec, Refraction.Normal);
            AstroVector vec_eqd = eqd.vec;
            RotationMatrix rot = Astronomy.Rotation_EQD_HOR(time, observer);
            AstroVector vec_hor = Astronomy.RotateVector(rot, vec_eqd);
            Spherical sphere = Astronomy.HorizonFromVector(vec_hor, Refraction.Normal);

            double diff_alt = abs(sphere.lat - hor.altitude);
            double diff_az = abs(sphere.lon - hor.azimuth);

            Debug("C# Test_EQD_HOR {0}: trusted alt={1}, az={2}; test alt={3}, az={4}; diff_alt={5}, diff_az={6}",
                body, hor.altitude, hor.azimuth, sphere.lat, sphere.lon, diff_alt, diff_az);

            if (diff_alt > 3.2e-14 || diff_az > 1.2e-13)
            {
                Console.WriteLine("C# Test_EQD_HOR: EXCESSIVE HORIZONTAL ERROR.");
                return 1;
            }

            // Confirm that we can convert back to horizontal vector.
            AstroVector check_hor = Astronomy.VectorFromHorizon(sphere, time, Refraction.Normal);
            double diff = VectorDiff(check_hor, vec_hor);
            Debug("C# Test_EQD_HOR {0}: horizontal recovery: diff = {1}", body, diff);
            if (diff > 3.0e-15)
            {
                Console.WriteLine("C# Test_EQD_HOR: EXCESSIVE ERROR IN HORIZONTAL RECOVERY.");
                return 1;
            }

            // Verify the inverse translation from horizontal vector to equatorial of-date vector.
            rot = Astronomy.Rotation_HOR_EQD(time, observer);
            AstroVector check_eqd = Astronomy.RotateVector(rot, vec_hor);
            diff = VectorDiff(check_eqd, vec_eqd);
            Debug("C# Test_EQD_HOR {0}: OFDATE inverse rotation diff = {1}", body, diff);
            if (diff > 2.1e-15)
            {
                Console.WriteLine("C# Test_EQD_HOR: EXCESSIVE OFDATE INVERSE HORIZONTAL ERROR.");
                return 1;
            }

            // Exercise HOR to EQJ translation.
            Equatorial eqj = Astronomy.Equator(body, time, observer, EquatorEpoch.J2000, Aberration.Corrected);
            AstroVector vec_eqj = eqj.vec;

            rot = Astronomy.Rotation_HOR_EQJ(time, observer);
            AstroVector check_eqj = Astronomy.RotateVector(rot, vec_hor);
            diff = VectorDiff(check_eqj, vec_eqj);
            Debug("C# Test_EQD_HOR {0}: J2000 inverse rotation diff = {1}", body, diff);
            if (diff > 6.0e-15)
            {
                Console.WriteLine("C# Test_EQD_HOR: EXCESSIVE J2000 INVERSE HORIZONTAL ERROR.");
                return 1;
            }

            // Verify the inverse translation: EQJ to HOR.
            rot = Astronomy.Rotation_EQJ_HOR(time, observer);
            check_hor = Astronomy.RotateVector(rot, vec_eqj);
            diff = VectorDiff(check_hor, vec_hor);
            Debug("C# Test_EQD_HOR {0}: EQJ inverse rotation diff = {1}", body, diff);
            if (diff > 3e-15)
            {
                Console.WriteLine("C# Test_EQD_HOR: EXCESSIVE EQJ INVERSE HORIZONTAL ERROR.");
                return 1;
            }

            Debug("C# Test_EQD_HOR: PASS");
            return 0;
        }


        static int Test_EQJ_GAL_NOVAS(string filename)
        {
            const double THRESHOLD_SECONDS = 8.8;
            RotationMatrix rot = Astronomy.Rotation_EQJ_GAL();
            RotationMatrix inv = Astronomy.Rotation_GAL_EQJ();
            if (0 != CheckInverse("EQJ_GAL", "GAL_EQJ", rot, inv)) return 1;

            var time = new AstroTime(0.0);  // placeholder time - value does not matter
            double max_diff = 0.0;

            using (StreamReader infile = File.OpenText(filename))
            {
                string line;
                int lnum = 0;
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    string[] token = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    if (token.Length != 4)
                    {
                        Console.WriteLine("C# Test_EQJ_GAL_NOVAS({0} line {1}): found {2} tokens instead of 4.", filename, lnum, token.Length);
                        return 1;
                    }

                    double ra = double.Parse(token[0]);
                    double dec = double.Parse(token[1]);
                    double glon = double.Parse(token[2]);
                    double glat = double.Parse(token[3]);

                    // Use Astronomy Engine to do the same EQJ/GAL conversion.
                    var eqj_sphere = new Spherical(dec, 15.0 * ra, 1.0);
                    AstroVector eqj_vec = Astronomy.VectorFromSphere(eqj_sphere, time);
                    AstroVector gal_vec = Astronomy.RotateVector(rot, eqj_vec);
                    Spherical gal_sphere = Astronomy.SphereFromVector(gal_vec);
                    double dlat = v(gal_sphere.lat - glat);
                    double dlon = cos(Astronomy.DEG2RAD * glat) * v(gal_sphere.lon - glon);
                    double diff = 3600.0 * sqrt(dlon*dlon + dlat*dlat);
                    if (diff > THRESHOLD_SECONDS)
                    {
                        Console.WriteLine("C# Test_EQJ_GAL_NOVAS({0} line {1}): EXCESSIVE ERROR = {2:F3} arcseconds.", filename, lnum, diff);
                        return 1;
                    }
                    if (diff > max_diff)
                        max_diff = diff;
                }
            }

            Debug("C# Test_EQJ_GAL_NOVAS: PASS. max_diff = {0:F3} arcseconds.", max_diff);
            return 0;
        }


        static int CheckInverse(string aname, string bname, RotationMatrix arot, RotationMatrix brot)
        {
            RotationMatrix crot = Astronomy.CombineRotation(arot, brot);

            var rot = new double[3,3];
            rot[0, 0] = 1; rot[1, 0] = 0; rot[2, 0] = 0;
            rot[0, 1] = 0; rot[1, 1] = 1; rot[2, 1] = 0;
            rot[0, 2] = 0; rot[1, 2] = 0; rot[2, 2] = 1;
            var identity = new RotationMatrix(rot);

            string caller = "CheckInverse(" + aname + ", " + bname + ")";
            return CompareMatrices(caller, crot, identity, 2.0e-15);
        }


        static int CheckCycle(
            string aname, string bname, string cname,
            RotationMatrix arot, RotationMatrix brot, RotationMatrix crot)
        {
            RotationMatrix xrot = Astronomy.CombineRotation(arot, brot);
            RotationMatrix irot = Astronomy.InverseRotation(xrot);
            string name = string.Format("CheckCycle({0}, {1}, {2})", aname, bname, cname);
            return CompareMatrices(name, crot, irot, 2.0e-15);
        }


        static int Test_RotRoundTrip()
        {
            AstroTime time = new AstroTime(2067, 5, 30, 14, 45, 0);
            Observer observer = new Observer(+28.0, -82.0, 0.0);

            // In each round trip, calculate a forward rotation and a backward rotation.
            // Verify the two are inverse matrices.

            // Round trip #1: EQJ <==> EQD.
            RotationMatrix eqj_eqd = Astronomy.Rotation_EQJ_EQD(time);
            RotationMatrix eqd_eqj = Astronomy.Rotation_EQD_EQJ(time);
            if (0 != CheckInverse(nameof(eqj_eqd), nameof(eqd_eqj), eqj_eqd, eqd_eqj)) return 1;

            // Round trip #2: EQJ <==> ECL.
            RotationMatrix eqj_ecl = Astronomy.Rotation_EQJ_ECL();
            RotationMatrix ecl_eqj = Astronomy.Rotation_ECL_EQJ();
            if (0 != CheckInverse(nameof(eqj_ecl), nameof(ecl_eqj), eqj_ecl, ecl_eqj)) return 1;

            // Round trip #3: EQJ <==> HOR.
            RotationMatrix eqj_hor = Astronomy.Rotation_EQJ_HOR(time, observer);
            RotationMatrix hor_eqj = Astronomy.Rotation_HOR_EQJ(time, observer);
            if (0 != CheckInverse(nameof(eqj_hor), nameof(hor_eqj), eqj_hor, hor_eqj)) return 1;

            // Round trip #4: EQD <==> HOR.
            RotationMatrix eqd_hor = Astronomy.Rotation_EQD_HOR(time, observer);
            RotationMatrix hor_eqd = Astronomy.Rotation_HOR_EQD(time, observer);
            if (0 != CheckInverse(nameof(eqd_hor), nameof(hor_eqd), eqd_hor, hor_eqd)) return 1;

            // Round trip #5: EQD <==> ECL.
            RotationMatrix eqd_ecl = Astronomy.Rotation_EQD_ECL(time);
            RotationMatrix ecl_eqd = Astronomy.Rotation_ECL_EQD(time);
            if (0 != CheckInverse(nameof(eqd_ecl), nameof(ecl_eqd), eqd_ecl, ecl_eqd)) return 1;

            // Round trip #6: HOR <==> ECL.
            RotationMatrix hor_ecl = Astronomy.Rotation_HOR_ECL(time, observer);
            RotationMatrix ecl_hor = Astronomy.Rotation_ECL_HOR(time, observer);
            if (0 != CheckInverse(nameof(hor_ecl), nameof(ecl_hor), hor_ecl, ecl_hor)) return 1;

            // Verify that combining different sequences of rotations result
            // in the expected combination.
            // For example, (EQJ ==> HOR ==> ECL) must be the same matrix as (EQJ ==> ECL).
            // Each of these is a "triangle" of relationships between 3 orientations.
            // There are 4 possible ways to pick 3 orientations from the 4 to form a triangle.
            // Because we have just proved that each transformation is reversible,
            // we only need to verify the triangle in one cyclic direction.
            if (0 != CheckCycle(nameof(eqj_ecl), nameof(ecl_eqd), nameof(eqd_eqj), eqj_ecl, ecl_eqd, eqd_eqj)) return 1;     // excluded corner = HOR
            if (0 != CheckCycle(nameof(eqj_hor), nameof(hor_ecl), nameof(ecl_eqj), eqj_hor, hor_ecl, ecl_eqj)) return 1;     // excluded corner = EQD
            if (0 != CheckCycle(nameof(eqj_hor), nameof(hor_eqd), nameof(eqd_eqj), eqj_hor, hor_eqd, eqd_eqj)) return 1;     // excluded corner = ECL
            if (0 != CheckCycle(nameof(ecl_eqd), nameof(eqd_hor), nameof(hor_ecl), ecl_eqd, eqd_hor, hor_ecl)) return 1;     // excluded corner = EQJ

            Debug("C# Test_RotRoundTrip: PASS");
            return 0;
        }

        static int Rotation_Pivot()
        {
            //astro_rotation_t a;
            //astro_vector_t v1, v2, ve;
            const double tolerance = 1.0e-15;

            // Test #1

            // Start with an identity matrix.
            RotationMatrix ident = Astronomy.IdentityMatrix();

            // Pivot 90 degrees counterclockwise around the z-axis.
            RotationMatrix r = Astronomy.Pivot(ident, 2, +90.0);

            // Put the expected answer in 'a'.
            var a = new RotationMatrix(new double[3,3]
            {
                {  0, +1,  0 },
                { -1,  0,  0 },
                {  0,  0, +1 },
            });

            // Compare actual 'r' with expected 'a'.
            if (0 != CompareMatrices("Rotation_Pivot #1", r, a, tolerance)) return 1;

            // Test #2.

            // Pivot again, -30 degrees around the x-axis.
            r = Astronomy.Pivot(r, 0, -30.0);

            // Pivot a third time, 180 degrees around the y-axis.
            r = Astronomy.Pivot(r, 1, +180.0);

            // Use the 'r' matrix to rotate a vector.
            var v1 = new AstroVector(1.0, 2.0, 3.0, new AstroTime(0.0));

            AstroVector v2 = Astronomy.RotateVector(r, v1);

            // Initialize the expected vector 've'.
            AstroVector ve = new AstroVector(+2.0, +2.3660254037844390, -2.0980762113533156, v1.t);

            if (0 != CompareVectors("Rotation_Pivot #2", v2, ve, tolerance)) return 1;

            Console.WriteLine("C# Rotation_Pivot: PASS");
            return 0;
        }

        static int RotationTest()
        {
            if (0 != Rotation_MatrixInverse()) return 1;
            if (0 != Rotation_MatrixMultiply()) return 1;
            if (0 != Rotation_Pivot()) return 1;
            if (0 != Test_EQJ_ECL()) return 1;
            if (0 != Test_EQJ_GAL_NOVAS("../../temp/galeqj.txt")) return 1;

            if (0 != Test_EQJ_EQD(Body.Mercury)) return 1;
            if (0 != Test_EQJ_EQD(Body.Venus)) return 1;
            if (0 != Test_EQJ_EQD(Body.Mars)) return 1;
            if (0 != Test_EQJ_EQD(Body.Jupiter)) return 1;
            if (0 != Test_EQJ_EQD(Body.Saturn)) return 1;

            if (0 != Test_EQD_HOR(Body.Mercury)) return 1;
            if (0 != Test_EQD_HOR(Body.Venus)) return 1;
            if (0 != Test_EQD_HOR(Body.Mars)) return 1;
            if (0 != Test_EQD_HOR(Body.Jupiter)) return 1;
            if (0 != Test_EQD_HOR(Body.Saturn)) return 1;

            if (0 != Test_RotRoundTrip()) return 1;

            Console.WriteLine("C# RotationTest: PASS");
            return 0;
        }

        static int RefractionTest()
        {
            for (double alt = -90.1; alt <= +90.1; alt += 0.001)
            {
                double refr = Astronomy.RefractionAngle(Refraction.Normal, alt);
                double corrected = alt + refr;
                double inv_refr = Astronomy.InverseRefractionAngle(Refraction.Normal, corrected);
                double check_alt = corrected + inv_refr;
                double diff = abs(check_alt - alt);
                if (diff > 2.0e-14)
                {
                    Console.WriteLine("C# ERROR(RefractionTest): alt={0}, refr={1}, diff={2}", alt, refr, diff);
                    return 1;
                }
            }

            Console.WriteLine("C# RefractionTest: PASS");
            return 0;
        }

        static int ConstellationTest()
        {
            const string inFileName = "../../constellation/test_input.txt";
            int failcount = 0;
            int lnum = 0;
            using (StreamReader infile = File.OpenText(inFileName))
            {
                string line;
                Regex reLine = new Regex(@"^\s*(\d+)\s+(\S+)\s+(\S+)\s+([A-Z][a-zA-Z]{2})\s*$");
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    Match m = reLine.Match(line);
                    if (!m.Success)
                    {
                        Console.WriteLine("C# ERROR(ConstellationTest): invalid line {0} in file {1}", lnum, inFileName);
                        return 1;
                    }
                    int id = int.Parse(m.Groups[1].Value);
                    double ra = double.Parse(m.Groups[2].Value);
                    double dec = double.Parse(m.Groups[3].Value);
                    string symbol = m.Groups[4].Value;
                    ConstellationInfo constel = Astronomy.Constellation(ra, dec);
                    if (constel.Symbol != symbol)
                    {
                        Console.WriteLine("Star {0,6}: expected {1}, found {2} at B1875 RA={3,10}, DEC={4,10}", id, symbol, constel.Symbol, constel.Ra1875, constel.Dec1875);
                        ++failcount;
                    }
                }
            }
            if (failcount > 0)
            {
                Console.WriteLine("C# ConstellationTest: {0} failures", failcount);
                return 1;
            }

            Console.WriteLine("C# ConstellationTest: PASS (verified {0})", lnum);
            return 0;
        }

        static int LunarEclipseIssue78()
        {
            LunarEclipseInfo eclipse = Astronomy.SearchLunarEclipse(new AstroTime(2020, 12, 19, 0, 0, 0));
            var expected_peak = new AstroTime(2021, 5, 26, 11, 18, 42);  // https://www.timeanddate.com/eclipse/lunar/2021-may-26
            double dt_seconds = SECONDS_PER_DAY * abs(expected_peak.tt - eclipse.peak.tt);
            if (dt_seconds > 40.0)
            {
                Console.WriteLine("C# LunarEclipseIssue78: Excessive prediction error = {0} seconds.", dt_seconds);
                return 1;
            }
            if (eclipse.kind != EclipseKind.Total)
            {
                Console.WriteLine("C# LunarEclipseIssue78: Expected total eclipse; found {0}", eclipse.kind);
                return 1;
            }
            Console.WriteLine("C# LunarEclipseIssue78: PASS");
            return 0;
        }

        static int LunarEclipseTest()
        {
            const string filename = "../../eclipse/lunar_eclipse.txt";
            const string statsFilename = "../../eclipse/cs_le_stats.csv";

            Astronomy.CalcMoonCount = 0;
            using (StreamReader infile = File.OpenText(filename))
            {
                using (StreamWriter outfile = File.CreateText(statsFilename))
                {
                    outfile.WriteLine("\"utc\",\"center\",\"partial\",\"total\"");
                    LunarEclipseInfo eclipse = Astronomy.SearchLunarEclipse(new AstroTime(1701, 1, 1, 0, 0, 0));
                    string line;
                    int lnum = 0;
                    int skip_count = 0;
                    int diff_count = 0;
                    double sum_diff_minutes = 0.0;
                    double max_diff_minutes = 0.0;
                    const double diff_limit = 2.0;
                    while (null != (line = infile.ReadLine()))
                    {
                        ++lnum;

                        // Make sure numeric data are finite numbers.
                        v(eclipse.obscuration);
                        v(eclipse.sd_partial);
                        v(eclipse.sd_penum);
                        v(eclipse.sd_total);

                        if (line.Length < 17)
                        {
                            Console.WriteLine("C# LunarEclipseTest({0} line {1}): line is too short.", filename, lnum);
                            return 1;
                        }
                        string time_text = line.Substring(0, 17);
                        AstroTime peak_time = ParseDate(time_text);
                        string[] token = Tokenize(line.Substring(17));
                        double partial_minutes, total_minutes;
                        if (token.Length != 2 || !double.TryParse(token[0], out partial_minutes) || !double.TryParse(token[1], out total_minutes))
                        {
                            Console.WriteLine("C# LunarEclipseTest({0} line {1}): invalid data format.", filename, lnum);
                            return 1;
                        }

                        // Verify that the calculated eclipse semi-durations are consistent with the kind.
                        // Verify that obscurations also make sense for the kind.
                        bool sd_valid = false;
                        bool frac_valid = false;
                        switch (eclipse.kind)
                        {
                        case EclipseKind.Penumbral:
                            sd_valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial == 0.0) && (eclipse.sd_total == 0.0);
                            frac_valid = (eclipse.obscuration == 0.0);
                            break;

                        case EclipseKind.Partial:
                            sd_valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial > 0.0) && (eclipse.sd_total == 0.0);
                            frac_valid = (eclipse.obscuration > 0.0) && (eclipse.obscuration < 1.0);
                            break;

                        case EclipseKind.Total:
                            sd_valid = (eclipse.sd_penum > 0.0) && (eclipse.sd_partial > 0.0) && (eclipse.sd_total > 0.0);
                            frac_valid = (eclipse.obscuration == 1.0);
                            break;

                        default:
                            Console.WriteLine("C# LunarEclipseTest({0} line {1}): invalid eclipse kind {2}.", filename, lnum, eclipse.kind);
                            return 1;
                        }

                        if (!sd_valid)
                        {
                            Console.WriteLine("C# LunarEclipseTest({0} line {1}): inalid semiduration(s) for kind {2}: penum={3}, partial={4}, total={5}",
                                filename, lnum, eclipse.kind, eclipse.sd_penum, eclipse.sd_partial, eclipse.sd_total);
                            return 1;
                        }

                        if (!frac_valid)
                        {
                            Console.WriteLine($"C# LunarEclipseTest({filename} line {lnum}): invalid obscuration {eclipse.obscuration} for kind {eclipse.kind}");
                            return 1;
                        }

                        // Check eclipse peak time.
                        double diff_days = eclipse.peak.ut - peak_time.ut;

                        // Tolerate missing penumbral eclipses - skip to next input line without calculating next eclipse.
                        if (partial_minutes == 0.0 && diff_days > 20.0)
                        {
                            ++skip_count;
                            continue;
                        }

                        outfile.WriteLine("\"{0}\",{1},{2},{3}",
                            time_text,
                            diff_days * (24.0 * 60.0),
                            eclipse.sd_partial - partial_minutes,
                            eclipse.sd_total - total_minutes
                        );

                        double diff_minutes = (24.0 * 60.0) * abs(diff_days);
                        sum_diff_minutes += diff_minutes;
                        ++diff_count;

                        if (diff_minutes > diff_limit)
                        {
                            Console.WriteLine("C# LunarEclipseTest expected peak: {0}", peak_time);
                            Console.WriteLine("C# LunarEclipseTest found    peak: {1}", eclipse.peak);
                            Console.WriteLine("C# LunarEclipseTest({0} line {1}): EXCESSIVE peak time error = {2} minutes ({3} days).", filename, lnum, diff_minutes, diff_days);
                            return 1;
                        }

                        if (diff_minutes > max_diff_minutes)
                            max_diff_minutes = diff_minutes;

                        // check partial eclipse duration

                        diff_minutes = abs(partial_minutes - eclipse.sd_partial);
                        sum_diff_minutes += diff_minutes;
                        ++diff_count;

                        if (diff_minutes > diff_limit)
                        {
                            Console.WriteLine("C# LunarEclipseTest({0} line {1}): EXCESSIVE partial eclipse semiduration error: {2} minutes", filename, lnum, diff_minutes);
                            return 1;
                        }

                        if (diff_minutes > max_diff_minutes)
                            max_diff_minutes = diff_minutes;

                        // check total eclipse duration

                        diff_minutes = abs(total_minutes - eclipse.sd_total);
                        sum_diff_minutes += diff_minutes;
                        ++diff_count;

                        if (diff_minutes > diff_limit)
                        {
                            Console.WriteLine("C# LunarEclipseTest({0} line {1}): EXCESSIVE total eclipse semiduration error: {2} minutes", filename, lnum, diff_minutes);
                            return 1;
                        }

                        if (diff_minutes > max_diff_minutes)
                            max_diff_minutes = diff_minutes;

                        // calculate for next iteration

                        eclipse = Astronomy.NextLunarEclipse(eclipse.peak);
                    }
                    Console.WriteLine("C# LunarEclipseTest: PASS (verified {0}, skipped {1}, max_diff_minutes = {2}, avg_diff_minutes = {3}, moon calcs = {4})", lnum, skip_count, max_diff_minutes, (sum_diff_minutes / diff_count), Astronomy.CalcMoonCount);
                }
            }
            return 0;
        }

        static int LunarFractionCase(int year, int month, int day, double obscuration)
        {
            // Search for the first lunar eclipse to occur after the given date.
            // It should always happen within 24 hours of the given date.
            AstroTime time = new AstroTime(year, month, day, 0, 0, 0.0);
            LunarEclipseInfo eclipse = Astronomy.SearchLunarEclipse(time);

            if (eclipse.kind != EclipseKind.Partial)
            {
                Console.WriteLine($"C# LunarFractionCase({year:0000}-{month:00}-{day:00}): expected partial eclipse, but found {eclipse.kind}.");
                return 1;
            }

            double dt = v(eclipse.peak.ut - time.ut);
            if (dt < 0.0 || dt > 1.0)
            {
                Console.WriteLine($"C# LunarFractionCase({year:0000}-{month:00}-{day:00}): eclipse occurs {dt:F4} days after predicted date.");
                return 1;
            }

            double diff = v(eclipse.obscuration - obscuration);
            if (abs(diff) > 0.00901)
            {
                Console.WriteLine($"C# LunarFractionCase({year:0000}-{month:00}-{day:00}) FAIL: obscuration error = {diff:F8}, expected = {obscuration:F3}, calculated = {eclipse.obscuration:F8}");
                return 1;
            }
            Debug($"C# LunarFractionCase({year:0000}-{month:00}-{day:00}) obscuration error = {diff:F8}");
            return 0;
        }

        static int LunarFractionTest()
        {
            // Verify calculation of the fraction of the Moon's disc covered by the Earth's umbra during a partial eclipse.
            // Data for this is more tedious to gather, because Espenak data does not contain it.
            // We already verify fraction=0.0 for penumbral eclipses and fraction=1.0 for total eclipses in LunarEclipseTest.

            if (0 != LunarFractionCase(2010,  6, 26, 0.506)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2010-june-26
            if (0 != LunarFractionCase(2012,  6,  4, 0.304)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2012-june-4
            if (0 != LunarFractionCase(2013,  4, 25, 0.003)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2013-april-25
            if (0 != LunarFractionCase(2017,  8,  7, 0.169)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2017-august-7
            if (0 != LunarFractionCase(2019,  7, 16, 0.654)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2019-july-16
            if (0 != LunarFractionCase(2021, 11, 19, 0.991)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2021-november-19
            if (0 != LunarFractionCase(2023, 10, 28, 0.060)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2023-october-28
            if (0 != LunarFractionCase(2024,  9, 18, 0.035)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2024-september-18
            if (0 != LunarFractionCase(2026,  8, 28, 0.962)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2026-august-28
            if (0 != LunarFractionCase(2028,  1, 12, 0.024)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2028-january-12
            if (0 != LunarFractionCase(2028,  7,  6, 0.325)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2028-july-6
            if (0 != LunarFractionCase(2030,  6, 15, 0.464)) return 1;  // https://www.timeanddate.com/eclipse/lunar/2030-june-15

            Console.WriteLine("C# LunarFractionTest: PASS");
            return 0;
        }

        static int GlobalSolarEclipseTest()
        {
            Astronomy.CalcMoonCount = 0;

            int lnum = 0;
            int skip_count = 0;
            double max_angle = 0.0;
            double max_minutes = 0.0;
            GlobalSolarEclipseInfo eclipse = Astronomy.SearchGlobalSolarEclipse(new AstroTime(1701, 1, 1, 0, 0, 0));

            const string inFileName = "../../eclipse/solar_eclipse.txt";
            using (StreamReader infile = File.OpenText(inFileName))
            {
                string line;
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    string[] token = Tokenize(line);
                    // 1889-12-22T12:54:15Z   -6 T   -12.7   -12.8
                    if (token.Length != 5)
                    {
                        Console.WriteLine("C# GlobalSolarEclipseTest({0} line {1}): wrong token count = {2}", inFileName, lnum, token.Length);
                        return 1;
                    }
                    AstroTime peak = ParseDate(token[0]);
                    string typeChar = token[2];
                    double lat = double.Parse(token[3]);
                    double lon = double.Parse(token[4]);
                    EclipseKind expected_kind;
                    switch (typeChar)
                    {
                        case "P": expected_kind = EclipseKind.Partial;  break;
                        case "A": expected_kind = EclipseKind.Annular;  break;
                        case "T": expected_kind = EclipseKind.Total;    break;
                        case "H": expected_kind = EclipseKind.Total;    break;
                        default:
                            Console.WriteLine("C# GlobalSolarEclipseTest({0} line {1}): invalid eclipse kind '{2}' in test data.", inFileName, lnum, typeChar);
                            return 1;
                    }

                    double diff_days = eclipse.peak.ut - peak.ut;

                    // Sometimes we find marginal eclipses that aren't listed in the test data.
                    // Ignore them if the distance between the Sun/Moon shadow axis and the Earth's center is large.
                    while (diff_days < -25.0 && eclipse.distance > 9000.0)
                    {
                        ++skip_count;
                        eclipse = Astronomy.NextGlobalSolarEclipse(eclipse.peak);
                        diff_days = eclipse.peak.ut - peak.ut;
                    }

                    // Validate the eclipse prediction.
                    double diff_minutes = (24 * 60) * abs(diff_days);
                    if (diff_minutes > 6.93)
                    {
                        Console.WriteLine("C# GlobalSolarEclipseTest({0} line {1}): EXCESSIVE TIME ERROR = {2} minutes", inFileName, lnum, diff_minutes);
                        return 1;
                    }

                    if (diff_minutes > max_minutes)
                        max_minutes = diff_minutes;

                    // Validate the eclipse kind, but only when it is not a "glancing" eclipse.
                    if ((eclipse.distance < 6360) && (eclipse.kind != expected_kind))
                    {
                        Console.WriteLine("C# GlobalSolarEclipseTest({0} line {1}): WRONG ECLIPSE KIND: expected {2}, found {3}", inFileName, lnum, expected_kind, eclipse.kind);
                        return 1;
                    }

                    if (eclipse.kind == EclipseKind.Total || eclipse.kind == EclipseKind.Annular)
                    {
                        // When the distance between the Moon's shadow ray and the Earth's center is beyond 6100 km,
                        // it creates a glancing blow whose geographic coordinates are excessively sensitive to
                        // slight changes in the ray. Therefore, it is unreasonable to count large errors there.
                        if (eclipse.distance < 6100.0)
                        {
                            double diff_angle = AngleDiff(lat, lon, eclipse.latitude, eclipse.longitude);
                            if (diff_angle > 0.247)
                            {
                                Console.WriteLine("C# GlobalSolarEclipseTest({0} line {1}): EXCESSIVE GEOGRAPHIC LOCATION ERROR = {2} degrees", inFileName, lnum, diff_angle);
                                return 1;
                            }
                            if (diff_angle > max_angle)
                                max_angle = diff_angle;
                        }
                    }

                    eclipse = Astronomy.NextGlobalSolarEclipse(eclipse.peak);
                }
            }

            const int expected_count = 1180;
            if (lnum != expected_count)
            {
                Console.WriteLine("C# GlobalSolarEclipseTest: WRONG LINE COUNT = {0}, expected {1}", lnum, expected_count);
                return 1;
            }

            if (skip_count > 2)
            {
                Console.WriteLine("C# GlobalSolarEclipseTest: EXCESSSIVE SKIP COUNT = {0}", skip_count);
                return 1;
            }

            Console.WriteLine("C# GlobalSolarEclipseTest: PASS ({0} verified, {1} skipped, {2} CalcMoons, max minutes = {3}, max angle = {4})", lnum, skip_count, Astronomy.CalcMoonCount, max_minutes, max_angle);
            return 0;
        }

        static void VectorFromAngles(double[] v, double lat, double lon)
        {
            double coslat = cos(Astronomy.DEG2RAD * lat);
            v[0] = cos(Astronomy.DEG2RAD * lon) * coslat;
            v[1] = sin(Astronomy.DEG2RAD * lon) * coslat;
            v[2] = sin(Astronomy.DEG2RAD * lat);
        }

        static double AngleDiff(double alat, double alon, double blat, double blon)
        {
            double[] a = new double[3];
            double[] b = new double[3];
            double dot;

            // Convert angles to vectors on a unit sphere.
            VectorFromAngles(a, alat, alon);
            VectorFromAngles(b, blat, blon);

            dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
            if (dot <= -1.0)
                return 180.0;

            if (dot >= +1.0)
                return 0.0;

            return v(Astronomy.RAD2DEG * Math.Acos(dot));
        }

        static int LocalSolarEclipseTest1()
        {
            // Re-use the test data for global solar eclipses, only feed the given coordinates
            // into the local solar eclipse predictor as the observer's location.
            // In each case, start the search 20 days before the expected eclipse.
            // Then verify that the peak time and eclipse type is correct in each case.
            Astronomy.CalcMoonCount = 0;

            int lnum = 0;
            int skip_count = 0;
            double max_minutes = 0.0;

            const string inFileName = "../../eclipse/solar_eclipse.txt";
            using (StreamReader infile = File.OpenText(inFileName))
            {
                string line;
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    string[] token = Tokenize(line);
                    // 1889-12-22T12:54:15Z   -6 T   -12.7   -12.8
                    if (token.Length != 5)
                    {
                        Console.WriteLine("C# LocalSolarEclipseTest1({0} line {1}): wrong token count = {2}", inFileName, lnum, token.Length);
                        return 1;
                    }
                    AstroTime peak = ParseDate(token[0]);
                    string typeChar = token[2];
                    double lat = double.Parse(token[3]);
                    double lon = double.Parse(token[4]);
                    var observer = new Observer(lat, lon, 0.0);

                    // Start the search 20 days before we know the eclipse should peak.
                    AstroTime search_start = peak.AddDays(-20.0);
                    LocalSolarEclipseInfo eclipse = Astronomy.SearchLocalSolarEclipse(search_start, observer);

                    // Validate the predicted eclipse peak time.
                    double diff_days = eclipse.peak.time.ut - peak.ut;
                    if (diff_days > 20.0)
                    {
                        ++skip_count;
                        continue;
                    }

                    double diff_minutes = (24 * 60) * abs(diff_days);
                    if (diff_minutes > 7.14)
                    {
                        Console.WriteLine("C LocalSolarEclipseTest1({0} line {1}): EXCESSIVE TIME ERROR = {2} minutes", inFileName, lnum, diff_minutes);
                        return 1;
                    }

                    if (diff_minutes > max_minutes)
                        max_minutes = diff_minutes;
                }
            }

            if (skip_count > 6)
            {
                Console.WriteLine("C# LocalSolarEclipseTest1: EXCESSSIVE SKIP COUNT = {0}", skip_count);
                return 1;
            }

            Console.WriteLine("C# LocalSolarEclipseTest1: PASS ({0} verified, {1} skipped, {2} CalcMoons, max minutes = {3})", lnum, skip_count, Astronomy.CalcMoonCount, max_minutes);
            return 0;
        }

        static string StripLine(string line)
        {
            int index = line.IndexOf('#');
            if (index >= 0)
                line = line.Substring(0, index);
            return line.Trim();
        }

        static int LocalSolarEclipseTest2()
        {
            // Test ability to calculate local solar eclipse conditions away from
            // the peak position on the Earth.

            Astronomy.CalcMoonCount = 0;

            const string inFileName = "../../eclipse/local_solar_eclipse.txt";
            double max_minutes=0.0, max_degrees=0.0;
            int verify_count = 0;

            using (StreamReader infile = File.OpenText(inFileName))
            {
                int lnum = 0;
                string line;
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    line = StripLine(line);
                    if (line.Length == 0) continue;
                    string[] token = Tokenize(line);
                    if (token.Length != 13)
                    {
                        Console.WriteLine("C# LocalSolarEclipseTest2({0} line {1}): Incorrect token count {2}", inFileName, lnum, token.Length);
                        return 1;
                    }
                    double latitude = double.Parse(token[0]);
                    double longitude = double.Parse(token[1]);
                    string typeCode = token[2];
                    AstroTime p1 = ParseDate(token[3]);
                    double p1alt = double.Parse(token[4]);
                    AstroTime t1 = OptionalParseDate(token[5]);
                    double t1alt = double.Parse(token[6]);
                    AstroTime peak = ParseDate(token[7]);
                    double peakalt = double.Parse(token[8]);
                    AstroTime t2 = OptionalParseDate(token[9]);
                    double t2alt = double.Parse(token[10]);
                    AstroTime p2 = ParseDate(token[11]);
                    double p2alt = double.Parse(token[12]);

                    EclipseKind expected_kind;
                    switch (typeCode)
                    {
                    case "P": expected_kind = EclipseKind.Partial; break;
                    case "A": expected_kind = EclipseKind.Annular; break;
                    case "T": expected_kind = EclipseKind.Total;   break;
                    default:
                        Console.WriteLine("C# LocalSolarEclipseTest2({0} line {1}): invalid eclipse type '{2}'", inFileName, lnum, typeCode);
                        return 1;
                    }

                    var observer = new Observer(latitude, longitude, 0.0);
                    AstroTime search_time = p1.AddDays(-20.0);
                    LocalSolarEclipseInfo eclipse = Astronomy.SearchLocalSolarEclipse(search_time, observer);
                    if (eclipse.kind != expected_kind)
                    {
                        Console.WriteLine("C# LocalSolarEclipseTest2({0} line {1}): expected {2}, found {3}", inFileName, lnum, expected_kind, eclipse.kind);
                        return 1;
                    }

                    if (CheckEvent(inFileName, lnum, "peak", peak, peakalt, eclipse.peak, ref max_minutes, ref max_degrees)) return 1;
                    if (CheckEvent(inFileName, lnum, "partial_begin", p1, p1alt, eclipse.partial_begin, ref max_minutes, ref max_degrees)) return 1;
                    if (CheckEvent(inFileName, lnum, "partial_end", p2, p2alt, eclipse.partial_end, ref max_minutes, ref max_degrees)) return 1;
                    if (typeCode != "P")
                    {
                        if (CheckEvent(inFileName, lnum, "total_begin", t1, t1alt, eclipse.total_begin, ref max_minutes, ref max_degrees)) return 1;
                        if (CheckEvent(inFileName, lnum, "total_end", t2, t2alt, eclipse.total_end, ref max_minutes, ref max_degrees)) return 1;
                    }

                    ++verify_count;
                }
            }

            Console.WriteLine("C# LocalSolarEclipseTest2: PASS ({0} verified, {1} CalcMoons, max minutes = {2}, max alt degrees = {3})",
                verify_count, Astronomy.CalcMoonCount, max_minutes, max_degrees);
            return 0;
        }

        static bool CheckEvent(
            string inFileName,
            int lnum,
            string name,
            AstroTime expected_time,
            double expected_altitude,
            EclipseEvent evt,
            ref double max_minutes,
            ref double max_degrees)
        {
            double diff_minutes = (24 * 60) * abs(expected_time.ut - evt.time.ut);
            if (diff_minutes > max_minutes)
                max_minutes = diff_minutes;

            if (diff_minutes > 1.0)
            {
                Console.WriteLine("CheckEvent({0} line {1}): EXCESSIVE TIME ERROR: {2} minutes", inFileName, lnum, diff_minutes);
                return true;
            }

            double diff_alt = abs(expected_altitude - evt.altitude);
            if (diff_alt > max_degrees) max_degrees = diff_alt;
            if (diff_alt > 0.5)
            {
                Console.WriteLine("CheckEvent({0} line {1}): EXCESSIVE ALTITUDE ERROR: {2} degrees", inFileName, lnum, diff_alt);
                return true;
            }

            return false;
        }

        static int LocalSolarEclipseTest()
        {
            if (0 != LocalSolarEclipseTest1())
                return 1;

            if (0 != LocalSolarEclipseTest2())
                return 1;

            return 0;
        }

        static int TransitFile(Body body, string filename, double limit_minutes, double limit_sep)
        {
            using (StreamReader infile = File.OpenText(filename))
            {
                string line;
                int lnum = 0;
                double max_minutes = 0.0;
                double max_sep = 0.0;
                TransitInfo transit = Astronomy.SearchTransit(body, new AstroTime(1600, 1, 1, 0, 0, 0));
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    // 22:17 1881-11-08T00:57Z 03:38  3.8633
                    string[] token = Tokenize(line);
                    if (token.Length != 4)
                    {
                        Console.WriteLine("C# TransitFile({0} line {1}): expected 4 tokens, found {2}", filename, lnum, token.Length);
                        return 1;
                    }

                    string textp = token[1];
                    string text1 = textp.Substring(0, 11) + token[0] + "Z";
                    string text2 = textp.Substring(0, 11) + token[2] + "Z";
                    AstroTime timep = ParseDate(textp);
                    AstroTime time1 = ParseDate(text1);
                    AstroTime time2 = ParseDate(text2);
                    double separation = double.Parse(token[3]);

                    // If the start time is after the peak time, it really starts on the previous day.
                    if (time1.ut > timep.ut)
                        time1 = time1.AddDays(-1.0);

                    // If the finish time is before the peak time, it really starts on the following day.
                    if (time2.ut < timep.ut)
                        time2 = time2.AddDays(+1.0);

                    double diff_start  = (24.0 * 60.0) * abs(time1.ut - transit.start.ut );
                    double diff_peak   = (24.0 * 60.0) * abs(timep.ut - transit.peak.ut  );
                    double diff_finish = (24.0 * 60.0) * abs(time2.ut - transit.finish.ut);
                    double diff_sep = abs(separation - transit.separation);

                    max_minutes = max(max_minutes, diff_start);
                    max_minutes = max(max_minutes, diff_peak);
                    max_minutes = max(max_minutes, diff_finish);
                    if (max_minutes > limit_minutes)
                    {
                        Console.WriteLine("C# TransitFile({0} line {1}): EXCESSIVE TIME ERROR = {2} minutes.", filename, lnum, max_minutes);
                        return 1;
                    }

                    max_sep = max(max_sep, diff_sep);
                    if (max_sep > limit_sep)
                    {
                        Console.WriteLine("C# TransitFile({0} line {1}): EXCESSIVE SEPARATION ERROR = {2} arcminutes.", filename, lnum, max_sep);
                        return 1;
                    }

                    transit = Astronomy.NextTransit(body, transit.finish);
                }
                Console.WriteLine("C# TransitFile({0}): PASS - verified {1}, max minutes = {2}, max sep arcmin = {3}", filename, lnum, max_minutes, max_sep);
                return 0;
            }
        }

        static int TransitTest()
        {
            if (0 != TransitFile(Body.Mercury, "../../eclipse/mercury.txt", 10.710, 0.2121))
                return 1;

            if (0 != TransitFile(Body.Venus, "../../eclipse/venus.txt", 9.109, 0.6772))
                return 1;

            return 0;
        }

        static int PlutoCheckDate(double ut, double arcmin_tolerance, double x, double y, double z)
        {
            var time = new AstroTime(ut);
            string timeText;

            try
            {
                timeText = time.ToString();
            }
            catch (ArgumentOutOfRangeException)
            {
                // One of the dates we pass in is before the year 0000.
                // The DateTime class can't handle this.
                timeText = "???";
            }

            Debug("C# PlutoCheck: {0} = {1} UT = {2} TT", timeText, time.ut, time.tt);

            AstroVector vector = Astronomy.HelioVector(Body.Pluto, time);
            double dx = v(vector.x - x);
            double dy = v(vector.y - y);
            double dz = v(vector.z - z);
            double diff = sqrt(dx*dx + dy*dy + dz*dz);
            double dist = (sqrt(x*x + y*y + z*z) - 1.0);       // worst-case distance between Pluto and Earth
            double arcmin = (diff / dist) * (180.0 * 60.0 / Math.PI);
            Debug("C# PlutoCheck: calc pos = [{0}, {1}, {2}]", vector.x, vector.y, vector.z);
            Debug("C# PlutoCheck: ref  pos = [{0}, {1}, {2}]", x, y, z);
            Debug("C# PlutoCheck: del  pos = [{0}, {1}, {2}]", vector.x - x, vector.y - y, vector.z - z);
            Debug("C# PlutoCheck: diff = {0} AU, {1} arcmin", diff, arcmin);
            Debug("");

            if (v(arcmin) > arcmin_tolerance)
            {
                Console.WriteLine("C# PlutoCheck: EXCESSIVE ERROR");
                return 1;
            }
            return 0;
        }

        static int PlutoCheck()
        {
            if (0 != PlutoCheckDate(  +18250.0,  0.089, +37.4377303523676090, -10.2466292454075898, -14.4773101310875809)) return 1;
            if (0 != PlutoCheckDate( -856493.0,  4.067, +23.4292113199166252, +42.1452685817740829,  +6.0580908436642940)) return 1;
            if (0 != PlutoCheckDate( +435633.0,  0.016, -27.3178902095231813, +18.5887022581070305, +14.0493896259306936)) return 1;
            if (0 != PlutoCheckDate(       0.0,   8e-9, -9.8753673425269000,  -27.9789270580402771,  -5.7537127596369588)) return 1;
            if (0 != PlutoCheckDate( +800916.0,  2.286, -29.5266052645301365, +12.0554287322176474, +12.6878484911631091)) return 1;
            Console.WriteLine("C# PlutoCheck: PASS");
            return 0;
        }

        static int GeoidTestCase(AstroTime time, Observer observer, EquatorEpoch epoch)
        {
            Equatorial topo_moon = Astronomy.Equator(Body.Moon, time, observer, epoch, Aberration.None);
            AstroVector surface = Astronomy.ObserverVector(time, observer, epoch);
            AstroVector geo_moon = Astronomy.GeoVector(Body.Moon, time, Aberration.None);

            if (epoch == EquatorEpoch.OfDate)
            {
                // Astronomy.GeoVector() returns J2000 coordinates. Convert to equator-of-date coordinates.
                RotationMatrix rot = Astronomy.Rotation_EQJ_EQD(time);
                geo_moon = Astronomy.RotateVector(rot, geo_moon);
            }

            double dx = Astronomy.KM_PER_AU * v((geo_moon.x - surface.x) - topo_moon.vec.x);
            double dy = Astronomy.KM_PER_AU * v((geo_moon.y - surface.y) - topo_moon.vec.y);
            double dz = Astronomy.KM_PER_AU * v((geo_moon.z - surface.z) - topo_moon.vec.z);
            double diff = sqrt(dx*dx + dy*dy + dz*dz);
            Debug("C# GeoidTestCase: epoch={0}, time={1}, lat={2}, lon={3}, ht={4}, surface=({5}, {6}, {7}), diff = {8} km",
                epoch,
                time,
                observer.latitude,
                observer.longitude,
                observer.height,
                Astronomy.KM_PER_AU * surface.x,
                Astronomy.KM_PER_AU * surface.y,
                Astronomy.KM_PER_AU * surface.z,
                diff);

            // Require 1 millimeter accuracy! (one millionth of a kilometer).
            if (diff > 1.0e-6)
            {
                Console.WriteLine("C# GeoidTestCase: EXCESSIVE POSITION ERROR.");
                return 1;
            }

            // Verify that we can convert the surface vector back to an observer.
            Observer vobs = Astronomy.VectorObserver(surface, epoch);
            double lat_diff = abs(vobs.latitude - observer.latitude);

            // Longitude is meaningless at the poles, so don't bother checking it there.
            double lon_diff;
            if (-89.99 <= observer.latitude && observer.latitude <= +89.99)
            {
                lon_diff = abs(vobs.longitude - observer.longitude);
                if (lon_diff > 180.0)
                    lon_diff = 360.0 - lon_diff;
                lon_diff = abs(lon_diff * cos(Astronomy.DEG2RAD * observer.latitude));
                if (lon_diff > 1.0e-6)
                {
                    Console.WriteLine("C# GeoidTestCase: EXCESSIVE longitude check error = {0}", lon_diff);
                    return 1;
                }
            }
            else
            {
                lon_diff = 0.0;
            }

            double h_diff = abs(vobs.height - observer.height);
            Debug("C# GeoidTestCase: vobs=(lat={0}, lon={1}, height={2}), lat_diff={3}, lon_diff={4}, h_diff={5}",
                vobs.latitude, vobs.longitude, vobs.height, lat_diff, lon_diff, h_diff);

            if (lat_diff > 1.0e-6)
            {
                Console.WriteLine("C# GeoidTestCase: EXCESSIVE latitude check error = {0}", lat_diff);
                return 1;
            }

            if (h_diff > 0.001)
            {
                Console.WriteLine("C# GeoidTestCase: EXCESSIVE height check error = {0}", h_diff);
                return 1;
            }

            return 0;
        }

        static int GeoidTest()
        {
            var time_list = new AstroTime[]
            {
                new AstroTime(1066,  9, 27, 18,  0,  0),
                new AstroTime(1970, 12, 13, 15, 42,  0),
                new AstroTime(1970, 12, 13, 15, 43,  0),
                new AstroTime(2015,  3,  5,  2, 15, 45)
            };

            var observer_list = new Observer[]
            {
                new Observer( +1.5,   +2.7,    7.4),
                new Observer(-53.7, +141.7, +100.0),
                new Observer(+30.0,  -85.2,  -50.0),
                new Observer(+90.0,  +45.0,  -50.0),
                new Observer(-90.0, -180.0,    0.0)
            };

            // Test a variety of times and locations, in both supported orientation systems.

            foreach (Observer observer in observer_list)
            {
                foreach (AstroTime time in time_list)
                {
                    if (0 != GeoidTestCase(time, observer, EquatorEpoch.J2000))
                        return 1;
                    if (0 != GeoidTestCase(time, observer, EquatorEpoch.OfDate))
                        return 1;
                }
            }

            var fixed_time = new AstroTime(2021, 6, 20, 15, 8, 0);
            for (int lat = -90; lat <= +90; lat += 1)
            {
                for (int lon = -175; lon <= +180; lon += 5)
                {
                    var observer = new Observer(lat, lon, 0.0);
                    if (0 != GeoidTestCase(fixed_time, observer, EquatorEpoch.OfDate))
                        return 1;
                }
            }

            Console.WriteLine("C# GeoidTest: PASS");
            return 0;
        }

        const int NUM_JUPITER_MOONS = 4;

        static StateVector SelectJupiterMoon(JupiterMoonsInfo jm, int mindex)
        {
            switch (mindex)
            {
                case 0: return jm.io;
                case 1: return jm.europa;
                case 2: return jm.ganymede;
                case 3: return jm.callisto;
                default: throw new ArgumentOutOfRangeException($"Invalid mindex = {mindex}");
            }
        }

        static int JupiterMoons_CheckJpl(int mindex, double tt, double[] pos, double[] vel)
        {
            const double pos_tolerance = 9.0e-4;
            const double vel_tolerance = 9.0e-4;
            AstroTime time = AstroTime.FromTerrestrialTime(tt);
            JupiterMoonsInfo jm = Astronomy.JupiterMoons(time);
            StateVector moon = SelectJupiterMoon(jm, mindex);

            double dx = v(pos[0] - moon.x);
            double dy = v(pos[1] - moon.y);
            double dz = v(pos[2] - moon.z);
            double mag = Math.Sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
            double pos_diff = Math.Sqrt(dx*dx + dy*dy + dz*dz) / mag;
            if (pos_diff > pos_tolerance)
            {
                Console.WriteLine($"C# JupiterMoons_CheckJpl(mindex={mindex}, tt={tt}): excessive position error {pos_diff}");
                return 1;
            }

            dx = v(vel[0] - moon.vx);
            dy = v(vel[1] - moon.vy);
            dz = v(vel[2] - moon.vz);
            mag = Math.Sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
            double vel_diff = Math.Sqrt(dx*dx + dy*dy + dz*dz) / mag;
            if (vel_diff > vel_tolerance)
            {
                Console.WriteLine($"C# JupiterMoons_CheckJpl(mindex={mindex}, tt={tt}): excessive velocity error {vel_diff}");
                return 1;
            }

            Debug($"C# JupiterMoons_CheckJpl: mindex={mindex}, tt={tt}, pos_diff={pos_diff}, vel_diff={vel_diff}");
            return 0;
        }

        static int JupiterMoonsTest()
        {
            const int expected_count = 5001;
            for (int mindex = 0; mindex < NUM_JUPITER_MOONS; ++mindex)
            {
                string filename = $"../../jupiter_moons/horizons/jm{mindex}.txt";
                using (StreamReader input = File.OpenText(filename))
                {
                    int lnum = 0;
                    bool found = false;
                    int part = -1;
                    int count = 0;
                    string line;
                    double tt = 1.0e+99;
                    var pos = new double[3];
                    var vel = new double[3];
                    while (null != (line = input.ReadLine()))
                    {
                        ++lnum;
                        if (!found)
                        {
                            if (line == "$$SOE")
                            {
                                found = true;
                                part = 0;
                            }
                            else if (line.StartsWith("Revised:"))
                            {
                                if (line.Length != 79)
                                {
                                    Console.WriteLine($"C# JupiterMoonsTest({filename} line {lnum}): unexpected line length.");
                                    return 1;
                                }
                                int check_mindex = int.Parse(line.Substring(76)) - 501;
                                if (mindex != check_mindex)
                                {
                                    Console.WriteLine($"C# JupiterMoonsTest({filename} line {lnum}): moon index does not match: check={check_mindex}, mindex={mindex}.");
                                    return 1;
                                }
                            }
                        }
                        else if (line == "$$EOE")
                        {
                            break;
                        }
                        else
                        {
                            Match match;
                            switch (part)
                            {
                                case 0:
                                    // 2446545.000000000 = A.D. 1986-Apr-24 12:00:00.0000 TDB
                                    tt = double.Parse(line.Split()[0]) - 2451545.0;    // convert JD to J2000 TT
                                    break;

                                case 1:
                                    // X = 1.134408131605554E-03 Y =-2.590904586750408E-03 Z =-7.490427225904720E-05
                                    match = Regex.Match(line, @"\s*X =\s*(\S+) Y =\s*(\S+) Z =\s*(\S+)");
                                    if (!match.Success)
                                    {
                                        Console.WriteLine($"C# JupiterMoonsTest({filename} line {lnum}): cannot parse position vector.");
                                        return 1;
                                    }
                                    pos[0] = double.Parse(match.Groups[1].Value);
                                    pos[1] = double.Parse(match.Groups[2].Value);
                                    pos[2] = double.Parse(match.Groups[3].Value);
                                    break;

                                case 2:
                                    // VX= 9.148038778472862E-03 VY= 3.973823407182510E-03 VZ= 2.765660368640458E-04
                                    match = Regex.Match(line, @"\s*VX=\s*(\S+) VY=\s*(\S+) VZ=\s*(\S+)");
                                    if (!match.Success)
                                    {
                                        Console.WriteLine($"C# JupiterMoonsTest({filename} line {lnum}): cannot parse velocity vector.");
                                        return 1;
                                    }
                                    vel[0] = double.Parse(match.Groups[1].Value);
                                    vel[1] = double.Parse(match.Groups[2].Value);
                                    vel[2] = double.Parse(match.Groups[3].Value);
                                    if (0 != JupiterMoons_CheckJpl(mindex, tt, pos, vel))
                                    {
                                        Console.WriteLine($"C# JupiterMoonsTest({filename} line {lnum}): FAILED VERIFICATION.");
                                        return 1;
                                    }
                                    ++count;
                                    break;

                                default:
                                    Console.WriteLine($"C# JupiterMoonsTest({filename} line {lnum}): unexpected part = {part}.");
                                    return 1;
                            }
                            part = (part + 1) % 3;
                        }
                    }
                    if (count != expected_count)
                    {
                        Console.WriteLine($"C# JupiterMoonsTest({filename}): expected {expected_count} test cases, found {count}.");
                        return 1;
                    }
                }
            }
            return 0;
        }

        static double StateVectorDiff(bool relative, double[] vec, double x, double y, double z)
        {
            double dx = vec[0] - x;
            double dy = vec[1] - y;
            double dz = vec[2] - z;
            double diff_squared = dx*dx + dy*dy + dz*dz;
            if (relative)
                diff_squared /= (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
            return sqrt(diff_squared);
        }

        interface IStateVectorFunc
        {
            StateVector Eval(AstroTime time);
        }

        static int VerifyState(
            IStateVectorFunc func,
            ref double max_rdiff,
            ref double max_vdiff,
            string filename,
            int lnum,
            AstroTime time,
            double[] pos,
            double[] vel,
            double r_thresh,
            double v_thresh)
        {
            StateVector state = func.Eval(time);

            double rdiff = StateVectorDiff((r_thresh > 0.0), pos, state.x, state.y, state.z);
            if (rdiff > max_rdiff)
                max_rdiff = rdiff;

            double vdiff = StateVectorDiff((v_thresh > 0.0), vel, state.vx, state.vy, state.vz);
            if (vdiff > max_vdiff)
                max_vdiff = vdiff;

            if (rdiff > Math.Abs(r_thresh))
            {
                Console.WriteLine($"C# VerifyState({filename} line {lnum}): EXCESSIVE position error = {rdiff:E3}");
                return 1;
            }

            if (vdiff > Math.Abs(v_thresh))
            {
                Console.WriteLine($"C# VerifyState({filename} line {lnum}): EXCESSIVE velocity error = {vdiff:E3}");
                return 1;
            }

            return 0;
        }

        static int VerifyStateBody(
            IStateVectorFunc func,
            string filename,
            double r_thresh,
            double v_thresh)
        {
            double max_rdiff = 0.0, max_vdiff = 0.0;
            var pos = new double[3];
            var vel = new double[3];
            int count = 0;
            foreach (JplStateRecord rec in JplHorizonsStateVectors(filename))
            {
                pos[0] = rec.state.x;
                pos[1] = rec.state.y;
                pos[2] = rec.state.z;
                vel[0] = rec.state.vx;
                vel[1] = rec.state.vy;
                vel[2] = rec.state.vz;
                if (0 != VerifyState(func, ref max_rdiff, ref max_vdiff, filename, rec.lnum, rec.state.t, pos, vel, r_thresh, v_thresh))
                    return 1;
                ++count;
            }
            Debug($"C# VerifyStateBody({filename}): PASS - Tested {count} cases. max rdiff={max_rdiff:E3}, vdiff={max_vdiff:E3}");
            return 0;
        }

        // Constants for use inside unit tests only; they doesn't make sense for public consumption.
        const Body Body_GeoMoon = (Body)(-100);
        const Body Body_Geo_EMB = (Body)(-101);

        class BaryStateFunc : IStateVectorFunc
        {
            private Body body;

            public BaryStateFunc(Body body)
            {
                this.body = body;
            }

            public StateVector Eval(AstroTime time)
            {
                if (body == Body_GeoMoon)
                    return Astronomy.GeoMoonState(time);

                if (body == Body_Geo_EMB)
                    return  Astronomy.GeoEmbState(time);

                return Astronomy.BaryState(body, time);
            }
        }

        static int BaryStateTest()
        {
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Sun),     "../../barystate/Sun.txt",     -1.224e-05, -1.134e-07)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Mercury), "../../barystate/Mercury.txt",  1.672e-04,  2.698e-04)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Venus),   "../../barystate/Venus.txt",    4.123e-05,  4.308e-05)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Earth),   "../../barystate/Earth.txt",    2.296e-05,  6.359e-05)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Mars),    "../../barystate/Mars.txt",     3.107e-05,  5.550e-05)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Jupiter), "../../barystate/Jupiter.txt",  7.389e-05,  2.471e-04)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Saturn),  "../../barystate/Saturn.txt",   1.067e-04,  3.220e-04)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Uranus),  "../../barystate/Uranus.txt",   9.035e-05,  2.519e-04)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Neptune), "../../barystate/Neptune.txt",  9.838e-05,  4.446e-04)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Pluto),   "../../barystate/Pluto.txt",    4.259e-05,  7.827e-05)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.Moon),    "../../barystate/Moon.txt",     2.354e-05,  6.604e-05)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body.EMB),     "../../barystate/EMB.txt",      2.353e-05,  6.511e-05)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body_GeoMoon), "../../barystate/GeoMoon.txt",  4.086e-05,  5.347e-05)) return 1;
            if (0 != VerifyStateBody(new BaryStateFunc(Body_Geo_EMB), "../../barystate/GeoEMB.txt",   4.076e-05,  5.335e-05)) return 1;
            Console.WriteLine("C# BaryStateTest: PASS");
            return 0;
        }

        class HelioStateFunc : IStateVectorFunc
        {
            private Body body;

            public HelioStateFunc(Body body)
            {
                this.body = body;
            }

            public StateVector Eval(AstroTime time)
            {
                return Astronomy.HelioState(body, time);
            }
        }

        static int HelioStateTest()
        {
            if (0 != VerifyStateBody(new HelioStateFunc(Body.SSB),     "../../heliostate/SSB.txt",     -1.209e-05, -1.125e-07)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Mercury), "../../heliostate/Mercury.txt",  1.481e-04,  2.756e-04)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Venus),   "../../heliostate/Venus.txt",    3.528e-05,  4.485e-05)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Earth),   "../../heliostate/Earth.txt",    1.476e-05,  6.105e-05)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Mars),    "../../heliostate/Mars.txt",     3.154e-05,  5.603e-05)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Jupiter), "../../heliostate/Jupiter.txt",  7.455e-05,  2.562e-04)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Saturn),  "../../heliostate/Saturn.txt",   1.066e-04,  3.150e-04)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Uranus),  "../../heliostate/Uranus.txt",   9.034e-05,  2.712e-04)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Neptune), "../../heliostate/Neptune.txt",  9.834e-05,  4.534e-04)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Pluto),   "../../heliostate/Pluto.txt",    4.271e-05,  1.198e-04)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.Moon),    "../../heliostate/Moon.txt",     1.477e-05,  6.195e-05)) return 1;
            if (0 != VerifyStateBody(new HelioStateFunc(Body.EMB),     "../../heliostate/EMB.txt",      1.476e-05,  6.106e-05)) return 1;
            Console.WriteLine("C# HelioStateTest: PASS");
            return 0;
        }

        class TopoStateFunc : IStateVectorFunc
        {
            private Body body;

            public TopoStateFunc(Body body)
            {
                this.body = body;
            }

            public StateVector Eval(AstroTime time)
            {
                var observer = new Observer(30.0, -80.0, 1000.0);

                StateVector observer_state = Astronomy.ObserverState(time, observer, EquatorEpoch.J2000);
                StateVector state;
                if (body == Body_Geo_EMB)
                {
                    state = Astronomy.GeoEmbState(time);
                }
                else if (body == Body.Earth)
                {
                    state = new StateVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time);
                }
                else
                {
                    throw new ArgumentException($"C# TopoStateFunction: unsupported body {body}");
                }

                state.x  -= observer_state.x;
                state.y  -= observer_state.y;
                state.z  -= observer_state.z;
                state.vx -= observer_state.vx;
                state.vy -= observer_state.vy;
                state.vz -= observer_state.vz;

                return state;
            }
        }

        static int TopoStateTest()
        {
            if (0 != VerifyStateBody(new TopoStateFunc(Body.Earth),   "../../topostate/Earth_N30_W80_1000m.txt",  2.108e-04, 2.430e-04)) return 1;
            if (0 != VerifyStateBody(new TopoStateFunc(Body_Geo_EMB), "../../topostate/EMB_N30_W80_1000m.txt",    7.195e-04, 2.497e-04)) return 1;
            Console.WriteLine("C# TopoStateTest: PASS");
            return 0;
        }

        static int AberrationTest()
        {
            const string filename = "../../equatorial/Mars_j2000_ofdate_aberration.txt";
            const double THRESHOLD_SECONDS = 0.4;

            using (StreamReader infile = File.OpenText(filename))
            {
                int lnum = 0;
                int count = 0;
                string line;
                bool found_begin = false;
                double max_diff_seconds = 0.0;
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    if (!found_begin)
                    {
                        if (line == "$$SOE")
                            found_begin = true;
                    }
                    else if (line == "$$EOE")
                    {
                        break;
                    }
                    else
                    {
                        // 2459371.500000000 *   118.566080210  22.210647456 118.874086738  22.155784122
                        string[] token = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                        if (token.Length < 5)
                        {
                            Console.WriteLine($"C# AberrationTest({filename} line {lnum}): not enough tokens");
                            return 1;
                        }

                        double jd = double.Parse(token[0]);
                        double jra = double.Parse(token[token.Length-4]);
                        double jdec = double.Parse(token[token.Length-3]);
                        double dra = double.Parse(token[token.Length-2]);
                        double ddec = double.Parse(token[token.Length-1]);

                        // Convert julian day value to AstroTime.
                        var time = new AstroTime(jd - 2451545.0);

                        // Convert EQJ angular coordinates (jra, jdec) to an EQJ vector.
                        // Make the maginitude of the vector the speed of light,
                        // to prepare for aberration correction.
                        var eqj_sphere = new Spherical(jdec, jra, Astronomy.C_AUDAY);
                        var eqj_vec = Astronomy.VectorFromSphere(eqj_sphere, time);

                        // Aberration correction: calculate the Earth's barycentric
                        // velocity vector in EQJ coordinates.
                        StateVector eqj_earth = Astronomy.BaryState(Body.Earth, time);

                        // Use non-relativistic approximation: add light vector to Earth velocity vector.
                        // This gives aberration-corrected apparent position of the start in EQJ.
                        eqj_vec.x += eqj_earth.vx;
                        eqj_vec.y += eqj_earth.vy;
                        eqj_vec.z += eqj_earth.vz;

                        // Calculate the rotation matrix that converts J2000 coordinates to of-date coordinates.
                        RotationMatrix rot = Astronomy.Rotation_EQJ_EQD(time);

                        // Use the rotation matrix to re-orient the EQJ vector to an EQD vector.
                        AstroVector eqd_vec = Astronomy.RotateVector(rot, eqj_vec);

                        // Convert the EQD vector back to spherical angular coordinates.
                        Spherical eqd_sphere = Astronomy.SphereFromVector(eqd_vec);

                        // Calculate the differences in RA and DEC between expected and calculated values.
                        double factor = cos(eqd_sphere.lat * Astronomy.DEG2RAD);    // RA errors are less important toward the poles.
                        double xra = factor * abs(eqd_sphere.lon - dra);
                        double xdec = abs(eqd_sphere.lat - ddec);
                        double diff_seconds = 3600.0 * sqrt(xra*xra + xdec*xdec);
                        Debug($"C# AberrationTest({filename} line {lnum}): xra={xra:F6} deg, xdec={xdec:F6} deg, diff_seconds={diff_seconds:F3}.");
                        if (diff_seconds > THRESHOLD_SECONDS)
                        {
                            Console.WriteLine($"C# AberrationTest({filename} line {lnum}): EXCESSIVE ANGULAR ERROR = {diff_seconds:F3} seconds.");
                            return 1;
                        }

                        if (diff_seconds > max_diff_seconds)
                            max_diff_seconds = diff_seconds;

                        // We have completed one more test case.
                        ++count;
                    }
                }
                Console.WriteLine($"C# AberrationTest({filename}): PASS - Tested {count} cases. max_diff_seconds = {max_diff_seconds:F3}");
            }

            return 0;
        }

        static int TwilightTest()
        {
            const string filename = "../../riseset/twilight.txt";
            const double tolerance_seconds = 60.0;
            double max_diff = 0.0;
            int lnum = 0;
            using (StreamReader infile = File.OpenText(filename))
            {
                string line;
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    string[] tokens = line.Split();
                    if (tokens.Length != 9)
                    {
                        Console.WriteLine($"C# TwilightTest({filename} line {lnum}): invalid number of tokens = {tokens.Length}");
                        return 1;
                    }

                    double lat = double.Parse(tokens[0]);
                    double lon = double.Parse(tokens[1]);
                    var observer = new Observer(lat, lon, 0.0);
                    AstroTime searchDate = ParseDate(tokens[2]);
                    AstroTime[] correctTimes = tokens.Skip(3).Select(s => ParseDate(s)).ToArray();
                    var calcTimes = new AstroTime[]
                    {
                        Astronomy.SearchAltitude(Body.Sun, observer, Direction.Rise, searchDate, 1.0, -18.0),  // astronomical dawn
                        Astronomy.SearchAltitude(Body.Sun, observer, Direction.Rise, searchDate, 1.0, -12.0),  // nautical dawn
                        Astronomy.SearchAltitude(Body.Sun, observer, Direction.Rise, searchDate, 1.0,  -6.0),  // civil dawn
                        Astronomy.SearchAltitude(Body.Sun, observer, Direction.Set,  searchDate, 1.0,  -6.0),  // civil dawn
                        Astronomy.SearchAltitude(Body.Sun, observer, Direction.Set,  searchDate, 1.0, -12.0),  // nautical dawn
                        Astronomy.SearchAltitude(Body.Sun, observer, Direction.Set,  searchDate, 1.0, -18.0),  // astronomical dawn
                    };

                    for (int i = 0; i < 6; ++i)
                    {
                        AstroTime correct = correctTimes[i];
                        AstroTime calc = calcTimes[i];
                        double diff = SECONDS_PER_DAY * abs(calc.ut - correct.ut);
                        if (diff > tolerance_seconds)
                        {
                            Console.WriteLine($"C# TwilightTest({filename} line {lnum}): EXCESSIVE ERROR = {diff} seconds in test {i}.");
                            return 1;
                        }
                        if (diff > max_diff)
                            max_diff = diff;
                    }
                }
            }
            Console.WriteLine($"C# TwilightTest: PASS ({lnum} test cases, max error = {max_diff} seconds)");
            return 0;
        }

        static int Libration(string filename)
        {
            using StreamReader infile = File.OpenText(filename);
            int lnum = 0;
            int count = 0;
            double max_diff_elon = 0.0;
            double max_diff_elat = 0.0;
            double max_diff_distance = 0.0;
            double max_diff_diam = 0.0;
            double max_eclip_lon = -900.0;
            string line;
            while (null != (line = infile.ReadLine()))
            {
                ++lnum;
                if (lnum == 1)
                {
                    //              0..2       3..4       5      6      7        8       9        10        11       12       13        14     15
                    if (line != "   Date       Time    Phase    Age    Diam    Dist     RA        Dec      Slon      Slat     Elon     Elat   AxisA")
                    {
                        Console.WriteLine($"C# Libration({filename} line {lnum}): unexpected header line.");
                        return 1;
                    }
                }
                else
                {
                    //  0   1   2    3    4    5       6      7       8        9         10       11         12      13       14      15
                    // 01 Jan 2020 00:00 UT  29.95   5.783  1774.5  403898  23.2609  -10.0824   114.557   -0.045   0.773    6.360  336.353
                    string[] token = Tokenize(line);
                    if (token.Length != 16)
                    {
                        Console.WriteLine($"C# Libration({filename} line {lnum}): expected 16 tokens, found {token.Length}.");
                        return 1;
                    }

                    int day = int.Parse(token[0]);
                    int month = MonthNumber(token[1]);
                    int year = int.Parse(token[2]);
                    string[] hmtoken = token[3].Split(':');
                    if (hmtoken.Length != 2)
                    {
                        Console.WriteLine($"C# Libration({filename} line {lnum}): expected hh:mm but found '{token[3]}'");
                        return 1;
                    }
                    int hour = int.Parse(hmtoken[0]);
                    int minute = int.Parse(hmtoken[1]);
                    var time = new AstroTime(year, month, day, hour, minute, 0);

                    double diam = double.Parse(token[7]) / 3600.0;
                    double dist = double.Parse(token[8]);
                    double elon = double.Parse(token[13]);
                    double elat = double.Parse(token[14]);

                    LibrationInfo lib = Astronomy.Libration(time);

                    double diff_elon = 60.0 * abs(lib.elon - elon);
                    if (diff_elon > max_diff_elon)
                        max_diff_elon = diff_elon;

                    double diff_elat = 60.0 * abs(lib.elat - elat);
                    if (diff_elat > max_diff_elat)
                        max_diff_elat = diff_elat;

                    double diff_distance = abs(lib.dist_km - dist);
                    if (diff_distance > max_diff_distance)
                        max_diff_distance = diff_distance;

                    double diff_diam = abs(lib.diam_deg - diam);
                    if (diff_diam > max_diff_diam)
                        max_diff_diam = diff_diam;

                    if (diff_elon > 0.1304)
                    {
                        Console.WriteLine($"C# Libration({filename} line {lnum}): EXCESSIVE diff_elon = {diff_elon} arcmin");
                        return 1;
                    }

                    if (diff_elat > 1.6476)
                    {
                        Console.WriteLine($"C# Libration({filename} line {lnum}): EXCESSIVE diff_elat = {diff_elat} arcmin");
                        return 1;
                    }

                    if (diff_distance > 54.377)
                    {
                        Console.WriteLine($"C# Libration({filename} line {lnum}): EXCESSIVE diff_distance = {diff_distance} km");
                        return 1;
                    }

                    if (lib.mlon > max_eclip_lon)
                        max_eclip_lon = lib.mlon;

                    ++count;
                }
            }

            if (max_eclip_lon < 359.0 || max_eclip_lon > 360.0)
            {
                Console.WriteLine($"C# Libration({filename}): INVALID max ecliptic longitude {max_eclip_lon:F3} degrees.");
                return 1;
            }

            Console.WriteLine($"C# Libration({filename}): PASS ({count} test cases, max_diff_elon = {max_diff_elon} arcmin, max_diff_elat = {max_diff_elat} arcmin, max_diff_distance = {max_diff_distance} km, max_diff_diam = {max_diff_diam} deg)");
            return 0;
        }

        static int LibrationTest()
        {
            if (0 != Libration("../../libration/mooninfo_2020.txt")) return 1;
            if (0 != Libration("../../libration/mooninfo_2021.txt")) return 1;
            if (0 != Libration("../../libration/mooninfo_2022.txt")) return 1;
            return 0;
        }

        static int AxisTest()
        {
            if (0 != AxisTestBody(Body.Sun,      "../../axis/Sun.txt",       0.0))        return 1;
            if (0 != AxisTestBody(Body.Mercury,  "../../axis/Mercury.txt",   0.074340))   return 1;
            if (0 != AxisTestBody(Body.Venus,    "../../axis/Venus.txt",     0.0))        return 1;
            if (0 != AxisTestBody(Body.Earth,    "../../axis/Earth.txt",     0.000591))   return 1;
            if (0 != AxisTestBody(Body.Moon,     "../../axis/Moon.txt",      0.264845))   return 1;
            if (0 != AxisTestBody(Body.Mars,     "../../axis/Mars.txt",      0.075323))   return 1;
            if (0 != AxisTestBody(Body.Jupiter,  "../../axis/Jupiter.txt",   0.000324))   return 1;
            if (0 != AxisTestBody(Body.Saturn,   "../../axis/Saturn.txt",    0.000304))   return 1;
            if (0 != AxisTestBody(Body.Uranus,   "../../axis/Uranus.txt",    0.0))        return 1;
            if (0 != AxisTestBody(Body.Neptune,  "../../axis/Neptune.txt",   0.000462))   return 1;
            if (0 != AxisTestBody(Body.Pluto,    "../../axis/Pluto.txt",     0.0))        return 1;
            Console.WriteLine("C# AxisTest: PASS");
            return 0;
        }

        static int AxisTestBody(Body body, string filename, double arcmin_tolerance)
        {
            double max_arcmin = 0.0;
            int count = 0;
            using (StreamReader infile = File.OpenText(filename))
            {
                bool found_data = false;
                int lnum = 0;
                string line;
                while (null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    if (!found_data)
                    {
                        if (line == "$$SOE")
                            found_data = true;
                    }
                    else
                    {
                        if (line == "$$EOE")
                            break;

                        if (line.Length < 61)
                        {
                            Console.WriteLine($"C# AxisBodyTest({filename} line {lnum}): line is too short.");
                            return 1;
                        }

                        string[] token = Tokenize(line.Substring(19));
                        if (token.Length != 3)
                        {
                            Console.WriteLine($"C# AxisBodyTest({filename} line {lnum}): expected 3 tokens but found {token.Length}.");
                            return 1;
                        }
                        if (!double.TryParse(token[0], out double jd) ||
                            !double.TryParse(token[1], out double ra) ||
                            !double.TryParse(token[2], out double dec))
                        {
                            Console.WriteLine($"C# AxisBodyTest({filename} line {lnum}): error parsing floating point numbers.");
                            return 1;
                        }

                        var time = new AstroTime(jd - 2451545.0);
                        AxisInfo axis = Astronomy.RotationAxis(body, time);

                        // Convert the reference angles to a reference north pole vector.
                        // tricky: `ra` is in degrees, not sidereal hours; so don't multiply by 15.
                        var sphere = new Spherical(dec, ra, 1.0);
                        AstroVector north = Astronomy.VectorFromSphere(sphere, time);

                        // Find angle between two versions of the north pole. Use that as the measure of error.
                        double arcmin = 60.0 * Astronomy.AngleBetween(north, axis.north);
                        if (arcmin > max_arcmin)
                            max_arcmin = arcmin;

                        ++count;
                    }
                }
            }
            Debug($"C# AxisTestBody({body}): {count} test cases, max arcmin error = {max_arcmin}.");
            if (max_arcmin > arcmin_tolerance)
            {
                Console.WriteLine($"C AxisTestBody({body}): EXCESSIVE ERROR = {max_arcmin} arcmin.");
                return 1;
            }
            return 0;
        }

        //-----------------------------------------------------------------------------------------

        internal struct JplStateRecord
        {
            public int lnum;            // the line number where the state vector ends in the JPL Horizons text file.
            public StateVector state;   // the state vector itself: position, velocity, and time.
        }

        static IEnumerable<JplStateRecord> JplHorizonsStateVectors(string filename)
        {
            using (StreamReader infile = File.OpenText(filename))
            {
                int lnum = 0;
                string line;
                bool found_begin = false;
                bool found_end = false;
                int part = 0;
                AstroTime time = null;
                var pos = new double[3];
                var vel = new double[3];
                while (!found_end && null != (line = infile.ReadLine()))
                {
                    ++lnum;
                    if (!found_begin)
                    {
                        if (line == "$$SOE")
                            found_begin = true;
                    }
                    else
                    {
                        // Input comes in triplets of lines:
                        //
                        // 2444249.500000000 = A.D. 1980-Jan-11 00:00:00.0000 TDB
                        // X =-3.314860345089456E-01 Y = 8.463418210972562E-01 Z = 3.667227830514760E-01
                        // VX=-1.642704711077836E-02 VY=-5.494770742558920E-03 VZ=-2.383170237527642E-03
                        //
                        // Track which of these 3 cases we are in using the 'part' variable...
                        Match match;
                        switch (part)
                        {
                            case 0:
                                if (line == "$$EOE")
                                {
                                    found_end = true;
                                }
                                else
                                {
                                    // 2444249.500000000 = A.D. 1980-Jan-11 00:00:00.0000 TDB
                                    // Convert JD to J2000 TT.
                                    double tt = double.Parse(line.Split()[0]) - 2451545.0;
                                    time = AstroTime.FromTerrestrialTime(tt);
                                }
                                break;

                            case 1:
                                // X = 1.134408131605554E-03 Y =-2.590904586750408E-03 Z =-7.490427225904720E-05
                                match = Regex.Match(line, @"\s*X =\s*(\S+) Y =\s*(\S+) Z =\s*(\S+)");
                                if (!match.Success)
                                    throw new Exception($"C# JplHorizonsStateVectors({filename} line {lnum}): cannot parse position vector.");
                                pos[0] = double.Parse(match.Groups[1].Value);
                                pos[1] = double.Parse(match.Groups[2].Value);
                                pos[2] = double.Parse(match.Groups[3].Value);
                                break;

                            case 2:
                                // VX= 9.148038778472862E-03 VY= 3.973823407182510E-03 VZ= 2.765660368640458E-04
                                match = Regex.Match(line, @"\s*VX=\s*(\S+) VY=\s*(\S+) VZ=\s*(\S+)");
                                if (!match.Success)
                                    throw new Exception($"C# JplHorizonsStateVectors({filename} line {lnum}): cannot parse velocity vector.");
                                vel[0] = double.Parse(match.Groups[1].Value);
                                vel[1] = double.Parse(match.Groups[2].Value);
                                vel[2] = double.Parse(match.Groups[3].Value);
                                var state = new StateVector(
                                    pos[0], pos[1], pos[2],
                                    vel[0], vel[1], vel[2],
                                    time
                                );
                                yield return new JplStateRecord { lnum = lnum, state = state };
                                break;

                            default:
                                throw new Exception($"C# JplHorizonsStateVectors({filename} line {lnum}): unexpected part = {part}.");
                        }
                        part = (part + 1) % 3;
                    }
                }
                yield break;
            }
        }

        //-----------------------------------------------------------------------------------------

        class LagrangeFunc : IStateVectorFunc
        {
            private int point;
            private Body major_body;
            private Body minor_body;

            public LagrangeFunc(int point, Body major_body, Body minor_body)
            {
                this.point = point;
                this.major_body = major_body;
                this.minor_body = minor_body;
            }

            public StateVector Eval(AstroTime time)
            {
                return Astronomy.LagrangePoint(point, time, major_body, minor_body);
            }
        }

        static int VerifyStateLagrange(
            Body major_body,
            Body minor_body,
            int point,
            string filename,
            double r_thresh,
            double v_thresh)
        {
            var func = new LagrangeFunc(point, major_body, minor_body);
            return VerifyStateBody(func, filename, r_thresh, v_thresh);
        }

        static int LagrangeTest()
        {
            // Test Sun/EMB Lagrange points.
            if (0 != VerifyStateLagrange(Body.Sun, Body.EMB, 1, "../../lagrange/semb_L1.txt",   1.33e-5, 6.13e-5)) return 1;
            if (0 != VerifyStateLagrange(Body.Sun, Body.EMB, 2, "../../lagrange/semb_L2.txt",   1.33e-5, 6.13e-5)) return 1;
            if (0 != VerifyStateLagrange(Body.Sun, Body.EMB, 4, "../../lagrange/semb_L4.txt",   3.75e-5, 5.28e-5)) return 1;
            if (0 != VerifyStateLagrange(Body.Sun, Body.EMB, 5, "../../lagrange/semb_L5.txt",   3.75e-5, 5.28e-5)) return 1;

            // Test Earth/Moon Lagrange points.
            if (0 != VerifyStateLagrange(Body.Earth, Body.Moon, 1, "../../lagrange/em_L1.txt",  3.79e-5, 5.06e-5)) return 1;
            if (0 != VerifyStateLagrange(Body.Earth, Body.Moon, 2, "../../lagrange/em_L2.txt",  3.79e-5, 5.06e-5)) return 1;
            if (0 != VerifyStateLagrange(Body.Earth, Body.Moon, 4, "../../lagrange/em_L4.txt",  3.79e-5, 1.59e-3)) return 1;
            if (0 != VerifyStateLagrange(Body.Earth, Body.Moon, 5, "../../lagrange/em_L5.txt",  3.79e-5, 1.59e-3)) return 1;

            Console.WriteLine("C# LagrangeTest: PASS");
            return 0;   // not yet implemented
        }

        //-----------------------------------------------------------------------------------------

        static int SiderealTimeTest()
        {
            const double correct = 9.398368460418821;
            var time = new AstroTime(2022, 3, 15, 21, 50, 0);
            double gast = Astronomy.SiderealTime(time);
            double diff = abs(gast - correct);
            Console.WriteLine($"C# SiderealTimeTest: gast={gast:F10}, correct={correct:F10}, diff={diff:E3}.");
            if (diff > 1.0e-15)
            {
                Console.WriteLine("C# SiderealTimeTest: EXCESSIVE ERROR");
                return 1;
            }
            Console.WriteLine("C# SiderealTimeTest: PASS");
            return 0;
        }

        //-----------------------------------------------------------------------------------------

        static int GravitySimulatorTest()
        {
            Debug("");

            if (0 != GravSimEmpty("barystate/Sun.txt",      Body.SSB, Body.Sun,      0.0269, 1.9635)) return 1;
            if (0 != GravSimEmpty("barystate/Mercury.txt",  Body.SSB, Body.Mercury,  0.5725, 0.9332)) return 1;
            if (0 != GravSimEmpty("barystate/Venus.txt",    Body.SSB, Body.Venus,    0.1433, 0.1458)) return 1;
            if (0 != GravSimEmpty("barystate/Earth.txt",    Body.SSB, Body.Earth,    0.0651, 0.2098)) return 1;
            if (0 != GravSimEmpty("barystate/Mars.txt",     Body.SSB, Body.Mars,     0.1150, 0.1896)) return 1;
            if (0 != GravSimEmpty("barystate/Jupiter.txt",  Body.SSB, Body.Jupiter,  0.2546, 0.8831)) return 1;
            if (0 != GravSimEmpty("barystate/Saturn.txt",   Body.SSB, Body.Saturn,   0.3660, 1.0818)) return 1;
            if (0 != GravSimEmpty("barystate/Uranus.txt",   Body.SSB, Body.Uranus,   0.3107, 0.9321)) return 1;
            if (0 != GravSimEmpty("barystate/Neptune.txt",  Body.SSB, Body.Neptune,  0.3382, 1.5586)) return 1;

            if (0 != GravSimEmpty("heliostate/Mercury.txt", Body.Sun, Body.Mercury,  0.5087, 0.9473)) return 1;
            if (0 != GravSimEmpty("heliostate/Venus.txt",   Body.Sun, Body.Venus,    0.1214, 0.1543)) return 1;
            if (0 != GravSimEmpty("heliostate/Earth.txt",   Body.Sun, Body.Earth,    0.0508, 0.2099)) return 1;
            if (0 != GravSimEmpty("heliostate/Mars.txt",    Body.Sun, Body.Mars,     0.1085, 0.1927)) return 1;
            if (0 != GravSimEmpty("heliostate/Jupiter.txt", Body.Sun, Body.Jupiter,  0.2564, 0.8805)) return 1;
            if (0 != GravSimEmpty("heliostate/Saturn.txt",  Body.Sun, Body.Saturn,   0.3664, 1.0826)) return 1;
            if (0 != GravSimEmpty("heliostate/Uranus.txt",  Body.Sun, Body.Uranus,   0.3106, 0.9322)) return 1;
            if (0 != GravSimEmpty("heliostate/Neptune.txt", Body.Sun, Body.Neptune,  0.3381, 1.5584)) return 1;

            Debug("");
            const int nsteps = 20;

            if (0 != GravSimFile("barystate/Ceres.txt",    Body.SSB,   nsteps, 0.6640, 0.6226)) return 1;
            if (0 != GravSimFile("barystate/Pallas.txt",   Body.SSB,   nsteps, 0.4687, 0.3474)) return 1;
            if (0 != GravSimFile("barystate/Vesta.txt",    Body.SSB,   nsteps, 0.5806, 0.5462)) return 1;
            if (0 != GravSimFile("barystate/Juno.txt",     Body.SSB,   nsteps, 0.6760, 0.5750)) return 1;
            if (0 != GravSimFile("barystate/Bennu.txt",    Body.SSB,   nsteps, 3.7444, 2.6581)) return 1;
            if (0 != GravSimFile("barystate/Halley.txt",   Body.SSB,   nsteps, 0.0539, 0.0825)) return 1;

            if (0 != GravSimFile("heliostate/Ceres.txt",   Body.Sun,   nsteps, 0.0445, 0.0355)) return 1;
            if (0 != GravSimFile("heliostate/Pallas.txt",  Body.Sun,   nsteps, 0.1062, 0.0854)) return 1;
            if (0 != GravSimFile("heliostate/Vesta.txt",   Body.Sun,   nsteps, 0.1432, 0.1308)) return 1;
            if (0 != GravSimFile("heliostate/Juno.txt",    Body.Sun,   nsteps, 0.1554, 0.1328)) return 1;

            if (0 != GravSimFile("geostate/Ceres.txt",     Body.Earth, nsteps, 6.5689, 6.4797)) return 1;
            if (0 != GravSimFile("geostate/Pallas.txt",    Body.Earth, nsteps, 9.3288, 7.3533)) return 1;
            if (0 != GravSimFile("geostate/Vesta.txt",     Body.Earth, nsteps, 3.2980, 3.8863)) return 1;
            if (0 != GravSimFile("geostate/Juno.txt",      Body.Earth, nsteps, 6.0962, 7.7147)) return 1;

            Debug("");
            Console.WriteLine("C# GravitySimulatorTest: PASS");

            return 0;
        }

        static int GravSimFile(string fileNameSuffix, Body originBody, int nsteps, double rthresh, double vthresh)
        {
            string filename = "../../" + fileNameSuffix;
            GravitySimulator sim = null;    // can't create until we see the first state vector.
            JplStateRecord prev = new JplStateRecord();
            AstroTime time = null;
            var smallBodyArray = new StateVector[1];
            double max_rdiff = 0.0;
            double max_vdiff = 0.0;
            foreach (JplStateRecord rec in JplHorizonsStateVectors(filename))
            {
                if (sim == null)
                {
                    sim = new GravitySimulator(originBody, rec.state.t, new StateVector[] { rec.state });
                    time = rec.state.t;
                }
                else
                {
                    double tt1 = prev.state.t.tt;
                    double tt2 = rec.state.t.tt;
                    double dt = (tt2 - tt1) / nsteps;
                    for (int k = 1; k <= nsteps; ++k)
                    {
                        time = AstroTime.FromTerrestrialTime(tt1 + k*dt);
                        sim.Update(time, null);     // confirm null works to indicate no output state vectors are desired.
                        if (time.tt != sim.Time.tt)
                        {
                            Console.WriteLine($"C# GravSimFile({filename} line {rec.lnum}): expected time {time} but simulator reports {sim.Time}.");
                            return 1;
                        }
                    }
                    // Confirm we can set to the same time and request output parameters.
                    sim.Update(time, smallBodyArray);

                    double rdiff = ArcminPosError(rec.state, smallBodyArray[0]);
                    if (rdiff > rthresh)
                    {
                        Console.WriteLine($"C# GravSimFile({filename} line {rec.lnum}): excessive position error = {rdiff} arcmin.");
                        return 1;
                    }
                    if (rdiff > max_rdiff)
                        max_rdiff = rdiff;

                    double vdiff = ArcminVelError(rec.state, smallBodyArray[0]);
                    if (vdiff > vthresh)
                    {
                        Console.WriteLine($"C# GravSimFile({filename} line {rec.lnum}): excessive velocity error = {vdiff} arcmin.");
                        return 1;
                    }
                    if (vdiff > max_vdiff)
                        max_vdiff = vdiff;
                }
                prev = rec;
            }

            Debug($"C# GravSimFile ({filename,-28}): PASS - max pos error = {max_rdiff:F4} arcmin, max vel error = {max_vdiff:F4} arcmin.");
            return 0;
        }

        static int GravSimEmpty(string fileNameSuffix, Body origin, Body body, double rthresh, double vthresh)
        {
            string filename = "../../" + fileNameSuffix;
            GravitySimulator sim = null;
            double max_rdiff = 0.0;
            double max_vdiff = 0.0;
            foreach (JplStateRecord rec in JplHorizonsStateVectors(filename))
            {
                if (sim == null)
                    sim = new GravitySimulator(origin, rec.state.t, new StateVector[0]);

                sim.Update(rec.state.t, null);
                StateVector calc = sim.SolarSystemBodyState(body);

                double rdiff = (
                    (origin==Body.SSB && body==Body.Sun)
                    ? SsbArcminPosError(rec.state, calc)
                    : ArcminPosError(rec.state, calc)
                );

                if (rdiff > rthresh)
                {
                    Console.WriteLine($"C# GravSimEmpty({filename} line {rec.lnum}): excessive position error = {rdiff} arcmin.");
                    return 1;
                }
                if (rdiff > max_rdiff)
                    max_rdiff = rdiff;

                double vdiff = ArcminVelError(rec.state, calc);
                if (vdiff > vthresh)
                {
                    Console.WriteLine($"C# GravSimEmpty({filename} line {rec.lnum}): excessive velocity error = {vdiff} arcmin.");
                    return 1;
                }
                if (vdiff > max_vdiff)
                    max_vdiff = vdiff;
            }

            Debug($"C# GravSimEmpty({filename,-28}): PASS - max pos error = {max_rdiff:F4} arcmin, max vel error = {max_vdiff:F4} arcmin.");
            return 0;
        }

        static double SsbArcminPosError(StateVector correct, StateVector calc)
        {
            // Scale the SSB based on 1 AU, not on its absolute magnitude, which can become very close to zero.
            double dx = calc.x - correct.x;
            double dy = calc.y - correct.y;
            double dz = calc.z - correct.z;
            double diffSquared = dx*dx + dy*dy + dz*dz;
            double radians = sqrt(diffSquared);
            return (60.0 * Astronomy.RAD2DEG) * radians;
        }

        static double ArcminPosError(StateVector correct, StateVector calc)
        {
            double dx = calc.x - correct.x;
            double dy = calc.y - correct.y;
            double dz = calc.z - correct.z;
            double diffSquared = dx*dx + dy*dy + dz*dz;
            double magSquared = correct.x*correct.x + correct.y*correct.y + correct.z*correct.z;
            double radians = sqrt(diffSquared / magSquared);
            return (60.0 * Astronomy.RAD2DEG) * radians;
        }

        static double ArcminVelError(StateVector correct, StateVector calc)
        {
            double dx = calc.vx - correct.vx;
            double dy = calc.vy - correct.vy;
            double dz = calc.vz - correct.vz;
            double diffSquared = dx*dx + dy*dy + dz*dz;
            double magSquared = correct.vx*correct.vx + correct.vy*correct.vy + correct.vz*correct.vz;
            double radians = sqrt(diffSquared / magSquared);
            return (60.0 * Astronomy.RAD2DEG) * radians;
        }

        //-----------------------------------------------------------------------------------------
    }
}
