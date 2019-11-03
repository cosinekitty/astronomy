using System;
using System.IO;
using System.Text.RegularExpressions;

using CosineKitty;

namespace csharp_test
{
    class Program
    {
        static int Main(string[] args)
        {
            try
            {
                Console.WriteLine("csharp_test: starting");
                if (TestTime() != 0) return 1;
                if (MoonTest() != 0) return 1;
                if (RiseSetTest("../../riseset/riseset.txt") != 0) return 1;
                if (AstroCheck() != 0) return 1;
                if (SeasonsTest("../../seasons/seasons.txt") != 0) return 1;
                if (MoonPhaseTest("../../moonphase/moonphases.txt") != 0) return 1;
                Console.WriteLine("csharp_test: PASS");
                return 0;
            }
            catch (Exception ex)
            {
                Console.WriteLine("charp_test: EXCEPTION: {0}", ex);
                return 1;
            }
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
            Console.WriteLine("TestTime: text={0}, ut={1}, tt={2}", time.ToString(), time.ut.ToString("F6"), time.tt.ToString("F6"));

            const double expected_ut = 6910.270978506945;
            double diff = time.ut - expected_ut;
            if (Math.Abs(diff) > 1.0e-12)
            {
                Console.WriteLine("TestTime: ERROR - excessive UT error {0}", diff);
                return 1;
            }

            const double expected_tt = 6910.271779431480;
            diff = time.tt - expected_tt;
            if (Math.Abs(diff) > 1.0e-12)
            {
                Console.WriteLine("TestTime: ERROR - excessive TT error {0}", diff);
                return 1;
            }

            DateTime utc = time.ToUtcDateTime();
            if (utc.Year != year || utc.Month != month || utc.Day != day || utc.Hour != hour || utc.Minute != minute || utc.Second != second || utc.Millisecond != milli)
            {
                Console.WriteLine("TestTime: ERROR - Expected {0:o}, found {1:o}", d, utc);
                return 1;
            }

            return 0;
        }

        static int MoonTest()
        {
            var time = new AstroTime(2019, 6, 24, 15, 45, 37);
            AstroVector vec = Astronomy.GeoVector(Body.Moon, time, Aberration.None);
            Console.WriteLine("MoonTest: {0} {1} {2}", vec.x.ToString("f17"), vec.y.ToString("f17"), vec.z.ToString("f17"));

            double dx = vec.x - (+0.002674036155459549);
            double dy = vec.y - (-0.0001531716308218381);
            double dz = vec.z - (-0.0003150201604895409);
            double diff = Math.Sqrt(dx*dx + dy*dy + dz*dz);
            Console.WriteLine("MoonTest: diff = {0}", diff.ToString("g5"));
            if (diff > 4.34e-19)
            {
                Console.WriteLine("MoonTest: EXCESSIVE ERROR");
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
                    Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto
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
                        outfile.WriteLine("v {0} {1} {2} {3} {4}", body, pos.t.tt.ToString("G17"), pos.x.ToString("G17"), pos.y.ToString("G17"), pos.z.ToString("G17"));
                        if (body != Body.Earth)
                        {
                            j2000 = Astronomy.Equator(body, time, observer, EquatorEpoch.J2000, Aberration.None);
                            ofdate = Astronomy.Equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
                            hor = Astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None);
                            outfile.WriteLine("s {0} {1} {2} {3} {4} {5} {6} {7}",
                                body,
                                time.tt.ToString("G17"), time.ut.ToString("G17"),
                                j2000.ra.ToString("G17"), j2000.dec.ToString("G17"), j2000.dist.ToString("G17"),
                                hor.azimuth.ToString("G17"), hor.altitude.ToString("G17"));
                        }
                    }

                    pos = Astronomy.GeoVector(Body.Moon, time, Aberration.None);
                    outfile.WriteLine("v GM {0} {1} {2} {3}", pos.t.tt.ToString("G17"), pos.x.ToString("G17"), pos.y.ToString("G17"), pos.z.ToString("G17"));
                    j2000 = Astronomy.Equator(Body.Moon, time, observer, EquatorEpoch.J2000, Aberration.None);
                    ofdate = Astronomy.Equator(Body.Moon, time, observer, EquatorEpoch.OfDate, Aberration.Corrected);
                    hor = Astronomy.Horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None);
                    outfile.WriteLine("s GM {0} {1} {2} {3} {4} {5} {6}",
                        time.tt.ToString("G17"), time.ut.ToString("G17"),
                        j2000.ra.ToString("G17"), j2000.dec.ToString("G17"), j2000.dist.ToString("G17"),
                        hor.azimuth.ToString("G17"), hor.altitude.ToString("G17"));

                    time = time.AddDays(10.0 + Math.PI/100.0);
                }
            }
            Console.WriteLine("AstroCheck: finished");
            return 0;
        }

        static int SeasonsTest(string filename)
        {
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
                    /*
                        2019-01-03T05:20Z Perihelion
                        2019-03-20T21:58Z Equinox
                        2019-06-21T15:54Z Solstice
                        2019-07-04T22:11Z Aphelion
                        2019-09-23T07:50Z Equinox
                        2019-12-22T04:19Z Solstice
                    */
                    Match m = re.Match(line);
                    if (!m.Success)
                    {
                        Console.WriteLine("SeasonsTest: ERROR {0} line {1}: cannot parse", filename, lnum);
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
                                Console.WriteLine("SeasonsTest: {0} line {1}: Invalid equinox date in test data.", filename, lnum);
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
                                Console.WriteLine("SeasonsTest: {0} line {1}: Invalid solstice date in test data.", filename, lnum);
                                return 1;
                        }
                    }
                    else if (name == "Aphelion")
                    {
                        /* not yet calculated */
                        continue;
                    }
                    else if (name == "Perihelion")
                    {
                        /* not yet calculated */
                        continue;
                    }
                    else
                    {
                        Console.WriteLine("SeasonsTest: {0} line {1}: unknown event type {2}", filename, lnum, name);
                        return 1;
                    }

                    /* Verify that the calculated time matches the correct time for this event. */
                    double diff_minutes = (24.0 * 60.0) * Math.Abs(calc_time.tt - correct_time.tt);
                    if (diff_minutes > max_minutes)
                        max_minutes = diff_minutes;

                    if (diff_minutes > 1.7)
                    {
                        Console.WriteLine("SeasonsTest: %s line %d: excessive error (%s): %lf minutes.\n", filename, lnum, name, diff_minutes);
                        return 1;
                    }
                }
                Console.WriteLine("SeasonsTest: verified {0} lines from file {1} : max error minutes = {2:0.000}", lnum, filename, max_minutes);
                Console.WriteLine("SeasonsTest: Event counts: mar={0}, jun={1}, sep={2}, dec={3}", mar_count, jun_count, sep_count, dec_count);
                return 0;
            }
        }

        static int MoonPhaseTest(string filename)
        {
            using (StreamReader infile = File.OpenText(filename))
            {
                const double threshold_seconds = 120.0;
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
                    /*
                        0 1800-01-25T03:21:00.000Z
                        1 1800-02-01T20:40:00.000Z
                        2 1800-02-09T17:26:00.000Z
                        3 1800-02-16T15:49:00.000Z
                    */
                    Match m = re.Match(line);
                    if (!m.Success)
                    {
                        Console.WriteLine("MoonPhaseTest: ERROR {0} line {1}: cannot parse", filename, lnum);
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
                    double degree_error = Math.Abs(calc_elong - expected_elong);
                    if (degree_error > 180.0)
                        degree_error = 360.0 - degree_error;
                    double arcmin = 60.0 * degree_error;
                    if (arcmin > 1.0)
                    {
                        Console.WriteLine("MoonPhaseTest({0} line {1}): EXCESSIVE ANGULAR ERROR: {2} arcmin", filename, lnum, arcmin);
                        return 1;
                    }
                    if (arcmin > max_arcmin)
                        max_arcmin = arcmin;

                    if (year != prev_year)
                    {
                        prev_year = year;
                        /* The test data contains a single year's worth of data for every 10 years. */
                        /* Every time we see the year value change, it breaks continuity of the phases. */
                        /* Start the search over again. */
                        AstroTime start_time = new AstroTime(year, 1, 1, 0, 0, 0);
                        mq = Astronomy.SearchMoonQuarter(start_time);
                        expected_quarter = -1;  /* we have no idea what the quarter should be */
                    }
                    else
                    {
                        /* Yet another lunar quarter in the same year. */
                        expected_quarter = (1 + mq.quarter) % 4;
                        mq = Astronomy.NextMoonQuarter(mq);

                        /* Make sure we find the next expected quarter. */
                        if (expected_quarter != mq.quarter)
                        {
                            Console.WriteLine("MoonPhaseTest({0} line {1}): SearchMoonQuarter returned quarter {2}, but expected {3}", filename, lnum, mq.quarter, expected_quarter);
                            return 1;
                        }
                    }
                    ++quarter_count;
                    /* Make sure the time matches what we expect. */
                    double diff_seconds = Math.Abs(mq.time.tt - expected_time.tt) * (24.0 * 3600.0);
                    if (diff_seconds > threshold_seconds)
                    {
                        Console.WriteLine("MoonPhaseTest({0} line {1}): excessive time error {2:0.000} seconds", filename, lnum, diff_seconds);
                        return 1;
                    }

                    if (diff_seconds > maxdiff)
                        maxdiff = diff_seconds;
                }

                Console.WriteLine("MoonPhaseTest: passed {0} lines for file {1} : max_arcmin = {2:0.000000}, maxdiff = {3:0.000} seconds, {4} quarters",
                    lnum, filename, max_arcmin, maxdiff, quarter_count);

                return 0;
            }
        }

        static int RiseSetTest(string filename)
        {
            using (StreamReader infile = File.OpenText(filename))
            {
                int lnum = 0;
                string line;
                var re = new Regex(@"^([A-Za-z]+)\s+([\-\+]?\d+\.?\d*)\s+([\-\+]?\d+\.?\d*)\s+(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z\s+([rs])\s*$");
                Body current_body = Body.Invalid;
                Observer observer = null;
                AstroTime r_search_date = null, s_search_date = null;
                AstroTime r_evt = null, s_evt = null;     /* rise event, set event: search results */
                AstroTime a_evt = null, b_evt = null;     /* chronologically first and second events */
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
                        Console.WriteLine("RiseSetTest({0} line {1}): invalid input format", filename, lnum);
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

                    /* Every time we see a new geographic location or body, start a new iteration */
                    /* of finding all rise/set times for that UTC calendar year. */
                    if (observer == null || observer.latitude != latitude || observer.longitude != longitude || current_body != body)
                    {
                        current_body = body;
                        observer = new Observer(latitude, longitude, 0.0);
                        r_search_date = s_search_date = new AstroTime(year, 1, 1, 0, 0, 0);
                        b_evt = null;
                        Console.WriteLine("RiseSetTest: {0} lat={1} lon={2}", body, latitude, longitude);
                    }

                    if (b_evt != null)
                    {
                        a_evt = b_evt;
                        a_dir = b_dir;
                        b_evt = null;
                    }
                    else
                    {
                        r_evt = Astronomy.SearchRiseSet(body, observer, Direction.Rise, r_search_date, 366.0);
                        if (r_evt == null)
                        {
                            Console.WriteLine("RiseSetTest({0} line {1}): Did not find {2} rise event.", filename, lnum, body);
                            return 1;
                        }

                        s_evt = Astronomy.SearchRiseSet(body, observer, Direction.Set, s_search_date, 366.0);
                        if (s_evt == null)
                        {
                            Console.WriteLine("RiseSetTest({0} line {1}): Did not find {2} rise event.", filename, lnum, body);
                            return 1;
                        }

                        /* Expect the current event to match the earlier of the found dates. */
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

                        /* Nudge the event times forward a tiny amount. */
                        r_search_date = r_evt.AddDays(nudge_days);
                        s_search_date = s_evt.AddDays(nudge_days);
                    }

                    if (a_dir != direction)
                    {
                        Console.WriteLine("RiseSetTest({0} line {1}): expected dir={2} but found {3}", filename, lnum, a_dir, direction);
                        return 1;
                    }
                    double error_minutes = (24.0 * 60.0) * Math.Abs(a_evt.tt - correct_date.tt);
                    sum_minutes += error_minutes * error_minutes;
                    if (error_minutes > max_minutes)
                        max_minutes = error_minutes;

                    if (error_minutes > 0.56)
                    {
                        Console.WriteLine("RiseSetTest({0} line {1}): excessive prediction time error = {2} minutes.\n", filename, lnum, error_minutes);
                        return 1;
                    }
                }

                double rms_minutes = Math.Sqrt(sum_minutes / lnum);
                Console.WriteLine("RiseSetTest: passed {0} lines: time errors in minutes: rms={1}, max={2}", lnum, rms_minutes, max_minutes);
                return 0;
            }
        }
    }
}
