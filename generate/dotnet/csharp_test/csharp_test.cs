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
                if (AstroCheck() != 0) return 1;
                if (SeasonsTest("../../seasons/seasons.txt") != 0) return 1;
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
    }
}
