using System;
using CosineKitty;
using demo_helper;

namespace horizon
{
    class Program
    {
        static int Main(string[] args)
        {
            Observer observer;
            AstroTime time;
            DemoHelper.ParseArgs("horizon", args, out observer, out time);
            return FindEclipticCrossings(observer, time);
        }

        static int HorizontalCoords(
            out Spherical hor,
            double ecliptic_longitude,
            AstroTime time,
            RotationMatrix rot_ecl_hor)
        {
            var eclip = new Spherical(
                0.0,         /* being "on the ecliptic plane" means ecliptic latitude is zero. */
                ecliptic_longitude,
                1.0         /* any positive distance value will work fine. */
            );

            /* Convert ecliptic angular coordinates to ecliptic vector. */
            AstroVector ecl_vec = Astronomy.VectorFromSphere(eclip, time);

            /* Use the rotation matrix to convert ecliptic vector to horizontal vector. */
            AstroVector hor_vec = Astronomy.RotateVector(rot_ecl_hor, ecl_vec);

            /* Find horizontal angular coordinates, correcting for atmospheric refraction. */
            hor = Astronomy.HorizonFromVector(hor_vec, Refraction.Normal);

            return 0;   /* success */
        }

        static int Search(
            out double ecliptic_longitude_crossing,
            out Spherical hor_crossing,
            AstroTime time,
            RotationMatrix rot_ecl_hor,
            double e1, double e2)
        {
            int error;
            double e3;
            Spherical h3;
            const double tolerance = 1.0e-6;        /* one-millionth of a degree is close enough! */

            /*
                Binary search: find the ecliptic longitude such that the horizontal altitude
                ascends through a zero value. The caller must pass e1, e2 such that the altitudes
                bound zero in ascending order.
            */

            ecliptic_longitude_crossing = 1.0e+99;      // initialize with impossible value
            hor_crossing = new Spherical();

            for(;;)
            {
                e3 = (e1 + e2) / 2.0;
                error = HorizontalCoords(out h3, e3, time, rot_ecl_hor);
                if (error != 0)
                    return error;

                if (Math.Abs(e2-e1) < tolerance)
                {
                    /* We have found the horizon crossing within tolerable limits. */
                    ecliptic_longitude_crossing = e3;
                    hor_crossing = h3;
                    return 0;
                }

                if (h3.lat < 0.0)
                    e1 = e3;
                else
                    e2 = e3;
            }
        }

        const int NUM_SAMPLES = 4;
        static double ECLIPLON(int i) => ((360.0 * i) / NUM_SAMPLES);

        static int FindEclipticCrossings(Observer observer, AstroTime time)
        {
            int i;
            var hor = new Spherical[NUM_SAMPLES];

            /*
                The ecliptic is a celestial circle that describes the mean plane of
                the Earth's orbit around the Sun. We use J2000 ecliptic coordinates,
                meaning the x-axis is defined to where the plane of the Earth's
                equator on January 1, 2000 at noon UTC intersects the ecliptic plane.
                The positive x-axis points toward the March equinox.
                Calculate a rotation matrix that converts J2000 ecliptic vectors
                to horizontal vectors for this observer and time.
            */
            RotationMatrix rot = Astronomy.Rotation_ECL_HOR(time, observer);

            /*
                Sample several points around the ecliptic.
                Remember the horizontal coordinates for each sample.
            */
            for (i=0; i<NUM_SAMPLES; ++i)
                if (0 != HorizontalCoords(out hor[i], ECLIPLON(i), time, rot))
                    return 1;   /* failure */

            for (i=0; i < NUM_SAMPLES; ++i)
            {
                double a1 = hor[i].lat;
                double a2 = hor[(i+1) % NUM_SAMPLES].lat;
                double e1 = ECLIPLON(i);
                double e2 = ECLIPLON(i+1);
                double ex;
                Spherical hx;
                int error;

                if (a1*a2 <= 0.0)
                {
                    /* Looks like a horizon crossing. Is altitude going up with longitude or down? */
                    if (a2 > a1)
                    {
                        /* Search for the ecliptic longitude and azimuth where altitude ascends through zero. */
                        error = Search(out ex, out hx, time, rot, e1, e2);
                    }
                    else
                    {
                        /* Search for the ecliptic longitude and azimuth where altitude descends through zero. */
                        error = Search(out ex, out hx, time, rot, e2, e1);
                    }

                    if (error != 0)
                        return error;

                    string direction;
                    if (hx.lon > 0.0 && hx.lon < 180.0)
                        direction = "ascends ";     /* azimuth is more toward the east than the west */
                    else
                        direction = "descends";     /* azimuth is more toward the west than the east */

                    Console.WriteLine("Ecliptic longitude {0,9:0.0000} {1} through horizon az {2,9:0.0000}, alt {3,12:0.000e+00}", ex, direction, hx.lon, hx.lat);
                }
            }
            return 0;
        }
    }
}
