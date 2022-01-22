/*
    horizon.c  -  Don Cross  -  2019-12-11

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This is a more advanced example. It shows how to use coordinate
    transforms and a binary search to find the two azimuths where the
    ecliptic intersects with an observer's horizon at a given date and time.
*/

#include <stdio.h>
#include <math.h>
#include "astro_demo_common.h"

static int FindEclipticCrossings(astro_observer_t observer, astro_time_t time);


int main(int argc, const char *argv[])
{
    int error;
    astro_observer_t observer;
    astro_time_t time;

    if (argc != 4)
    {
        fprintf(stderr, "USAGE: horizon latitude longitude yyyy-mm-ddThh:mm:ssZ\n");
        return 1;
    }

    error = ParseArgs(argc, argv, &observer, &time);
    if (error)
        return error;

    return FindEclipticCrossings(observer, time);
}

static int HorizontalCoords(
    astro_spherical_t *hor,
    double ecliptic_longitude,
    astro_time_t time,
    astro_rotation_t rot_ecl_hor)
{
    astro_vector_t ecl_vec;
    astro_vector_t hor_vec;
    astro_spherical_t eclip;

    eclip.status = ASTRO_SUCCESS;
    eclip.lat = 0.0;        /* being "on the ecliptic plane" means ecliptic latitude is zero. */
    eclip.lon = ecliptic_longitude;
    eclip.dist = 1.0;       /* any positive distance value will work fine. */

    /* Convert ecliptic angular coordinates to ecliptic vector. */
    ecl_vec = Astronomy_VectorFromSphere(eclip, time);
    if (ecl_vec.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Altitude: error %d converting ecliptic angles to vector.\n", ecl_vec.status);
        return 1;
    }

    /* Use the rotation matrix to convert ecliptic vector to horizontal vector. */
    hor_vec = Astronomy_RotateVector(rot_ecl_hor, ecl_vec);
    if (hor_vec.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Altitude: error %d rotating ecliptic vector to horizontal vector.\n", hor_vec.status);
        return 1;
    }

    /* Find horizontal angular coordinates, correcting for atmospheric refraction. */
    *hor = Astronomy_HorizonFromVector(hor_vec, REFRACTION_NORMAL);
    if (hor->status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Altitude: error %d converting horizontal vector to angles.\n", hor->status);
        return 1;
    }
    return 0;   /* success */
}

static int Search(
    double *ecliptic_longitude_crossing,
    astro_spherical_t *hor_crossing,
    astro_time_t time,
    astro_rotation_t rot_ecl_hor,
    double e1, double e2)
{
    int error;
    double e3;
    astro_spherical_t h3;
    const double tolerance = 1.0e-6;        /* one-millionth of a degree is close enough! */

    *ecliptic_longitude_crossing = NAN;
    hor_crossing->status = ASTRO_NOT_INITIALIZED;
    hor_crossing->lat = hor_crossing->lon = hor_crossing->dist = NAN;

    /*
        Binary search: find the ecliptic longitude such that the horizontal altitude
        ascends through a zero value. The caller must pass e1, e2 such that the altitudes
        bound zero in ascending order.
    */

    for(;;)
    {
        e3 = (e1 + e2) / 2.0;
        error = HorizontalCoords(&h3, e3, time, rot_ecl_hor);
        if (error)
            return error;

        if (fabs(e2-e1) < tolerance)
        {
            /* We have found the horizon crossing within tolerable limits. */
            *ecliptic_longitude_crossing = e3;
            *hor_crossing = h3;
            return 0;
        }

        if (h3.lat < 0.0)
            e1 = e3;
        else
            e2 = e3;
    }
}

#define NUM_SAMPLES 4
#define ECLIPLON(i)     ((360.0 * (i)) / NUM_SAMPLES)

static int FindEclipticCrossings(astro_observer_t observer, astro_time_t time)
{
    int i;
    astro_rotation_t rot;
    astro_spherical_t hor[NUM_SAMPLES];

    /*
        The ecliptic is a celestial circle that describes the mean plane of
        the Earth's orbit around the Sun. We use J2000 ecliptic coordinates,
        meaning the x-axis is defined to where the plane of the Earth's
        equator on January 1, 2000 at noon UTC intersects the ecliptic plane.
        The positive x-axis points toward the March equinox.
        Calculate a rotation matrix that converts J2000 ecliptic vectors
        to horizontal vectors for this observer and time.
    */
    rot = Astronomy_Rotation_ECL_HOR(&time, observer);
    if (rot.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "FindEclipticCrossings: error %d trying to find rotation matrix.\n", rot.status);
        return 1;
    }

    /*
        Sample several points around the ecliptic.
        Remember the horizontal coordinates for each sample.
    */
    for (i=0; i<NUM_SAMPLES; ++i)
        if (HorizontalCoords(&hor[i], ECLIPLON(i), time, rot))
            return 1;   /* failure */

    for (i=0; i < NUM_SAMPLES; ++i)
    {
        double a1 = hor[i].lat;
        double a2 = hor[(i+1) % NUM_SAMPLES].lat;
        double e1 = ECLIPLON(i);
        double e2 = ECLIPLON(i+1);
        double ex = NAN;
        astro_spherical_t hx;
        int error;

        if (a1*a2 <= 0.0)
        {
            const char *direction;

            /* Looks like a horizon crossing. Is altitude going up with longitude or down? */
            if (a2 > a1)
            {
                /* Search for the ecliptic longitude and azimuth where altitude ascends through zero. */
                error = Search(&ex, &hx, time, rot, e1, e2);
            }
            else
            {
                /* Search for the ecliptic longitude and azimuth where altitude descends through zero. */
                error = Search(&ex, &hx, time, rot, e2, e1);
            }

            if (error)
                return error;

            if (hx.lon > 0.0 && hx.lon < 180.0)
                direction = "ascends ";     /* azimuth is more toward the east than the west */
            else
                direction = "descends";     /* azimuth is more toward the west than the east */

            printf("Ecliptic longitude %9.4lf %s through horizon az %9.4lf, alt %12lg\n", ex, direction, hx.lon, hx.lat);
        }
    }
    return 0;
}
