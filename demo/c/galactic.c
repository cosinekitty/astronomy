/*
    galactic.c  -  Don Cross  -  2021-06-10

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program illustrates how to convert a location
    in the sky expressed in IAU 1958 galactic coordinates
    into the local altitude and azimuth of someone wanting
    to aim a radio dish at it.
*/

#include <stdio.h>
#include <math.h>
#include "astro_demo_common.h"


int GalaticToHorizontal(
    astro_time_t time,
    astro_observer_t observer,
    double glat,
    double glon,
    double *altitude,
    double *azimuth)
{
    astro_rotation_t    rot, adjust_rot;
    astro_spherical_t   gsphere, hsphere;
    astro_vector_t      gvec, hvec;

    /*
        Calculate a rotation matrix that converts
        galactic coordinates to J2000 equatorial coordinates.
    */
    rot = Astronomy_Rotation_GAL_EQJ();

    /*
        Adjust the rotation matrix to convert galatic to horizontal (HOR).
    */
    adjust_rot = Astronomy_Rotation_EQJ_HOR(&time, observer);
    rot = Astronomy_CombineRotation(rot, adjust_rot);

    /*
        Convert the galactic coordinates from angles to a unit vector.
    */
    gsphere.status = ASTRO_SUCCESS;
    gsphere.lat = glat;
    gsphere.lon = glon;
    gsphere.dist = 1.0;
    gvec = Astronomy_VectorFromSphere(gsphere, time);
    if (gvec.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Astronomy_VectorFromSphere returned error %d\n", gvec.status);
        return 1;
    }

    /*
        Use the rotation matrix to convert the galactic vector to a horizontal vector.
    */
    hvec = Astronomy_RotateVector(rot, gvec);
    if (hvec.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Astronomy_RotateVector returned error %d\n", hvec.status);
        return 1;
    }

    /*
        Convert the horizontal vector back to angular coordinates: altitude and azimuth.
        Assuming this is a radio source (not optical), do not correct for refraction.
    */
    hsphere = Astronomy_HorizonFromVector(hvec, REFRACTION_NONE);
    if (hsphere.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Astronomy_HorizonFromVector returned error %d\n", hsphere.status);
        return 1;
    }

    *altitude = hsphere.lat;
    *azimuth  = hsphere.lon;
    return 0;
}


int main(int argc, const char *argv[])
{
    astro_observer_t observer;
    astro_time_t time;
    double glat, glon;
    double azimuth, altitude;

    if (argc < 5 || argc > 6)
    {
        fprintf(stderr,
            "\n"
            "USAGE: galactic olat olon glat glon [yyyy-mm-ddThh:mm:ssZ]\n"
            "\n"
            "where\n"
            "\n"
            "    olat = observer's latitude on the Earth\n"
            "    olon = observer's longitude on the Earth\n"
            "    glat = IAU 1958 galatic latitude of the target\n"
            "    glon = IAU 1958 galatic longitude of the target\n"
            "    yyyy-mm-ddThh:mm:ssZ = optional UTC date/time\n"
            "\n"
            "Given the galactic coordinates of a point source in the sky,\n"
            "this program calculates horizontal aiming coordinates for an\n"
            "observer on or near the Earth's surface.\n"
            "\n"
            "If the date/time is given on the command line, it is used.\n"
            "Otherwise, the computer's current date/time is used.\n"
            "\n"
        );
        return 1;
    }

    observer.height = 0.0;

    if (1 != sscanf(argv[1], "%lf", &observer.latitude) ||
        observer.latitude < -90.0 ||
        observer.latitude > +90.0)
    {
        fprintf(stderr, "ERROR: Invalid observer latitude '%s' on command line\n", argv[1]);
        return 1;
    }

    if (1 != sscanf(argv[2], "%lf", &observer.longitude) ||
        observer.longitude < -180.0 ||
        observer.longitude > +180.0)
    {
        fprintf(stderr, "ERROR: Invalid observer longitude '%s' on command line\n", argv[2]);
        return 1;
    }

    if (1 != sscanf(argv[3], "%lf", &glat) || glat < -90.0 || glat > +90.0)
    {
        fprintf(stderr, "ERROR: Invalid galatic latitude '%s' on command line\n", argv[3]);
        return 1;
    }

    if (1 != sscanf(argv[4], "%lf", &glon) || glon <= -360.0 || glon >= +360.0)
    {
        fprintf(stderr, "ERROR: Invalid galatic longitude '%s' on command line\n", argv[3]);
        return 1;
    }

    if (argc > 5)
    {
        /* Time is present on command line, so use it. */
        if (ParseTime(argv[5], &time))
            return 1;
    }
    else
    {
        /* Time is absent on command line, so use current time. */
        time = Astronomy_CurrentTime();
    }

    if (GalaticToHorizontal(time, observer, glat, glon, &altitude, &azimuth))
        return 1;

    printf("altitude = %10.3lf, azimuth = %10.3lf\n", altitude, azimuth);
    return 0;
}
