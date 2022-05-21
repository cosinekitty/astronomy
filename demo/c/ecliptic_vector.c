/*
    ecliptic_vector.c  -  by Don Cross - 2022-05-20

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    An example of converting equatorial coordinates to ecliptic coordinates.
*/

#include <stdio.h>
#include "astro_demo_common.h"

int PrintUsage(void)
{
    fprintf(stderr,
        "\n"
        "USAGE: ecliptic_vector [yyyy-mm-ddThh:mm:ssZ]\n"
        "\n"
        "Displays cartesian ecliptic vectors for the Sun, Moon, and planets\n"
        "at the time specified on the command line (if present), or at the\n"
        "current time (if absent).\n"
        "\n"
    );
    return 1;
}

int main(int argc, const char *argv[])
{
    static const astro_body_t body[] = {
        BODY_SUN, BODY_MERCURY, BODY_VENUS, BODY_EARTH, BODY_MOON, BODY_MARS,
        BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE, BODY_PLUTO
    };

    int i;
    astro_time_t time;
    int num_bodies = sizeof(body) / sizeof(body[0]);
    astro_rotation_t rot;

    if (argc == 1)
    {
        /* Use the computer's current time. */
        time = Astronomy_CurrentTime();
    }
    else if (argc == 2)
    {
        /* Read the time from the command line. */
        if (0 != ParseTime(argv[1], &time))
            return 1;
    }
    else
    {
        return PrintUsage();
    }

    printf("UTC date = ");
    PrintTime(time);
    printf("\n");

    /* Get a rotation matrix that converts equatorial J2000 (EQJ) vectors to ecliptic vectors (ECL). */
    rot = Astronomy_Rotation_EQJ_ECL();

    printf("BODY               X           Y           Z\n");
    for (i=0; i < num_bodies; ++i)
    {
        astro_vector_t eqj = Astronomy_HelioVector(body[i], time);
        astro_vector_t ecl = Astronomy_RotateVector(rot, eqj);
        const char *name = Astronomy_BodyName(body[i]);
        if (ecl.status != ASTRO_SUCCESS)
        {
            fprintf(stderr, "ERROR %d calculating vector for %s.\n", (int)ecl.status, name);
            return 1;
        }
        printf("%-8s %11.6lf %11.6lf %11.6lf\n", name, ecl.x, ecl.y, ecl.z);
    }

    return 0;
}
