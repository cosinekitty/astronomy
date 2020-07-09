/*
    positions.c  -  by Don Cross - 2019-06-11

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Given an observer's geographic latitude and longitude,
    and an optional date and time, this program displays the
    equatorial and horizontal coordinates of the Sun, Moon, and planets.
    If the date and time is omitted from the command line, the
    program uses the computer's current date and time.
*/

#include <stdio.h>
#include "astro_demo_common.h"

int main(int argc, const char *argv[])
{
    static const astro_body_t body[] = {
        BODY_SUN, BODY_MOON, BODY_MERCURY, BODY_VENUS, BODY_MARS,
        BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE, BODY_PLUTO
    };

    int error;
    astro_observer_t observer;
    astro_time_t time;
    astro_equatorial_t equ_2000, equ_ofdate;
    astro_horizon_t hor;
    int i;
    int num_bodies = sizeof(body) / sizeof(body[0]);

    error = ParseArgs(argc, argv, &observer, &time);
    if (error)
        return error;

    printf("UTC date = ");
    PrintTime(time);
    printf("\n");

    printf("BODY           RA      DEC       AZ      ALT\n");
    for (i=0; i < num_bodies; ++i)
    {
        equ_2000 = Astronomy_Equator(body[i], &time, observer, EQUATOR_J2000, ABERRATION);
        if (equ_2000.status != ASTRO_SUCCESS)
        {
            fprintf(stderr, "ERROR: Astronomy_Equator returned status %d trying to get J2000 coordinates.\n", equ_2000.status);
            return 1;
        }

        equ_ofdate = Astronomy_Equator(body[i], &time, observer, EQUATOR_OF_DATE, ABERRATION);
        if (equ_ofdate.status != ASTRO_SUCCESS)
        {
            fprintf(stderr, "ERROR: Astronomy_Equator returned status %d trying to get coordinates of date.\n", equ_ofdate.status);
            return 1;
        }

        hor = Astronomy_Horizon(&time, observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NORMAL);
        printf("%-8s %8.2lf %8.2lf %8.2lf %8.2lf\n", Astronomy_BodyName(body[i]), equ_2000.ra, equ_2000.dec, hor.azimuth, hor.altitude);
    }

    return 0;
}
