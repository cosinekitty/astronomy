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


#include "astro_demo_common.h"

char printBuffer[80];

void setup() {
  Serial.begin(9600);
}

void loop() {
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

    const char *argvs[] = {"22.2","11.8","1988-02-02T12:55:33Z"};
    int args = 3;
    error = ParseArgs(args, argvs, &observer, &time);
    if (error)
        Serial.print("Error.....");

    Serial.println("UTC date = ");
    PrintTime(time);
    Serial.println();

    Serial.println("BODY           RA      DEC       AZ      ALT\n");
    for (i=0; i < num_bodies; ++i)
    {
        equ_2000 = Astronomy_Equator(body[i], &time, observer, EQUATOR_J2000, ABERRATION);
        if (equ_2000.status != ASTRO_SUCCESS)
        {
            snprintf(printBuffer, sizeof printBuffer, "ERROR: Astronomy_Equator returned status %d trying to get J2000 coordinates.\n", equ_2000.status);
            Serial.println(printBuffer);

        }

        equ_ofdate = Astronomy_Equator(body[i], &time, observer, EQUATOR_OF_DATE, ABERRATION);
        if (equ_ofdate.status != ASTRO_SUCCESS)
        {
            snprintf(printBuffer, sizeof printBuffer, "ERROR: Astronomy_Equator returned status %d trying to get coordinates of date.\n", equ_ofdate.status);
            Serial.println(printBuffer);

        }

        hor = Astronomy_Horizon(&time, observer, equ_ofdate.ra, equ_ofdate.dec, REFRACTION_NORMAL);
        snprintf(printBuffer, sizeof printBuffer,"%-8s %8.2lf %8.2lf %8.2lf %8.2lf\n", Astronomy_BodyName(body[i]), equ_2000.ra, equ_2000.dec, hor.azimuth, hor.altitude);
        Serial.println(printBuffer);
    }
  delay(2000);

}
