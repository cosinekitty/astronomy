/*
    moonphase.c  -  by Don Cross - 2019-05-25

    Example C program for Astronomy Engine:
    https://cosinekitty.github.io/astronomy/
*/

#include <stdio.h>
#include "astronomy.h"

static const char *QuarterName(int quarter)
{
    switch (quarter)
    {
    case 0:   return "new moon";
    case 1:   return "first quarter";
    case 2:   return "full moon";
    case 3:   return "third quarter";
    default:  return "INVALID QUARTER";
    }
}

static void PrintTime(astro_time_t time)
{
    astro_utc_t utc = Astronomy_UtcFromTime(time);
    printf("%04d-%02d-%02d %02d:%02d:%02.0lf UTC", utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second);
}

int main()
{
    astro_time_t now;
    astro_angle_result_t phase;
    astro_moon_quarter_t mq;
    int i;

    /* 
        Calculate the Moon's current phase angle,
        which ranges from 0 to 360 degrees.

          0 = new moon,
         90 = first quarter,
        180 = full moon,
        270 = third quarter.
    */

    now = Astronomy_CurrentTime();
    printf("Current date/time  = ");
    PrintTime(now);
    printf("\n");

    phase = Astronomy_MoonPhase(now);
    if (phase.status != ASTRO_SUCCESS)
    {
        printf("Astronomy_MoonPhase error %d\n", phase.status);
        return 1;
    }

    printf("Moon's phase angle = %0.6lf degrees.\n", phase.angle);

    /* Find the next 10 lunar quarter phases. */
    printf("\nThe next 10 lunar quarters are:\n");
    for (i=0; i < 10; ++i)
    {
        if (i == 0)
            mq = Astronomy_SearchMoonQuarter(now);
        else
            mq = Astronomy_NextMoonQuarter(mq);

        if (mq.status != ASTRO_SUCCESS)
        {
            printf("Error %d trying to find moon quarter.\n", mq.status);
            return 1;
        }

        PrintTime(mq.time);
        printf(" : %s\n", QuarterName(mq.quarter));
    }

    return 0;
}
