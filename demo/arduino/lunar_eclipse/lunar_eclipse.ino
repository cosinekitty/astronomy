/*
    lunar_eclipse.c  -  by Don Cross - 2020-05-17

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Searches for the next 10 partial/total lunar eclipses after
    the current date, or a date specified on the command line.
*/

#include <stdio.h>
#include "astro_demo_common.h"
char printBuffer[80];

void PrintEclipse(astro_lunar_eclipse_t eclipse)
{
    /*
        Calculate beginning/ending of different phases
        of an eclipse by subtracting/adding the peak time
        with the number of minutes indicated by the "semi-duration"
        fields sd_partial and sd_total.
    */
    const double MINUTES_PER_DAY = 24 * 60;

    PrintTime(Astronomy_AddDays(eclipse.peak, -eclipse.sd_partial / MINUTES_PER_DAY));
    Serial.println(" - Partial eclipse begins.\n");

    if (eclipse.sd_total > 0.0)
    {
        PrintTime(Astronomy_AddDays(eclipse.peak, -eclipse.sd_total / MINUTES_PER_DAY));
        Serial.println(" - Total eclipse begins.");
    }

    PrintTime(eclipse.peak);
    snprintf(printBuffer, sizeof printBuffer," - Peak of %s eclipse.\n", (eclipse.kind == ECLIPSE_TOTAL) ? "total" : "partial");
    Serial.println(printBuffer);
    if (eclipse.sd_total > 0.0)
    {
        PrintTime(Astronomy_AddDays(eclipse.peak, +eclipse.sd_total / MINUTES_PER_DAY));
        Serial.println(" - Total eclipse ends.\n");
    }

    PrintTime(Astronomy_AddDays(eclipse.peak, +eclipse.sd_partial / MINUTES_PER_DAY));
    Serial.println(" - Partial eclipse ends.\n\n");
}


void setup() {
  Serial.begin(9600);
}

void loop() {
    astro_time_t time;
    int count;
    astro_lunar_eclipse_t eclipse;

    const char *argvs[] = {"0","1988-02-02T12:55:33Z"};

    ParseTime(argvs[1], &time);




    eclipse = Astronomy_SearchLunarEclipse(time);
    count = 0;
    for(;;)
    {
        if (eclipse.status != ASTRO_SUCCESS)
        {
            snprintf(printBuffer, sizeof printBuffer, "Error %d: could not find next lunar eclipse.\n", eclipse.status);
            Serial.println(printBuffer);
            delay(1000);
        }
        if (eclipse.kind != ECLIPSE_PENUMBRAL)
        {
            PrintEclipse(eclipse);
            if (++count == 10)
                delay(5000);
        }
        eclipse = Astronomy_NextLunarEclipse(eclipse.peak);
    }
    delay(1000);
}
