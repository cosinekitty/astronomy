/*
    lunar_eclipse.c  -  by Don Cross - 2020-05-17

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Searches for the next 10 partial/total lunar eclipses after
    the current date, or a date specified on the command line.
*/

#include <stdio.h>
#include "astro_demo_common.h"


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
    printf(" - Partial eclipse begins.\n");

    if (eclipse.sd_total > 0.0)
    {
        PrintTime(Astronomy_AddDays(eclipse.peak, -eclipse.sd_total / MINUTES_PER_DAY));
        printf(" - Total eclipse begins.\n");
    }

    PrintTime(eclipse.peak);
    printf(" - Peak of %s eclipse.\n", (eclipse.kind == ECLIPSE_TOTAL) ? "total" : "partial");

    if (eclipse.sd_total > 0.0)
    {
        PrintTime(Astronomy_AddDays(eclipse.peak, +eclipse.sd_total / MINUTES_PER_DAY));
        printf(" - Total eclipse ends.\n");
    }

    PrintTime(Astronomy_AddDays(eclipse.peak, +eclipse.sd_partial / MINUTES_PER_DAY));
    printf(" - Partial eclipse ends.\n\n");
}


int main(int argc, const char *argv[])
{
    astro_time_t time;
    int count;
    astro_lunar_eclipse_t eclipse;

    switch (argc)
    {
    case 1:
        time = Astronomy_CurrentTime();
        break;

    case 2:
        if (ParseTime(argv[1], &time))
            return 1;
        break;

    default:
        fprintf(stderr, "USAGE: lunar_eclipse [date]\n");
        return 1;
    }

    eclipse = Astronomy_SearchLunarEclipse(time);
    count = 0;
    for(;;)
    {
        if (eclipse.status != ASTRO_SUCCESS)
        {
            fprintf(stderr, "Error %d: could not find next lunar eclipse.\n", eclipse.status);
            return 1;
        }
        if (eclipse.kind != ECLIPSE_PENUMBRAL)
        {
            PrintEclipse(eclipse);
            if (++count == 10)
                return 0;
        }
        eclipse = Astronomy_NextLunarEclipse(eclipse.peak);
    }
}
