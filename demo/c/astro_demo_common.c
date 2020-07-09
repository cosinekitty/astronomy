#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astro_demo_common.h"

int ParseArgs(int argc, const char *argv[], astro_observer_t *observer, astro_time_t *time)
{
    if (argc == 3 || argc == 4)
    {
        observer->height = 0.0;

        if (1 != sscanf(argv[1], "%lf", &observer->latitude) ||
            observer->latitude < -90.0 ||
            observer->latitude > +90.0)
        {
            fprintf(stderr, "ERROR: Invalid latitude '%s' on command line\n", argv[1]);
            return 1;
        }

        if (1 != sscanf(argv[2], "%lf", &observer->longitude) ||
            observer->longitude < -180.0 ||
            observer->longitude > +180.0)
        {
            fprintf(stderr, "ERROR: Invalid longitude '%s' on command line\n", argv[2]);
            return 1;
        }

        if (argc == 4)
        {
            /* Time is present on command line, so use it. */
            return ParseTime(argv[3], time);
        }

        /* Time is absent on command line, so use current time. */
        *time = Astronomy_CurrentTime();
        return 0;
    }

    fprintf(stderr, "USAGE: %s latitude longitude [yyyy-mm-ddThh:mm:ssZ]\n", argv[0]);
    return 1;
}

int ParseTime(const char *text, astro_time_t *time)
{
    astro_utc_t utc;
    int nscanned = sscanf(text, "%d-%d-%dT%d:%d:%lfZ",
        &utc.year, &utc.month, &utc.day,
        &utc.hour, &utc.minute, &utc.second);

    if (nscanned != 6)
    {
        fprintf(stderr, "ERROR: Invalid date/time format in '%s'\n", text);
        return 1;
    }

    *time = Astronomy_TimeFromUtc(utc);
    return 0;
}

void PrintTime(astro_time_t time)
{
    astro_status_t status;
    char text[TIME_TEXT_BYTES];

    status = Astronomy_FormatTime(time, TIME_FORMAT_SECOND, text, sizeof(text));
    if (status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "\nFATAL(PrintTime): status %d\n", status);
        exit(1);
    }
    printf("%s", text);
}

