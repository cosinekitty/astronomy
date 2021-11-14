/*
    astro_demo_common.h  -  by Don Cross <cosinekitty@gmail.com>
    https://github.com/cosinekitty/astronomy

    Helper code for Astronomy Engine demo programs.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "astro_demo_common.h"

static const char ObserverVarName[] = "ASTRONOMY_ENGINE_OBSERVER";

static int IsValid(double x, double lo, double hi)
{
    return isfinite(x) && (x >= lo) && (x <= hi);
}


int ParseArgs(int argc, const char *argv[], astro_observer_t *observer, astro_time_t *time)
{
    int time_arg = 0;

    observer->height = 0.0;

    if ((argc==2 || argc==3) && !strcmp(argv[1], "-e"))
    {
        /* Use an environment variable to specify the observer's location. */
        const char *obstext = getenv(ObserverVarName);
        if (obstext == NULL)
        {
            fprintf(stderr, "ERROR: The -e option was specified but environment variable %s is not defined.\n", ObserverVarName);
            return 1;
        }

        if (2 != sscanf(obstext, "%lf %lf", &observer->latitude, &observer->longitude) ||
            !IsValid(observer->latitude, -90.0, +90.0) ||
            !IsValid(observer->longitude, -180.0, +180.0))
        {
            fprintf(stderr, "ERROR: Invalid geographic coordinates in environment variable %s.\n", ObserverVarName);
            return 1;
        }

        time_arg = 2;
    }
    else if (argc==3 || argc==4)
    {
        /* The observer's location must appear in the command line. */

        if (1 != sscanf(argv[1], "%lf", &observer->latitude) ||
            !IsValid(observer->latitude, -90.0, +90.0))
        {
            fprintf(stderr, "ERROR: Invalid latitude '%s' on command line\n", argv[1]);
            return 1;
        }

        if (1 != sscanf(argv[2], "%lf", &observer->longitude) ||
            !IsValid(observer->longitude, -180.0, +180.0))
        {
            fprintf(stderr, "ERROR: Invalid longitude '%s' on command line\n", argv[2]);
            return 1;
        }

        time_arg = 3;
    }

    if (time_arg && argc == 1+time_arg)
    {
        /* Time is present on command line, so use it. */
        return ParseTime(argv[time_arg], time);
    }

    if (time_arg && argc == time_arg)
    {
        /* Time is absent on command line, so use current time. */
        *time = Astronomy_CurrentTime();
        return 0;
    }

    fprintf(stderr, "USAGE: %s [-e | latitude longitude] [yyyy-mm-ddThh:mm:ssZ]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "If '-e' is specified, an environment variable must be set as follows:\n");
    fprintf(stderr, "%s='latitude longitude'\n", ObserverVarName);
    fprintf(stderr, "Otherwise, the latitude and longitude must appear in the command line.\n");
    fprintf(stderr, "\n");
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

