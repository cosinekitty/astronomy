/*
    astro_demo_common.h  -  by Don Cross <cosinekitty@gmail.com>
    https://github.com/cosinekitty/astronomy

    Helper code for Astronomy Engine demo programs.
*/

// #include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "astro_demo_common.h"
#include "Arduino.h"
char printBufferCommon[80];

static const char ObserverVarName[] = "ASTRONOMY_ENGINE_OBSERVER";

static int IsValid(double x, double lo, double hi)
{
    return isfinite(x) && (x >= lo) && (x <= hi);
}

int ParseArgs(int argc, const char *argv[], astro_observer_t *observer, astro_time_t *time)
{

    int time_arg = 2;

    observer->height = 0.0;
    // fprintf(stderr, "USAGE: %s [0 first]\n", argv[0]);
    // fprintf(stderr, "USAGE: %s [2 doo]\n", argv[1]);
    //     fprintf(stderr, "USAGE: %s [3 se]\n", argv[2]);
    // fprintf(stderr, "%s='asdasd asdasd'\n", argv[2]);
    observer->latitude = atof(argv[0]);
    observer->longitude = atof(argv[1]);
    // time_arg = 2;

    if (!IsValid(observer->latitude, -90.0, +90.0) || !IsValid(observer->longitude, -180.0, +180.0))
    {
        fprintf(stderr, "ERROR: Invalid geographic coordinates in environment variable %s.\n", ObserverVarName);
        return 1;
    }

    time_arg = 2;
    return ParseTime(argv[time_arg], time);
    return 0;
}

int ParseTime(const char *text, astro_time_t *time)
{
    astro_utc_t utc;
    int nscanned = sscanf(text, "%d-%d-%dT%d:%d:%lfZ",
                          &utc.year, &utc.month, &utc.day,
                          &utc.hour, &utc.minute, &utc.second);

    if (nscanned != 6)
    {
        snprintf(printBufferCommon, sizeof printBufferCommon, "ERROR: Invalid date/time format in '%s'\n", text);
        Serial.println(printBufferCommon);
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
        snprintf(printBufferCommon, sizeof printBufferCommon, "\nFATAL(PrintTime): status %d\n", status);
        Serial.println(printBufferCommon);
        // exit(1);
    }
        snprintf(printBufferCommon, sizeof printBufferCommon, "%s", text);
        Serial.print(printBufferCommon);

}
