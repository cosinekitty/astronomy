/*
    Astronomy library for C/C++.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astronomy.h"

#ifdef __cplusplus
extern "C" {
#endif

#define T0                  2451545.0
#define MJD_BASIS           2400000.5
#define Y2000_IN_MJD        (T0 - MJD_BASIS)

static void FatalError(const char *message)
{
    fprintf(stderr, "FATAL: %s\n", message);
    exit(1);
}

typedef struct
{
    double mjd;
    double dt;
}
deltat_entry_t;


static const deltat_entry_t DT[] = $ASTRO_DELTA_T();

#define DT_LENGTH     (sizeof(DT) / sizeof(DT[0]))

static double DeltaT(double mjd)
{
    int lo, hi, c;
    double frac;

    if (mjd <= DT[0].mjd)
        return DT[0].dt;

    if (mjd >= DT[DT_LENGTH-1].mjd)
        return DT[DT_LENGTH-1].dt;

    /* Do a binary search to find the pair of indexes this mjd lies between. */

    lo = 0;
    hi = DT_LENGTH-2;   /* make sure there is always an array element after the one we are looking at. */
    for(;;)
    {
        if (lo > hi)
        {
            /* This should never happen unless there is a bug in the binary search. */
            FatalError("DeltaT: could not find delta-t value");
        }

        c = (lo + hi) / 2;
        if (mjd < DT[c].mjd)
            hi = c-1;
        else if (mjd > DT[c+1].mjd)
            lo = c+1;
        else
        {
            frac = (mjd - DT[c].mjd) / (DT[c+1].mjd - DT[c].mjd);
            return DT[c].dt + frac*(DT[c+1].dt - DT[c].dt);
        }
    }
}

static double TerrestrialTime(double ut)
{
    return ut + DeltaT(ut + Y2000_IN_MJD)/86400.0;
}

astro_time_t Astronomy_MakeTime(int year, int month, int day, int hour, int minute, double second)
{
    astro_time_t time;
    long int jd12h;
    long int y2000;

    /* This formula is adapted from NOVAS C 3.1 function julian_date() */
    jd12h = (long) day - 32075L + 1461L * ((long) year + 4800L
        + ((long) month - 14L) / 12L) / 4L
        + 367L * ((long) month - 2L - ((long) month - 14L) / 12L * 12L)
        / 12L - 3L * (((long) year + 4900L + ((long) month - 14L) / 12L)
        / 100L) / 4L;    

    y2000 = jd12h - 2451545L;

    time.ut = (double)y2000 - 0.5 + (hour / 24.0) + (minute / (24.0 * 60.0)) + (second / (24.0 * 3600.0));

    time.tt = TerrestrialTime(time.ut);

    return time;
}

astro_time_t Astronomy_AddDays(astro_time_t time, double days)
{
    /* 
        This is slightly wrong, but the error is tiny.
        We really should be adding to TT, not to UT.
        But using TT would require creating an inverse function for DeltaT,
        which would be quite a bit of extra calculation.
        I estimate the error is in practice on the order of 10^(-7)
        times the value of 'days'.
        This is based on a typical drift of 1 second per year between UT and TT.
    */

    astro_time_t sum;

    sum.ut = time.ut + days;
    sum.tt = TerrestrialTime(sum.ut);

    return sum;   
}

#ifdef __cplusplus
}
#endif
