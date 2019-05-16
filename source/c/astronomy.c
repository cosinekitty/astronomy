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


static const deltat_entry_t DT[] = {
{ -72638.0, 38 },
{ -65333.0, 26 },
{ -58028.0, 21 },
{ -50724.0, 21.1 },
{ -43419.0, 13.5 },
{ -39766.0, 13.7 },
{ -36114.0, 14.8 },
{ -32461.0, 15.7 },
{ -28809.0, 15.6 },
{ -25156.0, 13.3 },
{ -21504.0, 12.6 },
{ -17852.0, 11.2 },
{ -14200.0, 11.13 },
{ -10547.0, 7.95 },
{ -6895.0, 6.22 },
{ -3242.0, 6.55 },
{ -1416.0, 7.26 },
{ 410.0, 7.35 },
{ 2237.0, 5.92 },
{ 4063.0, 1.04 },
{ 5889.0, -3.19 },
{ 7715.0, -5.36 },
{ 9542.0, -5.74 },
{ 11368.0, -5.86 },
{ 13194.0, -6.41 },
{ 15020.0, -2.70 },
{ 16846.0, 3.92 },
{ 18672.0, 10.38 },
{ 20498.0, 17.19 },
{ 22324.0, 21.41 },
{ 24151.0, 23.63 },
{ 25977.0, 24.02 },
{ 27803.0, 23.91 },
{ 29629.0, 24.35 },
{ 31456.0, 26.76 },
{ 33282.0, 29.15 },
{ 35108.0, 31.07 },
{ 36934.0, 33.150 },
{ 38761.0, 35.738 },
{ 40587.0, 40.182 },
{ 42413.0, 45.477 },
{ 44239.0, 50.540 },
{ 44605.0, 51.3808 },
{ 44970.0, 52.1668 },
{ 45335.0, 52.9565 },
{ 45700.0, 53.7882 },
{ 46066.0, 54.3427 },
{ 46431.0, 54.8712 },
{ 46796.0, 55.3222 },
{ 47161.0, 55.8197 },
{ 47527.0, 56.3000 },
{ 47892.0, 56.8553 },
{ 48257.0, 57.5653 },
{ 48622.0, 58.3092 },
{ 48988.0, 59.1218 },
{ 49353.0, 59.9845 },
{ 49718.0, 60.7853 },
{ 50083.0, 61.6287 },
{ 50449.0, 62.2950 },
{ 50814.0, 62.9659 },
{ 51179.0, 63.4673 },
{ 51544.0, 63.8285 },
{ 51910.0, 64.0908 },
{ 52275.0, 64.2998 },
{ 52640.0, 64.4734 },
{ 53005.0, 64.5736 },
{ 53371.0, 64.6876 },
{ 53736.0, 64.8452 },
{ 54101.0, 65.1464 },
{ 54466.0, 65.4573 },
{ 54832.0, 65.7768 },
{ 55197.0, 66.0699 },
{ 55562.0, 66.3246 },
{ 55927.0, 66.6030 },
{ 56293.0, 66.9069 },
{ 56658.0, 67.2810 },
{ 57023.0, 67.6439 },
{ 57388.0, 68.1024 },
{ 57754.0, 68.5927 },
{ 58119.0, 68.9676 },
{ 58484.0, 69.2201 },
{ 58849.0, 69.87 },
{ 59214.0, 70.39 },
{ 59580.0, 70.91 },
{ 59945.0, 71.40 },
{ 60310.0, 71.88 },
{ 60675.0, 72.36 },
{ 61041.0, 72.83 },
{ 61406.0, 73.32 },
{ 61680.0, 73.66 }
};

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
