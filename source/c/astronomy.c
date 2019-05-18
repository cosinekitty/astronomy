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

#define ARRAYSIZE(x)    (sizeof(x) / sizeof(x[0]))

static const double T0        = 2451545.0;
static const double MJD_BASIS = 2400000.5;
#define Y2000_IN_MJD        (T0 - MJD_BASIS)
#define PI      3.14159265358979323846
static const double DEG2RAD = 0.017453292519943296;
static const double RAD2DEG = 57.295779513082321;
static const double ASEC360 = 1296000.0;
static const double ASEC2RAD = 4.848136811095359935899141e-6;
static const double PI2 = 2.0 * PI;
static const double ARC = 3600.0 * 180.0 / PI;          /* arcseconds per radian */
static const double C_AUDAY = 173.1446326846693;        /* speed of light in AU/day */
static const double ERAD = 6378136.6;                   /* mean earth radius in meters */
static const double AU = 1.4959787069098932e+11;        /* astronomical unit in meters */
static const double KM_PER_AU = 1.4959787069098932e+8;
static const double ANGVEL = 7.2921150e-5;

double Astronomy_VectorLength(astro_vector_t vector)
{
    return sqrt(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);
}

const char *Astronomy_BodyName(astro_body_t body)
{
    switch (body)
    {
    case BODY_MERCURY:  return "Mercury";
    case BODY_VENUS:    return "Venus";
    case BODY_EARTH:    return "Earth";
    case BODY_MARS:     return "Mars";
    case BODY_JUPITER:  return "Jupiter";
    case BODY_SATURN:   return "Saturn";
    case BODY_URANUS:   return "Uranus";
    case BODY_NEPTUNE:  return "Neptune";
    case BODY_PLUTO:    return "Pluto";
    case BODY_SUN:      return "Sun";
    case BODY_MOON:     return "Moon";
    default:            return "";
    }
}

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

astro_observer_t Astronomy_MakeObserver(double latitude, double longitude, double height)
{
    astro_observer_t observer;

    observer.latitude = latitude;
    observer.longitude = longitude;
    observer.height = height;

    return observer;
}

void iau2000b(astro_time_t time, double *dpsi, double *deps)
{
    /* Adapted from the NOVAS C 3.1 function of the same name. */

    static const short int nals_t[77][5] = 
    {
        { 0,    0,    0,    0,    1},
        { 0,    0,    2,   -2,    2},
        { 0,    0,    2,    0,    2},
        { 0,    0,    0,    0,    2},
        { 0,    1,    0,    0,    0},
        { 0,    1,    2,   -2,    2},
        { 1,    0,    0,    0,    0},
        { 0,    0,    2,    0,    1},
        { 1,    0,    2,    0,    2},
        { 0,   -1,    2,   -2,    2},
        { 0,    0,    2,   -2,    1},
        {-1,    0,    2,    0,    2},
        {-1,    0,    0,    2,    0},
        { 1,    0,    0,    0,    1},
        {-1,    0,    0,    0,    1},
        {-1,    0,    2,    2,    2},
        { 1,    0,    2,    0,    1},
        {-2,    0,    2,    0,    1},
        { 0,    0,    0,    2,    0},
        { 0,    0,    2,    2,    2},
        { 0,   -2,    2,   -2,    2},
        {-2,    0,    0,    2,    0},
        { 2,    0,    2,    0,    2},
        { 1,    0,    2,   -2,    2},
        {-1,    0,    2,    0,    1},
        { 2,    0,    0,    0,    0},
        { 0,    0,    2,    0,    0},
        { 0,    1,    0,    0,    1},
        {-1,    0,    0,    2,    1},
        { 0,    2,    2,   -2,    2},
        { 0,    0,   -2,    2,    0},
        { 1,    0,    0,   -2,    1},
        { 0,   -1,    0,    0,    1},
        {-1,    0,    2,    2,    1},
        { 0,    2,    0,    0,    0},
        { 1,    0,    2,    2,    2},
        {-2,    0,    2,    0,    0},
        { 0,    1,    2,    0,    2},
        { 0,    0,    2,    2,    1},
        { 0,   -1,    2,    0,    2},
        { 0,    0,    0,    2,    1},
        { 1,    0,    2,   -2,    1},
        { 2,    0,    2,   -2,    2},
        {-2,    0,    0,    2,    1},
        { 2,    0,    2,    0,    1},
        { 0,   -1,    2,   -2,    1},
        { 0,    0,    0,   -2,    1},
        {-1,   -1,    0,    2,    0},
        { 2,    0,    0,   -2,    1},
        { 1,    0,    0,    2,    0},
        { 0,    1,    2,   -2,    1},
        { 1,   -1,    0,    0,    0},
        {-2,    0,    2,    0,    2},
        { 3,    0,    2,    0,    2},
        { 0,   -1,    0,    2,    0},
        { 1,   -1,    2,    0,    2},
        { 0,    0,    0,    1,    0},
        {-1,   -1,    2,    2,    2},
        {-1,    0,    2,    0,    0},
        { 0,   -1,    2,    2,    2},
        {-2,    0,    0,    0,    1},
        { 1,    1,    2,    0,    2},
        { 2,    0,    0,    0,    1},
        {-1,    1,    0,    1,    0},
        { 1,    1,    0,    0,    0},
        { 1,    0,    2,    0,    0},
        {-1,    0,    2,   -2,    1},
        { 1,    0,    0,    0,    2},
        {-1,    0,    0,    1,    0},
        { 0,    0,    2,    1,    2},
        {-1,    0,    2,    4,    2},
        {-1,    1,    0,    1,    1},
        { 0,   -2,    2,   -2,    1},
        { 1,    0,    2,    2,    1},
        {-2,    0,    2,    2,    2},
        {-1,    0,    0,    0,    2},
        { 1,    1,    2,   -2,    2}
    };

   static const double cls_t[77][6] = 
   {
        {-172064161.0, -174666.0,  33386.0, 92052331.0,  9086.0, 15377.0},
        { -13170906.0,   -1675.0, -13696.0,  5730336.0, -3015.0, -4587.0},
        {  -2276413.0,    -234.0,   2796.0,   978459.0,  -485.0,  1374.0},
        {   2074554.0,     207.0,   -698.0,  -897492.0,   470.0,  -291.0},
        {   1475877.0,   -3633.0,  11817.0,    73871.0,  -184.0, -1924.0},
        {   -516821.0,    1226.0,   -524.0,   224386.0,  -677.0,  -174.0},
        {    711159.0,      73.0,   -872.0,    -6750.0,     0.0,   358.0},
        {   -387298.0,    -367.0,    380.0,   200728.0,    18.0,   318.0},
        {   -301461.0,     -36.0,    816.0,   129025.0,   -63.0,   367.0},
        {    215829.0,    -494.0,    111.0,   -95929.0,   299.0,   132.0},
        {    128227.0,     137.0,    181.0,   -68982.0,    -9.0,    39.0},
        {    123457.0,      11.0,     19.0,   -53311.0,    32.0,    -4.0},
        {    156994.0,      10.0,   -168.0,    -1235.0,     0.0,    82.0},
        {     63110.0,      63.0,     27.0,   -33228.0,     0.0,    -9.0},
        {    -57976.0,     -63.0,   -189.0,    31429.0,     0.0,   -75.0},
        {    -59641.0,     -11.0,    149.0,    25543.0,   -11.0,    66.0},
        {    -51613.0,     -42.0,    129.0,    26366.0,     0.0,    78.0},
        {     45893.0,      50.0,     31.0,   -24236.0,   -10.0,    20.0},
        {     63384.0,      11.0,   -150.0,    -1220.0,     0.0,    29.0},
        {    -38571.0,      -1.0,    158.0,    16452.0,   -11.0,    68.0},
        {     32481.0,       0.0,      0.0,   -13870.0,     0.0,     0.0},
        {    -47722.0,       0.0,    -18.0,      477.0,     0.0,   -25.0},
        {    -31046.0,      -1.0,    131.0,    13238.0,   -11.0,    59.0},
        {     28593.0,       0.0,     -1.0,   -12338.0,    10.0,    -3.0},
        {     20441.0,      21.0,     10.0,   -10758.0,     0.0,    -3.0},
        {     29243.0,       0.0,    -74.0,     -609.0,     0.0,    13.0},
        {     25887.0,       0.0,    -66.0,     -550.0,     0.0,    11.0},
        {    -14053.0,     -25.0,     79.0,     8551.0,    -2.0,   -45.0},
        {     15164.0,      10.0,     11.0,    -8001.0,     0.0,    -1.0},
        {    -15794.0,      72.0,    -16.0,     6850.0,   -42.0,    -5.0},
        {     21783.0,       0.0,     13.0,     -167.0,     0.0,    13.0},
        {    -12873.0,     -10.0,    -37.0,     6953.0,     0.0,   -14.0},
        {    -12654.0,      11.0,     63.0,     6415.0,     0.0,    26.0},
        {    -10204.0,       0.0,     25.0,     5222.0,     0.0,    15.0},
        {     16707.0,     -85.0,    -10.0,      168.0,    -1.0,    10.0},
        {     -7691.0,       0.0,     44.0,     3268.0,     0.0,    19.0},
        {    -11024.0,       0.0,    -14.0,      104.0,     0.0,     2.0},
        {      7566.0,     -21.0,    -11.0,    -3250.0,     0.0,    -5.0},
        {     -6637.0,     -11.0,     25.0,     3353.0,     0.0,    14.0},
        {     -7141.0,      21.0,      8.0,     3070.0,     0.0,     4.0},
        {     -6302.0,     -11.0,      2.0,     3272.0,     0.0,     4.0},
        {      5800.0,      10.0,      2.0,    -3045.0,     0.0,    -1.0},
        {      6443.0,       0.0,     -7.0,    -2768.0,     0.0,    -4.0},
        {     -5774.0,     -11.0,    -15.0,     3041.0,     0.0,    -5.0},
        {     -5350.0,       0.0,     21.0,     2695.0,     0.0,    12.0},
        {     -4752.0,     -11.0,     -3.0,     2719.0,     0.0,    -3.0},
        {     -4940.0,     -11.0,    -21.0,     2720.0,     0.0,    -9.0},
        {      7350.0,       0.0,     -8.0,      -51.0,     0.0,     4.0},
        {      4065.0,       0.0,      6.0,    -2206.0,     0.0,     1.0},
        {      6579.0,       0.0,    -24.0,     -199.0,     0.0,     2.0},
        {      3579.0,       0.0,      5.0,    -1900.0,     0.0,     1.0},
        {      4725.0,       0.0,     -6.0,      -41.0,     0.0,     3.0},
        {     -3075.0,       0.0,     -2.0,     1313.0,     0.0,    -1.0},
        {     -2904.0,       0.0,     15.0,     1233.0,     0.0,     7.0},
        {      4348.0,       0.0,    -10.0,      -81.0,     0.0,     2.0},
        {     -2878.0,       0.0,      8.0,     1232.0,     0.0,     4.0},
        {     -4230.0,       0.0,      5.0,      -20.0,     0.0,    -2.0},
        {     -2819.0,       0.0,      7.0,     1207.0,     0.0,     3.0},
        {     -4056.0,       0.0,      5.0,       40.0,     0.0,    -2.0},
        {     -2647.0,       0.0,     11.0,     1129.0,     0.0,     5.0},
        {     -2294.0,       0.0,    -10.0,     1266.0,     0.0,    -4.0},
        {      2481.0,       0.0,     -7.0,    -1062.0,     0.0,    -3.0},
        {      2179.0,       0.0,     -2.0,    -1129.0,     0.0,    -2.0},
        {      3276.0,       0.0,      1.0,       -9.0,     0.0,     0.0},
        {     -3389.0,       0.0,      5.0,       35.0,     0.0,    -2.0},
        {      3339.0,       0.0,    -13.0,     -107.0,     0.0,     1.0},
        {     -1987.0,       0.0,     -6.0,     1073.0,     0.0,    -2.0},
        {     -1981.0,       0.0,      0.0,      854.0,     0.0,     0.0},
        {      4026.0,       0.0,   -353.0,     -553.0,     0.0,  -139.0},
        {      1660.0,       0.0,     -5.0,     -710.0,     0.0,    -2.0},
        {     -1521.0,       0.0,      9.0,      647.0,     0.0,     4.0},
        {      1314.0,       0.0,      0.0,     -700.0,     0.0,     0.0},
        {     -1283.0,       0.0,      0.0,      672.0,     0.0,     0.0},
        {     -1331.0,       0.0,      8.0,      663.0,     0.0,     4.0},
        {      1383.0,       0.0,     -2.0,     -594.0,     0.0,    -2.0},
        {      1405.0,       0.0,      4.0,     -610.0,     0.0,     2.0},
        {      1290.0,       0.0,      0.0,     -556.0,     0.0,     0.0}
    };

    double t, el, elp, f, d, om, arg, dp, de, sarg, carg;
    int i;

    t = time.tt / 36525;
    el  = fmod(485868.249036 + t * 1717915923.2178, ASEC360) * ASEC2RAD;
    elp = fmod(1287104.79305 + t * 129596581.0481,  ASEC360) * ASEC2RAD;
    f   = fmod(335779.526232 + t * 1739527262.8478, ASEC360) * ASEC2RAD;
    d   = fmod(1072260.70369 + t * 1602961601.2090, ASEC360) * ASEC2RAD;
    om  = fmod(450160.398036 - t * 6962890.5431,    ASEC360) * ASEC2RAD;
    dp = 0;
    de = 0;
    for (i=76; i >= 0; --i) 
    {
        arg = fmod((nals_t[i][0]*el + nals_t[i][1]*elp + nals_t[i][2]*f + nals_t[i][3]*d + nals_t[i][4]*om), PI2);
        sarg = sin(arg);
        carg = cos(arg);
        dp += (cls_t[i][0] + cls_t[i][1] * t) * sarg + cls_t[i][2] * carg;
        de += (cls_t[i][3] + cls_t[i][4] * t) * carg + cls_t[i][5] * sarg;
    }

    *dpsi = -0.000135 + (dp * 1.0e-7);
    *deps = +0.000388 + (de * 1.0e-7);
}

static double mean_obliq(double tt)
{
    double t = tt / 36525.0;
    double asec = 
        (((( -  0.0000000434   * t
             -  0.000000576  ) * t
             +  0.00200340   ) * t
             -  0.0001831    ) * t
             - 46.836769     ) * t + 84381.406;

    return asec / 3600.0;
}

typedef struct 
{
    double tt;
    double dpsi;
    double deps;
    double ee;
    double mobl;
    double tobl;
}
earth_tilt_t;

static earth_tilt_t e_tilt(astro_time_t time)
{
    earth_tilt_t et;

    iau2000b(time, &et.dpsi, &et.deps);
    et.mobl = mean_obliq(time.tt);
    et.tobl = et.mobl + (et.deps / 3600.0);
    et.tt = time.tt;
    et.ee = et.dpsi * cos(et.mobl * DEG2RAD) / 15.0;

    return et;
}

static void ecl2equ_vec(astro_time_t time, const double ecl[3], double equ[3])
{
    double obl = mean_obliq(time.tt) * DEG2RAD;
    double cos_obl = cos(obl);
    double sin_obl = sin(obl);

    equ[0] = ecl[0];
    equ[1] = ecl[1]*cos_obl - ecl[2]*sin_obl;
    equ[2] = ecl[1]*sin_obl + ecl[2]*cos_obl;
}

static void precession(double tt1, const double pos1[3], double tt2, double pos2[3])
{
    double xx, yx, zx, xy, yy, zy, xz, yz, zz;
    double t, psia, omegaa, chia, sa, ca, sb, cb, sc, cc, sd, cd;
    double eps0 = 84381.406;

    if ((tt1 != 0.0) && (tt2 != 0.0))
        FatalError("precession: one of (tt1, tt2) must be zero.");
    
    t = (tt2 - tt1) / 36525;
    if (tt2 == 0)
        t = -t;

    psia   = (((((-    0.0000000951  * t
                 +    0.000132851 ) * t
                 -    0.00114045  ) * t
                 -    1.0790069   ) * t
                 + 5038.481507    ) * t);

    omegaa = (((((+    0.0000003337  * t
                 -    0.000000467 ) * t
                 -    0.00772503  ) * t
                 +    0.0512623   ) * t
                 -    0.025754    ) * t + eps0);

    chia   = (((((-    0.0000000560  * t
                 +    0.000170663 ) * t
                 -    0.00121197  ) * t
                 -    2.3814292   ) * t
                 +   10.556403    ) * t);

    eps0 = eps0 * ASEC2RAD;
    psia = psia * ASEC2RAD;
    omegaa = omegaa * ASEC2RAD;
    chia = chia * ASEC2RAD;

    sa = sin(eps0);
    ca = cos(eps0);
    sb = sin(-psia);
    cb = cos(-psia);
    sc = sin(-omegaa);
    cc = cos(-omegaa);
    sd = sin(chia);
    cd = cos(chia);

    xx =  cd * cb - sb * sd * cc;
    yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc;
    zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc;
    xy = -sd * cb - sb * cd * cc;
    yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc;
    zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc;
    xz =  sb * sc;
    yz = -sc * cb * ca - sa * cc;
    zz = -sc * cb * sa + cc * ca;

    if (tt2 == 0.0)
    { 
        /* Perform rotation from other epoch to J2000.0. */
        pos2[0] = xx * pos1[0] + xy * pos1[1] + xz * pos1[2];
        pos2[1] = yx * pos1[0] + yy * pos1[1] + yz * pos1[2];
        pos2[2] = zx * pos1[0] + zy * pos1[1] + zz * pos1[2];
    }
    else
    {
        /* Perform rotation from J2000.0 to other epoch. */
        pos2[0] = xx * pos1[0] + yx * pos1[1] + zx * pos1[2];
        pos2[1] = xy * pos1[0] + yy * pos1[1] + zy * pos1[2];
        pos2[2] = xz * pos1[0] + yz * pos1[1] + zz * pos1[2];
    }
}

static astro_equatorial_t vector2radec(const double pos[3])
{
    astro_equatorial_t equ;
    double xyproj;

    xyproj = pos[0]*pos[0] + pos[1]*pos[1];
    equ.dist = sqrt(xyproj + pos[2]*pos[2]);
    if (xyproj == 0.0)
    {
        if (pos[2] == 0.0)
        {
            /* Indeterminate coordinates; pos vector has zero length. */
            equ.ra = NAN;
            equ.dec = NAN;
        }
        else if (pos[2] < 0)
        {
            equ.ra = 0.0;
            equ.dec = -90.0;
        }
        else
        {
            equ.ra = 0.0;
            equ.dec = +90.0;
        }        
    }
    else
    {
        equ.ra = atan2(pos[1], pos[0]) / (DEG2RAD * 15.0);
        if (equ.ra < 0)
            equ.ra += 24.0;

        equ.dec = RAD2DEG * atan2(pos[2], sqrt(xyproj));
    }    

    return equ;
}

static void nutation(astro_time_t time, int direction, const double inpos[3], double outpos[3])
{
    earth_tilt_t tilt = e_tilt(time);
    double oblm = tilt.mobl * DEG2RAD;
    double oblt = tilt.tobl * DEG2RAD;
    double psi = tilt.dpsi * ASEC2RAD;
    double cobm = cos(oblm);
    double sobm = sin(oblm);
    double cobt = cos(oblt);
    double sobt = sin(oblt);
    double cpsi = cos(psi);
    double spsi = sin(psi);

    double xx = cpsi;
    double yx = -spsi * cobm;
    double zx = -spsi * sobm;
    double xy = spsi * cobt;
    double yy = cpsi * cobm * cobt + sobm * sobt;
    double zy = cpsi * sobm * cobt - cobm * sobt;
    double xz = spsi * sobt;
    double yz = cpsi * cobm * sobt - sobm * cobt;
    double zz = cpsi * sobm * sobt + cobm * cobt; 

    if (direction == 0) 
    {
        /* forward rotation */
        outpos[0] = xx * inpos[0] + yx * inpos[1] + zx * inpos[2];
        outpos[1] = xy * inpos[0] + yy * inpos[1] + zy * inpos[2];
        outpos[2] = xz * inpos[0] + yz * inpos[1] + zz * inpos[2];
    }
    else
    {
        /* inverse rotation */
        outpos[0] = xx * inpos[0] + xy * inpos[1] + xz * inpos[2];
        outpos[1] = yx * inpos[0] + yy * inpos[1] + yz * inpos[2];
        outpos[2] = zx * inpos[0] + zy * inpos[1] + zz * inpos[2];
    }
}

static double era(astro_time_t time)        /* Earth Rotation Angle */
{
    double thet1 = 0.7790572732640 + 0.00273781191135448 * time.ut;
    double thet3 = fmod(time.ut, 1.0);
    double theta = 360.0 * fmod(thet1 + thet3, 1.0);
    if (theta < 0.0)
        theta += 360.0;

    return theta;
}

static double sidereal_time(astro_time_t time)
{
    double t = time.tt / 36525.0;
    double eqeq = 15.0 * e_tilt(time).ee;    /* Replace with eqeq=0 to get GMST instead of GAST (if we ever need it) */
    double theta = era(time);
    double st = (eqeq + 0.014506 +
        (((( -    0.0000000368   * t
            -    0.000029956  ) * t
            -    0.00000044   ) * t
            +    1.3915817    ) * t
            + 4612.156534     ) * t);

    double gst = fmod(st/3600.0 + theta, 360.0) / 15.0;
    if (gst < 0.0) 
        gst += 24.0;

    return gst;
}

static void terra(astro_observer_t observer, double st, double pos[3], double vel[3]) 
{
    double erad_km = ERAD / 1000.0;
    double df = 1.0 - 0.003352819697896;    /* flattening of the Earth */
    double df2 = df * df;
    double phi = observer.latitude * DEG2RAD;
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double c = 1.0 / sqrt(cosphi*cosphi + df2*sinphi*sinphi);
    double s = df2 * c;
    double ht_km = observer.height / 1000.0;
    double ach = erad_km*c + ht_km;
    double ash = erad_km*s + ht_km;
    double stlocl = (15.0*st + observer.longitude) * DEG2RAD;
    double sinst = sin(stlocl);
    double cosst = cos(stlocl);

    pos[0] = ach * cosphi * cosst / KM_PER_AU;
    pos[1] = ach * cosphi * sinst / KM_PER_AU;
    pos[2] = ash * sinphi / KM_PER_AU;

    vel[0] = -ANGVEL * ach * cosphi * sinst * 86400.0;
    vel[1] = +ANGVEL * ach * cosphi * cosst * 86400.0;
    vel[2] = 0.0;
}

static void geo_pos(astro_time_t time, astro_observer_t observer, double outpos[3])
{
    double gast, vel[3], pos1[3], pos2[3];

    gast = sidereal_time(time);
    terra(observer, gast, pos1, vel);
    nutation(time, -1, pos1, pos2);
    precession(time.tt, pos2, 0.0, outpos);
}

static void spin(double angle, const double pos1[3], double vec2[3]) 
{
    double angr = angle * DEG2RAD;
    double cosang = cos(angr);
    double sinang = sin(angr);
    double xx = cosang;
    double yx = sinang;
    double zx = 0;
    double xy = -sinang;
    double yy = cosang;
    double zy = 0;
    double xz = 0;
    double yz = 0;
    double zz = 1;

    vec2[0] = xx*pos1[0] + yx*pos1[1] + zx*pos1[2];
    vec2[1] = xy*pos1[0] + yy*pos1[1] + zy*pos1[2];
    vec2[2] = xz*pos1[0] + yz*pos1[1] + zz*pos1[2];
}

static void ter2cel(astro_time_t time, const double vec1[3], double vec2[3])
{
    double gast = sidereal_time(time);
    spin(-15.0 * gast, vec1, vec2);
}

/*------------------ CalcMoon ------------------*/

#define DECLARE_PASCAL_ARRAY_1(elemtype,name,xmin,xmax) \
    elemtype name[(xmax)-(xmin)+1]

#define DECLARE_PASCAL_ARRAY_2(elemtype,name,xmin,xmax,ymin,ymax) \
    elemtype name[(xmax)-(xmin)+1][(ymax)-(ymin)+1]

#define ACCESS_PASCAL_ARRAY_1(name,xmin,x) \
    ((name)[(x)-(xmin)])

#define ACCESS_PASCAL_ARRAY_2(name,xmin,ymin,x,y) \
    ((name)[(x)-(xmin)][(y)-(ymin)])

typedef struct 
{
    double t;
    double dgam;
    double dlam, n, gam1c, sinpi;
    double l0, l, ls, f, d, s;
    double dl0, dl, dls, df, dd, ds;
    DECLARE_PASCAL_ARRAY_2(double,co,-6,6,1,4);   /* ARRAY[-6..6,1..4] OF REAL */
    DECLARE_PASCAL_ARRAY_2(double,si,-6,6,1,4);   /* ARRAY[-6..6,1..4] OF REAL */
}
MoonContext;

#define T           (ctx->t)
#define DGAM        (ctx->dgam)
#define DLAM        (ctx->dlam)
#define N           (ctx->n)
#define GAM1C       (ctx->gam1c)
#define SINPI       (ctx->sinpi)
#define L0          (ctx->l0)
#define L           (ctx->l)
#define LS          (ctx->ls)
#define F           (ctx->f)
#define D           (ctx->d)
#define S           (ctx->s)
#define DL0         (ctx->dl0)
#define DL          (ctx->dl)
#define DLS         (ctx->dls)
#define DF          (ctx->df)
#define DD          (ctx->dd)
#define DS          (ctx->ds)
#define CO(x,y)     ACCESS_PASCAL_ARRAY_2(ctx->co,-6,1,x,y)
#define SI(x,y)     ACCESS_PASCAL_ARRAY_2(ctx->si,-6,1,x,y)

static double Frac(double x)
{
    return x - floor(x);
}

static void AddThe(
    double c1, double s1, double c2, double s2,
    double *c, double *s)
{
    *c = c1*c2 - s1*s2;
    *s = s1*c2 + c1*s2;
}

static double Sine(double phi)
{
    /* sine, of phi in revolutions, not radians */
    return sin(PI2 * phi);
}

static void LongPeriodic(MoonContext *ctx)
{
    double S1 = Sine(0.19833+0.05611*T); 
    double S2 = Sine(0.27869+0.04508*T);
    double S3 = Sine(0.16827-0.36903*T); 
    double S4 = Sine(0.34734-5.37261*T);
    double S5 = Sine(0.10498-5.37899*T); 
    double S6 = Sine(0.42681-0.41855*T);
    double S7 = Sine(0.14943-5.37511*T);

    DL0 = 0.84*S1+0.31*S2+14.27*S3+ 7.26*S4+ 0.28*S5+0.24*S6;
    DL  = 2.94*S1+0.31*S2+14.27*S3+ 9.34*S4+ 1.12*S5+0.83*S6;
    DLS =-6.40*S1                                   -1.89*S6;
    DF  = 0.21*S1+0.31*S2+14.27*S3-88.70*S4-15.30*S5+0.24*S6-1.86*S7;
    DD  = DL0-DLS;
    DGAM  = -3332E-9 * Sine(0.59734-5.37261*T)
             -539E-9 * Sine(0.35498-5.37899*T)
              -64E-9 * Sine(0.39943-5.37511*T);
}

static void Init(MoonContext *ctx)
{
    int I, J, MAX;
    double T2, ARG, FAC;

    T2 = T*T;
    DLAM = 0; 
    DS = 0; 
    GAM1C = 0; 
    SINPI = 3422.7000;
    LongPeriodic(ctx);
    L0 = PI2*Frac(0.60643382+1336.85522467*T-0.00000313*T2) + DL0/ARC;
    L  = PI2*Frac(0.37489701+1325.55240982*T+0.00002565*T2) + DL /ARC;
    LS = PI2*Frac(0.99312619+  99.99735956*T-0.00000044*T2) + DLS/ARC;
    F  = PI2*Frac(0.25909118+1342.22782980*T-0.00000892*T2) + DF /ARC;
    D  = PI2*Frac(0.82736186+1236.85308708*T-0.00000397*T2) + DD /ARC;
    for (I=1; I<=4; ++I)
    {
        switch(I)
        {
            case 1: ARG=L;  MAX=4; FAC=1.000002208;               break;
            case 2: ARG=LS; MAX=3; FAC=0.997504612-0.002495388*T; break;
            case 3: ARG=F;  MAX=4; FAC=1.000002708+139.978*DGAM;  break;
            case 4: ARG=D;  MAX=6; FAC=1.0;                       break;
        }
        CO(0,1) = 1.0; 
        CO(1,I) = cos(ARG)*FAC;
        SI(0,I) = 0.0; 
        SI(1,I) = sin(ARG)*FAC;
        for (J=2; J<=MAX; ++J)
            AddThe(CO(J-1,I), SI(J-1,I), CO(1,I), SI(1,I), &CO(J,I), &SI(J,I));

        for (J=1; J<=MAX; ++J)
        {
            CO(-J,I) =  CO(J,I);
            SI(-J,I) = -SI(J,I);
        }
    }
}

static void Term(MoonContext *ctx, int p, int q, int r, int s, double *x, double *y)
{
    int k;
    DECLARE_PASCAL_ARRAY_1(int, i, 1, 4);
    #define I(n) ACCESS_PASCAL_ARRAY_1(i,1,n)

    I(1) = p;
    I(2) = q;
    I(3) = r;
    I(4) = s;
    *x = 1.0;
    *y = 0.0;

    for (k=1; k<=4; ++k)
        if (I(k) != 0.0)
            AddThe(*x, *y, CO(I(k), k), SI(I(k), k), x, y);

    #undef I
}

static void AddSol(
    MoonContext *ctx, 
    double coeffl,
    double coeffs,
    double coeffg,
    double coeffp,
    int p,
    int q,
    int r,
    int s)
{
    double x, y;
    Term(ctx, p, q, r, s, &x, &y);
    DLAM += coeffl*y;
    DS += coeffs*y;
    GAM1C += coeffg*x;
    SINPI += coeffp*x;
}

static void Solar1(MoonContext *ctx)
{
    AddSol(ctx,    13.902,   14.06,-0.001,   0.2607,0, 0, 0, 4);
    AddSol(ctx,     0.403,   -4.01,+0.394,   0.0023,0, 0, 0, 3);
    AddSol(ctx,  2369.912, 2373.36,+0.601,  28.2333,0, 0, 0, 2);
    AddSol(ctx,  -125.154, -112.79,-0.725,  -0.9781,0, 0, 0, 1);
    AddSol(ctx,     1.979,    6.98,-0.445,   0.0433,1, 0, 0, 4);
    AddSol(ctx,   191.953,  192.72,+0.029,   3.0861,1, 0, 0, 2);
    AddSol(ctx,    -8.466,  -13.51,+0.455,  -0.1093,1, 0, 0, 1);
    AddSol(ctx, 22639.500,22609.07,+0.079, 186.5398,1, 0, 0, 0);
    AddSol(ctx,    18.609,    3.59,-0.094,   0.0118,1, 0, 0,-1);
    AddSol(ctx, -4586.465,-4578.13,-0.077,  34.3117,1, 0, 0,-2);
    AddSol(ctx,    +3.215,    5.44,+0.192,  -0.0386,1, 0, 0,-3);
    AddSol(ctx,   -38.428,  -38.64,+0.001,   0.6008,1, 0, 0,-4);
    AddSol(ctx,    -0.393,   -1.43,-0.092,   0.0086,1, 0, 0,-6);
    AddSol(ctx,    -0.289,   -1.59,+0.123,  -0.0053,0, 1, 0, 4);
    AddSol(ctx,   -24.420,  -25.10,+0.040,  -0.3000,0, 1, 0, 2);
    AddSol(ctx,    18.023,   17.93,+0.007,   0.1494,0, 1, 0, 1);
    AddSol(ctx,  -668.146, -126.98,-1.302,  -0.3997,0, 1, 0, 0);
    AddSol(ctx,     0.560,    0.32,-0.001,  -0.0037,0, 1, 0,-1);
    AddSol(ctx,  -165.145, -165.06,+0.054,   1.9178,0, 1, 0,-2);
    AddSol(ctx,    -1.877,   -6.46,-0.416,   0.0339,0, 1, 0,-4);
    AddSol(ctx,     0.213,    1.02,-0.074,   0.0054,2, 0, 0, 4);
    AddSol(ctx,    14.387,   14.78,-0.017,   0.2833,2, 0, 0, 2);
    AddSol(ctx,    -0.586,   -1.20,+0.054,  -0.0100,2, 0, 0, 1);
    AddSol(ctx,   769.016,  767.96,+0.107,  10.1657,2, 0, 0, 0);
    AddSol(ctx,    +1.750,    2.01,-0.018,   0.0155,2, 0, 0,-1);
    AddSol(ctx,  -211.656, -152.53,+5.679,  -0.3039,2, 0, 0,-2);
    AddSol(ctx,    +1.225,    0.91,-0.030,  -0.0088,2, 0, 0,-3);
    AddSol(ctx,   -30.773,  -34.07,-0.308,   0.3722,2, 0, 0,-4);
    AddSol(ctx,    -0.570,   -1.40,-0.074,   0.0109,2, 0, 0,-6);
    AddSol(ctx,    -2.921,  -11.75,+0.787,  -0.0484,1, 1, 0, 2);
    AddSol(ctx,    +1.267,    1.52,-0.022,   0.0164,1, 1, 0, 1);
    AddSol(ctx,  -109.673, -115.18,+0.461,  -0.9490,1, 1, 0, 0);
    AddSol(ctx,  -205.962, -182.36,+2.056,  +1.4437,1, 1, 0,-2);
    AddSol(ctx,     0.233,    0.36, 0.012,  -0.0025,1, 1, 0,-3);
    AddSol(ctx,    -4.391,   -9.66,-0.471,   0.0673,1, 1, 0,-4);
}

static void Solar2(MoonContext *ctx)
{
    AddSol(ctx,     0.283,    1.53,-0.111,  +0.0060,1,-1, 0,+4);
    AddSol(ctx,    14.577,   31.70,-1.540,  +0.2302,1,-1, 0, 2);
    AddSol(ctx,   147.687,  138.76,+0.679,  +1.1528,1,-1, 0, 0);
    AddSol(ctx,    -1.089,    0.55,+0.021,   0.0   ,1,-1, 0,-1);
    AddSol(ctx,    28.475,   23.59,-0.443,  -0.2257,1,-1, 0,-2);
    AddSol(ctx,    -0.276,   -0.38,-0.006,  -0.0036,1,-1, 0,-3);
    AddSol(ctx,     0.636,    2.27,+0.146,  -0.0102,1,-1, 0,-4);
    AddSol(ctx,    -0.189,   -1.68,+0.131,  -0.0028,0, 2, 0, 2);
    AddSol(ctx,    -7.486,   -0.66,-0.037,  -0.0086,0, 2, 0, 0);
    AddSol(ctx,    -8.096,  -16.35,-0.740,   0.0918,0, 2, 0,-2);
    AddSol(ctx,    -5.741,   -0.04, 0.0  ,  -0.0009,0, 0, 2, 2);
    AddSol(ctx,     0.255,    0.0 , 0.0  ,   0.0   ,0, 0, 2, 1);
    AddSol(ctx,  -411.608,   -0.20, 0.0  ,  -0.0124,0, 0, 2, 0);
    AddSol(ctx,     0.584,    0.84, 0.0  ,  +0.0071,0, 0, 2,-1);
    AddSol(ctx,   -55.173,  -52.14, 0.0  ,  -0.1052,0, 0, 2,-2);
    AddSol(ctx,     0.254,    0.25, 0.0  ,  -0.0017,0, 0, 2,-3);
    AddSol(ctx,    +0.025,   -1.67, 0.0  ,  +0.0031,0, 0, 2,-4);
    AddSol(ctx,     1.060,    2.96,-0.166,   0.0243,3, 0, 0,+2);
    AddSol(ctx,    36.124,   50.64,-1.300,   0.6215,3, 0, 0, 0);
    AddSol(ctx,   -13.193,  -16.40,+0.258,  -0.1187,3, 0, 0,-2);
    AddSol(ctx,    -1.187,   -0.74,+0.042,   0.0074,3, 0, 0,-4);
    AddSol(ctx,    -0.293,   -0.31,-0.002,   0.0046,3, 0, 0,-6);
    AddSol(ctx,    -0.290,   -1.45,+0.116,  -0.0051,2, 1, 0, 2);
    AddSol(ctx,    -7.649,  -10.56,+0.259,  -0.1038,2, 1, 0, 0);
    AddSol(ctx,    -8.627,   -7.59,+0.078,  -0.0192,2, 1, 0,-2);
    AddSol(ctx,    -2.740,   -2.54,+0.022,   0.0324,2, 1, 0,-4);
    AddSol(ctx,     1.181,    3.32,-0.212,   0.0213,2,-1, 0,+2);
    AddSol(ctx,     9.703,   11.67,-0.151,   0.1268,2,-1, 0, 0);
    AddSol(ctx,    -0.352,   -0.37,+0.001,  -0.0028,2,-1, 0,-1);
    AddSol(ctx,    -2.494,   -1.17,-0.003,  -0.0017,2,-1, 0,-2);
    AddSol(ctx,     0.360,    0.20,-0.012,  -0.0043,2,-1, 0,-4);
    AddSol(ctx,    -1.167,   -1.25,+0.008,  -0.0106,1, 2, 0, 0);
    AddSol(ctx,    -7.412,   -6.12,+0.117,   0.0484,1, 2, 0,-2);
    AddSol(ctx,    -0.311,   -0.65,-0.032,   0.0044,1, 2, 0,-4);
    AddSol(ctx,    +0.757,    1.82,-0.105,   0.0112,1,-2, 0, 2);
    AddSol(ctx,    +2.580,    2.32,+0.027,   0.0196,1,-2, 0, 0);
    AddSol(ctx,    +2.533,    2.40,-0.014,  -0.0212,1,-2, 0,-2);
    AddSol(ctx,    -0.344,   -0.57,-0.025,  +0.0036,0, 3, 0,-2);
    AddSol(ctx,    -0.992,   -0.02, 0.0  ,   0.0   ,1, 0, 2, 2);
    AddSol(ctx,   -45.099,   -0.02, 0.0  ,  -0.0010,1, 0, 2, 0);
    AddSol(ctx,    -0.179,   -9.52, 0.0  ,  -0.0833,1, 0, 2,-2);
    AddSol(ctx,    -0.301,   -0.33, 0.0  ,   0.0014,1, 0, 2,-4);
    AddSol(ctx,    -6.382,   -3.37, 0.0  ,  -0.0481,1, 0,-2, 2);
    AddSol(ctx,    39.528,   85.13, 0.0  ,  -0.7136,1, 0,-2, 0);
    AddSol(ctx,     9.366,    0.71, 0.0  ,  -0.0112,1, 0,-2,-2);
    AddSol(ctx,     0.202,    0.02, 0.0  ,   0.0   ,1, 0,-2,-4);
}

static void Solar3(MoonContext *ctx)
{
    AddSol(ctx,     0.415,    0.10, 0.0  ,  0.0013,0, 1, 2, 0);
    AddSol(ctx,    -2.152,   -2.26, 0.0  , -0.0066,0, 1, 2,-2);
    AddSol(ctx,    -1.440,   -1.30, 0.0  , +0.0014,0, 1,-2, 2);
    AddSol(ctx,     0.384,   -0.04, 0.0  ,  0.0   ,0, 1,-2,-2);
    AddSol(ctx,    +1.938,   +3.60,-0.145, +0.0401,4, 0, 0, 0);
    AddSol(ctx,    -0.952,   -1.58,+0.052, -0.0130,4, 0, 0,-2);
    AddSol(ctx,    -0.551,   -0.94,+0.032, -0.0097,3, 1, 0, 0);
    AddSol(ctx,    -0.482,   -0.57,+0.005, -0.0045,3, 1, 0,-2);
    AddSol(ctx,     0.681,    0.96,-0.026,  0.0115,3,-1, 0, 0);
    AddSol(ctx,    -0.297,   -0.27, 0.002, -0.0009,2, 2, 0,-2);
    AddSol(ctx,     0.254,   +0.21,-0.003,  0.0   ,2,-2, 0,-2);
    AddSol(ctx,    -0.250,   -0.22, 0.004,  0.0014,1, 3, 0,-2);
    AddSol(ctx,    -3.996,    0.0 , 0.0  , +0.0004,2, 0, 2, 0);
    AddSol(ctx,     0.557,   -0.75, 0.0  , -0.0090,2, 0, 2,-2);
    AddSol(ctx,    -0.459,   -0.38, 0.0  , -0.0053,2, 0,-2, 2);
    AddSol(ctx,    -1.298,    0.74, 0.0  , +0.0004,2, 0,-2, 0);
    AddSol(ctx,     0.538,    1.14, 0.0  , -0.0141,2, 0,-2,-2);
    AddSol(ctx,     0.263,    0.02, 0.0  ,  0.0   ,1, 1, 2, 0);
    AddSol(ctx,     0.426,   +0.07, 0.0  , -0.0006,1, 1,-2,-2);
    AddSol(ctx,    -0.304,   +0.03, 0.0  , +0.0003,1,-1, 2, 0);
    AddSol(ctx,    -0.372,   -0.19, 0.0  , -0.0027,1,-1,-2, 2);
    AddSol(ctx,    +0.418,    0.0 , 0.0  ,  0.0   ,0, 0, 4, 0);
    AddSol(ctx,    -0.330,   -0.04, 0.0  ,  0.0   ,3, 0, 2, 0);
}

#define ADDN(coeffn,p,q,r,s)    ( Term(ctx, (p),(q),(r),(s),&x,&y), (N += (coeffn)*y) )

static void SolarN(MoonContext *ctx)
{
    double x, y;

    N = 0.0;
    ADDN(-526.069, 0, 0,1,-2); 
    ADDN(  -3.352, 0, 0,1,-4);
    ADDN( +44.297,+1, 0,1,-2); 
    ADDN(  -6.000,+1, 0,1,-4);
    ADDN( +20.599,-1, 0,1, 0); 
    ADDN( -30.598,-1, 0,1,-2);
    ADDN( -24.649,-2, 0,1, 0); 
    ADDN(  -2.000,-2, 0,1,-2);
    ADDN( -22.571, 0,+1,1,-2); 
    ADDN( +10.985, 0,-1,1,-2);
}

static void Planetary(MoonContext *ctx)
{
    DLAM +=
        +0.82*Sine(0.7736  -62.5512*T)+0.31*Sine(0.0466 -125.1025*T)
        +0.35*Sine(0.5785  -25.1042*T)+0.66*Sine(0.4591+1335.8075*T)
        +0.64*Sine(0.3130  -91.5680*T)+1.14*Sine(0.1480+1331.2898*T)
        +0.21*Sine(0.5918+1056.5859*T)+0.44*Sine(0.5784+1322.8595*T)
        +0.24*Sine(0.2275   -5.7374*T)+0.28*Sine(0.2965   +2.6929*T)
        +0.33*Sine(0.3132   +6.3368*T);
}

static void CalcMoon(
    double centuries_since_j2000,
    double *geo_eclip_lon,      /* (LAMBDA) equinox of date */
    double *geo_eclip_lat,      /* (BETA)   equinox of date */
    double *distance_au)        /* (R) */
{
    double lat_seconds;
    MoonContext context;
    MoonContext *ctx = &context;    /* goofy, but makes macros work inside this function */

    context.t = centuries_since_j2000;
    Init(ctx);
    Solar1(ctx);
    Solar2(ctx);
    Solar3(ctx);
    SolarN(ctx);
    Planetary(ctx);
    S = F + DS/ARC;

    lat_seconds = (1.000002708 + 139.978*DGAM)*(18518.511+1.189+GAM1C)*sin(S)-6.24*sin(3*S) + N;

    *geo_eclip_lon = PI2 * Frac((L0+DLAM/ARC) / PI2);
    *geo_eclip_lat = lat_seconds * (DEG2RAD / 3600.0);
    *distance_au = (ARC * (ERAD / AU)) / (0.999953253 * SINPI);
}

#undef T
#undef DGAM
#undef DLAM
#undef N
#undef GAM1C
#undef SINPI
#undef L0
#undef L
#undef LS
#undef F
#undef D
#undef S
#undef DL0
#undef DL
#undef DLS
#undef DF
#undef DD
#undef DS
#undef CO
#undef SI

astro_vector_t Astronomy_GeoMoon(astro_time_t time)
{
    double geo_eclip_lon, geo_eclip_lat, distance_au;
    double dist_cos_lat;
    astro_vector_t vector;
    double gepos[3];
    double mpos1[3];
    double mpos2[3];

    CalcMoon(time.tt / 36525.0, &geo_eclip_lon, &geo_eclip_lat, &distance_au);

    /* Convert geocentric ecliptic spherical coordinates to Cartesian coordinates. */
    dist_cos_lat = distance_au * cos(geo_eclip_lat);
    gepos[0] = dist_cos_lat * cos(geo_eclip_lon);
    gepos[1] = dist_cos_lat * sin(geo_eclip_lon);
    gepos[2] = distance_au * sin(geo_eclip_lat);

    /* Convert ecliptic coordinates to equatorial coordinates, both in mean equinox of date. */
    ecl2equ_vec(time, gepos, mpos1);

    /* Convert from mean equinox of date to J2000. */
    precession(time.tt, mpos1, 0, mpos2);

    vector.status = ASTRO_SUCCESS;
    vector.x = mpos2[0];
    vector.y = mpos2[1];
    vector.z = mpos2[2];
    vector.t = time;
    return vector;
}

/*------------------ VSOP ------------------*/

typedef struct 
{
    double amplitude;
    double phase;
    double frequency;
}
vsop_term_t;

typedef struct 
{
    int nterms;
    const vsop_term_t *term; 
}
vsop_series_t;

typedef struct
{
    int nseries;
    const vsop_series_t *series;
}
vsop_formula_t;

typedef struct 
{
    const vsop_formula_t formula[3];
}
vsop_model_t;

static const vsop_term_t vsop_lat_Mercury_0[] = 
{
    { 4.40250710144, 0.00000000000, 0.00000000000 },
    { 0.40989414977, 1.48302034195, 26087.90314157420 },
    { 0.05046294200, 4.47785489551, 52175.80628314840 },
    { 0.00855346844, 1.16520322459, 78263.70942472259 },
    { 0.00165590362, 4.11969163423, 104351.61256629678 },
    { 0.00034561897, 0.77930768443, 130439.51570787099 },
    { 0.00007583476, 3.71348404924, 156527.41884944518 }
};

static const vsop_term_t vsop_lat_Mercury_1[] = 
{
    { 26087.90313685529, 0.00000000000, 0.00000000000 },
    { 0.01131199811, 6.21874197797, 26087.90314157420 },
    { 0.00292242298, 3.04449355541, 52175.80628314840 },
    { 0.00075775081, 6.08568821653, 78263.70942472259 },
    { 0.00019676525, 2.80965111777, 104351.61256629678 }
};

static const vsop_series_t vsop_lat_Mercury[] = 
{
    { 7, vsop_lat_Mercury_0 },
    { 5, vsop_lat_Mercury_1 }
};

static const vsop_term_t vsop_lon_Mercury_0[] = 
{
    { 0.11737528961, 1.98357498767, 26087.90314157420 },
    { 0.02388076996, 5.03738959686, 52175.80628314840 },
    { 0.01222839532, 3.14159265359, 0.00000000000 },
    { 0.00543251810, 1.79644363964, 78263.70942472259 },
    { 0.00129778770, 4.83232503958, 104351.61256629678 },
    { 0.00031866927, 1.58088495658, 130439.51570787099 },
    { 0.00007963301, 4.60972126127, 156527.41884944518 }
};

static const vsop_term_t vsop_lon_Mercury_1[] = 
{
    { 0.00274646065, 3.95008450011, 26087.90314157420 },
    { 0.00099737713, 3.14159265359, 0.00000000000 }
};

static const vsop_series_t vsop_lon_Mercury[] = 
{
    { 7, vsop_lon_Mercury_0 },
    { 2, vsop_lon_Mercury_1 }
};

static const vsop_term_t vsop_rad_Mercury_0[] = 
{
    { 0.39528271651, 0.00000000000, 0.00000000000 },
    { 0.07834131818, 6.19233722598, 26087.90314157420 },
    { 0.00795525558, 2.95989690104, 52175.80628314840 },
    { 0.00121281764, 6.01064153797, 78263.70942472259 },
    { 0.00021921969, 2.77820093972, 104351.61256629678 },
    { 0.00004354065, 5.82894543774, 130439.51570787099 }
};

static const vsop_term_t vsop_rad_Mercury_1[] = 
{
    { 0.00217347740, 4.65617158665, 26087.90314157420 },
    { 0.00044141826, 1.42385544001, 52175.80628314840 }
};

static const vsop_series_t vsop_rad_Mercury[] = 
{
    { 6, vsop_rad_Mercury_0 },
    { 2, vsop_rad_Mercury_1 }
};

;
static const vsop_term_t vsop_lat_Venus_0[] = 
{
    { 3.17614666774, 0.00000000000, 0.00000000000 },
    { 0.01353968419, 5.59313319619, 10213.28554621100 },
    { 0.00089891645, 5.30650047764, 20426.57109242200 },
    { 0.00005477194, 4.41630661466, 7860.41939243920 },
    { 0.00003455741, 2.69964447820, 11790.62908865880 },
    { 0.00002372061, 2.99377542079, 3930.20969621960 },
    { 0.00001317168, 5.18668228402, 26.29831979980 },
    { 0.00001664146, 4.25018630147, 1577.34354244780 },
    { 0.00001438387, 4.15745084182, 9683.59458111640 },
    { 0.00001200521, 6.15357116043, 30639.85663863300 }
};

static const vsop_term_t vsop_lat_Venus_1[] = 
{
    { 10213.28554621638, 0.00000000000, 0.00000000000 },
    { 0.00095617813, 2.46406511110, 10213.28554621100 },
    { 0.00007787201, 0.62478482220, 20426.57109242200 }
};

static const vsop_series_t vsop_lat_Venus[] = 
{
    { 10, vsop_lat_Venus_0 },
    { 3, vsop_lat_Venus_1 }
};

static const vsop_term_t vsop_lon_Venus_0[] = 
{
    { 0.05923638472, 0.26702775812, 10213.28554621100 },
    { 0.00040107978, 1.14737178112, 20426.57109242200 },
    { 0.00032814918, 3.14159265359, 0.00000000000 }
};

static const vsop_term_t vsop_lon_Venus_1[] = 
{
    { 0.00287821243, 1.88964962838, 10213.28554621100 }
};

static const vsop_series_t vsop_lon_Venus[] = 
{
    { 3, vsop_lon_Venus_0 },
    { 1, vsop_lon_Venus_1 }
};

static const vsop_term_t vsop_rad_Venus_0[] = 
{
    { 0.72334820891, 0.00000000000, 0.00000000000 },
    { 0.00489824182, 4.02151831717, 10213.28554621100 },
    { 0.00001658058, 4.90206728031, 20426.57109242200 }
};

static const vsop_term_t vsop_rad_Venus_1[] = 
{
    { 0.00034551041, 0.89198706276, 10213.28554621100 }
};

static const vsop_series_t vsop_rad_Venus[] = 
{
    { 3, vsop_rad_Venus_0 },
    { 1, vsop_rad_Venus_1 }
};

;
static const vsop_term_t vsop_lat_Earth_0[] = 
{
    { 1.75347045673, 0.00000000000, 0.00000000000 },
    { 0.03341656453, 4.66925680415, 6283.07584999140 },
    { 0.00034894275, 4.62610242189, 12566.15169998280 },
    { 0.00003417572, 2.82886579754, 3.52311834900 },
    { 0.00003497056, 2.74411783405, 5753.38488489680 },
    { 0.00003135899, 3.62767041756, 77713.77146812050 },
    { 0.00002676218, 4.41808345438, 7860.41939243920 },
    { 0.00002342691, 6.13516214446, 3930.20969621960 },
    { 0.00001273165, 2.03709657878, 529.69096509460 },
    { 0.00001324294, 0.74246341673, 11506.76976979360 },
    { 0.00000901854, 2.04505446477, 26.29831979980 },
    { 0.00001199167, 1.10962946234, 1577.34354244780 },
    { 0.00000857223, 3.50849152283, 398.14900340820 },
    { 0.00000779786, 1.17882681962, 5223.69391980220 },
    { 0.00000990250, 5.23268072088, 5884.92684658320 },
    { 0.00000753141, 2.53339052847, 5507.55323866740 },
    { 0.00000505267, 4.58292599973, 18849.22754997420 },
    { 0.00000492392, 4.20505711826, 775.52261132400 },
    { 0.00000356672, 2.91954114478, 0.06731030280 },
    { 0.00000284125, 1.89869240932, 796.29800681640 },
    { 0.00000242879, 0.34481445893, 5486.77784317500 },
    { 0.00000317087, 5.84901948512, 11790.62908865880 },
    { 0.00000271112, 0.31486255375, 10977.07880469900 },
    { 0.00000206217, 4.80646631478, 2544.31441988340 },
    { 0.00000205478, 1.86953770281, 5573.14280143310 },
    { 0.00000202318, 2.45767790232, 6069.77675455340 },
    { 0.00000126225, 1.08295459501, 20.77539549240 },
    { 0.00000155516, 0.83306084617, 213.29909543800 }
};

static const vsop_term_t vsop_lat_Earth_1[] = 
{
    { 6283.07584999140, 0.00000000000, 0.00000000000 },
    { 0.00206058863, 2.67823455808, 6283.07584999140 },
    { 0.00004303419, 2.63512233481, 12566.15169998280 }
};

static const vsop_term_t vsop_lat_Earth_2[] = 
{
    { 0.00008721859, 1.07253635559, 6283.07584999140 }
};

static const vsop_series_t vsop_lat_Earth[] = 
{
    { 28, vsop_lat_Earth_0 },
    { 3, vsop_lat_Earth_1 },
    { 1, vsop_lat_Earth_2 }
};

static const vsop_term_t vsop_lon_Earth_0[] = 
{
};

static const vsop_term_t vsop_lon_Earth_1[] = 
{
    { 0.00227777722, 3.41376620530, 6283.07584999140 },
    { 0.00003805678, 3.37063423795, 12566.15169998280 }
};

static const vsop_series_t vsop_lon_Earth[] = 
{
    { 0, vsop_lon_Earth_0 },
    { 2, vsop_lon_Earth_1 }
};

static const vsop_term_t vsop_rad_Earth_0[] = 
{
    { 1.00013988784, 0.00000000000, 0.00000000000 },
    { 0.01670699632, 3.09846350258, 6283.07584999140 },
    { 0.00013956024, 3.05524609456, 12566.15169998280 },
    { 0.00003083720, 5.19846674381, 77713.77146812050 },
    { 0.00001628463, 1.17387558054, 5753.38488489680 },
    { 0.00001575572, 2.84685214877, 7860.41939243920 },
    { 0.00000924799, 5.45292236722, 11506.76976979360 },
    { 0.00000542439, 4.56409151453, 3930.20969621960 },
    { 0.00000472110, 3.66100022149, 5884.92684658320 }
};

static const vsop_term_t vsop_rad_Earth_1[] = 
{
    { 0.00103018607, 1.10748968172, 6283.07584999140 },
    { 0.00001721238, 1.06442300386, 12566.15169998280 }
};

static const vsop_term_t vsop_rad_Earth_2[] = 
{
    { 0.00004359385, 5.78455133808, 6283.07584999140 }
};

static const vsop_series_t vsop_rad_Earth[] = 
{
    { 9, vsop_rad_Earth_0 },
    { 2, vsop_rad_Earth_1 },
    { 1, vsop_rad_Earth_2 }
};

;
static const vsop_term_t vsop_lat_Mars_0[] = 
{
    { 6.20347711581, 0.00000000000, 0.00000000000 },
    { 0.18656368093, 5.05037100270, 3340.61242669980 },
    { 0.01108216816, 5.40099836344, 6681.22485339960 },
    { 0.00091798406, 5.75478744667, 10021.83728009940 },
    { 0.00027744987, 5.97049513147, 3.52311834900 },
    { 0.00010610235, 2.93958560338, 2281.23049651060 },
    { 0.00012315897, 0.84956094002, 2810.92146160520 },
    { 0.00008926784, 4.15697846427, 0.01725365220 },
    { 0.00008715691, 6.11005153139, 13362.44970679920 },
    { 0.00006797556, 0.36462229657, 398.14900340820 },
    { 0.00007774872, 3.33968761376, 5621.84292321040 },
    { 0.00003575078, 1.66186505710, 2544.31441988340 },
    { 0.00004161108, 0.22814971327, 2942.46342329160 },
    { 0.00003075252, 0.85696614132, 191.44826611160 },
    { 0.00002628117, 0.64806124465, 3337.08930835080 },
    { 0.00002937546, 6.07893711402, 0.06731030280 },
    { 0.00002389414, 5.03896442664, 796.29800681640 },
    { 0.00002579844, 0.02996736156, 3344.13554504880 },
    { 0.00001528141, 1.14979301996, 6151.53388830500 },
    { 0.00001798806, 0.65634057445, 529.69096509460 },
    { 0.00001264357, 3.62275122593, 5092.15195811580 },
    { 0.00001286228, 3.06796065034, 2146.16541647520 },
    { 0.00001546404, 2.91579701718, 1751.53953141600 },
    { 0.00001024902, 3.69334099279, 8962.45534991020 },
    { 0.00000891566, 0.18293837498, 16703.06213349900 },
    { 0.00000858759, 2.40093811940, 2914.01423582380 },
    { 0.00000832715, 2.46418619474, 3340.59517304760 },
    { 0.00000832720, 4.49495782139, 3340.62968035200 },
    { 0.00000712902, 3.66335473479, 1059.38193018920 },
    { 0.00000748723, 3.82248614017, 155.42039943420 },
    { 0.00000723861, 0.67497311481, 3738.76143010800 },
    { 0.00000635548, 2.92182225127, 8432.76438481560 },
    { 0.00000655162, 0.48864064125, 3127.31333126180 },
    { 0.00000550474, 3.81001042328, 0.98032106820 },
    { 0.00000552750, 4.47479317037, 1748.01641306700 },
    { 0.00000425966, 0.55364317304, 6283.07584999140 },
    { 0.00000415131, 0.49662285038, 213.29909543800 },
    { 0.00000472167, 3.62547124025, 1194.44701022460 },
    { 0.00000306551, 0.38052848348, 6684.74797174860 },
    { 0.00000312141, 0.99853944405, 6677.70173505060 },
    { 0.00000293198, 4.22131299634, 20.77539549240 },
    { 0.00000302375, 4.48618007156, 3532.06069281140 },
    { 0.00000274027, 0.54222167059, 3340.54511639700 },
    { 0.00000281079, 5.88163521788, 1349.86740965880 },
    { 0.00000231183, 1.28242156993, 3870.30339179440 },
    { 0.00000283602, 5.76885434940, 3149.16416058820 },
    { 0.00000236117, 5.75503217933, 3333.49887969900 },
    { 0.00000274033, 0.13372524985, 3340.67973700260 },
    { 0.00000299395, 2.78323740866, 6254.62666252360 }
};

static const vsop_term_t vsop_lat_Mars_1[] = 
{
    { 3340.61242700512, 0.00000000000, 0.00000000000 },
    { 0.01457554523, 3.60433733236, 3340.61242669980 },
    { 0.00168414711, 3.92318567804, 6681.22485339960 },
    { 0.00020622975, 4.26108844583, 10021.83728009940 },
    { 0.00003452392, 4.73210393190, 3.52311834900 },
    { 0.00002586332, 4.60670058555, 13362.44970679920 },
    { 0.00000841535, 4.45864030426, 2281.23049651060 }
};

static const vsop_term_t vsop_lat_Mars_2[] = 
{
    { 0.00058152577, 2.04961712429, 3340.61242669980 },
    { 0.00013459579, 2.45738706163, 6681.22485339960 }
};

static const vsop_series_t vsop_lat_Mars[] = 
{
    { 49, vsop_lat_Mars_0 },
    { 7, vsop_lat_Mars_1 },
    { 2, vsop_lat_Mars_2 }
};

static const vsop_term_t vsop_lon_Mars_0[] = 
{
    { 0.03197134986, 3.76832042431, 3340.61242669980 },
    { 0.00298033234, 4.10616996305, 6681.22485339960 },
    { 0.00289104742, 0.00000000000, 0.00000000000 },
    { 0.00031365539, 4.44651053090, 10021.83728009940 },
    { 0.00003484100, 4.78812549260, 13362.44970679920 }
};

static const vsop_term_t vsop_lon_Mars_1[] = 
{
    { 0.00217310991, 6.04472194776, 3340.61242669980 },
    { 0.00020976948, 3.14159265359, 0.00000000000 },
    { 0.00012834709, 1.60810667915, 6681.22485339960 }
};

static const vsop_series_t vsop_lon_Mars[] = 
{
    { 5, vsop_lon_Mars_0 },
    { 3, vsop_lon_Mars_1 }
};

static const vsop_term_t vsop_rad_Mars_0[] = 
{
    { 1.53033488271, 0.00000000000, 0.00000000000 },
    { 0.14184953160, 3.47971283528, 3340.61242669980 },
    { 0.00660776362, 3.81783443019, 6681.22485339960 },
    { 0.00046179117, 4.15595316782, 10021.83728009940 },
    { 0.00008109733, 5.55958416318, 2810.92146160520 },
    { 0.00007485318, 1.77239078402, 5621.84292321040 },
    { 0.00005523191, 1.36436303770, 2281.23049651060 },
    { 0.00003825160, 4.49407183687, 13362.44970679920 },
    { 0.00002306537, 0.09081579001, 2544.31441988340 },
    { 0.00001999396, 5.36059617709, 3337.08930835080 },
    { 0.00002484394, 4.92545639920, 2942.46342329160 },
    { 0.00001960195, 4.74249437639, 3344.13554504880 },
    { 0.00001167119, 2.11260868341, 5092.15195811580 },
    { 0.00001102816, 5.00908403998, 398.14900340820 },
    { 0.00000899066, 4.40791133207, 529.69096509460 },
    { 0.00000992252, 5.83861961952, 6151.53388830500 },
    { 0.00000807354, 2.10217065501, 1059.38193018920 },
    { 0.00000797915, 3.44839203899, 796.29800681640 },
    { 0.00000740975, 1.49906336885, 2146.16541647520 }
};

static const vsop_term_t vsop_rad_Mars_1[] = 
{
    { 0.01107433345, 2.03250524857, 3340.61242669980 },
    { 0.00103175887, 2.37071847807, 6681.22485339960 },
    { 0.00012877200, 0.00000000000, 0.00000000000 },
    { 0.00010815880, 2.70888095665, 10021.83728009940 }
};

static const vsop_term_t vsop_rad_Mars_2[] = 
{
    { 0.00044242249, 0.47930604954, 3340.61242669980 },
    { 0.00008138042, 0.86998389204, 6681.22485339960 }
};

static const vsop_series_t vsop_rad_Mars[] = 
{
    { 19, vsop_rad_Mars_0 },
    { 4, vsop_rad_Mars_1 },
    { 2, vsop_rad_Mars_2 }
};

;
static const vsop_term_t vsop_lat_Jupiter_0[] = 
{
    { 0.59954691494, 0.00000000000, 0.00000000000 },
    { 0.09695898719, 5.06191793158, 529.69096509460 },
    { 0.00573610142, 1.44406205629, 7.11354700080 },
    { 0.00306389205, 5.41734730184, 1059.38193018920 },
    { 0.00097178296, 4.14264726552, 632.78373931320 },
    { 0.00072903078, 3.64042916389, 522.57741809380 },
    { 0.00064263975, 3.41145165351, 103.09277421860 },
    { 0.00039806064, 2.29376740788, 419.48464387520 },
    { 0.00038857767, 1.27231755835, 316.39186965660 },
    { 0.00027964629, 1.78454591820, 536.80451209540 },
    { 0.00013589730, 5.77481040790, 1589.07289528380 },
    { 0.00008246349, 3.58227925840, 206.18554843720 },
    { 0.00008768704, 3.63000308199, 949.17560896980 },
    { 0.00007368042, 5.08101194270, 735.87651353180 },
    { 0.00006263150, 0.02497628807, 213.29909543800 },
    { 0.00006114062, 4.51319998626, 1162.47470440780 },
    { 0.00004905396, 1.32084470588, 110.20632121940 },
    { 0.00005305285, 1.30671216791, 14.22709400160 },
    { 0.00005305441, 4.18625634012, 1052.26838318840 },
    { 0.00004647248, 4.69958103684, 3.93215326310 },
    { 0.00003045023, 4.31676431084, 426.59819087600 },
    { 0.00002609999, 1.56667394063, 846.08283475120 },
    { 0.00002028191, 1.06376530715, 3.18139373770 },
    { 0.00001764763, 2.14148655117, 1066.49547719000 },
    { 0.00001722972, 3.88036268267, 1265.56747862640 },
    { 0.00001920945, 0.97168196472, 639.89728631400 },
    { 0.00001633223, 3.58201833555, 515.46387109300 },
    { 0.00001431999, 4.29685556046, 625.67019231240 },
    { 0.00000973272, 4.09764549134, 95.97922721780 }
};

static const vsop_term_t vsop_lat_Jupiter_1[] = 
{
    { 529.69096508814, 0.00000000000, 0.00000000000 },
    { 0.00489503243, 4.22082939470, 529.69096509460 },
    { 0.00228917222, 6.02646855621, 7.11354700080 },
    { 0.00030099479, 4.54540782858, 1059.38193018920 },
    { 0.00020720920, 5.45943156902, 522.57741809380 },
    { 0.00012103653, 0.16994816098, 536.80451209540 },
    { 0.00006067987, 4.42422292017, 103.09277421860 },
    { 0.00005433968, 3.98480737746, 419.48464387520 },
    { 0.00004237744, 5.89008707199, 14.22709400160 }
};

static const vsop_term_t vsop_lat_Jupiter_2[] = 
{
    { 0.00047233601, 4.32148536482, 7.11354700080 },
    { 0.00030649436, 2.92977788700, 529.69096509460 },
    { 0.00014837605, 3.14159265359, 0.00000000000 }
};

static const vsop_series_t vsop_lat_Jupiter[] = 
{
    { 29, vsop_lat_Jupiter_0 },
    { 9, vsop_lat_Jupiter_1 },
    { 3, vsop_lat_Jupiter_2 }
};

static const vsop_term_t vsop_lon_Jupiter_0[] = 
{
    { 0.02268615702, 3.55852606721, 529.69096509460 },
    { 0.00109971634, 3.90809347197, 1059.38193018920 },
    { 0.00110090358, 0.00000000000, 0.00000000000 },
    { 0.00008101428, 3.60509572885, 522.57741809380 },
    { 0.00006043996, 4.25883108339, 1589.07289528380 },
    { 0.00006437782, 0.30627119215, 536.80451209540 }
};

static const vsop_term_t vsop_lon_Jupiter_1[] = 
{
    { 0.00078203446, 1.52377859742, 529.69096509460 }
};

static const vsop_series_t vsop_lon_Jupiter[] = 
{
    { 6, vsop_lon_Jupiter_0 },
    { 1, vsop_lon_Jupiter_1 }
};

static const vsop_term_t vsop_rad_Jupiter_0[] = 
{
    { 5.20887429326, 0.00000000000, 0.00000000000 },
    { 0.25209327119, 3.49108639871, 529.69096509460 },
    { 0.00610599976, 3.84115365948, 1059.38193018920 },
    { 0.00282029458, 2.57419881293, 632.78373931320 },
    { 0.00187647346, 2.07590383214, 522.57741809380 },
    { 0.00086792905, 0.71001145545, 419.48464387520 },
    { 0.00072062974, 0.21465724607, 536.80451209540 },
    { 0.00065517248, 5.97995884790, 316.39186965660 },
    { 0.00029134542, 1.67759379655, 103.09277421860 },
    { 0.00030135335, 2.16132003734, 949.17560896980 },
    { 0.00023453271, 3.54023522184, 735.87651353180 },
    { 0.00022283743, 4.19362594399, 1589.07289528380 },
    { 0.00023947298, 0.27458037480, 7.11354700080 },
    { 0.00013032614, 2.96042965363, 1162.47470440780 },
    { 0.00009703360, 1.90669633585, 206.18554843720 },
    { 0.00012749023, 2.71550286592, 1052.26838318840 }
};

static const vsop_term_t vsop_rad_Jupiter_1[] = 
{
    { 0.01271801520, 2.64937512894, 529.69096509460 },
    { 0.00061661816, 3.00076460387, 1059.38193018920 },
    { 0.00053443713, 3.89717383175, 522.57741809380 },
    { 0.00031185171, 4.88276958012, 536.80451209540 },
    { 0.00041390269, 0.00000000000, 0.00000000000 }
};

static const vsop_series_t vsop_rad_Jupiter[] = 
{
    { 16, vsop_rad_Jupiter_0 },
    { 5, vsop_rad_Jupiter_1 }
};

;
static const vsop_term_t vsop_lat_Saturn_0[] = 
{
    { 0.87401354025, 0.00000000000, 0.00000000000 },
    { 0.11107659762, 3.96205090159, 213.29909543800 },
    { 0.01414150957, 4.58581516874, 7.11354700080 },
    { 0.00398379389, 0.52112032699, 206.18554843720 },
    { 0.00350769243, 3.30329907896, 426.59819087600 },
    { 0.00206816305, 0.24658372002, 103.09277421860 },
    { 0.00079271300, 3.84007056878, 220.41264243880 },
    { 0.00023990355, 4.66976924553, 110.20632121940 },
    { 0.00016573588, 0.43719228296, 419.48464387520 },
    { 0.00014906995, 5.76903183869, 316.39186965660 },
    { 0.00015820290, 0.93809155235, 632.78373931320 },
    { 0.00014609559, 1.56518472000, 3.93215326310 },
    { 0.00013160301, 4.44891291899, 14.22709400160 },
    { 0.00015053543, 2.71669915667, 639.89728631400 },
    { 0.00013005299, 5.98119023644, 11.04570026390 },
    { 0.00010725067, 3.12939523827, 202.25339517410 },
    { 0.00005863206, 0.23656938524, 529.69096509460 },
    { 0.00005227757, 4.20783365759, 3.18139373770 },
    { 0.00006126317, 1.76328667907, 277.03499374140 },
    { 0.00005019687, 3.17787728405, 433.71173787680 },
    { 0.00004592550, 0.61977744975, 199.07200143640 },
    { 0.00004005867, 2.24479718502, 63.73589830340 },
    { 0.00002953796, 0.98280366998, 95.97922721780 },
    { 0.00003873670, 3.22283226966, 138.51749687070 },
    { 0.00002461186, 2.03163875071, 735.87651353180 },
    { 0.00003269484, 0.77492638211, 949.17560896980 },
    { 0.00001758145, 3.26580109940, 522.57741809380 },
    { 0.00001640172, 5.50504453050, 846.08283475120 },
    { 0.00001391327, 4.02333150505, 323.50541665740 },
    { 0.00001580648, 4.37265307169, 309.27832265580 },
    { 0.00001123498, 2.83726798446, 415.55249061210 },
    { 0.00001017275, 3.71700135395, 227.52618943960 },
    { 0.00000848642, 3.19150170830, 209.36694217490 }
};

static const vsop_term_t vsop_lat_Saturn_1[] = 
{
    { 213.29909521690, 0.00000000000, 0.00000000000 },
    { 0.01297370862, 1.82834923978, 213.29909543800 },
    { 0.00564345393, 2.88499717272, 7.11354700080 },
    { 0.00093734369, 1.06311793502, 426.59819087600 },
    { 0.00107674962, 2.27769131009, 206.18554843720 },
    { 0.00040244455, 2.04108104671, 220.41264243880 },
    { 0.00019941774, 1.27954390470, 103.09277421860 },
    { 0.00010511678, 2.74880342130, 14.22709400160 },
    { 0.00006416106, 0.38238295041, 639.89728631400 },
    { 0.00004848994, 2.43037610229, 419.48464387520 },
    { 0.00004056892, 2.92133209468, 110.20632121940 },
    { 0.00003768635, 3.64965330780, 3.93215326310 }
};

static const vsop_term_t vsop_lat_Saturn_2[] = 
{
    { 0.00116441330, 1.17988132879, 7.11354700080 },
    { 0.00091841837, 0.07325195840, 213.29909543800 },
    { 0.00036661728, 0.00000000000, 0.00000000000 },
    { 0.00015274496, 4.06493179167, 206.18554843720 }
};

static const vsop_series_t vsop_lat_Saturn[] = 
{
    { 33, vsop_lat_Saturn_0 },
    { 12, vsop_lat_Saturn_1 },
    { 4, vsop_lat_Saturn_2 }
};

static const vsop_term_t vsop_lon_Saturn_0[] = 
{
    { 0.04330678039, 3.60284428399, 213.29909543800 },
    { 0.00240348302, 2.85238489373, 426.59819087600 },
    { 0.00084745939, 0.00000000000, 0.00000000000 },
    { 0.00030863357, 3.48441504555, 220.41264243880 },
    { 0.00034116062, 0.57297307557, 206.18554843720 },
    { 0.00014734070, 2.11846596715, 639.89728631400 },
    { 0.00009916667, 5.79003188904, 419.48464387520 },
    { 0.00006993564, 4.73604689720, 7.11354700080 },
    { 0.00004807588, 5.43305312061, 316.39186965660 }
};

static const vsop_term_t vsop_lon_Saturn_1[] = 
{
    { 0.00198927992, 4.93901017903, 213.29909543800 },
    { 0.00036947916, 3.14159265359, 0.00000000000 },
    { 0.00017966989, 0.51979431110, 426.59819087600 }
};

static const vsop_series_t vsop_lon_Saturn[] = 
{
    { 9, vsop_lon_Saturn_0 },
    { 3, vsop_lon_Saturn_1 }
};

static const vsop_term_t vsop_rad_Saturn_0[] = 
{
    { 9.55758135486, 0.00000000000, 0.00000000000 },
    { 0.52921382865, 2.39226219573, 213.29909543800 },
    { 0.01873679867, 5.23549604660, 206.18554843720 },
    { 0.01464663929, 1.64763042902, 426.59819087600 },
    { 0.00821891141, 5.93520042303, 316.39186965660 },
    { 0.00547506923, 5.01532618980, 103.09277421860 },
    { 0.00371684650, 2.27114821115, 220.41264243880 },
    { 0.00361778765, 3.13904301847, 7.11354700080 },
    { 0.00140617506, 5.70406606781, 632.78373931320 },
    { 0.00108974848, 3.29313390175, 110.20632121940 },
    { 0.00069006962, 5.94099540992, 419.48464387520 },
    { 0.00061053367, 0.94037691801, 639.89728631400 },
    { 0.00048913294, 1.55733638681, 202.25339517410 },
    { 0.00034143772, 0.19519102597, 277.03499374140 },
    { 0.00032401773, 5.47084567016, 949.17560896980 },
    { 0.00020936596, 0.46349251129, 735.87651353180 }
};

static const vsop_term_t vsop_rad_Saturn_1[] = 
{
    { 0.06182981340, 0.25843511480, 213.29909543800 },
    { 0.00506577242, 0.71114625261, 206.18554843720 },
    { 0.00341394029, 5.79635741658, 426.59819087600 },
    { 0.00188491195, 0.47215589652, 220.41264243880 },
    { 0.00186261486, 3.14159265359, 0.00000000000 },
    { 0.00143891146, 1.40744822888, 7.11354700080 }
};

static const vsop_term_t vsop_rad_Saturn_2[] = 
{
    { 0.00436902572, 4.78671677509, 213.29909543800 }
};

static const vsop_series_t vsop_rad_Saturn[] = 
{
    { 16, vsop_rad_Saturn_0 },
    { 6, vsop_rad_Saturn_1 },
    { 1, vsop_rad_Saturn_2 }
};

;
static const vsop_term_t vsop_lat_Uranus_0[] = 
{
    { 5.48129294297, 0.00000000000, 0.00000000000 },
    { 0.09260408234, 0.89106421507, 74.78159856730 },
    { 0.01504247898, 3.62719260920, 1.48447270830 },
    { 0.00365981674, 1.89962179044, 73.29712585900 },
    { 0.00272328168, 3.35823706307, 149.56319713460 },
    { 0.00070328461, 5.39254450063, 63.73589830340 },
    { 0.00068892678, 6.09292483287, 76.26607127560 },
    { 0.00061998615, 2.26952066061, 2.96894541660 },
    { 0.00061950719, 2.85098872691, 11.04570026390 },
    { 0.00026468770, 3.14152083966, 71.81265315070 },
    { 0.00025710476, 6.11379840493, 454.90936652730 },
    { 0.00021078850, 4.36059339067, 148.07872442630 },
    { 0.00017818647, 1.74436930289, 36.64856292950 },
    { 0.00014613507, 4.73732166022, 3.93215326310 },
    { 0.00011162509, 5.82681796350, 224.34479570190 },
    { 0.00010997910, 0.48865004018, 138.51749687070 },
    { 0.00009527478, 2.95516862826, 35.16409022120 },
    { 0.00007545601, 5.23626582400, 109.94568878850 },
    { 0.00004220241, 3.23328220918, 70.84944530420 },
    { 0.00004051900, 2.27755017300, 151.04766984290 },
    { 0.00003354596, 1.06549007380, 4.45341812490 },
    { 0.00002926718, 4.62903718891, 9.56122755560 },
    { 0.00003490340, 5.48306144511, 146.59425171800 },
    { 0.00003144069, 4.75199570434, 77.75054398390 },
    { 0.00002922333, 5.35235361027, 85.82729883120 },
    { 0.00002272788, 4.36600400036, 70.32818044240 },
    { 0.00002051219, 1.51773566586, 0.11187458460 },
    { 0.00002148602, 0.60745949945, 38.13303563780 },
    { 0.00001991643, 4.92437588682, 277.03499374140 },
    { 0.00001376226, 2.04283539351, 65.22037101170 },
    { 0.00001666902, 3.62744066769, 380.12776796000 },
    { 0.00001284107, 3.11347961505, 202.25339517410 },
    { 0.00001150429, 0.93343589092, 3.18139373770 },
    { 0.00001533221, 2.58594681212, 52.69019803950 },
    { 0.00001281604, 0.54271272721, 222.86032299360 },
    { 0.00001372139, 4.19641530878, 111.43016149680 },
    { 0.00001221029, 0.19900650030, 108.46121608020 },
    { 0.00000946181, 1.19253165736, 127.47179660680 },
    { 0.00001150989, 4.17898916639, 33.67961751290 }
};

static const vsop_term_t vsop_lat_Uranus_1[] = 
{
    { 74.78159860910, 0.00000000000, 0.00000000000 },
    { 0.00154332863, 5.24158770553, 74.78159856730 },
    { 0.00024456474, 1.71260334156, 1.48447270830 },
    { 0.00009258442, 0.42829732350, 11.04570026390 },
    { 0.00008265977, 1.50218091379, 63.73589830340 },
    { 0.00009150160, 1.41213765216, 149.56319713460 }
};

static const vsop_series_t vsop_lat_Uranus[] = 
{
    { 39, vsop_lat_Uranus_0 },
    { 6, vsop_lat_Uranus_1 }
};

static const vsop_term_t vsop_lon_Uranus_0[] = 
{
    { 0.01346277648, 2.61877810547, 74.78159856730 },
    { 0.00062341400, 5.08111189648, 149.56319713460 },
    { 0.00061601196, 3.14159265359, 0.00000000000 },
    { 0.00009963722, 1.61603805646, 76.26607127560 },
    { 0.00009926160, 0.57630380333, 73.29712585900 }
};

static const vsop_term_t vsop_lon_Uranus_1[] = 
{
    { 0.00034101978, 0.01321929936, 74.78159856730 }
};

static const vsop_series_t vsop_lon_Uranus[] = 
{
    { 5, vsop_lon_Uranus_0 },
    { 1, vsop_lon_Uranus_1 }
};

static const vsop_term_t vsop_rad_Uranus_0[] = 
{
    { 19.21264847206, 0.00000000000, 0.00000000000 },
    { 0.88784984413, 5.60377527014, 74.78159856730 },
    { 0.03440836062, 0.32836099706, 73.29712585900 },
    { 0.02055653860, 1.78295159330, 149.56319713460 },
    { 0.00649322410, 4.52247285911, 76.26607127560 },
    { 0.00602247865, 3.86003823674, 63.73589830340 },
    { 0.00496404167, 1.40139935333, 454.90936652730 },
    { 0.00338525369, 1.58002770318, 138.51749687070 },
    { 0.00243509114, 1.57086606044, 71.81265315070 },
    { 0.00190522303, 1.99809394714, 1.48447270830 },
    { 0.00161858838, 2.79137786799, 148.07872442630 },
    { 0.00143706183, 1.38368544947, 11.04570026390 },
    { 0.00093192405, 0.17437220467, 36.64856292950 },
    { 0.00071424548, 4.24509236074, 224.34479570190 },
    { 0.00089806014, 3.66105364565, 109.94568878850 },
    { 0.00039009723, 1.66971401684, 70.84944530420 },
    { 0.00046677296, 1.39976401694, 35.16409022120 },
    { 0.00039025624, 3.36234773834, 277.03499374140 },
    { 0.00036755274, 3.88649278513, 146.59425171800 },
    { 0.00030348723, 0.70100838798, 151.04766984290 },
    { 0.00029156413, 3.18056336700, 77.75054398390 }
};

static const vsop_term_t vsop_rad_Uranus_1[] = 
{
    { 0.01479896629, 3.67205697578, 74.78159856730 }
};

static const vsop_series_t vsop_rad_Uranus[] = 
{
    { 21, vsop_rad_Uranus_0 },
    { 1, vsop_rad_Uranus_1 }
};

;
static const vsop_term_t vsop_lat_Neptune_0[] = 
{
    { 5.31188633046, 0.00000000000, 0.00000000000 },
    { 0.01798475530, 2.90101273890, 38.13303563780 },
    { 0.01019727652, 0.48580922867, 1.48447270830 },
    { 0.00124531845, 4.83008090676, 36.64856292950 },
    { 0.00042064466, 5.41054993053, 2.96894541660 },
    { 0.00037714584, 6.09221808686, 35.16409022120 },
    { 0.00033784738, 1.24488874087, 76.26607127560 },
    { 0.00016482741, 0.00007727998, 491.55792945680 },
    { 0.00009198584, 4.93747051954, 39.61750834610 },
    { 0.00008994250, 0.27462171806, 175.16605980020 }
};

static const vsop_term_t vsop_lat_Neptune_1[] = 
{
    { 38.13303563957, 0.00000000000, 0.00000000000 },
    { 0.00016604172, 4.86323329249, 1.48447270830 },
    { 0.00015744045, 2.27887427527, 38.13303563780 }
};

static const vsop_series_t vsop_lat_Neptune[] = 
{
    { 10, vsop_lat_Neptune_0 },
    { 3, vsop_lat_Neptune_1 }
};

static const vsop_term_t vsop_lon_Neptune_0[] = 
{
    { 0.03088622933, 1.44104372644, 38.13303563780 },
    { 0.00027780087, 5.91271884599, 76.26607127560 },
    { 0.00027623609, 0.00000000000, 0.00000000000 },
    { 0.00015355489, 2.52123799551, 36.64856292950 },
    { 0.00015448133, 3.50877079215, 39.61750834610 }
};

static const vsop_series_t vsop_lon_Neptune[] = 
{
    { 5, vsop_lon_Neptune_0 }
};

static const vsop_term_t vsop_rad_Neptune_0[] = 
{
    { 30.07013205828, 0.00000000000, 0.00000000000 },
    { 0.27062259632, 1.32999459377, 38.13303563780 },
    { 0.01691764014, 3.25186135653, 36.64856292950 },
    { 0.00807830553, 5.18592878704, 1.48447270830 },
    { 0.00537760510, 4.52113935896, 35.16409022120 },
    { 0.00495725141, 1.57105641650, 491.55792945680 },
    { 0.00274571975, 1.84552258866, 175.16605980020 }
};

static const vsop_series_t vsop_rad_Neptune[] = 
{
    { 7, vsop_rad_Neptune_0 }
};

;

#define VSOPFORMULA(x)    { ARRAYSIZE(x), x }

static const vsop_model_t vsop[] = 
{
    { { VSOPFORMULA(vsop_lat_Mercury),  VSOPFORMULA(vsop_lon_Mercury),  VSOPFORMULA(vsop_rad_Mercury) } },
    { { VSOPFORMULA(vsop_lat_Venus),    VSOPFORMULA(vsop_lon_Venus),    VSOPFORMULA(vsop_rad_Venus)   } },
    { { VSOPFORMULA(vsop_lat_Earth),    VSOPFORMULA(vsop_lon_Earth),    VSOPFORMULA(vsop_rad_Earth)   } },
    { { VSOPFORMULA(vsop_lat_Mars),     VSOPFORMULA(vsop_lon_Mars),     VSOPFORMULA(vsop_rad_Mars)    } },
    { { VSOPFORMULA(vsop_lat_Jupiter),  VSOPFORMULA(vsop_lon_Jupiter),  VSOPFORMULA(vsop_rad_Jupiter) } },
    { { VSOPFORMULA(vsop_lat_Saturn),   VSOPFORMULA(vsop_lon_Saturn),   VSOPFORMULA(vsop_rad_Saturn)  } },
    { { VSOPFORMULA(vsop_lat_Uranus),   VSOPFORMULA(vsop_lon_Uranus),   VSOPFORMULA(vsop_rad_Uranus)  } },
    { { VSOPFORMULA(vsop_lat_Neptune),  VSOPFORMULA(vsop_lon_Neptune),  VSOPFORMULA(vsop_rad_Neptune) } }
};

#define CalcEarth(time)     CalcVsop(&vsop[BODY_EARTH], (time))

static astro_vector_t CalcVsop(const vsop_model_t *model, astro_time_t time)
{
    int k, s, i;
    double t = time.tt / 365250;    /* millennia since 2000 */
    double sphere[3];
    double r_coslat;
    double eclip[3];
    astro_vector_t vector;

    for (k=0; k < 3; ++k)
    {
        double tpower = 1.0;
        const vsop_formula_t *formula = &model->formula[k];
        sphere[k] = 0.0;
        for (s=0; s < formula->nseries; ++s)
        {
            double sum = 0.0;
            const vsop_series_t *series = &formula->series[s];
            for (i=0; i < series->nterms; ++i)
            {
                const vsop_term_t *term = &series->term[i];
                sum += term->amplitude * cos(term->phase + (t * term->frequency));
            }
            sphere[k] += tpower * sum;
            tpower *= t;
        }
    }

    /* Convert ecliptic spherical coordinates to ecliptic Cartesian coordinates. */
    r_coslat = sphere[2] * cos(sphere[1]);
    eclip[0] = r_coslat * cos(sphere[0]);
    eclip[1] = r_coslat * sin(sphere[0]);
    eclip[2] = sphere[2] * sin(sphere[1]);

    /* Convert ecliptic Cartesian coordinates to equatorial Cartesian coordinates. */
    vector.status = ASTRO_SUCCESS;
    vector.x = eclip[0] + 0.000000440360*eclip[1] - 0.000000190919*eclip[2];
    vector.y = -0.000000479966*eclip[0] + 0.917482137087*eclip[1] - 0.397776982902*eclip[2];
    vector.z = 0.397776982902*eclip[1] + 0.917482137087*eclip[2];
    vector.t = time;

    return vector;
}

astro_vector_t Astronomy_HelioVector(astro_body_t body, astro_time_t time)
{
    astro_vector_t vector;

    switch (body)
    {
    case BODY_SUN:
        vector.status = ASTRO_SUCCESS;
        vector.x = 0.0;
        vector.y = 0.0;
        vector.z = 0.0;
        vector.t = time;
        return vector;

    case BODY_MERCURY:
    case BODY_VENUS:
    case BODY_EARTH:
    case BODY_MARS:
    case BODY_JUPITER:
    case BODY_SATURN:
    case BODY_URANUS:
    case BODY_NEPTUNE:
        return CalcVsop(&vsop[body], time);

    case BODY_MOON:     /* FIXFIXFIX: GeoMoon not yet implemented. */
    case BODY_PLUTO:    /* FIXFIXFIX: Chebyshev models not yet implemented */
    default:
        vector.status = ASTRO_INVALID_BODY;
        vector.x = vector.y = vector.z = 1.0e+99;   /* Invalid coordinates */
        vector.t = time;
        return vector;
    }
}

static astro_vector_t RawGeoVector(astro_body_t body, astro_time_t time)
{
    astro_vector_t earth;
    astro_vector_t helio;
    astro_vector_t geo;

    earth = CalcEarth(time);
    if (earth.status != ASTRO_SUCCESS)
        return earth;

    helio = Astronomy_HelioVector(body, time);
    if (helio.status != ASTRO_SUCCESS)
        return helio;

    geo.status = ASTRO_SUCCESS;
    geo.x = helio.x - earth.x;
    geo.y = helio.y - earth.y;
    geo.z = helio.z - earth.z;
    geo.t = time;
    return geo;
}

astro_vector_t Astronomy_GeoVector(astro_body_t body, astro_time_t time)
{
    astro_vector_t vector;
    astro_time_t ltime;
    astro_time_t ltime2;
    double dt;
    int iter;

    switch (body)
    {
    case BODY_EARTH:
        /* The Earth's geocentric coordinates are always (0,0,0). */
        vector.status = ASTRO_SUCCESS;
        vector.x = 0.0;
        vector.y = 0.0;
        vector.z = 0.0;
        break;

    case BODY_SUN:
        /* The Sun's heliocentric coordinates are always (0,0,0). No need for light travel correction. */
        vector = CalcEarth(time);
        vector.x *= -1.0;
        vector.y *= -1.0;
        vector.z *= -1.0;
        break;

    case BODY_MOON:
        vector.status = ASTRO_INVALID_BODY;     /* FIXFIXFIX: GeoMoon not yet implemented. */
        vector.x = vector.y = vector.z = 1.0e+99;
        break;

    default:
        /* For all other bodies, apply light travel time correction. */
        ltime = time;
        for (iter=0; iter < 10; ++iter)
        {            
            vector = RawGeoVector(body, ltime);
            if (vector.status != ASTRO_SUCCESS)
                return vector;

            ltime2 = Astronomy_AddDays(time, -Astronomy_VectorLength(vector) / C_AUDAY);
            dt = fabs(ltime2.tt - ltime.tt);
            if (dt < 1.0e-9)
                goto finished;  /* Tricky: ensures we patch 'vector.t' with current time, not ante-dated time. */

            ltime = ltime2;
        }
        vector.status = ASTRO_NO_CONVERGE;  /* light travel time solver did not converge */
        break;
    }

finished:
    vector.t = time;
    return vector;
}

astro_sky_t Astronomy_SkyPos(astro_vector_t gc_vector, astro_observer_t observer)
{
    astro_sky_t sky;
    double gc_observer[3];
    double j2000_vector[3];
    double ofdate_vector[3];
    double pos7[3];

    geo_pos(gc_vector.t, observer, gc_observer);

    j2000_vector[0] = gc_vector.x - gc_observer[0];
    j2000_vector[1] = gc_vector.y - gc_observer[1];
    j2000_vector[2] = gc_vector.z - gc_observer[2];

    sky.j2000 = vector2radec(j2000_vector);
    precession(0.0, j2000_vector, gc_vector.t.tt, pos7);
    nutation(gc_vector.t, 0.0, pos7, ofdate_vector);
    sky.ofdate = vector2radec(ofdate_vector);
    sky.t = gc_vector.t;

    return sky;
}

astro_horizon_t Astronomy_Horizon(
    astro_time_t time, astro_observer_t observer, double ra, double dec, astro_refraction_t refraction)
{
    astro_horizon_t hor;
    double uze[3], une[3], uwe[3];
    double uz[3], un[3], uw[3];
    double p[3], pz, pn, pw, proj;
    double az, zd;

    double sinlat = sin(observer.latitude * DEG2RAD);
    double coslat = cos(observer.latitude * DEG2RAD);
    double sinlon = sin(observer.longitude * DEG2RAD);
    double coslon = cos(observer.longitude * DEG2RAD);
    double sindc = sin(dec * DEG2RAD);
    double cosdc = cos(dec * DEG2RAD);
    double sinra = sin(ra * 15 * DEG2RAD);
    double cosra = cos(ra * 15 * DEG2RAD);

    uze[0] = coslat * coslon;
    uze[1] = coslat * sinlon;
    uze[2] = sinlat;

    une[0] = -sinlat * coslon;
    une[1] = -sinlat * sinlon;
    une[2] = coslat;

    uwe[0] = sinlon;
    uwe[1] = -coslon;
    uwe[2] = 0.0;

    ter2cel(time, uze, uz);
    ter2cel(time, une, un);
    ter2cel(time, uwe, uw);

    p[0] = cosdc * cosra;
    p[1] = cosdc * sinra;
    p[2] = sindc;

    pz = p[0]*uz[0] + p[1]*uz[1] + p[2]*uz[2];
    pn = p[0]*un[0] + p[1]*un[1] + p[2]*un[2];
    pw = p[0]*uw[0] + p[1]*uw[1] + p[2]*uw[2];

    proj = sqrt(pn*pn + pw*pw);
    az = 0.0;
    if (proj > 0.0) {
        az = -atan2(pw, pn) * RAD2DEG;
        if (az < 0) 
            az += 360;
        if (az >= 360) 
            az -= 360;
    }
    zd = atan2(proj, pz) * RAD2DEG;
    hor.ra = ra;
    hor.dec = dec;

    if (refraction == REFRACTION_NORMAL || refraction == REFRACTION_JPLHOR) 
    {
        double zd0, refr, hd;
        int j;

        zd0 = zd;

        // http://extras.springer.com/1999/978-1-4471-0555-8/chap4/horizons/horizons.pdf
        // JPL Horizons says it uses refraction algorithm from 
        // Meeus "Astronomical Algorithms", 1991, p. 101-102.
        // I found the following Go implementation:
        // https://github.com/soniakeys/meeus/blob/master/v3/refraction/refract.go
        // This is a translation from the function "Saemundsson" there.
        // I found experimentally that JPL Horizons clamps the angle to 1 degree below the horizon.
        // This is important because the 'refr' formula below goes crazy near hd = -5.11.
        hd = 90.0 - zd;
        if (hd < -1.0)
            hd = -1.0;

        refr = (1.02 / tan((hd+10.3/(hd+5.11))*DEG2RAD)) / 60.0;

        if (refraction == REFRACTION_NORMAL && zd > 91.0) 
        {
            // In "normal" mode we gradually reduce refraction toward the nadir
            // so that we never get an altitude angle less than -90 degrees.
            // When horizon angle is -1 degrees, zd = 91, and the factor is exactly 1.
            // As zd approaches 180 (the nadir), the fraction approaches 0 linearly.
            refr *= (180.0 - zd) / 89.0;
        }

        zd -= refr;

        if (refr > 0.0 && zd > 3.0e-4) 
        {
            double sinzd = sin(zd * DEG2RAD);
            double coszd = cos(zd * DEG2RAD);
            double sinzd0 = sin(zd0 * DEG2RAD);
            double coszd0 = cos(zd0 * DEG2RAD);
            double pr[3];

            for (j=0; j<3; ++j)
                pr[j] = ((p[j] - coszd0 * uz[j]) / sinzd0)*sinzd + uz[j]*coszd;

            proj = sqrt(pr[0]*pr[0] + pr[1]*pr[1]);
            if (proj > 0) 
            {
                hor.ra = atan2(pr[1], pr[0]) * RAD2DEG / 15;
                if (hor.ra < 0)
                    hor.ra += 24;
                if (hor.ra >= 24)
                    hor.ra -= 24;
            } 
            else 
            {
                hor.ra = 0;
            }
            hor.dec = atan2(pr[2], proj) * RAD2DEG;
        }
    }

    hor.azimuth = az;
    hor.altitude = 90.0 - zd;
    return hor;
}


#ifdef __cplusplus
}
#endif
