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

$ASTRO_C_VSOP(Mercury);
$ASTRO_C_VSOP(Venus);
$ASTRO_C_VSOP(Earth);
$ASTRO_C_VSOP(Mars);
$ASTRO_C_VSOP(Jupiter);
$ASTRO_C_VSOP(Saturn);
$ASTRO_C_VSOP(Uranus);
$ASTRO_C_VSOP(Neptune);

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

    vector.t = time;

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
                return vector;

            ltime = ltime2;
        }
        vector.status = ASTRO_NO_CONVERGE;  /* light travel time solver did not converge */
        break;
    }

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


#ifdef __cplusplus
}
#endif
