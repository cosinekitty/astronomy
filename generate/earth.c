/*
    earth.c    
    Special-purpose logic for calculating Earth-Moon Barycenter.
    Based on sun_eph() from NOVAS C 3.1 file solsys3.c.
*/

#include <math.h>
#include "novas.h"
#include "eph_manager.h"
#include "earth.h"

const double EarthMoonMassRatio = 81.30056;

static void sun_eph (double jd,
              double *ra, double *dec, double *dis)
/*
------------------------------------------------------------------------

   PURPOSE:
      To compute equatorial spherical coordinates of Sun referred to
      the mean equator and equinox of date.

   REFERENCES:
      Bretagnon, P. and Simon, J.L. (1986).  Planetary Programs and
         Tables from -4000 to + 2800. (Richmond, VA: Willmann-Bell).
      Kaplan, G.H. (2005). US Naval Observatory Circular 179.

   INPUT
   ARGUMENTS:
      jd (double)
         Julian date on TDT or ET time scale.

   OUTPUT
   ARGUMENTS:
      ra (double)
         Right ascension referred to mean equator and equinox of date
         (hours).
      dec (double)
         Declination referred to mean equator and equinox of date
         (degrees).
      dis (double)
         Geocentric distance (AU).

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      T0, TWOPI, ASEC2RAD

   FUNCTIONS
   CALLED:
      sin           math.h
      cos           math.h
      asin          math.h
      atan2         math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/08-94/JAB (USNO/AA)
      V1.1/05-96/JAB (USNO/AA): Compute mean coordinates instead of
                                apparent.
      V1.2/01-07/JAB (USNO/AA): Use 'ASEC2RAD' instead of 'RAD2SEC'.
      V1.3/04-09/JAB (USNO/AA): Update the equation for mean
                                obliquity of the ecliptic, and correct
                                longitude based on a linear fit to DE405
                                in the interval 1900-2100 (see notes).

   NOTES:
      1. Quoted accuracy is 2.0 + 0.03 * T^2 arcsec, where T is
      measured in units of 1000 years from J2000.0.  See reference.
      2. The obliquity equation is updated to equation 5.12 of the
      second reference.
      3. The linear fit to DE405 primarily corrects for the
      difference between "old" (Lieske) and "new" (IAU 2006)
      precession.  The difference, new - old, is -0.3004 arcsec/cy.

------------------------------------------------------------------------
*/
{
   short int i;

   double sum_lon = 0.0;
   double sum_r = 0.0;
   const double factor = 1.0e-07;
   double u, arg, lon, t, emean, sin_lon;

   struct sun_con
   {
   double l;
   double r;
   double alpha;
   double nu;
   };

   static const struct sun_con con[50] =
      {{403406.0,      0.0, 4.721964,     1.621043},
       {195207.0, -97597.0, 5.937458, 62830.348067},
       {119433.0, -59715.0, 1.115589, 62830.821524},
       {112392.0, -56188.0, 5.781616, 62829.634302},
       {  3891.0,  -1556.0, 5.5474  , 125660.5691 },
       {  2819.0,  -1126.0, 1.5120  , 125660.9845 },
       {  1721.0,   -861.0, 4.1897  ,  62832.4766 },
       {     0.0,    941.0, 1.163   ,      0.813  },
       {   660.0,   -264.0, 5.415   , 125659.310  },
       {   350.0,   -163.0, 4.315   ,  57533.850  },
       {   334.0,      0.0, 4.553   ,    -33.931  },
       {   314.0,    309.0, 5.198   , 777137.715  },
       {   268.0,   -158.0, 5.989   ,  78604.191  },
       {   242.0,      0.0, 2.911   ,      5.412  },
       {   234.0,    -54.0, 1.423   ,  39302.098  },
       {   158.0,      0.0, 0.061   ,    -34.861  },
       {   132.0,    -93.0, 2.317   , 115067.698  },
       {   129.0,    -20.0, 3.193   ,  15774.337  },
       {   114.0,      0.0, 2.828   ,   5296.670  },
       {    99.0,    -47.0, 0.52    ,  58849.27   },
       {    93.0,      0.0, 4.65    ,   5296.11   },
       {    86.0,      0.0, 4.35    ,  -3980.70   },
       {    78.0,    -33.0, 2.75    ,  52237.69   },
       {    72.0,    -32.0, 4.50    ,  55076.47   },
       {    68.0,      0.0, 3.23    ,    261.08   },
       {    64.0,    -10.0, 1.22    ,  15773.85   },
       {    46.0,    -16.0, 0.14    ,  188491.03  },
       {    38.0,      0.0, 3.44    ,   -7756.55  },
       {    37.0,      0.0, 4.37    ,     264.89  },
       {    32.0,    -24.0, 1.14    ,  117906.27  },
       {    29.0,    -13.0, 2.84    ,   55075.75  },
       {    28.0,      0.0, 5.96    ,   -7961.39  },
       {    27.0,     -9.0, 5.09    ,  188489.81  },
       {    27.0,      0.0, 1.72    ,    2132.19  },
       {    25.0,    -17.0, 2.56    ,  109771.03  },
       {    24.0,    -11.0, 1.92    ,   54868.56  },
       {    21.0,      0.0, 0.09    ,   25443.93  },
       {    21.0,     31.0, 5.98    ,  -55731.43  },
       {    20.0,    -10.0, 4.03    ,   60697.74  },
       {    18.0,      0.0, 4.27    ,    2132.79  },
       {    17.0,    -12.0, 0.79    ,  109771.63  },
       {    14.0,      0.0, 4.24    ,   -7752.82  },
       {    13.0,     -5.0, 2.01    ,  188491.91  },
       {    13.0,      0.0, 2.65    ,     207.81  },
       {    13.0,      0.0, 4.98    ,   29424.63  },
       {    12.0,      0.0, 0.93    ,      -7.99  },
       {    10.0,      0.0, 2.21    ,   46941.14  },
       {    10.0,      0.0, 3.59    ,     -68.29  },
       {    10.0,      0.0, 1.50    ,   21463.25  },
       {    10.0,     -9.0, 2.55    ,  157208.40  }};

/*
   Define the time units 'u', measured in units of 10000 Julian years
   from J2000.0, and 't', measured in Julian centuries from J2000.0.
*/

   u = (jd - T0) / 3652500.0;
   t = u * 100.0;

/*
   Compute longitude and distance terms from the series.
*/

   for (i = 0; i < 50; i++)
   {
      arg = con[i].alpha + con[i].nu * u;
      sum_lon += con[i].l * sin (arg);
      sum_r += con[i].r * cos (arg);
   }

/*
   Compute longitude, latitude, and distance referred to mean equinox
   and ecliptic of date.  Apply correction to longitude based on a
   linear fit to DE405 in the interval 1900-2100.
*/

   lon = 4.9353929 + 62833.1961680 * u + factor * sum_lon;
   lon += ((-0.1371679461 - 0.2918293271 * t) * ASEC2RAD);

   lon = fmod (lon, TWOPI);
   if (lon < 0.0)
      lon += TWOPI;

   *dis = 1.0001026 + factor * sum_r;

/*
   Compute mean obliquity of the ecliptic.
*/

   emean = (84381.406 + (-46.836769 +
      (-0.0001831 + 0.00200340 * t) * t) * t)  * ASEC2RAD;

/*
   Compute equatorial spherical coordinates referred to the mean equator
   and equinox of date.
*/

   sin_lon = sin (lon);
   *ra = atan2 ((cos (emean) * sin_lon), cos (lon)) * RAD2DEG;
   *ra = fmod (*ra, 360.0);
   if (*ra < 0.0)
      *ra += 360.0;
   *ra = *ra / 15.0;

   *dec = asin (sin (emean) * sin_lon) * RAD2DEG;
}

void CalcEarth(double jd, double pos[3])
{
    double ra, dec, dist;
    double date_pos[3];
    int error;

    sun_eph(jd, &ra, &dec, &dist);
    radec2vector(ra, dec, dist, date_pos);
    error = precession(jd, date_pos, T0, pos);
    if (error)
    {
        fprintf(stderr, "CalcEarth(FATAL): precession returned %d\n", error);
        exit(1);
    }

    /* Convert Sun seen from Earth, to Earth seen from Sun. */
    pos[0] *= -1.0;
    pos[1] *= -1.0;
    pos[2] *= -1.0;
}

int NovasEarth(double jd, double pos[3])
{
    int error;
    int k;
    double jed[2];
    double emb_pos[3], sun_pos[3], moon_pos[3], vel[3];

    jed[0] = jd;
    jed[1] = 0.0;

    error = state(jed, 2, emb_pos, vel);    /* Earth/Moon Barycenter with respect to Solar System Barycenter */
    if (error) return error;

    error = state(jed, 10, sun_pos, vel);   /* Position of the Sun with respect to the Solar System Barycenter */
    if (error) return error;

    error = state(jed, 9, moon_pos, vel);   /* Geocentric Moon */
    if (error) return error;

    /* Calculate heliocentric Earth. */
    for (k=0; k<3; ++k)
        pos[k] = (emb_pos[k] - sun_pos[k]) - moon_pos[k]/(1.0 + EarthMoonMassRatio);

    return 0;
}
