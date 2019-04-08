/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1
 
  checkout-stars-full.c: Checkout program for use with solsys1 in full-accuracy mode 
 
  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC 
  http://www.usno.navy.mil/USNO/astronomical-applications
*/

#include <stdio.h>
#include <stdlib.h>
#include "eph_manager.h"
#include "novas.h"

#define N_STARS 3
#define N_TIMES 4

int main (void)
{
/*
   Main function to check out many parts of NOVAS-C by calling
   function 'topo_star' with version 1 of function 'solarsystem'.

   For use with NOVAS-C Version 3.1.
*/

   short int error = 0;
   short int accuracy = 0;
   short int i, j, de_num;

/*
   'deltat' is the difference in time scales, TT - UT1.

    The array 'tjd' contains four selected Julian dates at which the
    star positions will be evaluated.
*/

   double deltat = 60.0;
   double tjd[N_TIMES] = {2450203.5, 2450203.5, 2450417.5, 2450300.5};
   double jd_beg, jd_end, ra, dec;

/*
   Hipparcos (ICRS) catalog data for three selected stars.
*/

   cat_entry stars[N_STARS] = {
      {"POLARIS", "HIP",   0,  2.530301028,  89.264109444,
               44.22, -11.75,  7.56, -17.4},
      {"Delta ORI", "HIP", 1,  5.533444639,  -0.299091944,
                1.67,   0.56,  3.56,  16.0},
      {"Theta CAR", "HIP", 2, 10.715944806, -64.394450000,
               -18.87, 12.06,  7.43,  24.0}};

/*
   The observer's terrestrial coordinates (latitude, longitude, height).
*/

   on_surface geo_loc = {45.0, -75.0, 0.0, 10.0, 1010.0};

/*
   Open the JPL ephemeris file.
*/

   if ((error = ephem_open ("JPLEPH", &jd_beg,&jd_end,&de_num)) != 0)
   {
      printf ("Error %d from ephem_open\n", error);
      return (error);
   }
    else
   {
      printf ("JPL Ephemeris DE%d open. jd_beg = %10.2f  jd_end = %10.2f\n",
         de_num, jd_beg, jd_end);
      printf ("\n");
   }

/*
   Compute the topocentric places of the three stars at the four
   selected Julian dates.
*/

   for (i = 0; i < N_TIMES; i++)
   {
      for (j = 0; j < N_STARS; j++)
      {
         if ((error = topo_star (tjd[i],deltat,&stars[j],&geo_loc,
            accuracy, &ra,&dec)) != 0)
         {
            printf ("Error %d from topo_star. Star %d  Time %d\n",
               error, j, i);
            return (error);
         }
          else
         {
            printf ("JD = %f  Star = %s\n", tjd[i], stars[j].starname);
            printf ("RA = %12.9f  Dec = %12.8f\n", ra, dec);
            printf ("\n");
         }
      }
      printf ("\n");
   }

   ephem_close ();
   return (0);
}
