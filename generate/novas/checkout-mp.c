/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1
 
  checkout-mp.c: Checkout program for use with minor planet ephemerides
 
  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC 
  http://www.usno.navy.mil/USNO/astronomical-applications
*/

#include <stdio.h>
#include <stdlib.h>
#include "eph_manager.h"
#include "novas.h"

#define N_TIMES 4

int main (void)
{
/*
   Main function to check out many parts of NOVAS-C by computing the
   places of a minor planet and the Moon.  The USNO/AE98 minor planet
   ephemerides and DE405 are used along with version 1 of function
   'solarsystem'.

   For use with NOVAS-C Version 3.1.
*/

   short int error = 0;
   short int accuracy = 0;
   short int i, de_num;

/*
   'deltat' is the difference in time scales, TT - UT1.

    The array 'tjd' contains four selected Julian dates at which the
    minor planet positions will be evaluated.
*/

   double deltat = 60.0;
   double tjd_tt[N_TIMES] = {2450203.5, 2450203.5, 2450417.5,
      2450300.5};
   double jd_beg, jd_end, ra, dec, dis;

   on_surface geo_loc;

   cat_entry dummy_star;

   object pallas, moon;

/*
   The observer's terrestrial coordinates (latitude, longitude, height).
*/

    make_on_surface (45.0, -75.0, 0.0, 10.0, 1010.0, &geo_loc);

/*
   Set up the structure containing the body designation for Pallas.
   First, create a dummy star-catalog entry.
*/

   make_cat_entry ("dummy","   ",0,0.0,0.0,0.0,0.0,0.0,0.0,
      &dummy_star);

   if ((error = make_object (1,2,"Pallas",&dummy_star, &pallas)))
   {
      printf ("Error %d from set_body.\n", error);
      exit (1);
   }

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
   Compute the apparent place of Pallas at the four selected Julian
   dates.
*/

   printf ("\nPallas:\n\n");

   printf ("Apparent Places:\n");
   for (i = 0; i < N_TIMES; i++)
   {
      if ((error = app_planet (tjd_tt[i],&pallas,accuracy,
         &ra,&dec,&dis)) != 0)
      {
         printf ("Error %d from app_planet.\n", error);
         exit (1);
      }
       else
      {
         printf ("JD = %f  RA = %12.9f  Dec = %12.8f  Dis = %12.10f\n",
            tjd_tt[i], ra, dec, dis);
      }
   }

/*
   Compute the topocentric place of Pallas at the four selected Julian
   dates.
*/

   printf ("\nTopocentric Places:\n");
   for (i = 0; i < N_TIMES; i++)
   {
      if ((error = topo_planet (tjd_tt[i],&pallas,deltat,&geo_loc,
         accuracy,   &ra,&dec,&dis)) != 0)
      {
         printf ("Error %d from topo_planet.\n", error);
         return (error);
      }
       else
      {
         printf ("JD = %f  RA = %12.9f  Dec = %12.8f  Dis = %12.10f\n",
            tjd_tt[i], ra, dec, dis);
      }
   }

/*
   Make an 'object' structure for the Moon.
*/

   if ((error = make_object (0,11,"Moon",&dummy_star, &moon)) != 0)
   {
      printf ("Error %d from make_object\n", error);
      return (error);
   }

/*
   Compute the topocentric place of the Moon at the four selected Julian
   dates.
*/

   printf ("\nMoon:\n\n");

   printf ("Apparent Places:\n");
   for (i = 0; i < N_TIMES; i++)
   {
      if ((error = app_planet (tjd_tt[i],&moon,accuracy,
         &ra,&dec,&dis)) != 0)
      {
         printf ("Error %d from app_planet. Time %d\n", error, i);
         return (error);
      }
       else
      {
         printf ("JD = %f  RA = %12.9f  Dec = %12.8f  Dis = %12.10f\n",
            tjd_tt[i], ra, dec, dis);
      }
   }

   ephem_close ();
   return (0);
}
