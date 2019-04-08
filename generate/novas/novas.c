/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1

  novas.c: Main library

  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC
  http://www.usno.navy.mil/USNO/astronomical-applications
*/

#ifndef _NOVAS_
   #include "novas.h"
#endif

#include <math.h>

/*
   Global variables.

   'PSI_COR' and 'EPS_COR' are celestial pole offsets for high-
   precision applications.  See function 'cel_pole' for more details.
*/

static double PSI_COR = 0.0;
static double EPS_COR = 0.0;



/********app_star */

short int app_star (double jd_tt, cat_entry *star, short int accuracy,

                    double *ra, double *dec)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the apparent place of a star at date 'jd_tt', given its
      catalog mean place, proper motion, parallax, and radial velocity.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for apparent place.
      *star (struct cat_entry)
         Pointer to catalog entry structure containing catalog data for
         the object in the ICRS (defined in novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Apparent right ascension in hours, referred to true equator
         and equinox of date 'jd_tt'.
      *dec (double)
         Apparent declination in degrees, referred to true equator
         and equinox of date 'jd_tt'.

   RETURNED
   VALUE:
      (short int)
          =  0 ... Everything OK.
          > 10 ... Error code from function 'make_object'.
          > 20 ... Error code from function 'place'.

   GLOBALS
   USED:
      SIZE_OF_OBJ_NAME   novas.h

   FUNCTIONS
   CALLED:
      make_object        novas.c
      place              novas.c
      strcpy             string.h

   VER./DATE/
   PROGRAMMER:
      V1.0/12-05/JAB (USNO/AA)
      V1.1/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.2/10-08/JAB (USNO/AA): Add 'accuracy' option to input.
      V1.3/05-10/JAB (USNO/AA): Fix bug in set-up of 'obj_name'.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'apstar'.
      2. SIZE_OF_OBJ_NAME is defined in novas.h

------------------------------------------------------------------------
*/
{
   char obj_name[SIZE_OF_OBJ_NAME];

   short int error = 0;
   short int type, number, coord_sys;

   double delta_t = 0.0;

   object cel_obj;

   observer location;

   sky_pos output;

/*
   Set 'obj_name' equal to 'starname' in the 'star' structure.  Length
   will be checked in 'make_object'.
*/

   strcpy (obj_name, star->starname);

/*
   Set up a structure of type 'object' containing the star data.
*/

   type = 2;      /* Object located outside the solar system */
   number = 0;    /* Set to zero */

   if ((error = make_object (type,number,obj_name,star, &cel_obj)) != 0)
   {
      error += 10;
      return (error);
   }

/*
   Compute the apparent place with a call to function 'place'.
*/

   location.where = 0;   /* Geocenter */
   coord_sys = 1;        /* True equator and equinox of date */

   if ((error = place (jd_tt,&cel_obj,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      error += 20;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
   }

   return (error);
}

/********virtual_star */

short int virtual_star (double jd_tt, cat_entry *star,
                        short int accuracy,

                        double *ra, double *dec)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the virtual place of a star at date 'jd_tt', given its
      catalog mean place, proper motion, parallax, and radial velocity.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for virtual place.
      *star (struct cat_entry)
         Pointer to catalog entry structure containing catalog data for
         the object in the ICRS (defined in novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Virtual right ascension in hours, referred to the GCRS.
      *dec (double)
         Virtual declination in degrees, referred to the GCRS.

   RETURNED
   VALUE:
      (short int)
          =  0 ... Everything OK.
          > 10 ... Error code from function 'make_object'.
          > 20 ... Error code from function 'place'.

   GLOBALS
   USED:
      SIZE_OF_OBJ_NAME   novas.h


   FUNCTIONS
   CALLED:
      make_object        novas.c
      place              novas.c
      strcpy             string.h

   VER./DATE/
   PROGRAMMER:
      V1.0/12-05/JAB (USNO/AA)
      V1.1/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.2/10-08/JAB (USNO/AA): Add 'accuracy' option to input.
      V1.3/05-10/JAB (USNO/AA): Fix bug in set-up of 'obj_name'.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'vpstar'.
      2. SIZE_OF_OBJ_NAME is defined in novas.h

------------------------------------------------------------------------
*/
{
   char obj_name[SIZE_OF_OBJ_NAME];

   short int error = 0;
   short int type, number, coord_sys;

   double delta_t = 0.0;

   object cel_obj;

   observer location;

   sky_pos output;

/*
   Set 'obj_name' equal to 'starname' in the 'star' structure.  Length
   will be checked in 'make_object'.
*/

   strcpy (obj_name, star->starname);

/*
   Set up a structure of type 'object' containing the star data.
*/

   type = 2;      /* Object located outside the solar system */
   number = 0;    /* Set to zero */

   if ((error = make_object (type,number,obj_name,star, &cel_obj)) != 0)
   {
      error += 10;
      return (error);
   }

/*
   Compute the virtual place with a call to function 'place'.
*/

   location.where = 0;   /* Geocenter */
   coord_sys = 0;        /* GCRS */

   if ((error = place (jd_tt,&cel_obj,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      error += 20;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
   }

   return (error);
}

/********astro_star */

short int astro_star (double jd_tt, cat_entry *star, short int accuracy,

                      double *ra, double *dec)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the astrometric place of a star at date 'jd_tt', given
      its catalog mean place, proper motion, parallax, and radial
      velocity.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for astrometric place.
      *star (struct cat_entry)
         Pointer to catalog entry structure containing catalog data for
         the object in the ICRS (defined in novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Astrometric right ascension in hours (referred to the ICRS,
         without light deflection or aberration).
      *dec (double)
         Astrometric declination in degrees (referred to the ICRS,
         without light deflection or aberration).

   RETURNED
   VALUE:
      (short int)
          =  0 ... Everything OK.
          > 10 ... Error code from function 'make_object'.
          > 20 ... Error code from function 'place'.

   GLOBALS
   USED:
      SIZE_OF_OBJ_NAME   novas.h

   FUNCTIONS
   CALLED:
      make_object        novas.c
      place              novas.c
      strcpy             string.h

   VER./DATE/
   PROGRAMMER:
      V1.0/12-05/JAB (USNO/AA)
      V1.1/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.2/10-08/JAB (USNO/AA): Add 'accuracy' option to input.
      V1.3/05-10/JAB (USNO/AA): Fix bug in set-up of 'obj_name'.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'asstar'.
      2. SIZE_OF_OBJ_NAME is defined in novas.h

------------------------------------------------------------------------
*/
{
   char obj_name[SIZE_OF_OBJ_NAME];

   short int error = 0;
   short int type, number, coord_sys;

   double delta_t = 0.0;

   object cel_obj;

   observer location;

   sky_pos output;

/*
   Set 'obj_name' equal to 'starname' in the 'star' structure.  Length
   will be checked in 'make_object'.
*/

   strcpy (obj_name, star->starname);

/*
   Set up a structure of type 'object' containing the star data.
*/

   type = 2;      /* Object located outside the solar system */
   number = 0;    /* Set to zero */

   if ((error = make_object (type,number,obj_name,star, &cel_obj)) != 0)
   {
      error += 10;
      return (error);
   }

/*
   Compute the astrometric place with a call to function 'place'.
*/

   location.where = 0;   /* Geocenter */
   coord_sys = 3;        /* ICRS astrometric coordinates */

   if ((error = place (jd_tt,&cel_obj,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      error += 20;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
   }

   return (error);
}

/********app_planet */

short int app_planet (double jd_tt, object *ss_body, short int accuracy,

                      double *ra, double *dec, double *dis)
/*
------------------------------------------------------------------------

   PURPOSE:
      Compute the apparent place of a solar system body.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for apparent place.
      *ss_body (struct object)
         Pointer to structure containing the body designation for the
         solar system body (defined in novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Apparent right ascension in hours, referred to true equator
         and equinox of date.
      *dec (double)
         Apparent declination in degrees, referred to true equator
         and equinox of date.
      *dis (double)
         True distance from Earth to the body at 'jd_tt' in AU.

   RETURNED
   VALUE:
      (short int)
          =  0 ... Everything OK.
          =  1 ... Invalid value of 'type' in structure 'ss_body'.
          > 10 ... Error code from function 'place'.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      place              novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/12-05/JAB (USNO/AA)
      V1.1/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.2/10-08/JAB (USNO/AA): Add 'accuracy' option to input.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'applan'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int coord_sys;

   double delta_t = 0.0;

   observer location;

   sky_pos output;

/*
   Check for a valid value of 'type' in structure 'ss_body'.
*/

   if ((ss_body->type < 0) || (ss_body->type > 1))
   {
      error = 1;
      return (error);
   }

/*
   Compute the apparent place with a call to function 'place'.
*/

   location.where = 0;   /* Geocenter */
   coord_sys = 1;        /* True equator and equinox of date */

   if ((error = place (jd_tt,ss_body,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      *dis = 0.0;
      error += 10;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
      *dis = output.dis;
   }

   return (error);
}

/********virtual_planet */

short int virtual_planet (double jd_tt, object *ss_body,
                          short int accuracy,

                          double *ra, double *dec, double *dis)
/*
------------------------------------------------------------------------

   PURPOSE:
      Compute the virtual place of a solar system body.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for virtual place.
      *ss_body (struct object)
         Pointer to structure containing the body designation for the
         solar system body (defined in novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Virtual right ascension in hours, referred to the GCRS.
      *dec (double)
         Virtual declination in degrees, referred to the GCRS.
      *dis (double)
         True distance from Earth to the body in AU.

   RETURNED
   VALUE:
      (short int)
          =  0 ... Everything OK.
          =  1 ... Invalid value of 'type' in structure 'ss_body'.
          > 10 ... Error code from function 'place'.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      place              novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/12-05/JAB (USNO/AA)
      V1.1/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.2/10-08/JAB (USNO/AA): Add 'accuracy' option to input.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'vpplan'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int coord_sys;

   double delta_t = 0.0;

   observer location;

   sky_pos output;

/*
   Check for a valid value of 'type' in structure 'ss_body'.
*/

   if ((ss_body->type < 0) || (ss_body->type > 1))
   {
      error = 1;
      return (error);
   }

/*
   Compute the virtual place with a call to function 'place'.
*/

   location.where = 0;   /* Geocenter */
   coord_sys = 0;        /* GCRS */

   if ((error = place (jd_tt,ss_body,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      *dis = 0.0;
      error += 10;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
      *dis = output.dis;
   }

   return (error);
}

/********astro_planet */

short int astro_planet (double jd_tt, object *ss_body,
                        short int accuracy,

                        double *ra, double *dec, double *dis)
/*
------------------------------------------------------------------------

   PURPOSE:
      Compute the astrometric place of a solar system body.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for astrometric place.
      *ss_body (struct object)
         Pointer to structure containing the body designation for the
         solar system body (defined in novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Astrometric right ascension in hours (referred to the ICRS,
         without light deflection or aberration).
      *dec (double)
         Astrometric declination in degrees (referred to the ICRS,
         without light deflection or aberration).
      *dis (double)
         True distance from Earth to the body in AU.

   RETURNED
   VALUE:
      (short int)
          =  0 ... Everything OK.
          =  1 ... Invalid value of 'type' in structure 'ss_body'.
          > 10 ... Error code from function 'place'.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      place              novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/12-05/JAB (USNO/AA)
      V1.1/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.2/10-08/JAB (USNO/AA): Add 'accuracy' option to input.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'asplan'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int coord_sys;

   double delta_t = 0.0;

   observer location;

   sky_pos output;

/*
   Check for a valid value of 'type' in structure 'ss_body'.
*/

   if ((ss_body->type < 0) || (ss_body->type > 1))
   {
      error = 1;
      return (error);
   }

/*
   Compute the astrometric place with a call to function 'place'.
*/

   location.where = 0;   /* Geocenter */
   coord_sys = 3;        /* ICRS astrometric coordinates */

   if ((error = place (jd_tt,ss_body,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      *dis = 0.0;
      error += 10;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
      *dis = output.dis;
   }

   return (error);
}

/********topo_star */

short int topo_star (double jd_tt, double delta_t, cat_entry *star,
                     on_surface *position, short int accuracy,

                     double *ra, double *dec)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the topocentric place of a star at date 'jd_tt', given its
      catalog mean place, proper motion, parallax, and radial velocity.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for topocentric place.
      delta_t (double)
         Difference TT-UT1 at 'jd_tt', in seconds of time.
      *star (struct cat_entry)
         Pointer to catalog entry structure containing catalog data for
         the object in the ICRS (defined in novas.h).
      *position (struct on_surface)
         Specifies the position of the observer (structure defined in
         novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Topocentric right ascension in hours, referred to true equator
         and equinox of date 'jd_tt'.
      *dec (double)
         Topocentric declination in degrees, referred to true equator
         and equinox of date 'jd_tt'.

   RETURNED
   VALUE:
      (short int)
           =  0 ... Everything OK.
           =  1 ... Invalid value of 'where' in structure 'location'.
           > 10 ... Error code from function 'make_object'.
           > 20 ... Error code from function 'place'.

   GLOBALS
   USED:
      SIZE_OF_OBJ_NAME   novas.h

   FUNCTIONS
   CALLED:
      make_observer      novas.c
      make_object        novas.c
      place              novas.c
      strcpy             string.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-06/JAB (USNO/AA)
      V1.1/01-06/WKP (USNO/AA): Fixed minor syntax problem.
      V1.2/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.3/10-08/JAB (USNO/AA): Add 'accuracy' option to input.
      V1.4/05-10/JAB (USNO/AA): Fix bug in set-up of 'obj_name'.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'tpstar'.
      2. SIZE_OF_OBJ_NAME is defined in novas.h

------------------------------------------------------------------------
*/
{
   char obj_name[SIZE_OF_OBJ_NAME];

   short int error = 0;
   short int type, number, coord_sys;

   in_space dummy;

   observer location;

   object cel_obj;

   sky_pos output;

/*
   Initialize the 'dummy' structure.
*/

   dummy.sc_pos[0] = dummy.sc_pos[1] = dummy.sc_pos[2] = 0.0;
   dummy.sc_vel[0] = dummy.sc_vel[1] = dummy.sc_vel[2] = 0.0;

/*
   Set up a structure of type 'observer' containing the position
   of the observer.
*/

   if ((error = make_observer (1,position,&dummy, &location)) != 0)
   {
      error = 1;
      return (error);
   }

/*
   Set 'obj_name' equal to 'starname' in the 'star' structure.  Length
   will be checked in 'make_object'.
*/

   strcpy (obj_name, star->starname);

/*
   Set up a structure of type 'object' containing the star data.
*/

   type = 2;      /* Object located outside the solar system */
   number = 0;    /* Set to zero */

   if ((error = make_object (type,number,obj_name,star, &cel_obj)) != 0)
   {
      error += 10;
      return (error);
   }

/*
   Compute the topocentric place with a call to function 'place'.
*/

   coord_sys = 1;        /* True equator and equinox of date */

   if ((error = place (jd_tt,&cel_obj,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      error += 20;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
   }

   return (error);
}

/********local_star */

short int local_star (double jd_tt, double delta_t, cat_entry *star,
                      on_surface *position, short int accuracy,

                      double *ra, double *dec)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the local place of a star at date 'jd_tt', given its
      catalog mean place, proper motion, parallax, and radial velocity.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for local place.
      delta_t (double)
         Difference TT-UT1 at 'jd_tt', in seconds of time.
      *star (struct cat_entry)
         Pointer to catalog entry structure containing catalog data for
         the object in the ICRS (defined in novas.h).
      *position (struct on_surface)
         Specifies the position of the observer (structure defined in
         novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Local right ascension in hours, referred to the 'local GCRS'.
      *dec (double)
         Local declination in degrees, referred to the 'local GCRS'.

   RETURNED
   VALUE:
      (short int)
           =  0 ... Everything OK.
           =  1 ... Invalid value of 'where' in structure 'location'.
           > 10 ... Error code from function 'make_object'.
           > 20 ... Error code from function 'place'.

   GLOBALS
   USED:
      SIZE_OF_OBJ_NAME   novas.h

   FUNCTIONS
   CALLED:
      make_observer      novas.c
      make_object        novas.c
      place              novas.c
      strcpy             string.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-06/JAB (USNO/AA)
      V1.1/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.2/10-08/JAB (USNO/AA): Add 'accuracy' option to input.
      V1.3/05-10/JAB (USNO/AA): Fix bug in set-up of 'obj_name'.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'lpstar'.
      2. SIZE_OF_OBJ_NAME is defined in novas.h

------------------------------------------------------------------------
*/
{
   char obj_name[SIZE_OF_OBJ_NAME];

   short int error = 0;
   short int type, number, coord_sys;

   in_space dummy;

   observer location;

   object cel_obj;

   sky_pos output;

/*
   Initialize the 'dummy' structure.
*/

   dummy.sc_pos[0] = dummy.sc_pos[1] = dummy.sc_pos[2] = 0.0;
   dummy.sc_vel[0] = dummy.sc_vel[1] = dummy.sc_vel[2] = 0.0;

/*
   Set up a structure of type 'observer' containing the position
   of the observer.
*/

   if ((error = make_observer (1,position,&dummy, &location)) != 0)
   {
      error = 1;
      return (error);
   }

/*
   Set 'obj_name' equal to 'starname' in the 'star' structure.  Length
   will be checked in 'make_object'.
*/

   strcpy (obj_name, star->starname);

/*
   Set up a structure of type 'object' containing the star data.
*/

   type = 2;      /* Object located outside the solar system */
   number = 0;    /* Set to zero */

   if ((error = make_object (type,number,obj_name,star, &cel_obj)) != 0)
   {
      error += 10;
      return (error);
   }

/*
   Compute the local place with a call to function 'place'.
*/

   coord_sys = 0;        /* "Local GCRS" */

   if ((error = place (jd_tt,&cel_obj,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      error += 20;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
   }

   return (error);
}

/********topo_planet */

short int topo_planet (double jd_tt, object *ss_body, double delta_t,
                       on_surface *position, short int accuracy,

                       double *ra, double *dec, double *dis)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the topocentric place of a solar system body.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for topocentric place.
      *ss_body (struct object)
         Pointer to structure containing the body designation for the
         solar system body (defined in novas.h).
      delta_t (double)
         Difference TT-UT1 at 'jd_tt', in seconds of time.
      *position (struct on_surface)
         Specifies the position of the observer (structure defined in
         novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Topocentric right ascension in hours, referred to true equator
         and equinox of date.
      *dec (double)
         Topocentric declination in degrees, referred to true equator
         and equinox of date.
      *dis (double)
         True distance from Earth to the body at 'jd_tt' in AU.

   RETURNED
   VALUE:
      (short int)
           =  0 ... Everything OK.
           =  1 ... Invalid value of 'where' in structure 'location'.
           > 10 ... Error code from function 'place'.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      make_observer      novas.c
      place              novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/01-06/JAB (USNO/AA)
      V1.1/01-06/WKP (USNO/AA): Fixed minor syntax problem.
      V1.2/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.3/10-08/JAB (USNO/AA): Add 'accuracy' option to input.
      V1.4/07-10/JLB (USNO/AA): Fixed minor documentation error in prolog

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'tpplan'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int coord_sys;

   in_space dummy;

   observer location;

   sky_pos output;

/*
   Initialize the 'dummy' structure.
*/

   dummy.sc_pos[0] = dummy.sc_pos[1] = dummy.sc_pos[2] = 0.0;
   dummy.sc_vel[0] = dummy.sc_vel[1] = dummy.sc_vel[2] = 0.0;

/*
   Set up a structure of type 'observer' containing the position
   of the observer.
*/

   if ((error = make_observer (1,position,&dummy, &location)) != 0)
   {
      error = 1;
      return (error);
   }

/*
   Compute the topocentric place with a call to function 'place'.
*/

   coord_sys = 1;        /* True equator and equinox of date */

   if ((error = place (jd_tt,ss_body,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      *dis = 0.0;
      error += 10;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
      *dis = output.dis;
   }

   return (error);
}

/********local_planet */

short int local_planet (double jd_tt, object *ss_body,
                        double delta_t, on_surface *position,
                        short int accuracy,

                       double *ra, double *dec, double *dis)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the local place of a solar system body.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for local place.
      *ss_body (struct object)
         Pointer to structure containing the body designation for the
         solar system body (defined in novas.h).
      delta_t (double)
         Difference TT-UT1 at 'jd_tt', in seconds of time.
      *position (struct on_surface)
         Specifies the position of the observer (structure defined in
         novas.h).
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra (double)
         Local right ascension in hours, referred to the 'local GCRS'.
      *dec (double)
         Local declination in degrees, referred to the 'local GCRS'.
      *dis (double)
         True distance from Earth to the body in AU.

   RETURNED
   VALUE:
      (short int)
           =  0 ... Everything OK.
           =  1 ... Invalid value of 'where' in structure 'location'.
           > 10 ... Error code from function 'place'.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      make_observer      novas.c
      place              novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/01-06/JAB (USNO/AA)
      V1.1/01-06/WKP (USNO/AA): Fixed minor syntax problem.
      V1.2/10-06/JAB (USNO/AA): Incorporate 'output' structure.
      V1.3/10-08/JAB (USNO/AA): Add 'accuracy' option to input.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'lpplan'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int coord_sys;

   in_space dummy;

   observer location;

   sky_pos output;

/*
   Initialize the 'dummy' structure.
*/

   dummy.sc_pos[0] = dummy.sc_pos[1] = dummy.sc_pos[2] = 0.0;
   dummy.sc_vel[0] = dummy.sc_vel[1] = dummy.sc_vel[2] = 0.0;

/*
   Set up a structure of type 'observer' containing the position
   of the observer.
*/

   if ((error = make_observer (1,position,&dummy, &location)) != 0)
   {
      error = 1;
      return (error);
   }

/*
   Compute the local place with a call to function 'place'.
*/

   coord_sys = 0;        /* "Local GCRS" */

   if ((error = place (jd_tt,ss_body,&location,delta_t,coord_sys,
      accuracy, &output)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      *dis = 0.0;
      error += 10;
   }
    else
   {
      *ra = output.ra;
      *dec = output.dec;
      *dis = output.dis;
   }

   return (error);
}

/********mean_star */

short int mean_star (double jd_tt, double ra, double dec,
                     short int accuracy,

                     double *ira, double *idec)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the ICRS position of a star, given its apparent place
       at date 'jd_tt'.  Proper motion, parallax and radial velocity
      are assumed to be zero.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Explanatory Supplement to the Astronomical Almanac (1992),
         Chapter 3.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date of apparent place.
      ra (double)
         Apparent right ascension in hours, referred to true equator
         and equinox of date.
      dec (double)
         Apparent declination in degrees, referred to true equator
         and equinox of date.
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ira (double)
         ICRS right ascension in hours.
      *idec (double)
         ICRS declination in degrees.

   RETURNED
   VALUE:
      (short int)
          =  0 ... Everything OK.
          =  1 ... Iterative process did not converge after 30 iterations.
          =  2 ... length of dummy star name out of bounds.
          =  3 ... length of dummy catalog name out of bounds.
          > 10 ... Error from function 'vector2radec'.
          > 20 ... Error from function 'app_star'.

   GLOBALS
   USED:
      T0                 novascon.c

   FUNCTIONS
   CALLED:
      make_cat_entry     novas.c
      starvectors        novas.c
      precession         novas.c
      app_star           novas.c
      vector2radec       novas.c
      fabs               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/12-05/JAB (USNO/AA)
      V1.1/03-08/WKP (USNO/AA): Fixed subtle bug in exit condition of
                                iteration and updated variable names.
      V1.2/10-08/JAB (USNO/AA): Add 'accuracy' option to input.
      V1.3/05-10/JAB (USNO/AA): Use 'make_cat_entry' to set up
                                local structure 'tempstar'.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'mpstar'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int iter = 0;

   double pos[3], dum[3], pos2[3], newira, newidec, oldira, oldidec, ra2,
      dec2, deltara, deltadec;

   cat_entry tempstar;

/*
   Set up the 'tempstar' structure, then use it to create a position
   vector based on the apparent RA and declination of the star.
*/

   if ((error = make_cat_entry ("dummy","CAT",0,ra,dec,0.0,0.0,0.0,0.0,
      &tempstar)) != 0)
   {
      return (error + 1);
   }

   starvectors (&tempstar, pos,dum);

/*
   Get initial approximation by precessing star position at 'jd_tt'
   to its position at J2000.0.
*/

   precession (jd_tt,pos,T0, pos2);
   if ((error = vector2radec (pos2, &newira,&newidec)) != 0)
   {
      return (error + 10);
   }

/*
   Iteratively find ICRS coordinates that produce input apparent place
   of star at date 'jd_tt'.
*/

   do
   {
      oldira = newira;
      oldidec = newidec;
      tempstar.ra = oldira;
      tempstar.dec = oldidec;
      if ((error = app_star (jd_tt,&tempstar,accuracy, &ra2,&dec2)) != 0)
      {
         *ira = 0.0;
         *idec = 0.0;
         return (error + 20);
      }

      deltara = ra2 - oldira;
      deltadec = dec2 - oldidec;
      if (deltara < -12.0)
         deltara += 24.0;
      if (deltara > 12.0)
         deltara -= 24.0;
      newira = ra - deltara;
      newidec = dec - deltadec;

      if (iter >= 30)
      {
         *ira = 0.0;
         *idec = 0.0;
         return (error = 1);
      }
       else
      {
         iter++;
      }
   } while (!(fabs (newira - oldira) <= 1.0e-12) ||
            !(fabs (newidec - oldidec) <= 1.0e-11));

   *ira = newira;
   *idec = newidec;
   if (*ira < 0.0)
      *ira += 24.0;
   if (*ira >= 24.0)
      *ira -= 24.0;

   return (error);
}


/********place */

short int place (double jd_tt, object *cel_object,
                 observer *location, double delta_t,
                 short int coord_sys, short int accuracy,

                 sky_pos *output)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function computes the apparent direction of a star or solar
      system body at a specified time and in a specified coordinate
      system.

   REFERENCES:
      Kaplan, G. et al. (1989), Astronomical Journal 97, 1197-1210.
      Klioner, S. (2003), Astronomical Journal 125, 1580-1597.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date for place.
      *cel_object (struct object)
         Specifies the celestial object of interest (structure defined
         in novas.h).
      *location (struct observer)
         Specifies the location of the observer (structure defined in
         novas.h).
      delta_t (double)
         Difference TT-UT1 at 'jd_tt', in seconds of time.
      coord_sys (short int)
         Code specifying coordinate system of the output position.
            = 0 ... GCRS or "local GCRS"
            = 1 ... true equator and equinox of date
            = 2 ... true equator and CIO of date
            = 3 ... astrometric coordinates, i.e., without light
                    deflection or aberration.
      accuracy (short int)
         Code specifying the relative accuracy of the output position.
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *output (struct sky_pos)
         Output data specifying object's place on the sky at time
         'jd_tt', with respect to the specified output coordinate system
         (struct defined in novas.h).

   RETURNED
   VALUE:
      = 0         ... No problems.
      = 1         ... invalid value of 'coord_sys'
      = 2         ... invalid value of 'accuracy'
      = 3         ... Earth is the observed object, and the observer is
                      either at the geocenter or on the Earth's surface
                      (not permitted)
      > 10, < 40  ... 10 + error from function 'ephemeris'
      > 40, < 50  ... 40 + error from function 'geo_posvel'
      > 50, < 70  ... 50 + error from function 'light_time'
      > 70, < 80  ... 70 + error from function 'grav_def'
      > 80, < 90  ... 80 + error from function 'cio_location'
      > 90, < 100 ... 90 + error from function 'cio_basis'

   GLOBALS
   USED:
      T0, C_AUDAY        novascon.c

   FUNCTIONS
   CALLED:
      make_cat_entry     novas.c
      make_object        novas.c
      tdb2tt             novas.c
      ephemeris          novas.c
      geo_posvel         novas.c
      starvectors        novas.c
      d_light            novas.c
      proper_motion      novas.c
      bary2obs           novas.c
      light_time         novas.c
      limb_angle         novas.c
      grav_def           novas.c
      aberration         novas.c
      frame_tie          novas.c
      precession         novas.c
      nutation           novas.c
      cio_location       novas.c
      cio_basis          novas.c
      rad_vel            novas.c
      vector2radec       novas.c
      fabs               math.h
      sqrt               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-05/JAB (USNO/AA)
      V1.1/09-05/WKP (USNO/AA) Removed 'star' from the input argument
                               list since it is nested in the 'object'
                               structure.
      V1.2/01-06/WKP (USNO/AA) Updated error codes and replaced 'mode'
                               with 'accuracy' where appropriate.
      V1.3/01-06/JAB (USNO/AA) Added more static variables.
      V1.4/07-06/JAB (USNO/AA) Implement 'cio_location' construct.
      V1.5/10-06/JAB (USNO/AA) Add radial velocity calculation and
                               implement 'sky_pos' structure.
      V1.6/01-07/JAB (USNO/AA) Update to accommodate new radial velocity
                               algorithm.
      V1.7/10-08/JAB (USNO/AA) Modify calls to 'ephemeris' to support
                               two-part input Julian date.
      V1.8/07-10/JLB (USNO/AA) Corrected citation to Kaplan et al.

   NOTES:
      1. Values of 'location->where' and 'coord_sys' dictate the various
      standard kinds of place:
        location->where = 0 and coord_sys = 1: apparent place
        location->where = 1 and coord_sys = 1: topocentric place
        location->where = 0 and coord_sys = 0: virtual place
        location->where = 1 and coord_sys = 0: local place
        location->where = 0 and coord_sys = 3: astrometric place
        location->where = 1 and coord_sys = 3: topocentric astrometric
                                               place
      2. Input value of 'delta_t' is used only when 'location->where'
         equals 1 or 2 (observer is on surface of Earth or in a
         near-Earth satellite).
------------------------------------------------------------------------
*/
{
   static short int first_time = 1;
   short int error = 0;
   short int loc, rs, i;

   static double tlast1 = 0.0;
   static double tlast2 = 0.0;
   static double jd_tdb, peb[3], veb[3], psb[3], vsb[3], px[3], py[3],
      pz[3];
   double x, secdif, jd[2], pog[3], vog[3], pob[3], vob[3], pos1[3],
      vel1[3], dt, pos2[3], pos3[3], t_light, t_light0, pos4[3], frlimb,
      pos5[3], pos6[3], pos7[3], pos8[3], r_cio, d_obs_geo, d_obs_sun,
      d_obj_sun;

   cat_entry null_star;

   static object earth, sun;

/*
   Check for invalid value of 'coord_sys' or 'accuracy'.
*/

   if ((coord_sys < 0) || (coord_sys > 3))
      return (error = 1);

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 2);

/*
   Create a null star 'cat_entry' and  Earth and Sun 'object's.
*/

   if (first_time)
   {
      make_cat_entry ("NULL_STAR","   ",0L,0.0,0.0,0.0,0.0,0.0,0.0,
         &null_star);

      make_object (0,3,"Earth",&null_star, &earth);
      make_object (0,10,"Sun",&null_star, &sun);

      first_time = 0;
   }

/*
   ---------------------------------------------------------------------
   Check on Earth as an observed object.  Earth can only be an observed
   object when 'location' is a near-Earth satellite.
   ---------------------------------------------------------------------
*/

   if ((cel_object->type == 0) && (cel_object->number == 3) &&
      (location->where != 2))
      return (error = 3);

/*
   ---------------------------------------------------------------------
   Get position and velocity of Earth (geocenter) and Sun.
   ---------------------------------------------------------------------
*/

   if (fabs (jd_tt - tlast1) > 1.0e-8)
   {

/*
   Compute 'jd_tdb', the TDB Julian date corresponding to 'jd_tt'.
*/

      jd_tdb = jd_tt;
      tdb2tt (jd_tdb, &x,&secdif);
      jd_tdb = jd_tt + secdif / 86400.0;

/*
   Get position and velocity of Earth wrt barycenter of solar system,
   in ICRS.
*/

      jd[0] = jd_tdb;
      jd[1] = 0.0;

      if ((error = ephemeris (jd,&earth,0,accuracy, peb,veb)) != 0)
         return (error += 10);

/*
   Get position and velocity of Sun wrt barycenter of solar system,
   in ICRS.
*/

      if ((error = ephemeris (jd,&sun,0,accuracy, psb,vsb)) != 0)
         return (error += 10);

      tlast1 = jd_tt;
   }

/*
   ---------------------------------------------------------------------
   Get position and velocity of observer.
   ---------------------------------------------------------------------
*/

   if ((location->where == 1) || (location->where == 2))
   {

/*
   For topocentric place, get geocentric position and velocity vectors
   of observer (observer is on surface of Earth or in a near-Earth
   satellite).
*/

      if ((error = geo_posvel (jd_tt,delta_t,accuracy,location, pog,vog))
         != 0)
         return (error += 40);

      loc = 1;
   }
    else
   {

/*
   For geocentric place, there is nothing to do (observer is at
   geocenter).
*/

    for (i = 0; i < 3; i++)
    {
       pog[i] = 0.0;
       vog[i] = 0.0;
    }

    loc = 0;
   }

/*
   Compute position and velocity of observer wrt barycenter of
   solar system (Galilean transformation).
*/

    for (i = 0; i < 3; i++)
    {
       pob[i] = peb[i] + pog[i];
       vob[i] = veb[i] + vog[i];
    }

/*
   ---------------------------------------------------------------------
   Find geometric position of observed object.
   ---------------------------------------------------------------------
*/

   if (cel_object->type == 2)  /* Observed object is star. */
   {

/*
   Get position of star updated for its space motion.
*/

      starvectors (&cel_object->star, pos1,vel1);
      dt = d_light (pos1,pob);
      proper_motion (T0,pos1,vel1,(jd_tdb + dt), pos2);

/*
   Get position of star wrt observer (corrected for parallax).
*/

      bary2obs (pos2,pob, pos3,&t_light);
      output->dis = 0.0;
   }

    else   /* Observed object is solar system body. */
   {

/*
   Get position of body wrt barycenter of solar system.
*/

      jd[0] = jd_tdb;
      jd[1] = 0.0;

      if ((error = ephemeris (jd,cel_object,0,accuracy, pos1,vel1)) != 0)
         return (error += 10);

/*
   Get position of body wrt observer, and true (Euclidian) distance.
*/

      bary2obs (pos1,pob, pos2,&t_light0);
      output->dis = t_light0 * C_AUDAY;

/*
   Get position of body wrt observer, antedated for light-time.
*/

      if ((error = light_time (jd_tdb,cel_object,pob,t_light0,accuracy,
           pos3,&t_light)) != 0)
         return (error += 50);
   }

/*
   ---------------------------------------------------------------------
   Apply gravitational deflection of light and aberration.
   ---------------------------------------------------------------------
*/

   if (coord_sys == 3)
   {

/*
   These calculations are skipped for astrometric place.
*/

      for (i = 0; i < 3; i++)
      {
         pos5[i] = pos3[i];
      }
   }

    else
   {

/*
   Variable 'loc' determines whether Earth deflection is included.
*/

      if (loc == 1)
      {
         limb_angle (pos3,pog, &x,&frlimb);
         if (frlimb < 0.8)
            loc = 0;
      }

/*
   Compute gravitational deflection and aberration.
*/

      if ((error = grav_def (jd_tdb,loc,accuracy,pos3,pob, pos4)) != 0)
         return (error += 70);

      aberration (pos4,vob,t_light, pos5);
   }

/*
   ---------------------------------------------------------------------
   Transform, if necessary, to output coordinate system.
   ---------------------------------------------------------------------
*/

   switch (coord_sys)
   {
      case (1):    /* Transform to equator and equinox of date. */

         frame_tie (pos5,1, pos6);
         precession (T0,pos6,jd_tdb, pos7);
         nutation (jd_tdb,0,accuracy,pos7, pos8);
         break;

      case (2):    /* Transform to equator and CIO of date. */

         if (fabs (jd_tdb - tlast2) > 1.0e-8 )
         {

/*
   Obtain the basis vectors, in the GCRS, of the celestial intermediate
   system.
*/

            if ((error = cio_location (jd_tdb,accuracy, &r_cio,
               &rs)) != 0)
               return (error += 80);
            if ((error = cio_basis (jd_tdb,r_cio,rs,accuracy,
               px,py,pz)) != 0)
               return (error += 90);

            tlast2 = jd_tdb;
         }

/*
   Transform position vector to celestial intermediate system.
*/

         pos8[0] = px[0] * pos5[0] + px[1] * pos5[1] + px[2] * pos5[2];
         pos8[1] = py[0] * pos5[0] + py[1] * pos5[1] + py[2] * pos5[2];
         pos8[2] = pz[0] * pos5[0] + pz[1] * pos5[1] + pz[2] * pos5[2];
         break;

      default:     /* No transformation -- keep coordinates in GCRS, */
                   /* or ICRS for astrometric coordinates.           */

         for (i = 0; i < 3; i++)
         {
            pos8[i] = pos5[i];
         }
   }
/*
   ---------------------------------------------------------------------
   Compute radial velocity.
   ---------------------------------------------------------------------
*/

/*
   Compute distances: observer-geocenter, observer-Sun, object-Sun.
*/

   d_obs_geo = sqrt ((pob[0] - peb[0]) * (pob[0] - peb[0]) +
                    ( pob[1] - peb[1]) * (pob[1] - peb[1]) +
                    ( pob[2] - peb[2]) * (pob[2] - peb[2]));

   d_obs_sun = sqrt ((pob[0] - psb[0]) * (pob[0] - psb[0]) +
                    ( pob[1] - psb[1]) * (pob[1] - psb[1]) +
                    ( pob[2] - psb[2]) * (pob[2] - psb[2]));

   d_obj_sun = sqrt ((pos1[0] - psb[0]) * (pos1[0] - psb[0]) +
                    ( pos1[1] - psb[1]) * (pos1[1] - psb[1]) +
                    ( pos1[2] - psb[2]) * (pos1[2] - psb[2]));

   rad_vel (cel_object,pos3,vel1,vob,d_obs_geo,d_obs_sun,d_obj_sun,
      &output->rv);

/*
   ---------------------------------------------------------------------
   Finish up.
   ---------------------------------------------------------------------
*/

   vector2radec (pos8, &output->ra,&output->dec);

   x = sqrt (pos8[0] * pos8[0] + pos8[1] * pos8[1] + pos8[2] * pos8[2]);

   for (i = 0; i < 3; i++)
   {
      output->r_hat[i] = pos8[i] / x;
   }

   return (error);
}

/********equ2gal */

void equ2gal (double rai, double deci,

              double *glon, double *glat)
/*
------------------------------------------------------------------------

   PURPOSE:
      To convert ICRS right ascension and declination to galactic
      longitude and latitude.

   REFERENCES:
      Hipparcos and Tycho Catalogues, Vol. 1, Section 1.5.3.

   INPUT
   ARGUMENTS:
      rai (double)
         ICRS right ascension in hours.
      deci (double)
         ICRS declination in degrees.

   OUTPUT
   ARGUMENTS:
      *glon (double)
         Galactic longitude in degrees.
      *glat (double)
         Galactic latitude in degrees.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      DEG2RAD, RAD2DEG   novascon.c

   FUNCTIONS
   CALLED:
      sin                math.h
      cos                math.h
      sqrt               math.h
      atan2              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/11-03/JAB (USNO/AA)
      V1.1/03-06/JAB (USNO/AA): Fixed initialization of 'ag'.
      V1.2/03-11/WKP (USNO/AA): Added braces to 2-D array initialization
                                to quiet gcc warnings.

   NOTES:
      1. This function uses the coordinate transformation specified
      in the reference.
      2. This function is the C version of NOVAS Fortran routine
      'eqgal'.

------------------------------------------------------------------------
*/
{
   double r, d, pos1[3], pos2[3], xyproj, g;

/*
   Rotation matrix A_g from Hipparcos documentation eq. 1.5.11.
*/

   double ag[3][3] = {
      {-0.0548755604, +0.4941094279, -0.8676661490},
      {-0.8734370902, -0.4448296300, -0.1980763734},
      {-0.4838350155, +0.7469822445, +0.4559837762}};

/*
   Form position vector in equatorial system from input coordinates.
*/

   r = rai * 15.0 * DEG2RAD;
   d = deci * DEG2RAD;
   pos1[0] = cos (d) * cos (r);
   pos1[1] = cos (d) * sin (r);
   pos1[2] = sin (d);

/*
   Rotate position vector to galactic system, using Hipparcos
   documentation eq. 1.5.13.
*/

   pos2[0] = ag[0][0] * pos1[0] + ag[1][0] * pos1[1] +

             ag[2][0] * pos1[2];

   pos2[1] = ag[0][1] * pos1[0] + ag[1][1] * pos1[1] +

             ag[2][1] * pos1[2];

   pos2[2] = ag[0][2] * pos1[0] + ag[1][2] * pos1[1] +

             ag[2][2] * pos1[2];

/*
   Decompose galactic vector into longitude and latitude.
*/

   xyproj = sqrt (pos2[0] * pos2[0] + pos2[1] * pos2[1]);

   if (xyproj > 0.0)
      g = atan2 (pos2[1], pos2[0]);
    else
      g = 0.0;

   *glon = g * RAD2DEG;
   if (*glon < 0.0)
      *glon += 360.0;

   g = atan2 (pos2[2], xyproj);
   *glat = g * RAD2DEG;

   return;
}

/********equ2ecl */

short int equ2ecl (double jd_tt, short int coord_sys,
                   short int accuracy, double ra, double dec,

                   double *elon, double *elat)
/*
------------------------------------------------------------------------

   PURPOSE:
      To convert right ascension and declination to ecliptic longitude
      and latitude.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date of equator, equinox, and ecliptic used for
         coordinates.
      coord_sys (short int)
         Coordinate system selection.
            = 0 ... mean equator and equinox of date 'jd_tt'
            = 1 ... true equator and equinox of date 'jd_tt'
            = 2 ... ICRS
            (ecliptic is always the mean plane)
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      ra (double)
         Right ascension in hours, referred to specified equator and
         equinox of date.
      dec (double)
         Declination in degrees, referred to specified equator and
         equinox of date.

   OUTPUT
   ARGUMENTS:
      *elon (double)
         Ecliptic longitude in degrees, referred to specified ecliptic
         and equinox of date.
      *elat (double)
         Ecliptic latitude in degrees, referred to specified ecliptic
         and equinox of date.

   RETURNED
   VALUE:
      (short int)
         = 0 ... everything OK
         = 1 ... invalid value of 'coord_sys'

   GLOBALS
   USED:
      DEG2RAD, RAD2DEG   novascon.c

   FUNCTIONS
   CALLED:
      equ2ecl_vec        novas.c
      sin                math.h
      cos                math.h
      sqrt               math.h
      atan2              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/11-03/JAB (USNO/AA)
      V1.1/01-06/WKP (USNO/AA) Changed 'mode' to 'accuracy'.
      V1.2/05-06/JAB (USNO/AA) Use vector transformation function.
      V1.3/05-08/WKP (USNO/AA) Changed values of coord_sys to be
                               more consistent with gcrs2equ.

   NOTES:
      1. To convert ICRS RA and dec to ecliptic coordinates (mean
      ecliptic and equinox of J2000.0), set 'coord_sys' = 2; the value
      of 'jd_tt' can be set to anything, since J2000.0 is assumed.
      Except for the input to this case, all input coordinates
      are dynamical.
      2. This function is the C version of NOVAS Fortran routine
      'eqecl'.

------------------------------------------------------------------------
*/
{
   short int error = 0;

   double r, d, pos1[3], pos2[3], xyproj, e;

/*
   Form position vector in equatorial system from input coordinates.
*/

   r = ra * 15.0 * DEG2RAD;
   d = dec * DEG2RAD;
   pos1[0] = cos (d) * cos (r);
   pos1[1] = cos (d) * sin (r);
   pos1[2] = sin (d);

/*
   Convert the vector from equatorial to ecliptic system.
*/

   if ((error = equ2ecl_vec (jd_tt,coord_sys,accuracy,pos1, pos2)) != 0)
      return (error);

/*
   Decompose ecliptic vector into ecliptic longitude and latitude.
*/

   xyproj = sqrt (pos2[0] * pos2[0] + pos2[1] * pos2[1]);

   if (xyproj > 0.0)
      e = atan2 (pos2[1], pos2[0]);
    else
      e = 0.0;

   *elon = e * RAD2DEG;
   if (*elon < 0.0)
      *elon += 360.0;

   e = atan2 (pos2[2], xyproj);
   *elat = e * RAD2DEG;

   return (error);
}

/********equ2ecl_vec */

short int equ2ecl_vec (double jd_tt, short int coord_sys,
                       short int accuracy, double *pos1,

                       double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      Converts an equatorial position vector to an ecliptic position
      vector.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date of equator, equinox, and ecliptic used for
         coordinates.
      coord_sys (short int)
         Coordinate system selection.
            = 0 ... mean equator and equinox of date 'jd_tt'
            = 1 ... true equator and equinox of date 'jd_tt'
            = 2 ... ICRS
            (ecliptic is always the mean plane)
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      pos1[3] (double)
         Position vector, referred to specified equator and equinox of
         date.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, referred to specified ecliptic and equinox
         of date.

   RETURNED
   VALUE:
      (short int)
         = 0 ... everything OK
         = 1 ... invalid value of 'coord_sys'

   GLOBALS
   USED:
      T0, DEG2RAD        novascon.c

   FUNCTIONS
   CALLED:
      tdb2tt             novas.c
      e_tilt             novas.c
      frame_tie          novas.c
      fabs               math.h
      sin                math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/05-06/JAB (USNO/AA)
      V1.1/05-08/WKP (USNO/AA) Changed values of coord_sys to be
                               more consistent with gcrs2equ.

   NOTES:
      1. To convert an ICRS vector to an ecliptic vector (mean ecliptic
      and equinox of J2000.0 only), set 'coord_sys' = 2; the value
      of 'jd_tt' can be set to anything, since J2000.0 is assumed.
      Except for the input to this case, all vectors are assumed
      to be with respect to a dynamical system.
      2. This function is the C version of NOVAS Fortran routine
      'eqec'.

------------------------------------------------------------------------
*/
{
   short int error = 0;

   static double t_last = 0.0;
   static double ob2000 = 0.0;
   static double oblm, oblt;
   double t, secdiff, jd_tdb, pos0[3], w, x, y, z, obl;

/*
   'jd_tdb' is the TDB Julian date corresponding to 'jd_tt'.
*/

   tdb2tt (jd_tt, &t,&secdiff);
   jd_tdb = jd_tt + secdiff / 86400.0;

/*
   Get obliquity, depending upon the "system" of the input coordinates.
*/

   switch (coord_sys)
   {
      case 0:             /* Input: mean equator and equinox of date */
      case 1:             /* Input: true equator and equinox of date */
         pos0[0] = pos1[0];
         pos0[1] = pos1[1];
         pos0[2] = pos1[2];
         if (fabs (jd_tt - t_last) > 1.0e-8)
         {
            e_tilt (jd_tdb,accuracy, &oblm,&oblt,&x,&y,&z);
            t_last = jd_tt;
         }

         switch (coord_sys)
         {
            case 0:       /* Use mean obliquity of date */
               obl = oblm * DEG2RAD;
               break;
            case 1:       /* Use true obliquity of date */
               obl = oblt * DEG2RAD;
               break;
         }
         break;

      case 2:             /* Input: ICRS */
         frame_tie (pos1,1, pos0);

         if (ob2000 == 0.0)
         {
            e_tilt (T0,accuracy, &oblm,&w,&x,&y,&z);
            ob2000 = oblm;
         }
         obl = ob2000 * DEG2RAD;
         break;

      default:
         return (error = 1);
   }

/*
   Rotate position vector to ecliptic system.
*/

   pos2[0] =  pos0[0];
   pos2[1] =  pos0[1] * cos (obl) + pos0[2] * sin (obl);
   pos2[2] = -pos0[1] * sin (obl) + pos0[2] * cos (obl);

   return (error);
}

/********ecl2equ_vec */

short int ecl2equ_vec (double jd_tt, short int coord_sys,
                       short int accuracy, double *pos1,

                       double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      Converts an ecliptic position vector to an equatorial position
      vector.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date of equator, equinox, and ecliptic used for
         coordinates.
      coord_sys (short int)
         Coordinate system selection.
            = 0 ... mean equator and equinox of date
            = 1 ... true equator and equinox of date
            = 2 ... ICRS
            (ecliptic is always the mean plane)
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      pos1[3] (double)
         Position vector, referred to specified ecliptic and equinox of
         date.  If 'coord_sys' = 2, 'pos1' must be on mean ecliptic and
         equinox of J2000.0; see Note 1 below.

   OUTPUT
   ARGUMENTS:
      POS2[3] (double)
         Position vector, referred to specified equator and equinox
         of date.

   RETURNED
   VALUE:
      (short int)
         = 0 ... everything OK
         = 1 ... invalid value of 'coord_sys'

   GLOBALS
   USED:
      T0, DEG2RAD        novascon.c

   FUNCTIONS
   CALLED:
      tdb2tt             novas.c
      e_tilt             novas.c
      frame_tie          novas.c
      fabs               math.h
      sin                math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/05-06/JAB (USNO/AA)
      V1.1/05-08/WKP (USNO/AA) Changed values of coord_sys to be
                               more consistent with gcrs2equ.
      V1.2/09-10/WKP (USNO/AA) Initialized 'obl' variable to silence
                               compiler warning.

   NOTES:
      1. To convert an ecliptic vector (mean ecliptic and equinox of
      J2000.0 only) to an ICRS vector, set 'coord_sys' = 2; the value
      of 'jd_tt' can be set to anything, since J2000.0 is assumed.
      Except for the output from this case, all vectors are assumed to
      be with respect to a dynamical system.
      2. This function is the C version of NOVAS Fortran routine
      'eceq'.

------------------------------------------------------------------------
*/
{
   short int error = 0;

   static double t_last = 0.0;
   static double ob2000 = 0.0;
   static double oblm, oblt;
   double t, secdiff, jd_tdb, pos0[3], w, x, y, z, obl = 0.0;

/*
   'jd_tdb' is the TDB Julian date.
*/

   tdb2tt (jd_tt, &t,&secdiff);
   jd_tdb = jd_tt + secdiff / 86400.0;

/*
   Get obliquity, depending upon the "system" of the input coordinates.
*/

   switch (coord_sys)
   {
      case 0:             /* Output: mean equator and equinox of date */
      case 1:             /* Output: true equator and equinox of date */
         if (fabs (jd_tt - t_last) > 1.0e-8)
         {
            e_tilt (jd_tdb,accuracy, &oblm,&oblt,&x,&y,&z);
            t_last = jd_tt;
         }

         switch (coord_sys)
         {
            case 0:       /* Use mean obliquity of date */
               obl = oblm * DEG2RAD;
               break;
            case 1:       /* Use true obliquity of date */
               obl = oblt * DEG2RAD;
               break;
         }
         break;

      case 2:             /* Output: ICRS */
         if (ob2000 == 0.0)
         {
            e_tilt (T0,accuracy, &oblm,&w,&x,&y,&z);
            ob2000 = oblm;
         }
         obl = ob2000 * DEG2RAD;
         break;

      default:
         return (error = 1);
   }

/*
   Rotate position vector to ecliptic system.
*/

   pos2[0] = pos1[0];
   pos2[1] = pos1[1] * cos (obl) - pos1[2] * sin (obl);
   pos2[2] = pos1[1] * sin (obl) + pos1[2] * cos (obl);

/*
   Case where output vector is to be in ICRS, rotate from dynamical
   system to ICRS.
*/

   if (coord_sys == 2)
   {
      pos0[0] = pos2[0];
      pos0[1] = pos2[1];
      pos0[2] = pos2[2];
      frame_tie (pos0,-1, pos2);
   }

   return (error);
}

/********equ2hor */

void equ2hor (double jd_ut1, double delta_t, short int accuracy,
              double xp, double yp, on_surface *location, double ra,
              double dec, short int ref_option,

              double *zd, double *az, double *rar, double *decr)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function transforms topocentric right ascension and
      declination to zenith distance and azimuth.  It uses a method
      that properly accounts for polar motion, which is significant at
      the sub-arcsecond level.  This function can also adjust
      coordinates for atmospheric refraction.

   REFERENCES:
      Kaplan, G. (2008). USNO/AA Technical Note of 28 Apr 2008,
         "Refraction as a Vector."

   INPUT
   ARGUMENTS:
      jd_ut1 (double)
         UT1 Julian date.
      delta_t (double)
         Difference TT-UT1 at 'jd_ut1', in seconds.
      accuracy (short int)
         Selection for method and accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      xp (double)
         Conventionally-defined x coordinate of celestial intermediate
         pole with respect to ITRS reference pole, in arcseconds.
      yp (double)
         Conventionally-defined y coordinate of celestial intermediate
         pole with respect to ITRS reference pole, in arcseconds.
      *location (struct on_surface)
         Pointer to structure containing observer's location (defined
         in novas.h).
      ra (double)
         Topocentric right ascension of object of interest, in hours,
         referred to true equator and equinox of date.
      dec (double)
         Topocentric declination of object of interest, in degrees,
         referred to true equator and equinox of date.
      ref_option (short int)
         = 0 ... no refraction
         = 1 ... include refraction, using 'standard' atmospheric
                 conditions.
         = 2 ... include refraction, using atmospheric parameters
                 input in the 'location' structure.

   OUTPUT
   ARGUMENTS:
      *zd (double)
         Topocentric zenith distance in degrees, affected by
         refraction if 'ref_option' is non-zero.
      *az (double)
         Topocentric azimuth (measured east from north) in degrees.
      *rar (double)
         Topocentric right ascension of object of interest, in hours,
         referred to true equator and equinox of date, affected by
         refraction if 'ref_option' is non-zero.
      *decr (double)
         Topocentric declination of object of interest, in degrees,
         referred to true equator and equinox of date, affected by
         refraction if 'ref_option' is non-zero.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      DEG2RAD, RAD2DEG   novascon.c

   FUNCTIONS
   CALLED:
      ter2cel            novas.c
      refract            novas.c
      sin                math.h
      cos                math.h
      sqrt               math.h
      atan2              math.h
      fabs               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-98/JAB (USNO/AA)
      V2.0/01-06/JAB (USNO/AA)
      V2.1/01-06/WKP (USNO/AA): Changed 'mode' to 'accuracy'.
      V2.2/11-07/JAB (USNO/AA): Cosmetic changes.
      V2.2/05-08/WKP (USNO/AA): Updated the refraction algorithm.
      V2.3/06-08/WKP (USNO/AA): Tweaked convergence criteria.
      V2.4/03-09/JAB (USNO/AA): Conformed input variables to IERS
                                conventions.

   NOTES:
      1. 'xp' and 'yp' can be set to zero if sub-arcsecond accuracy is
      not needed.  'ra' and 'dec' can be obtained from functions
      'tpstar' or 'tpplan'.
      2. The directions 'zd'= 0 (zenith) and 'az'= 0 (north) are here
      considered fixed in the terrestrial system.  Specifically, the
      zenith is along the geodetic normal, and north is toward
      the ITRS pole.
      3. If 'ref_option'= 0, then 'rar'='ra' and 'decr'='dec'.
      4. This function is the C version of NOVAS Fortran routine
      'zdaz'.

------------------------------------------------------------------------
*/
{
   short int j;

   double sinlat, coslat, sinlon, coslon, sindc, cosdc, sinra, cosra,
      uze[3], une[3], uwe[3], uz[3], un[3], uw[3], p[3], pz, pn, pw,
      proj, zd0, zd1, refr, sinzd, coszd, sinzd0, coszd0, pr[3];

/*
   Preliminaries.
*/

   *rar = ra;
   *decr = dec;

   sinlat = sin (location->latitude * DEG2RAD);
   coslat = cos (location->latitude * DEG2RAD);
   sinlon = sin (location->longitude * DEG2RAD);
   coslon = cos (location->longitude * DEG2RAD);
   sindc = sin (dec * DEG2RAD);
   cosdc = cos (dec * DEG2RAD);
   sinra = sin (ra * 15.0 * DEG2RAD);
   cosra = cos (ra * 15.0 * DEG2RAD);

/*
   Set up orthonormal basis vectors in local Earth-fixed system.

   Define vector toward local zenith in Earth-fixed system (z axis).
*/
   uze[0] = coslat * coslon;
   uze[1] = coslat * sinlon;
   uze[2] = sinlat;

/*
   Define vector toward local north in Earth-fixed system (x axis).
*/

   une[0] = -sinlat * coslon;
   une[1] = -sinlat * sinlon;
   une[2] = coslat;

/*
   Define vector toward local west in Earth-fixed system (y axis).
*/

   uwe[0] = sinlon;
   uwe[1] = -coslon;
   uwe[2] = 0.0;

/*
   Obtain vectors in celestial system.

   Rotate Earth-fixed orthonormal basis vectors to celestial system
   (wrt equator and equinox of date).
*/

   ter2cel (jd_ut1,0.0,delta_t,1,accuracy,1,xp,yp,uze, uz);
   ter2cel (jd_ut1,0.0,delta_t,1,accuracy,1,xp,yp,une, un);
   ter2cel (jd_ut1,0.0,delta_t,1,accuracy,1,xp,yp,uwe, uw);

/*
   Define unit vector 'p' toward object in celestial system
   (wrt equator and equinox of date).
*/

   p[0] = cosdc * cosra;
   p[1] = cosdc * sinra;
   p[2] = sindc;

/*
   Compute coordinates of object wrt orthonormal basis.

   Compute components of 'p' - projections of 'p' onto rotated
   Earth-fixed basis vectors.
*/

   pz = p[0] * uz[0] + p[1] * uz[1] + p[2] * uz[2];
   pn = p[0] * un[0] + p[1] * un[1] + p[2] * un[2];
   pw = p[0] * uw[0] + p[1] * uw[1] + p[2] * uw[2];

/*
   Compute azimuth and zenith distance.
*/

   proj = sqrt (pn * pn + pw * pw);

   if (proj > 0.0)
      *az = -atan2 (pw, pn) * RAD2DEG;

   if (*az < 0.0)
      *az += 360.0;

   if (*az >= 360.0)
      *az -= 360.0;

   *zd = atan2 (proj, pz) * RAD2DEG;

/*
   Apply atmospheric refraction if requested.
*/

   if (ref_option != 0)
   {

/*
   Get refraction in zenith distance.

   Iterative process is required because refraction algorithms are
   always a function of observed (not computed) zenith distance.
   Require convergence to 0.1 arcsec (actual accuracy less).
*/

      zd0 = *zd;

      do
      {
         zd1 = *zd;
         refr = refract (location,ref_option,*zd);
         *zd = zd0 - refr;
      } while (fabs (*zd - zd1) > 3.0e-5);

/*
   Apply refraction to celestial coordinates of object.
*/

      if ((refr > 0.0) && (*zd > 3.0e-4))
      {

/*
   Shift position vector of object in celestial system to account
   for refraction (see USNO/AA Technical Note 1998-09).
*/

         sinzd = sin (*zd * DEG2RAD);
         coszd = cos (*zd * DEG2RAD);
         sinzd0 = sin (zd0 * DEG2RAD);
         coszd0 = cos (zd0 * DEG2RAD);

/*
   Compute refracted position vector.
*/

         for (j = 0; j < 3; j++)
            pr[j] = ((p[j] - coszd0 * uz[j]) / sinzd0) * sinzd + uz[j] *
               coszd;

/*
   Compute refracted right ascension and declination.
*/

         proj = sqrt (pr[0] * pr[0] + pr[1] * pr[1]);

         if (proj > 0.0)
           *rar = atan2 (pr[1],pr[0]) * RAD2DEG / 15.0;

         if (*rar < 0.0)
           *rar += 24.0;

         if (*rar >= 24.0)
           *rar -= 24.0;

         *decr = atan2 (pr[2],proj) * RAD2DEG;
      }
   }
   return;
}

/********gcrs2equ */

short int gcrs2equ (double jd_tt, short int coord_sys,
                    short int accuracy, double rag, double decg,

                    double *ra, double *dec)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function converts GCRS right ascension and declination
      to coordinates with respect to the equator of date (mean or true).
      For coordinates with respect to the true equator of date, the
      origin of right ascension can be either the true equinox or the
      celestial intermediate origin (CIO).

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date of equator to be used for output coordinates.
      coord_sys (short int)
         Coordinate system selection for output coordinates.
            = 0 ... mean equator and equinox of date
            = 1 ... true equator and equinox of date
            = 2 ... true equator and CIO of date
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      rag (double)
         GCRS right ascension in hours.
      decg (double)
         GCRS declination in degrees.

   OUTPUT
   ARGUMENTS:
      ra (double)
         Right ascension in hours, referred to specified equator and
         right ascension origin of date.
      dec (double)
         Declination in degrees, referred to specified equator of date.


   RETURNED
   VALUE:
      (short int)
         =  0 ... everything OK.
         <  0 ... error from function 'vector2radec'
         > 10 ... 10 + error from function 'cio_location'
         > 20 ... 20 + error from function 'cio_basis'

   GLOBALS
   USED:
      DEG2RAD, T0        novascon.c

   FUNCTIONS
   CALLED:
      tdb2tt             novas.c
      frame_tie          novas.c
      precession         novas.c
      nutation           novas.c
      cio_location       novas.c
      cio_basis          novas.c
      vector2radec       novas.c
      sin                math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/09-04/JAB (USNO/AA)
      V1.1/01-06/WKP (USNO/AA) Changed 'mode' to 'accuracy'.
      V1.2/07-06/JAB (USNO/AA) Implement 'cio_location' construct.
      V1.3/04-08/WKP (USNO/AA) Updated variable names.

   NOTES:
      1. Set input value of 'accuracy' equal to any short int if
      'coord_sys' equals 0 or 1.  It is not used in these cases.
      2. This function only supports the CIO-based method.
      3. This function is the C version of NOVAS Fortran routine
      'gcrseq'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int rs;

   double t1, t, secdiff, r, d, pos1[3], pos2[3], pos3[3], pos4[3],
      r_cio, x[3], y[3], z[3];

/*
   't1' is the TDB Julian date.
*/

   tdb2tt (jd_tt, &t,&secdiff);
   t1 = jd_tt + secdiff / 86400.0;

/*
   Form position vector in equatorial system from input coordinates.
*/

   r = rag * 15.0 * DEG2RAD;
   d = decg * DEG2RAD;

   pos1[0] = cos (d) * cos (r);
   pos1[1] = cos (d) * sin (r);
   pos1[2] = sin (d);

/*
   Transform the position vector based on the value of 'coord_sys'.
*/

   if (coord_sys <= 1 )
   {

/*
   Transform the position vector from GCRS to mean equator and equinox
   of date.
*/

      frame_tie (pos1,1, pos2);
      precession (T0,pos2,t1, pos3);

/*
   If requested, transform further to true equator and equinox of date.
*/

          if (coord_sys == 1)
          {
             nutation (t1,0,accuracy,pos3, pos4);
          }
           else
          {
             pos4[0] = pos3[0];
             pos4[1] = pos3[1];
             pos4[2] = pos3[2];
          }
    }
     else

    {

/*
   Obtain the basis vectors, in the GCRS, of the celestial intermediate
   system.
*/

      if ((error = cio_location (t1,accuracy, &r_cio,&rs)) != 0)
         return (error += 10);

      if ((error = cio_basis (t1,r_cio,rs,accuracy, x,y,z)) != 0)
         return (error += 20);

/*
   Transform position vector to the celestial intermediate system
   (which has the CIO as its origin of right ascension).
*/

      pos4[0] = x[0] * pos1[0] + x[1] * pos1[1] + x[2] * pos1[2];
      pos4[1] = y[0] * pos1[0] + y[1] * pos1[1] + y[2] * pos1[2];
      pos4[2] = z[0] * pos1[0] + z[1] * pos1[1] + z[2] * pos1[2];
   }

/*
   Convert the position vector to equatorial spherical coordinates.
*/

   if ((error = vector2radec (pos4, ra,dec)) != 0)
   {
      *ra = 0.0;
      *dec = 0.0;
      return (-1 * error);
   }

   return (error);
}

/********sidereal_time */

short int sidereal_time (double jd_high, double jd_low,
                         double delta_t,short int gst_type,
                         short int method, short int accuracy,

                         double *gst)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the Greenwich sidereal time, either mean or apparent, at
      Julian date 'jd_high' + 'jd_low'.

   REFERENCES:
      Kaplan, G. (2005), US Naval Observatory Circular 179.

   INPUT
   ARGUMENTS:
      jd_high (double)
         High-order part of UT1 Julian date.
      jd_low (double)
         Low-order part of UT1 Julian date.
      delta_t (double)
         Difference TT-UT1 at 'jd_high'+'jd_low', in seconds
         of time.
      gst_type (short int)
         = 0 ... compute Greenwich mean sidereal time
         = 1 ... compute Greenwich apparent sidereal time
      method (short int)
         Selection for method
            = 0 ... CIO-based method
            = 1 ... equinox-based method
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *gst (double)
         Greenwich (mean or apparent) sidereal time, in hours.

   RETURNED
   VALUE:
      (short int)
         = 0         ... everything OK
         = 1         ... invalid value of 'accuracy'
         = 2         ... invalid value of 'method'
         > 10, < 30  ... 10 + error from function 'cio_rai'

   GLOBALS
   USED:
      T0, RAD2DEG        novascon.c

   FUNCTIONS
   CALLED:
      tdb2tt             novas.c
      era                novas.c
      e_tilt             novas.c
      cio_location       novas.c
      cio_basis          novas.c
      nutation           novas.c
      precession         novas.c
      frame_tie          novas.c
      fabs               math.h
      atan2              math.h
      fmod               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-92/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C programing standards.
      V1.2/03-98/JAB (USNO/AA) Expand documentation.
      V1.3/08-98/JAB (USNO/AA) Match flow of the Fortran counterpart.
      V2.0/09-03/JAB (USNO/AA) Incorporate the 2nd-reference changes.
      V2.1/08-04/JAB (USNO/AA) Incorporate the 1st-reference changes.
      V2.2/12-05/WKP (USNO/AA) Updated error handling.
      V2.3/01-06/WKP (USNO/AA) Changed 'mode' to 'method' and 'accuracy'.
      V2.4/04-06/JAB (USNO/AA) Use precession-in-RA terms in mean
                               sidereal time from third reference.
      V2.5/07-06/JAB (USNO/AA) Implement 'cio_location' construct.
      V2.6/06-08/WKP (USNO/AA) Changed value of direction argument in
                               call to 'nutation' from 1 to -1 for
                               consistency.
      V2.7/03-11/WKP (USNO/AA) Updated prolog description to clarify
                               this function computes either mean or
                               apparent sidereal time, and removed
                               Note 1 for consistency with Fortran.

   NOTES:
      1. The Julian date may be split at any point, but for highest
      precision, set 'jd_high' to be the integral part of the Julian
      date, and set 'jd_low' to be the fractional part.
      2. This function is the C version of NOVAS Fortran routine
      'sidtim'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int ref_sys;

   static double ee;
   static double jd_last = -99.0;
   double unitx[3] = {1.0, 0.0, 0.0};
   double jd_ut, jd_tt, jd_tdb, tt_temp, t, theta, a, b, c, d,
      ra_cio, x[3], y[3], z[3], w1[3], w2[3], eq[3], ha_eq, st,
      secdiff, eqeq;

/*
   Invalid value of 'accuracy'.
*/

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 1);

/*
   Time argument for precession and nutation components of sidereal
   time is TDB.  First approximation is TDB = TT, then refine.
*/

   jd_ut = jd_high + jd_low;
   jd_tt = jd_ut + (delta_t / 86400.0);
   jd_tdb = jd_tt;
   tdb2tt (jd_tdb, &tt_temp,&secdiff);
   jd_tdb = jd_tt + (secdiff / 86400.0);

   t = (jd_tdb - T0) / 36525.0;

/*
   Compute the Earth Rotation Angle.  Time argument is UT1.
*/

   theta = era (jd_high, jd_low);

/*
   Compute the equation of the equinoxes if needed, depending upon the
   input values of 'gst_type' and 'method'.  If not needed, set to zero.
*/

   if (((gst_type == 0) && (method == 0)) ||       /* GMST; CIO-TIO */
       ((gst_type == 1) && (method == 1)))         /* GAST; equinox */
   {
      if (fabs (jd_tdb - jd_last) > 1.0e-8)
      {
         e_tilt (jd_tdb,accuracy, &a,&b,&ee,&c,&d);
         jd_last = jd_tdb;
      }
      eqeq = ee * 15.0;
   }
    else
   {
      eqeq = 0.0;
   }

/*
   Compute Greenwich sidereal time depending upon input values of
   'method' and 'gst_type'.
*/

   switch (method)
   {
      case (0):

/*
   Use 'CIO-TIO-theta' method.  See Circular 179, Section 6.5.4.
*/

/*
   Obtain the basis vectors, in the GCRS, of the celestial intermediate
   system.
*/

         if ((error = cio_location (jd_tdb,accuracy, &ra_cio,
            &ref_sys)) != 0)
         {
            *gst = 99.0;
            return (error += 10);
         }

         cio_basis (jd_tdb,ra_cio,ref_sys,accuracy, x,y,z);

/*
   Compute the direction of the true equinox in the GCRS.
*/

         nutation (jd_tdb,-1,accuracy,unitx, w1);
         precession (jd_tdb,w1,T0, w2);
         frame_tie (w2,-1, eq);

/*
   Compute the hour angle of the equinox wrt the TIO meridian
   (near Greenwich, but passes through the CIP and TIO).
*/

         ha_eq = theta - atan2 ((eq[0] * y[0] + eq[1] * y[1] +
            eq[2] * y[2]), (eq[0] * x[0] + eq[1] * x[1] +
            eq[2] * x[2])) * RAD2DEG;

/*
   For mean sidereal time, subtract the equation of the equinoxes.
*/

         ha_eq -= (eqeq / 240.0);

         ha_eq = fmod (ha_eq, 360.0) / 15.0;
         if (ha_eq < 0.0)
            ha_eq += 24.0;
         *gst = ha_eq;
         break;

      case (1):

/*
   Use equinox method.  See Circular 179, Section 2.6.2.
*/

/*
   Precession-in-RA terms in mean sidereal time taken from third
   reference, eq. (42), with coefficients in arcseconds.
*/

         st = eqeq + 0.014506 +
               (((( -    0.0000000368   * t
                    -    0.000029956  ) * t
                    -    0.00000044   ) * t
                    +    1.3915817    ) * t
                    + 4612.156534     ) * t;

/*
   Form the Greenwich sidereal time.
*/

         *gst = fmod ((st / 3600.0 + theta), 360.0) / 15.0;

         if (*gst < 0.0)
            *gst += 24.0;
         break;

/*
   Invalid value of 'method'.
*/

      default:
         *gst = 99.0;
         error = 2;
         break;
   }

   return (error);
}

/********era */

double era (double jd_high, double jd_low)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function returns the value of the Earth Rotation Angle
      (theta) for a given UT1 Julian date.  The expression used is
      taken from the note to IAU Resolution B1.8 of 2000.

   REFERENCES:
      IAU Resolution B1.8, adopted at the 2000 IAU General Assembly,
         Manchester, UK.
      Kaplan, G. (2005), US Naval Observatory Circular 179.

   INPUT
   ARGUMENTS:
      jd_high (double)
         High-order part of UT1 Julian date.
      jd_low (double)
         Low-order part of UT1 Julian date.

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (double)
         The Earth Rotation Angle (theta) in degrees.

   GLOBALS
   USED:
      T0                 novascon.c

   FUNCTIONS
   CALLED:
      fmod               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/09-03/JAB (USNO/AA)

   NOTES:
      1. The algorithm used here is equivalent to the canonical
      theta = 0.7790572732640 + 1.00273781191135448 * t,
      where t is the time in days from J2000 (t = jd_high +
      jd_low - T0), but it avoids many two-PI 'wraps' that decrease
      precision (adopted from SOFA Fortran routine iau_era00; see
      also expression at top of page 35 of IERS Conventions (1996)).
      2. This function is the C version of NOVAS Fortran routine
      'erot'.

------------------------------------------------------------------------
*/

{
   double theta, thet1, thet2, thet3;

   thet1 = 0.7790572732640 + 0.00273781191135448 * (jd_high - T0);
   thet2 = 0.00273781191135448 * jd_low;
   thet3 = fmod (jd_high, 1.0) + fmod (jd_low, 1.0);

   theta = fmod (thet1 + thet2 + thet3, 1.0) * 360.0;
   if (theta < 0.0)
      theta += 360.0;

   return theta;
}

/********ter2cel */

short int ter2cel (double jd_ut_high, double jd_ut_low, double delta_t,
                   short int method, short int accuracy, short int option,
                   double xp, double yp, double *vec1,

                   double *vec2)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function rotates a vector from the terrestrial to the
      celestial system.  Specifically, it transforms a vector in the
      ITRS (rotating earth-fixed system) to the GCRS (a local space-
      fixed system) by applying rotations for polar motion, Earth
      rotation, nutation, precession, and the dynamical-to-GCRS
      frame tie.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Kaplan, G. H. (2003), 'Another Look at Non-Rotating Origins',
         Proceedings of IAU XXV Joint Discussion 16.

   INPUT
   ARGUMENTS:
      jd_ut_high (double)
         High-order part of UT1 Julian date.
      jd_ut_low (double)
         Low-order part of UT1 Julian date.
      delta_t (double)
         Value of Delta T (= TT - UT1) at the input UT1 Julian date.
      method (short int)
         Selection for method
            = 0 ... CIO-based method
            = 1 ... equinox-based method
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      option (short int)
            = 0 ... The output vector is referred to GCRS axes.
            = 1 ... The output vector is produced with respect to the
                    equator and equinox of date.
                    See Note 2 below.
      xp (double)
         Conventionally-defined X coordinate of celestial intermediate
         pole with respect to ITRS pole, in arcseconds.
      yp (double)
         Conventionally-defined Y coordinate of celestial intermediate
         pole with respect to ITRS pole, in arcseconds.
      vec1[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to ITRS axes (terrestrial system) in the normal case
         where 'option' = 0.

   OUTPUT
   ARGUMENTS:
      vec2[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to GCRS axes (celestial system) or with respect to
         the equator and equinox of date, depending on 'option'.

   RETURNED
   VALUE:
      =  0  ... everything is ok.
      =  1  ... invalid value of 'accuracy'
      =  2  ... invalid value of 'method'
      > 10 ... 10 + error from function 'cio_location'
      > 20 ... 20 + error from function 'cio_basis'

   GLOBALS
   USED:
      T0                 novascon.c

   FUNCTIONS
   CALLED:
      tdb2tt             novas.c
      wobble             novas.c
      cio_location       novas.c
      cio_basis          novas.c
      era                novas.c
      spin               novas.c
      sidereal_time      novas.c
      nutation           novas.c
      precession         novas.c
      frame_tie          novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/08-93/WTH (USNO/AA) Translate Fortran.
      V2.0/10-03/JAB (USNO/AA) Formerly function 'pnsw.'  Update
                               for IAU 2000 resolutions.
      V2.1/01-05/JAB (USNO/AA) Update to include 'mode' flexibility.
      V2.2/12-05/WKP (USNO/AA) Check error codes.
      V2.3/01-06/WKP (USNO/AA) Changed 'mode' to 'method' and 'accuracy'.
      V2.4/07-06/JAB (USNO/AA) Implement 'cio_location' construct.
      V2.5/06-08/WKP (USNO/AA) Corrected sign of direction argument in
                               'frame_tie' call and changed value of
                               direction argument in call to 'nutation'
                               for consistency.
      V2.6/10-08/WKP (USNO/AA) Renamed input JD variables to include
                               time system.
      V2.7/03-09/JAB (USNO/AA) Conformed input variables to IERS
                               conventions.
      V2.8/12-10/JAB (USNO/AA) Cosmetic changes for consistency with
                               Fortran.

   NOTES:
      1. 'xp' = 'yp' = 0 means no polar motion transformation.
      2. The 'option' flag only works for the equinox-based method.
      3. This function is the C version of NOVAS Fortran routine
      'tercel'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int rs, j;

   double jd_ut1, jd_tt, dummy, secdiff, jd_tdb, gast, r_cio, theta,
   v1[3], v2[3], v3[3], v4[3], x[3], y[3], z[3];

/*
   Invalid value of 'accuracy'.
*/

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 1);

/*
   Compute the TT Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_ut1 = jd_ut_high + jd_ut_low;
   jd_tt = jd_ut1 + (delta_t / 86400.0);

/*
   Compute the TDB Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_tdb = jd_tt;
   tdb2tt (jd_tdb, &dummy,&secdiff);
   jd_tdb = jd_tt + secdiff / 86400.0;

   switch (method)
   {
      case (0):

/*
   'CIO-TIO-THETA' method.

   See second reference, eq. (3) and (4).

   Apply polar motion, transforming the vector to the terrestrial
   intermediate system.
*/

         if ((xp == 0.0) && (yp == 0.0))
         {
            v1[0] = vec1[0];
            v1[1] = vec1[1];
            v1[2] = vec1[2];
         }
          else
            wobble (jd_tdb,0,xp,yp,vec1, v1);

/*
   Obtain the basis vectors, in the GCRS, of the celestial intermediate
   system.
*/

      if ((error = cio_location (jd_tdb,accuracy, &r_cio,&rs)) != 0)
         return (error += 10);

      if ((error = cio_basis (jd_tdb,r_cio,rs,accuracy, x,y,z)) != 0)
         return (error += 20);

/*
   Compute and apply the Earth rotation angle, 'theta', transforming the
   vector to the celestial intermediate system.
*/

         theta = era (jd_ut_high,jd_ut_low);
         spin (-theta,v1, v2);

/*
   Transform the vector from the celestial intermediate system to the
   GCRS.
*/

         vec2[0] = x[0] * v2[0] + y[0] * v2[1] + z[0] * v2[2];
         vec2[1] = x[1] * v2[0] + y[1] * v2[1] + z[1] * v2[2];
         vec2[2] = x[2] * v2[0] + y[2] * v2[1] + z[2] * v2[2];
         break;

      case (1):

/*
   Equinox mode.

   Apply polar motion.
*/

         if ((xp == 0.0) && (yp == 0.0))
         {
            for (j = 0; j < 3; j++)
            {
               v1[j] = vec1[j];
            }
         }
          else
            wobble (jd_tdb,0,xp,yp,vec1, v1);

/*
   Apply Earth rotation.
*/

         sidereal_time (jd_ut_high,jd_ut_low,delta_t,1,1,accuracy, &gast);
         spin (-gast * 15.0,v1, v2);

/*
   'option' = 1 skips remaining transformations.
*/

         if (option == 1)
         {
            vec2[0] = v2[0];
            vec2[1] = v2[1];
            vec2[2] = v2[2];
         }
          else
         {

/*
   Apply precession, nutation, and frame tie.
*/

            nutation (jd_tdb,-1,accuracy,v2, v3);
            precession (jd_tdb,v3,T0, v4);
            frame_tie (v4,-1, vec2);
         }
         break;

/*
   Invalid value of 'method'.
*/

      default:
         error = 2;
         break;
   }

   return (error);
}

/********cel2ter */

short int cel2ter (double jd_ut_high, double jd_ut_low, double delta_t,
                   short int method, short int accuracy, short int option,
                   double xp, double yp, double *vec1,

                   double *vec2)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function rotates a vector from the celestial to the
      terrestrial system.  Specifically, it transforms a vector in the
      GCRS (a local space-fixed system) to the ITRS (a rotating earth-
      fixed system) by applying rotations for the GCRS-to-dynamical
      frame tie, precession, nutation, Earth rotation, and polar motion.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Kaplan, G. H. (2003), 'Another Look at Non-Rotating Origins',
         Proceedings of IAU XXV Joint Discussion 16.

   INPUT
   ARGUMENTS:
      jd_ut_high (double)
         High-order part of UT1 Julian date.
      jd_ut_low (double)
         Low-order part of UT1 Julian date.
      delta_t (double)
         Value of Delta T (= TT - UT1) at the input UT1 Julian date.
      method (short int)
         Selection for method
            = 0 ... CIO-based method
            = 1 ... equinox-based method
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      option (short int)
            = 0 ... The input vector is referred to GCRS axes.
            = 1 ... The input vector is produced with respect to the
                    equator and equinox of date.
                    See Note 2 below.
      xp (double)
         Conventionally-defined X coordinate of celestial intermediate
         pole with respect to ITRS pole, in arcseconds.
      yp (double)
         Conventionally-defined Y coordinate of celestial intermediate
         pole with respect to ITRS pole, in arcseconds.
      vec1[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to GCRS axes (celestial system) or with respect to
         the equator and equinox of date, depending on 'option'.

   OUTPUT
   ARGUMENTS:
      vec2[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to ITRS axes (terrestrial system).

   RETURNED
   VALUE:
      =  0  ... everything is ok.
      =  1  ... invalid value of 'accuracy'
      =  2  ... invalid value of 'method'
      > 10 ... 10 + error from function 'cio_location'
      > 20 ... 20 + error from function 'cio_basis'

   GLOBALS
   USED:
      T0                 novascon.c

   FUNCTIONS
   CALLED:
      tdb2tt             novas.c
      wobble             novas.c
      cio_location       novas.c
      cio_basis          novas.c
      era                novas.c
      spin               novas.c
      sidereal_time      novas.c
      nutation           novas.c
      precession         novas.c
      frame_tie          novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/03-09/JAB (USNO/AA)
      V1.1/12-10/JAB (USNO/AA) Fix confusion with input/output vectors;
                               cosmetic changes to better match Fortran.

   NOTES:
      1. 'xp' = 'yp' = 0 means no polar motion transformation.
      2. 'option' = 1 only works for the equinox-based method
      ('method' = 1).
      3. This function is the C version of NOVAS Fortran routine
      'celter'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int rs, j;

   double jd_ut1, jd_tt, dummy, secdiff, jd_tdb, gast, r_cio, theta,
   v1[3], v2[3], v3[3], v4[3], x[3], y[3], z[3];

/*
   Invalid value of 'accuracy'.
*/

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 1);

/*
   Compute the TT Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_ut1 = jd_ut_high + jd_ut_low;
   jd_tt = jd_ut1 + (delta_t / 86400.0);

/*
   Compute the TDB Julian date corresponding to the input UT1 Julian
   date.
*/

   jd_tdb = jd_tt;
   tdb2tt (jd_tdb, &dummy,&secdiff);
   jd_tdb = jd_tt + secdiff / 86400.0;

   switch (method)
   {
      case (0):

/*
   'CIO-TIO-THETA' method.

   See second reference, eq. (3) and (4).

   Obtain the basis vectors, in the GCRS, of the celestial intermediate
   system.
*/

         if ((error = cio_location (jd_tdb,accuracy, &r_cio,&rs)) != 0)
            return (error += 10);

         if ((error = cio_basis (jd_tdb,r_cio,rs,accuracy, x,y,z)) != 0)
            return (error += 20);

/*
   Transform the vector from the GCRS to the celestial intermediate
   system.
*/

         v1[0] = x[0] * vec1[0] + x[1] * vec1[1] + x[2] * vec1[2];
         v1[1] = y[0] * vec1[0] + y[1] * vec1[1] + y[2] * vec1[2];
         v1[2] = z[0] * vec1[0] + z[1] * vec1[1] + z[2] * vec1[2];

/*
   Compute and apply the Earth rotation angle, 'theta', transforming the
   vector to the terrestrial intermediate system.
*/

         theta = era (jd_ut_high,jd_ut_low);
         spin (theta,v1, v2);

/*
   Apply polar motion, transforming the vector to the ITRS.
*/

         if ((xp == 0.0) && (yp == 0.0))
         {
            vec2[0] = v2[0];
            vec2[1] = v2[1];
            vec2[2] = v2[2];
         }
          else
            wobble (jd_tdb,1,xp,yp,v2, vec2);

         break;

      case (1):

/*
   Equinox mode.

   'option' = 1 skips initial transformations.
*/

         if (option == 1)
         {
            v3[0] = vec1[0];
            v3[1] = vec1[1];
            v3[2] = vec1[2];
         }
          else
         {

/*
   Apply frame tie, nutation, and precession.
*/

            frame_tie (vec1,1, v1);
            precession (T0,v1,jd_tdb, v2);
            nutation (jd_tdb,0,accuracy,v2, v3);
         }

/*
   Apply Earth rotation.
*/

         sidereal_time (jd_ut_high,jd_ut_low,delta_t,1,1,accuracy, &gast);
         spin (gast * 15.0,v3, v4);

/*
   Apply polar motion.
*/

         if ((xp == 0.0) && (yp == 0.0))
         {
            for (j = 0; j < 3; j++)
            {
               vec2[j] = v4[j];
            }
         }
          else
            wobble (jd_tdb,1,xp,yp,v4, vec2);

         break;

/*
   Invalid value of 'method'.
*/

      default:
         error = 2;
         break;
   }

   return (error);
}

/********spin */

void spin (double angle, double *pos1,

           double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function transforms a vector from one coordinate system
      to another with same origin and axes rotated about the z-axis.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.

   INPUT
   ARGUMENTS:
      angle (double)
         Angle of coordinate system rotation, positive counterclockwise
         when viewed from +z, in degrees.
      pos1[3] (double)
         Position vector.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector expressed in new coordinate system rotated
         about z by 'angle'.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      DEG2RAD            novascon.c

   FUNCTIONS
   CALLED:
      sin                math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/08-93/WTH (USNO/AA) Translate Fortran.
      V2.0/10-03/JAB (USNO/AA) Update for IAU 2000 resolutions.
      V2.1/01-05/JAB (USNO/AA) Generalize the function.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine 'spin'.

------------------------------------------------------------------------
*/
{
   static double ang_last = -999.0;
   static double xx, yx, zx, xy, yy, zy, xz, yz, zz;
   double angr, cosang, sinang;

   if (fabs (angle - ang_last) >= 1.0e-12)
   {
      angr = angle * DEG2RAD;
      cosang = cos (angr);
      sinang = sin (angr);

/*
   Rotation matrix follows.
*/

      xx =  cosang;
      yx =  sinang;
      zx =  0.0;
      xy =  -sinang;
      yy =  cosang;
      zy =  0.0;
      xz =  0.0;
      yz =  0.0;
      zz =  1.0;

      ang_last = angle;
   }

/*
   Perform rotation.
*/

   pos2[0] = xx * pos1[0] + yx * pos1[1] + zx * pos1[2];
   pos2[1] = xy * pos1[0] + yy * pos1[1] + zy * pos1[2];
   pos2[2] = xz * pos1[0] + yz * pos1[1] + zz * pos1[2];

   return;
}

/********wobble */

void wobble (double tjd, short int direction, double xp, double yp,
             double *pos1,

             double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      Corrects a vector in the ITRS (rotating Earth-fixed system)
      for polar motion, and also corrects the longitude origin
      (by a tiny amount) to the Terrestrial Intermediate Origin
      (TIO).  The ITRS vector is thereby transformed to the terrestrial
      intermediate system, based on the true (rotational) equator and
      TIO.  Because the true equator is the plane orthogonal to the
      direction of the Celestial Intermediate Pole (CIP), the components
      of the output vector are referred to z and x axes toward the CIP
      and TIO, respectively.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.
      Lambert & Bizouard (2002), Astronomy and Astrophysics 394,
         317-321.
   INPUT
   ARGUMENTS:
      tjd (double)
         TT or UT1 Julian date.
      direction (short int)
         Flag determining 'direction' of transformation;
            direction  = 0 transformation applied, ITRS to terrestrial
                           intermediate system
            direction != 0 inverse transformation applied, terrestrial
                           intermediate system to ITRS
      xp (double)
         Conventionally-defined X coordinate of Celestial Intermediate
         Pole with respect to ITRS pole, in arcseconds.
      yp (double)
         Conventionally-defined Y coordinate of Celestial Intermediate
         Pole with respect to ITRS pole, in arcseconds.
      pos1[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to ITRS axes.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to true equator and TIO.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      T0, ASEC2RAD       novascon.c

   FUNCTIONS
   CALLED:
      sin                math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/08-93/WTH (USNO/AA) Translate Fortran.
      V2.0/10-03/JAB (USNO/AA) Update for IAU 2000 resolutions.
      V2.1/01-05/JAB (USNO/AA) Include higher-order terms in the
                               rotation.
      V2.2/06-05/JAB (USNO/AA) Corrected sign of first term in
                               expression for 'yy' (should be -).
      V2.3/03-09/JAB (USNO/AA) Add 'direction' flag and inverse
                               transformation; conformed input
                               variables to IERS conventions.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'wobble'.

------------------------------------------------------------------------
*/
{
   double xpole, ypole, t, sprime, tiolon, sinx, cosx, siny, cosy, sinl,
      cosl, xx, yx, zx, xy, yy, zy, xz, yz, zz;

   xpole = xp * ASEC2RAD;
   ypole = yp * ASEC2RAD;

/*
   Compute approximate longitude of TIO, using eq. (10) of the second
   reference.

   Note that 'tiolon,' the longitude correction, is negligible for
   most astronomical purposes.
*/

   t = (tjd - T0) / 36525.0;

   sprime = -47.0e-6 * t;
   tiolon = -sprime * ASEC2RAD;

/*
   Compute elements of rotation matrix.
   Equivalent to R3(-s')R2(x)R1(y) as per IERS Conventions (2003).
*/

   sinx = sin (xpole);
   cosx = cos (xpole);
   siny = sin (ypole);
   cosy = cos (ypole);
   sinl = sin (tiolon);
   cosl = cos (tiolon);

   xx =  cosx * cosl;
   yx =  sinx * siny * cosl + cosy * sinl;
   zx = -sinx * cosy * cosl + siny * sinl;
   xy = -cosx * sinl;
   yy = -sinx * siny * sinl + cosy * cosl;
   zy =  sinx * cosy * sinl + siny * cosl;
   xz =  sinx;
   yz = -cosx * siny;
   zz =  cosx * cosy;

/*
   Perform rotation from ITRS to terrestrial intermediate system.
*/

   if (direction == 0)
   {
      pos2[0] = xx * pos1[0] + yx * pos1[1] + zx * pos1[2];
      pos2[1] = xy * pos1[0] + yy * pos1[1] + zy * pos1[2];
      pos2[2] = xz * pos1[0] + yz * pos1[1] + zz * pos1[2];
   }
    else

/*
   Perform rotation from terrestrial intermediate system to ITRS.
*/

   {
      pos2[0] = xx * pos1[0] + xy * pos1[1] + xz * pos1[2];
      pos2[1] = yx * pos1[0] + yy * pos1[1] + yz * pos1[2];
      pos2[2] = zx * pos1[0] + zy * pos1[1] + zz * pos1[2];
   }

   return;
}

/********terra */

void terra (on_surface *location, double st,

            double *pos, double *vel)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the position and velocity vectors of a terrestrial
      observer with respect to the center of the Earth.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.

   INPUT
   ARGUMENTS:
      *location (struct on_surface)
         Pointer to structure containing observer's location (defined
         in novas.h).
      st (double)
         Local apparent sidereal time at reference meridian in hours.

   OUTPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector of observer with respect to center of Earth,
         equatorial rectangular coordinates, referred to true equator
         and equinox of date, components in AU.
      vel[3] (double)
         Velocity vector of observer with respect to center of Earth,
         equatorial rectangular coordinates, referred to true equator
         and equinox of date, components in AU/day.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      AU_KM, ERAD, F     novascon.c
      ANGVEL, DEG2RAD    novascon.c

   FUNCTIONS
   CALLED:
      sin                math.h
      cos                math.h
      sqrt               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/04-93/WTH (USNO/AA):  Translate Fortran.
      V1.1/06-98/JAB (USNO/AA):  Move constants 'f' and 'omega' to
                                 file 'novascon.c'.
      V1.2/10-03/JAB (USNO/AA):  Updates Notes; removed call to 'pow'.
      V1.3/12-04/JAB (USNO/AA):  Update to use 'on_surface" structure.
      V1.4/09-09/WKP (USNO/AA):  Moved ht_km calculation from first_entry
                                 block.

   NOTES:
      1. If reference meridian is Greenwich and st=0, 'pos' is
      effectively referred to equator and Greenwich.
      2. This function ignores polar motion, unless the
      observer's longitude and latitude have been corrected for it,
      and variation in the length of day (angular velocity of earth).
      3. The true equator and equinox of date do not form an
      inertial system.  Therefore, with respect to an inertial system,
      the very small velocity component (several meters/day) due to
      the precession and nutation of the Earth's axis is not accounted
      for here.
      4. This function is the C version of NOVAS Fortran routine
      'terra'.

------------------------------------------------------------------------
*/
{
   static short int first_entry = 1;
   short int j;

   static double erad_km, ht_km;
   double df, df2, phi, sinphi, cosphi, c, s, ach, ash, stlocl, sinst,
      cosst;

   if (first_entry)
   {
      erad_km = ERAD / 1000.0;
      first_entry = 0;
   }

/*
   Compute parameters relating to geodetic to geocentric conversion.
*/

   df = 1.0 - F;
   df2 = df * df;

   phi = location->latitude * DEG2RAD;
   sinphi = sin (phi);
   cosphi = cos (phi);
   c = 1.0 / sqrt (cosphi * cosphi + df2 * sinphi * sinphi);
   s = df2 * c;
   ht_km = location->height / 1000.0;
   ach = erad_km * c + ht_km;
   ash = erad_km * s + ht_km;

/*
   Compute local sidereal time factors at the observer's longitude.
*/

   stlocl = (st * 15.0 + location->longitude) * DEG2RAD;
   sinst = sin (stlocl);
   cosst = cos (stlocl);

/*
   Compute position vector components in kilometers.
*/

   pos[0] = ach * cosphi * cosst;
   pos[1] = ach * cosphi * sinst;
   pos[2] = ash * sinphi;

/*
   Compute velocity vector components in kilometers/sec.
*/

   vel[0] = -ANGVEL * ach * cosphi * sinst;
   vel[1] =  ANGVEL * ach * cosphi * cosst;
   vel[2] =  0.0;

/*
   Convert position and velocity components to AU and AU/DAY.
*/

   for (j = 0; j < 3; j++)
   {
      pos[j] /= AU_KM;
      vel[j] /= AU_KM;
      vel[j] *= 86400.0;
   }

   return;
}

/********e_tilt */

void e_tilt (double jd_tdb, short int accuracy,

             double *mobl, double *tobl, double *ee, double *dpsi,
             double *deps)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes quantities related to the orientation of the Earth's
      rotation axis at Julian date 'jd_tdb'.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd_tdb (double)
         TDB Julian Date.
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *mobl (double)
         Mean obliquity of the ecliptic in degrees at 'jd_tdb'.
      *tobl (double)
         True obliquity of the ecliptic in degrees at 'jd_tdb'.
      *ee (double)
         Equation of the equinoxes in seconds of time at 'jd_tdb'.
      *dpsi (double)
         Nutation in longitude in arcseconds at 'jd_tdb'.
      *deps (double)
         Nutation in obliquity in arcseconds at 'jd_tdb'.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      PSI_COR, EPS_COR   novas.c
      T0, ASEC2RAD       novascon.c
      DEG2RAD            novascon.c

   FUNCTIONS
   CALLED:
      nutation_angles    novas.c
      ee_ct              novas.c
      mean_obliq         novas.c
      fabs               math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/08-93/WTH (USNO/AA) Translate Fortran.
      V1.1/06-97/JAB (USNO/AA) Incorporate IAU (1994) and IERS (1996)
                               adjustment to the "equation of the
                               equinoxes".
      V1.2/10-97/JAB (USNO/AA) Implement function that computes
                               arguments of the nutation series.
      V1.3/07-98/JAB (USNO/AA) Use global variables 'PSI_COR' and
                               'EPS_COR' to apply celestial pole offsets
                               for high-precision applications.
      V2.0/10-03/JAB (USNO/AA) Update function for IAU 2000 resolutions.
      V2.1/12-04/JAB (USNO/AA) Add 'mode' argument.
      V2.2/01-06/WKP (USNO/AA) Changed 'mode' to 'accuracy'.

   NOTES:
      1. Values of the celestial pole offsets 'PSI_COR' and 'EPS_COR'
      are set using function 'cel_pole', if desired.  See the prolog
      of 'cel_pole' for details.
      2. This function is the C version of NOVAS Fortran routine
      'etilt'.

------------------------------------------------------------------------
*/
{
   static short int accuracy_last = 0;
   short int acc_diff;

   static double jd_last = 0.0;
   static double dp, de, c_terms;
   double t, d_psi, d_eps, mean_ob, true_ob, eq_eq;

/*
   Compute time in Julian centuries from epoch J2000.0.
*/

   t = (jd_tdb - T0) / 36525.0;

/*
   Check for difference in accuracy mode from last call.
*/

   acc_diff = accuracy - accuracy_last;

/*
   Compute the nutation angles (arcseconds) if the input Julian date
   is significantly different from the last Julian date, or the
   accuracy mode has changed from the last call.
*/

   if (((fabs (jd_tdb - jd_last)) > 1.0e-8) || (acc_diff != 0))
   {
      nutation_angles (t,accuracy, &dp,&de);

/*
   Obtain complementary terms for equation of the equinoxes in
   arcseconds.
*/

      c_terms = ee_ct (jd_tdb,0.0,accuracy) / ASEC2RAD;

/*
   Reset the values of the last Julian date and last mode.
*/

      jd_last = jd_tdb;
      accuracy_last = accuracy;
   }

/*
   Apply observed celestial pole offsets.
*/

   d_psi = dp + PSI_COR;
   d_eps = de + EPS_COR;

/*
   Compute mean obliquity of the ecliptic in arcseconds.
*/

   mean_ob = mean_obliq (jd_tdb);

/*
   Compute true obliquity of the ecliptic in arcseconds.
*/

   true_ob = mean_ob + d_eps;

/*
   Convert obliquity values to degrees.
*/

   mean_ob /= 3600.0;
   true_ob /= 3600.0;

/*
   Compute equation of the equinoxes in seconds of time.
*/

   eq_eq = d_psi * cos (mean_ob * DEG2RAD) + c_terms;
   eq_eq /= 15.0;

/*
   Set output values.
*/

   *dpsi = d_psi;
   *deps = d_eps;
   *ee   = eq_eq;
   *mobl = mean_ob;
   *tobl = true_ob;

   return;
}

/********cel_pole */

short int cel_pole (double tjd, short int type, double dpole1,
                    double dpole2)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function allows for the specification of celestial pole
      offsets for high-precision applications.  Each set of offsets is
      a correction to the modeled position of the pole for a specific
      date, derived from observations and published by the IERS.

   REFERENCES:
      Kaplan, G. (2005), US Naval Observatory Circular 179.
      Kaplan, G. (2003), USNO/AA Technical Note 2003-03.

   INPUT
   ARGUMENTS:
      tjd (double)
         TDB or TT Julian date for pole offsets.
      type (short int)
         Type of pole offset
            = 1 for corrections to angular coordinates of modeled pole
                referred to mean ecliptic of date, that is,
                delta-delta-psi and delta-delta-epsilon.
            = 2 for corrections to components of modeled pole unit
                vector referred to GCRS axes, that is, dx and dy.
      dpole1 (double)
         Value of celestial pole offset in first coordinate,
         (delta-delta-psi or dx) in milliarcseconds.
      dpole2 (double)
         Value of celestial pole offset in second coordinate,
         (delta-delta-epsilon or dy) in milliarcseconds.

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (short int)
          = 0 ... Everything OK.
          = 1 ... Invalid value of 'type'.

   GLOBALS
   USED:
      PSI_COR, EPS_COR   novas.c
      T0, ASEC2RAD       novascon.c

   FUNCTIONS
   CALLED:
      frame_tie          novas.c
      precession         novas.c
      mean_obliq         novas.c
      sin                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/12-04/JAB (USNO/AA)

   NOTES:
      1. This function sets the values of global variables 'PSI_COR'
      and 'EPS_COR' declared at the top of file 'novas.c'.  These
      global variables are used only in NOVAS function 'e_tilt'.
      2. This function, if used, should be called before any other
      NOVAS functions for a given date.  Values of the pole offsets
      specified via a call to this function will be used until
      explicitly changed.
      3. 'tjd' is used only for 'type' = 2, to transform dx and dy to
      the equivalent delta-delta-psi and delta-delta-epsilon values.
      4. For 'type' = 2, dx and dy are unit vector component
      corrections, but are expressed in milliarcseconds simply by
      multiplying by 206264806, the number of milliarcseconds in one
      radian.
      5. This function is the C version of NOVAS Fortran routine
      'celpol'.

------------------------------------------------------------------------
*/
{
   short int error = 0;

   double dx, dy, t, mean_ob, sin_e, x, dz, dp1[3], dp2[3], dp3[3];

   switch (type)
   {
      case (1):

/*
   Angular coordinates of modeled pole referred to mean ecliptic of
   date, that is,delta-delta-psi and delta-delta-epsilon.
*/

          PSI_COR = dpole1 * 1.0e-3;
          EPS_COR = dpole2 * 1.0e-3;
          break;

      case (2):

/*
   Components of modeled pole unit vector referred to GCRS axes, that
   is, dx and dy.
*/

          dx = dpole1;
          dy = dpole2;

          t = (tjd - T0) / 36525.0;

/*
   Compute sin_e of mean obliquity of date.
*/

          mean_ob = mean_obliq (tjd);
          sin_e = sin (mean_ob * ASEC2RAD);

/*
   The following algorithm, to transform dx and dy to
   delta-delta-psi and delta-delta-epsilon, is from eqs. (7)-(9) of the
   second reference.

   Trivial model of pole trajectory in GCRS allows computation of dz.
*/

          x = (2004.190 * t) * ASEC2RAD;
          dz = - (x + 0.5 * x * x * x) * dx;

/*
   Form pole offset vector (observed - modeled) in GCRS.
*/

          dp1[0] = dx * 1.0e-3 * ASEC2RAD;
          dp1[1] = dy * 1.0e-3 * ASEC2RAD;
          dp1[2] = dz * 1.0e-3 * ASEC2RAD;

/*
   Precess pole offset vector to mean equator and equinox of date.
*/

          frame_tie (dp1,1, dp2);
          precession (T0,dp2,tjd, dp3);

/*
   Compute delta-delta-psi and delta-delta-epsilon in arcseconds.
*/

          PSI_COR = (dp3[0] / sin_e) / ASEC2RAD;
          EPS_COR = dp3[1] / ASEC2RAD;
          break;

      default:

/*
   Invalid value of 'type'.
*/
         error = 1;
         break;
   }

   return (error);
}

/********ee_ct */

double ee_ct (double jd_high, double jd_low, short int accuracy)
/*
------------------------------------------------------------------------

   PURPOSE:
      To compute the "complementary terms" of the equation of the
      equinoxes.

   REFERENCES:
      Capitaine, N., Wallace, P.T., and McCarthy, D.D. (2003). Astron. &
         Astrophys. 406, p. 1135-1149. Table 3.
      IERS Conventions (2010), Chapter 5, p. 60, Table 5.2e.
         (Table 5.2e presented in the printed publication is a truncated
         series. The full series, which is used in NOVAS, is available
         on the IERS Conventions Center website in file tab5.2e.txt.)
         ftp://tai.bipm.org/iers/conv2010/chapter5/

   INPUT
   ARGUMENTS:
      jd_high (double)
         High-order part of TT Julian date.
      jd_low (double)
         Low-order part of TT Julian date.
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      None

   RETURNED
   VALUE:
      (double)
         Complementary terms, in radians.

   GLOBALS
   USED:
      T0, ASEC2RAD       novascon.c
      TWOPI              novascon.c

   FUNCTIONS
   CALLED:
      norm_ang           novas.c
      fund_args          novas.c
      fmod               math.h
      sin                math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/09-03/JAB (USNO/AA)
      V1.1/12-04/JAB (USNO/AA) Added low-accuracy formula.
      V1.2/01-06/WKP (USNO/AA) Changed 'mode' to 'accuracy'.
      V1.3/11-10/JAB (USNO/AA) Updated reference and notes.
      V1.4/03-11/WKP (USNO/AA) Added braces to 2-D array initialization
                               to quiet gcc warnings.

   NOTES:
      1. The series used in this function was derived from the first
      reference.  This same series was also adopted for use in the IAU's
      Standards of Fundamental Astronomy (SOFA) software (i.e.,
      subroutine eect00.for and function eect00.c).
      2. The low-accuracy series used in this function is a simple
      implementation derived from the first reference, in which terms
      smaller than 2 microarcseconds have been omitted.
      3. This function is based on NOVAS Fortran routine 'eect2000',
      with the low-accuracy formula taken from NOVAS Fortran routine
      'etilt'.

------------------------------------------------------------------------
*/
{
   short int i, j;

   double t, fa[14], fa2[5], s0, s1, a, c_terms;

/*
   Argument coefficients for t^0.
*/

   const short int ke0_t[33][14] = {
      {0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0},
      {0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1},
      {0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0},
      {1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0}};

/*
   Argument coefficients for t^1.
*/

   const short int ke1[14] =
      {0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0};

/*
   Sine and cosine coefficients for t^0.
*/

   const double se0_t[33][2] = {
      {+2640.96e-6,          -0.39e-6},
      {  +63.52e-6,          -0.02e-6},
      {  +11.75e-6,          +0.01e-6},
      {  +11.21e-6,          +0.01e-6},
      {   -4.55e-6,          +0.00e-6},
      {   +2.02e-6,          +0.00e-6},
      {   +1.98e-6,          +0.00e-6},
      {   -1.72e-6,          +0.00e-6},
      {   -1.41e-6,          -0.01e-6},
      {   -1.26e-6,          -0.01e-6},
      {   -0.63e-6,          +0.00e-6},
      {   -0.63e-6,          +0.00e-6},
      {   +0.46e-6,          +0.00e-6},
      {   +0.45e-6,          +0.00e-6},
      {   +0.36e-6,          +0.00e-6},
      {   -0.24e-6,          -0.12e-6},
      {   +0.32e-6,          +0.00e-6},
      {   +0.28e-6,          +0.00e-6},
      {   +0.27e-6,          +0.00e-6},
      {   +0.26e-6,          +0.00e-6},
      {   -0.21e-6,          +0.00e-6},
      {   +0.19e-6,          +0.00e-6},
      {   +0.18e-6,          +0.00e-6},
      {   -0.10e-6,          +0.05e-6},
      {   +0.15e-6,          +0.00e-6},
      {   -0.14e-6,          +0.00e-6},
      {   +0.14e-6,          +0.00e-6},
      {   -0.14e-6,          +0.00e-6},
      {   +0.14e-6,          +0.00e-6},
      {   +0.13e-6,          +0.00e-6},
      {   -0.11e-6,          +0.00e-6},
      {   +0.11e-6,          +0.00e-6},
      {   +0.11e-6,          +0.00e-6}};
/*
   Sine and cosine coefficients for t^1.
*/

   const double se1[2] =
      {   -0.87e-6,          +0.00e-6};

/*
   Interval between fundamental epoch J2000.0 and current date.
*/

      t = ((jd_high - T0) + jd_low) / 36525.0;

/*
   High accuracy mode.
*/

   if (accuracy == 0)
   {

/*
   Fundamental Arguments.

   Mean Anomaly of the Moon.
*/

      fa[0] = norm_ang ((485868.249036 +
                         (715923.2178 +
                         (    31.8792 +
                         (     0.051635 +
                         (    -0.00024470)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (1325.0*t, 1.0) * TWOPI);

/*
   Mean Anomaly of the Sun.
*/

      fa[1] = norm_ang ((1287104.793048 +
                         (1292581.0481 +
                         (     -0.5532 +
                         (     +0.000136 +
                         (     -0.00001149)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (99.0*t, 1.0) * TWOPI);

/*
   Mean Longitude of the Moon minus Mean Longitude of the Ascending
   Node of the Moon.
*/

      fa[2] = norm_ang (( 335779.526232 +
                         ( 295262.8478 +
                         (    -12.7512 +
                         (     -0.001037 +
                         (      0.00000417)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (1342.0*t, 1.0) * TWOPI);

/*
   Mean Elongation of the Moon from the Sun.
*/

      fa[3] = norm_ang ((1072260.703692 +
                         (1105601.2090 +
                         (     -6.3706 +
                         (      0.006593 +
                         (     -0.00003169)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (1236.0*t, 1.0) * TWOPI);

/*
   Mean Longitude of the Ascending Node of the Moon.
*/

      fa[4] = norm_ang (( 450160.398036 +
                         (-482890.5431 +
                         (      7.4722 +
                         (      0.007702 +
                         (     -0.00005939)
                         * t) * t) * t) * t) * ASEC2RAD
                         + fmod (-5.0*t, 1.0) * TWOPI);

      fa[ 5] = norm_ang (4.402608842 + 2608.7903141574 * t);
      fa[ 6] = norm_ang (3.176146697 + 1021.3285546211 * t);
      fa[ 7] = norm_ang (1.753470314 +  628.3075849991 * t);
      fa[ 8] = norm_ang (6.203480913 +  334.0612426700 * t);
      fa[ 9] = norm_ang (0.599546497 +   52.9690962641 * t);
      fa[10] = norm_ang (0.874016757 +   21.3299104960 * t);
      fa[11] = norm_ang (5.481293872 +    7.4781598567 * t);
      fa[12] = norm_ang (5.311886287 +    3.8133035638 * t);
      fa[13] =          (0.024381750 +    0.00000538691 * t) * t;

/*
   Evaluate the complementary terms.
*/

      s0 = 0.0;
      s1 = 0.0;

      for (i = 32; i >= 0; i--)
      {
         a = 0.0;

         for (j = 0; j < 14; j++)
         {
            a += (double) ke0_t[i][j] * fa[j];
         }

         s0 += (se0_t[i][0] * sin (a) + se0_t[i][1] * cos (a));
      }

      a = 0.0;

      for (j = 0; j < 14; j++)
      {
         a += (double) (ke1[j]) * fa[j];
      }

      s1 += (se1[0] * sin (a) + se1[1] * cos (a));

      c_terms = (s0 + s1 * t);
   }

    else

/*
   Low accuracy mode: Terms smaller than 2 microarcseconds omitted.
*/

   {
      fund_args (t, fa2);
      c_terms =
          2640.96e-6 * sin (fa2[4])
         +  63.52e-6 * sin (2.0 * fa2[4])
         +  11.75e-6 * sin (2.0 * fa2[2] - 2.0 * fa2[3] + 3.0 * fa2[4])
         +  11.21e-6 * sin (2.0 * fa2[2] - 2.0 * fa2[3] +       fa2[4])
         -   4.55e-6 * sin (2.0 * fa2[2] - 2.0 * fa2[3] + 2.0 * fa2[4])
         +   2.02e-6 * sin (2.0 * fa2[2]                + 3.0 * fa2[4])
         +   1.98e-6 * sin (2.0 * fa2[2]                +       fa2[4])
         -   1.72e-6 * sin (3.0 * fa2[4])
         -   0.87e-6 * t * sin (fa2[4]);
   }

   return (c_terms *= ASEC2RAD);
}

/********frame_tie */

void frame_tie (double *pos1, short int direction,

                double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      To transform a vector from the dynamical reference system to the
      International Celestial Reference System (ICRS), or vice versa.
      The dynamical reference system is based on the dynamical mean
      equator and equinox of J2000.0.  The ICRS is based on the space-
      fixed ICRS axes defined by the radio catalog positions of several
      hundred extragalactic objects.

   REFERENCES:
      Hilton, J. and Hohenkerk, C. (2004), Astronomy and Astrophysics
         413, 765-770, eq. (6) and (8).
      IERS (2003) Conventions, Chapter 5.

   INPUT
   ARGUMENTS:
      pos1[3] (double)
         Position vector, equatorial rectangular coordinates.
      direction (short int)
         Set 'direction' < 0 for dynamical to ICRS transformation.
         Set 'direction' >= 0 for ICRS to dynamical transformation.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, equatorial rectangular coordinates.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      ASEC2RAD           novascon.c

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/09-03/JAB (USNO/AA)
      V1.1/02-06/WKP (USNO/AA) Added second-order corrections to diagonal
                               elements.

   NOTES:
      1. For geocentric coordinates, the same transformation is used
      between the dynamical reference system and the GCRS.
      2. This function is the C version of NOVAS Fortran routine 'frame'.

------------------------------------------------------------------------
*/
{
   static short int compute_matrix = 1;

/*
   'xi0', 'eta0', and 'da0' are ICRS frame biases in arcseconds taken
   from IERS (2003) Conventions, Chapter 5.
*/

   const double xi0  = -0.0166170;
   const double eta0 = -0.0068192;
   const double da0  = -0.01460;
   static double xx, yx, zx, xy, yy, zy, xz, yz, zz;

/*
   Compute elements of rotation matrix to first order the first time
   this function is called.  Elements will be saved for future use and
   not recomputed.
*/

   if (compute_matrix == 1)
   {
      xx =  1.0;
      yx = -da0  * ASEC2RAD;
      zx =  xi0  * ASEC2RAD;
      xy =  da0  * ASEC2RAD;
      yy =  1.0;
      zy =  eta0 * ASEC2RAD;
      xz = -xi0  * ASEC2RAD;
      yz = -eta0 * ASEC2RAD;
      zz =  1.0;

/*
   Include second-order corrections to diagonal elements.
*/

      xx = 1.0 - 0.5 * (yx * yx + zx * zx);
      yy = 1.0 - 0.5 * (yx * yx + zy * zy);
      zz = 1.0 - 0.5 * (zy * zy + zx * zx);

      compute_matrix = 0;
   }

/*
   Perform the rotation in the sense specified by 'direction'.
*/

   if (direction < 0)
   {

/*
   Perform rotation from dynamical system to ICRS.
*/

      pos2[0] = xx * pos1[0] + yx * pos1[1] + zx * pos1[2];
      pos2[1] = xy * pos1[0] + yy * pos1[1] + zy * pos1[2];
      pos2[2] = xz * pos1[0] + yz * pos1[1] + zz * pos1[2];
   }
    else
   {

/*
   Perform rotation from ICRS to dynamical system.
*/

      pos2[0] = xx * pos1[0] + xy * pos1[1] + xz * pos1[2];
      pos2[1] = yx * pos1[0] + yy * pos1[1] + yz * pos1[2];
      pos2[2] = zx * pos1[0] + zy * pos1[1] + zz * pos1[2];
   }

   return;
}

/********proper_motion */

void proper_motion (double jd_tdb1, double *pos, double *vel,
                    double jd_tdb2,

                    double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      Applies proper motion, including foreshortening effects, to a
      star's position.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.

   INPUT
   ARGUMENTS:
      jd_tdb1 (double)
         TDB Julian date of first epoch.
      pos[3] (double)
         Position vector at first epoch.
      vel[3] (double)
         Velocity vector at first epoch.
      jd_tdb2 (double)
         TDB Julian date of second epoch.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector at second epoch.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Updated to C programming standards.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'propmo'.

------------------------------------------------------------------------
*/
{
    short int j;

    for (j = 0; j < 3; j++)
    {
       pos2[j] = pos[j] + (vel[j] * (jd_tdb2 - jd_tdb1));
    }

    return;
}

/********bary2obs */

void bary2obs (double *pos, double *pos_obs,

               double *pos2, double *lighttime)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function moves the origin of coordinates from the barycenter
      of the solar system to the observer (or the geocenter); i.e.,
      this function accounts for parallax (annual+geocentric or just
      annual).

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.

   INPUT
   ARGUMENTS:
      pos1[3] (double)
         Position vector, referred to origin at solar system barycenter,
         components in AU.
      pos_obs[3] (double)
         Position vector of observer (or the geocenter), with respect
         to origin at solar system barycenter, components in AU.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, referred to origin at center of mass of the
         Earth, components in AU.
      *lighttime (double)
         Light time from object to Earth in days.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      C_AUDAY            novascon.c

   FUNCTIONS
   CALLED:
      sqrt               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.
      V1.2/10-03/JAB (USNO/AA) Remove calls to 'pow'.
      V2.0/01-05/JAB (USNO/AA) Update documentation; change
                               'earthvector' to 'pos_obs'.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'geocen'.

------------------------------------------------------------------------
*/
{
   short int j;

/*
   Translate vector to geocentric coordinates.
*/

   for (j = 0; j < 3; j++)
   {
      pos2[j] = pos[j] - pos_obs[j];
   }

/*
   Calculate length of vector in terms of light time.
*/

   *lighttime = sqrt (pos2[0] * pos2[0] + pos2[1] * pos2[1] +
       pos2[2] * pos2[2]) / C_AUDAY;

   return;
}

/********geo_posvel */

short int geo_posvel (double jd_tt, double delta_t, short int accuracy,
                      observer *obs,

                      double *pos, double *vel)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function computes the geocentric position and velocity
      of an observer on the surface of the earth or on a near-earth
      spacecraft.  The final vectors are expressed in the GCRS.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date.
      delta_t (double)
         Value of Delta T (= TT - UT1) at 'jd_tt'.
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      *obs (struct observer)
         Data specifying the location of the observer (struct defined
         in novas.h).

   OUTPUT
   ARGUMENTS:
      *pos (double)
         Position vector of observer, with respect to origin at
         geocenter, referred to GCRS axes, components in AU.
      *vel (double)
         Velocity vector of observer, with respect to origin at
         geocenter, referred to GCRS axes, components in AU/day.

   RETURNED
   VALUE:
      (short int)
         = 0 ... everything OK.
         = 1 ... invalid value of 'accuracy'.

   GLOBALS
   USED:
      AU_KM, T0          novascon.c

   FUNCTIONS
   CALLED:
      tdb2tt             novas.c
      sidereal_time      novas.c
      e_tilt             novas.c
      terra              novas.c
      nutation           novas.c
      precession         novas.c
      frame_tie          novas.c
      fabs               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/12-04/JAB (USNO/AA)
      V2.2/12-05/WKP (USNO/AA) Updated error handling.
      V2.3/01-06/WKP (USNO/AA) Changed 'mode' to 'accuracy'.
      V2.4/02-06/WKP (USNO/AA) Changed time argument in sidereal_time
                               call from 'jd_tt' to 'jd_ut1'.
      V2.5/04-06/WKP (USNO/AA) Corrected calculation of 'fac' and
                               changed it to a static.
      V2.6/10-06/JAB (USNO/AA) Remove 'jd_tdb' input argument; TDB
                               is approximated by TT.
      V2.7/02-07/JAB (USNO/AA) Compute 'jd_tdb' corresponding to input
                               'jd_tt'.


   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'geopos'.  This C function differs from its Fortran counterpart
      in that it does not use the input TT time as a substitute for
      the TDB time.

------------------------------------------------------------------------
*/
{
   static double t_last = 0;
   static double gast, fac;
   static short int first_time = 1;

   double x, secdif, gmst, x1, x2, x3, x4, eqeq, pos1[3], vel1[3],
      pos2[3], vel2[3], pos3[3], vel3[3], jd_tdb, jd_ut1;
   short int error = 0;

/*
   Invalid value of 'accuracy'.
*/

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 1);

/*
   Compute 'jd_tdb', the TDB Julian date corresponding to 'jd_tt'.
*/

   jd_tdb = jd_tt;
   tdb2tt (jd_tdb, &x,&secdif);
   jd_tdb = jd_tt + secdif / 86400.0;

   switch (obs->where)
   {

/*
   Observer at geocenter.  Trivial case.
*/

      case (0):
         pos[0] = 0.0;
         pos[1] = 0.0;
         pos[2] = 0.0;
         vel[0] = 0.0;
         vel[1] = 0.0;
         vel[2] = 0.0;
         return (error = 0);
         break;

/*
   Other two cases: Get geocentric position and velocity vectors of
   observer wrt equator and equinox of date.

   Observer on surface of Earth.
*/

      case (1):

/*
   Compute UT1 and sidereal time.
*/

         jd_ut1 = jd_tt - (delta_t / 86400.0);
         if (fabs (jd_ut1 - t_last) > 1.0e-8 )
         {
            sidereal_time (jd_ut1,0.0,delta_t,0,1,accuracy, &gmst);
            e_tilt (jd_tdb,accuracy, &x1,&x2,&eqeq,&x3,&x4);
            gast = gmst + eqeq / 3600.0;
            t_last = jd_ut1;
         }

/*
   Function 'terra' does the hard work, given sidereal time.
*/

         terra (&obs->on_surf,gast, pos1,vel1);
         break;

/*
   Observer on near-earth spacecraft.
*/

      case (2):

/*
   Convert units to AU and AU/day.
*/

         if (first_time)
         {
            fac = AU_KM / 86400.0;
            first_time = 0;
         }

         pos1[0] = obs->near_earth.sc_pos[0] / AU_KM;
         pos1[1] = obs->near_earth.sc_pos[1] / AU_KM;
         pos1[2] = obs->near_earth.sc_pos[2] / AU_KM;

         vel1[0] = obs->near_earth.sc_vel[0] / fac;
         vel1[1] = obs->near_earth.sc_vel[1] / fac;
         vel1[2] = obs->near_earth.sc_vel[2] / fac;
         break;
   }

/*
   Transform geocentric position vector of observer to GCRS.
*/

   nutation (jd_tdb,-1,accuracy,pos1, pos2);
   precession (jd_tdb,pos2,T0, pos3);
   frame_tie (pos3,-1, pos);

/*
   Transform geocentric velocity vector of observer to GCRS.
*/

   nutation (jd_tdb,-1,accuracy,vel1, vel2);
   precession (jd_tdb,vel2,T0, vel3);
   frame_tie (vel3,-1, vel);

   return (error = 0);
}

/********light_time */

short int light_time (double jd_tdb, object *ss_object, double pos_obs[3],
                      double tlight0, short int accuracy,

                      double pos[3], double *tlight)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function computes the geocentric position of a solar system
      body, as antedated for light-time.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd_tdb (double)
         TDB Julian date of observation.
      *ss_object (struct object)
         Pointer to structure containing the designation for the
         solar system body (defined in novas.h).
      pos_obs[3] (double)
         Position vector of observer (or the geocenter), with respect to
         origin at solar system barycenter, referred to ICRS axes,
         components in AU.
      tlight0 (double)
         First approximation to light-time, in days (can be set to 0.0
         if unknown).
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector of body, with respect to origin at observer (or
         the geocenter), referred to ICRS axes, components in AU.
      tlight (double)
         Final light-time, in days.

   RETURNED
   VALUE:
      (short int)
         =  0 ... everything OK.
         =  1 ... algorithm failed to converge after 10 iterations.
         > 10 ... error is 10 + error from function 'solarsystem'.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      ephemeris        novas.c
      bary2obs         novas.c
      fabs             math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/09-04/JAB (USNO/AA)
      V1.1/06-05/JAB (USNO/AA): Changed the convergence tolerance from
                                1.0e-8 to 1.0e-10
      V1.2/04-06/JAB (USNO/AA): Changed the convergence tolerance from
                                1.0e-10 to 1.0e-9
      V2.0/12-06/JAB (USNO/AA): Support use of split Julian date to
                                improve accuracy when the full accuracy
                                option is selected.
      V2.1/11-07/JAB (USNO/AA): Changed high-accuracy convergence
                                tolerance from 1.0e-11 to 1.0e-12.
      V2.2/10-08/JAB (USNO/AA): Use function 'ephemeris' instead of
                                direct calls to the solar system
                                functions.
      V2.3/09-10/WKP (USNO/AA): Initialized 't3' variable to silence
                                compiler warning.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'littim'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int iter = 0;

   double tol, jd[2], t1, t2, t3 = 0.0, pos1[3], vel1[3];

/*
   Set light-time convergence tolerance.  If full-accuracy option has
   been selected, split the Julian date into whole days + fraction of
   day.
*/

   if (accuracy == 0)
   {
      tol = 1.0e-12;

      jd[0] = (double) ((long int) jd_tdb);
      t1 = jd_tdb - jd[0];
      t2 = t1 - tlight0;
   }
    else
   {
      tol = 1.0e-9;

      jd[0] = 0.0;
      jd[1] = 0.0;
      t1 = jd_tdb;
      t2 = jd_tdb - tlight0;
   }

/*
   Iterate to obtain correct light-time (usually converges rapidly).
*/

   do
   {
      if (iter > 10)
      {
         error = 1;
         *tlight = 0.0;
         break;
      }

      if (iter > 0)
      {
         t2 = t3;
      }

      jd[1] = t2;
      error = ephemeris (jd,ss_object,0,accuracy, pos1, vel1);

      if (error != 0)
      {
         error += 10;
         *tlight = 0.0;
         break;
      }

      bary2obs (pos1,pos_obs, pos,tlight);

      t3 = t1 - *tlight;
      iter++;
   } while (fabs (t3 - t2) > tol);

   return (error);
}

/********d_light */

double d_light (double *pos1, double *pos_obs)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function returns the difference in light-time, for a star,
      between the barycenter of the solar system and the observer (or
      the geocenter).

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      pos1[3] (double)
         Position vector of star, with respect to origin at solar
         system barycenter.
      pos_obs[3] (double)
         Position vector of observer (or the geocenter), with respect
         to origin at solar system barycenter, components in AU.

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (double)
         Difference in light time, in the sense star to barycenter minus
         star to earth, in days.

   GLOBALS
   USED:
      C_AUDAY            novascon.c

   FUNCTIONS
   CALLED:
      sqrt               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/12-04/JAB (USNO/AA)

   NOTES:
      1.  Alternatively, this function returns the light-time from the
      observer (or the geocenter) to a point on a light ray that is
      closest to a specific solar system body.  For this purpose, 'pos1'
      is the position vector toward observed object, with respect to
      origin at observer (or the geocenter); 'pos_obs' is the position
      vector of solar system body, with respect to origin at observer
      (or the geocenter), components in AU; and the returned value is
      the light time to point on line defined by 'pos1' that is closest
      to solar system body (positive if light passes body before hitting
      observer, i.e., if 'pos1' is within 90 degrees of 'pos_obs').
      2. This function is the C version of NOVAS Fortran routine
      'dlight'.

------------------------------------------------------------------------
*/
{
   double dis, u1[3], diflt;

/*
   From 'pos1', form unit vector 'u1' in direction of star or light
   source.
*/

   dis = sqrt (pos1[0] * pos1[0] + pos1[1] * pos1[1] +
      pos1[2] * pos1[2]);

   u1[0] = pos1[0] / dis;
   u1[1] = pos1[1] / dis;
   u1[2] = pos1[2] / dis;

/*
   Light-time returned is the projection of vector 'pos_obs' onto the
   unit vector 'u1' (formed from 'pos1'), divided by the speed of light.
*/

   diflt = (pos_obs[0] * u1[0] + pos_obs[1] * u1[1] +
      pos_obs[2] * u1[2]) / C_AUDAY;

   return (diflt);
}

/********grav_def */

short int grav_def (double jd_tdb, short int loc_code,
                    short int accuracy, double *pos1, double *pos_obs,

                    double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function computes the total gravitational deflection of
      light for the observed object due to the major gravitating bodies
      in the solar system.  This function valid for an observed body
      within the solar system as well as for a star.

   REFERENCES:
      Klioner, S. (2003), Astronomical Journal 125, 1580-1597,
         Section 6.

   INPUT
   ARGUMENTS:
      jd_tdb (double)
         TDB Julian date of observation.
      loc_code (short int)
         Code for location of observer, determining whether the
         gravitational deflection due to the earth itself is applied.
            = 0 ... No earth deflection (normally means observer
                          is at geocenter)
            = 1 ... Add in earth deflection (normally means
         observer is on or above surface of earth, including earth
         orbit)
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      *pos1 (double)
         Position vector of observed object, with respect to origin at
         observer (or the geocenter), referred to ICRS axes, components
         in AU.
      *pos_obs (double)
         Position vector of observer (or the geocenter), with respect to
         origin at solar system barycenter, referred to ICRS axes,
         components in AU.

   OUTPUT
   ARGUMENTS:
      *pos2 (double)
         Position vector of observed object, with respect to origin at
         observer (or the geocenter), referred to ICRS axes, corrected
         for gravitational deflection, components in AU.

   RETURNED
   VALUE:
      (short int)
         =  0 ... Everything OK.
         < 30 ... Error from function 'ephemeris'.
         > 30 ... Error from function 'make_object'.

   GLOBALS
   USED:
      C_AUDAY            novascon.c
      RMASS[12]          novascon.c

   FUNCTIONS
   CALLED:
      make_cat_entry     novas.c
      make_object        novas.c
      ephemeris          novas.c
      bary2obs           novas.c
      d_light            novas.c
      grav_vec           novas.c
      sqrt               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-05/JAB (USNO/AA)
      V1.1/04-06/WKP (USNO/AA): Added second note (superceded).
      V1.2/10-08/JAB (USNO/AA): Substituted calls to 'ephemeris' for
                                calls to 'solarsystem'; added
                                actions based on 'accuracy' input.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'grvdef'.
      2. If 'accuracy' is set to zero (full accuracy), three bodies
      (Sun, Jupiter, and Saturn) are used in the calculation.  If
      the reduced-accuracy option is set, only the Sun is used in the
      calculation.  In both cases, if the observer is not at the
      geocenter, the deflection due to the Earth is included.
      3. The number of bodies used at full and reduced accuracy can be
      set by making a change to the code in this function as indicated
      in the comments.

------------------------------------------------------------------------
*/
{
/*
   The following list of body numbers and corresponding names identifies
   which gravitating bodies (aside from the Earth) are potentially
   used -- list is taken from Klioner's table 1, the order based on area
   of sky affected (col 2).  Order is Sun, Jupiter, Saturn, Moon, Venus,
   Uranus, Neptune.
*/

   char body_name[7][8] = {"Sun", "Jupiter", "Saturn", "Moon",
      "Venus", "Uranus", "Neptune"};

   const short int body_num[7] = {10, 5, 6, 11, 2, 7, 8};

   static short int first_time = 1;
   static short int nbodies_last = 0;

   short int error = 0;
   short int nbodies, i;

   double jd[2], tlt, pbody[3], vbody[3], pbodyo[3], x, dlt, tclose;

   cat_entry dummy_star;

   static object body[7], earth;

   jd[1] = 0.0;

/*
   Initialize output vector of observed object to equal input vector.
*/

   for (i = 0; i < 3; i++)
   {
      pos2[i] = pos1[i];
   }

/*
   Set the number of bodies -- and hence the bodies used -- based on the
   value of the 'accuracy' flag.

   Change value of 'nbodies' to include or exclude gravitating bodies
   ('nbodies' <= 0 means no deflection calculated, 'nbodies' = 1 means
   Sun only, 'nbodies' = 2 means Sun + Jupiter, etc.)  Default is
   'nbodies' = 3: Sun + Jupiter + Saturn.
*/

   if (accuracy == 0)
   {
      nbodies = 3;
   }
    else
   {
      nbodies = 1;
   }

/*
   Set up the structures of type 'object' containing the body
   information.
*/

   if ((first_time == 1) || (nbodies != nbodies_last))
   {
      for (i = 0; i < nbodies; i++)
      {
         if (i == 0)
         {
            make_cat_entry ("dummy","   ",0,0.0,0.0,0.0,0.0,0.0,0.0,
               &dummy_star);

            make_object (0,3,"Earth",&dummy_star, &earth);
         }

         if ((error = make_object (0,body_num[i],body_name[i],
            &dummy_star, &body[i])) != 0)
         {
            return (error += 30);
         }
      }
      first_time = 0;
      nbodies_last = nbodies;
   }

/*
   Compute light-time to observed object.
*/

   tlt = sqrt (pos1[0] * pos1[0] + pos1[1] * pos1[1] +
   pos1[2] *pos1[2]) / C_AUDAY;

/*
   Cycle through gravitating bodies.
*/

   for (i = 0; i < nbodies; i++)
   {


/*
   Get position of gravitating body wrt ss barycenter at time 'jd_tdb'.
*/

      jd[0] = jd_tdb;
      if ((error = ephemeris (jd,&body[i],0,accuracy, pbody,vbody))
         != 0)
      {
         return (error);
      }

/*
   Get position of gravitating body wrt observer at time 'jd_tdb'.
*/

      bary2obs (pbody,pos_obs, pbodyo,&x);

/*
   Compute light-time from point on incoming light ray that is closest
   to gravitating body.
*/

      dlt = d_light (pos2,pbodyo);

/*
   Get position of gravitating body wrt ss barycenter at time when
   incoming photons were closest to it.
*/

      tclose = jd_tdb;

      if (dlt > 0.0)
         tclose = jd_tdb - dlt;

      if (tlt < dlt)
         tclose = jd_tdb - tlt;

      jd[0] = tclose;
      if ((error = ephemeris (jd,&body[i],0,accuracy, pbody,vbody))
         != 0)
      {
         return (error);
      }

/*
   Compute deflection due to gravitating body.
*/

      grav_vec (pos2,pos_obs,pbody,RMASS[body_num[i]], pos2);
   }

/*
   If observer is not at geocenter, add in deflection due to Earth.
*/

   if (loc_code != 0 )
   {

/*
   Get position of Earth wrt solar system barycenter at time 'jd_tdb'.
*/

      jd[0] = jd_tdb;
      if ((error = ephemeris (jd,&earth,0,accuracy, pbody,vbody))
         != 0)
      {
         return (error);
      }

/*
   Compute deflection due to Earth.
*/

      grav_vec (pos2,pos_obs,pbody,RMASS[3], pos2);
   }

   return (error);
}

/********grav_vec */

void grav_vec (double *pos1, double *pos_obs, double *pos_body,
               double rmass,

               double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function corrects position vector for the deflection
      of light in the gravitational field of an arbitrary body.  This
      function valid for an observed body within the solar system as
      well as for a star.

   REFERENCES:
      Murray, C.A. (1981) Mon. Notices Royal Ast. Society 195, 639-648.
      See also formulae in Section B of the Astronomical Almanac, or
         Kaplan, G. et al. (1989) Astronomical Journal 97, 1197-1210,
         section iii f.

   INPUT
   ARGUMENTS:
      pos1[3] (double)
         Position vector of observed object, with respect to origin at
         observer (or the geocenter), components in AU.
      pos_obs [3] (double)
         Position vector of observer (or the geocenter), with respect
         to origin at solar system barycenter, components in AU.
      pos_body[3] (double)
         Position vector of gravitating body, with respect to origin
         at solar system barycenter, components in AU.
      rmass (double)
         Reciprocal mass of gravitating body in solar mass units,
         that is, Sun mass / body mass.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector of observed object, with respect to origin at
         observer (or the geocenter), corrected for gravitational
         deflection, components in AU.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      C, AU, GS          novascon.c

   FUNCTIONS
   CALLED:
      sqrt               math.h
      fabs               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/12-04/JAB (USNO/AA)

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'grvd'.

------------------------------------------------------------------------
*/
{
   short int i;

   double pq[3], pe[3], pmag, emag, qmag, phat[3], ehat[3], qhat[3],
      pdotq, edotp, qdote, fac1, fac2, p2i;

/*
   Construct vector 'pq' from gravitating body to observed object and
   construct vector 'pe' from gravitating body to observer.
*/

   for (i = 0; i < 3; i++)
   {
      pq[i] = pos_obs[i] + pos1[i] - pos_body[i];
      pe[i] = pos_obs[i] - pos_body[i];
   }

/*
   Compute vector magnitudes and unit vectors.
*/

   pmag = sqrt (pos1[0] * pos1[0] + pos1[1] * pos1[1] +
      pos1[2] * pos1[2]);
   emag = sqrt (pe[0] * pe[0] + pe[1] * pe[1] + pe[2] * pe[2]);
   qmag = sqrt (pq[0] * pq[0] + pq[1] * pq[1] + pq[2] * pq[2]);

   for (i = 0; i < 3; i++)
   {
      phat[i] = pos1[i] / pmag;
      ehat[i] = pe[i] / emag;
      qhat[i] = pq[i] / qmag;
   }

/*
   Compute dot products of vectors.
*/

   pdotq = phat[0] * qhat[0] + phat[1] * qhat[1] + phat[2] * qhat[2];
   edotp = ehat[0] * phat[0] + ehat[1] * phat[1] + ehat[2] * phat[2];
   qdote = qhat[0] * ehat[0] + qhat[1] * ehat[1] + qhat[2] * ehat[2];

/*
   If gravitating body is observed object, or is on a straight line
   toward or away from observed object to within 1 arcsec, deflection
   is set to zero; set 'pos2' equal to 'pos1'.
*/

   if (fabs (edotp) > 0.99999999999)
   {
      for (i = 0; i < 3; i++)
      {
         pos2[i] = pos1[i];
      }
   }

    else
   {

/*
   Compute scalar factors.
*/

      fac1 = 2.0 * GS / (C * C * emag * AU * rmass);
      fac2 = 1.0 + qdote;

/*
   Construct corrected position vector 'pos2'.
*/

      for (i = 0; i < 3; i++)
      {
         p2i = phat[i] + fac1 * (pdotq * ehat[i] - edotp * qhat[i]) /
            fac2;
         pos2[i] = p2i * pmag;
      }
   }

   return;
}

/********aberration */

void aberration (double *pos, double *ve, double lighttime,

                 double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      Corrects position vector for aberration of light.  Algorithm
      includes relativistic terms.

   REFERENCES:
      Murray, C. A. (1981) Mon. Notices Royal Ast. Society 195, 639-648.
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.

   INPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector, referred to origin at center of mass of the
         Earth, components in AU.
      ve[3] (double)
         Velocity vector of center of mass of the Earth, referred to
         origin at solar system barycenter, components in AU/day.
      lighttime (double)
         Light time from object to Earth in days.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, referred to origin at center of mass of the
         Earth, corrected for aberration, components in AU

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      C_AUDAY            novascon.c

   FUNCTIONS
   CALLED:
      sqrt               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.
      V1.2/12-02/JAB (USNO/AA) Improve efficiency by removing calls to
                               'pow' and replacing final 'for' loop
                               with direct assignment.
      V1.3/11-03/JAB (USNO/AA) Remove returned value.
      V1.4/02-06/WKP (USNO/AA) Changed C to C_AUDAY.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'aberat'.
      2. If 'lighttime' = 0 on input, this function will compute it.

------------------------------------------------------------------------
*/
{
   double p1mag, vemag, beta, dot, cosd, gammai, p, q, r;

   if (lighttime == 0.0)
   {
      p1mag = sqrt (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] *
         pos[2]);
      lighttime = p1mag / C_AUDAY;
   }
    else
      p1mag = lighttime * C_AUDAY;

   vemag = sqrt (ve[0] * ve[0] + ve[1] * ve[1] + ve[2] * ve[2]);
   beta = vemag / C_AUDAY;
   dot = pos[0] * ve[0] + pos[1] * ve[1] + pos[2] * ve[2];

   cosd = dot / (p1mag * vemag);
   gammai = sqrt (1.0 - beta * beta);
   p = beta * cosd;
   q = (1.0 + p / (1.0 + gammai)) * lighttime;
   r = 1.0 + p;

   pos2[0] = (gammai * pos[0] + q * ve[0]) / r;
   pos2[1] = (gammai * pos[1] + q * ve[1]) / r;
   pos2[2] = (gammai * pos[2] + q * ve[2]) / r;

   return;
}

/********rad_vel */

void rad_vel (object *cel_object, double *pos, double *vel,
              double *vel_obs, double d_obs_geo, double d_obs_sun,
              double d_obj_sun,

              double *rv)
/*
------------------------------------------------------------------------

   PURPOSE:
      Predicts the radial velocity of the observed object
      as it would be measured by spectroscopic means.  Radial
      velocity is here defined as the radial velocity measure (z)
      times the speed of light.  For a solar system body, it applies
      to a fictitious emitter at the center of the observed object,
      assumed massless (no gravitational red shift), and does not
      in general apply to reflected light.  For stars, it includes
      all effects, such as gravitational red shift, contained
      in the catalog barycentric radial velocity measure, a scalar
      derived from spectroscopy.  Nearby stars with a known kinematic
      velocity vector (obtained independently of spectroscopy) can be
      treated like solar system objects.

   REFERENCES:
      Lindegren & Dravins (2003), Astronomy & Astrophysics 401,
         1185-1201.

   INPUT
   ARGUMENTS:
      *cel_object (struct object)
         Specifies the celestial object of interest (structure defined
         in novas.h).
      *pos (double)
         Geometric position vector of object with respect to observer,
         corrected for light-time, in AU.
      *vel (double)
         Velocity vector of object with respect to solar system
         barycenter, in AU/day.
      *vel_obs (double)
         Velocity vector of observer with respect to solar system
         barycenter, in AU/day.
      d_obs_geo (double)
         Distance from observer to geocenter, in AU.
      d_obs_sun (double)
         Distance from observer to Sun, in AU.
      d_obj_sun (double)
         Distance from object to Sun, in AU.

   OUTPUT
   ARGUMENTS:
      *rv (double)
         The observed radial velocity measure times the speed of light,
         in kilometers/second.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      AU, C, GS, GE      novascon.c
      DEG2RAD            novascon.c

   FUNCTIONS
   CALLED:
      sqrt               math.h
      sin                math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/09-06/JAB (USNO/AA)
      V2.0/01-07/JAB (USNO/AA): Implement new algorithm

   NOTES:
      1. All the input arguments are BCRS quantities, expressed
      with respect to the ICRS axes.  'vel' and 'vel_obs' are kinematic
      velocities - derived from geometry or dynamics, not spectroscopy.
      2. If the object is outside the solar system, the algorithm
      used will be consistent with the IAU definition of stellar
      radial velocity, specifically, the barycentric radial velocity
      measure, which is derived from spectroscopy.  In that case,
      the vector 'vel' can be very approximate -- or, for distant stars
      or galaxies, zero -- as it will be used only for a small geometric
      correction that is proportional to proper motion.
      3. Any of the distances (last three input arguments) can be set to
      zero (0.0) if the corresponding general relativistic gravitational
      potential term is not to be evaluated.  These terms generally
      are important only at the meter/second level.  If 'd_obs_geo' and
      'd_obs_sun' are both zero, an average value will be used for the
      relativistic term for the observer, appropriate for an observer on
      the surface of the Earth.  'd_obj_sun', if given, is used only for
      solar system objects.
      4. This function is the C version of NOVAS Fortran routine
      'radvl'.

------------------------------------------------------------------------
*/
{
   static short int first_call = 1;
   short int i;

   static double c2, toms, toms2;
   double v[3], ra, dec, radvel, posmag, uk[3], v2, vo2, r, phigeo,
      phisun, rel, rar, dcr, cosdec, du[3], zc, kv, zb1, kvobs, zobs1;

/*
   Set up local constants the first time this function is called.
*/

   if (first_call)
   {
      c2 = C * C;
      toms = AU / 86400.0;
      toms2 = toms * toms;

      first_call = 0;
   }

/*
   Initialize variables needed for radial velocity calculation.
*/

   for (i = 0; i < 3; i++)
   {
      v[i] = vel[i];
   }

   switch (cel_object->type)
   {
      case 2:     /* Data for objects outside the solar system. */
         ra = cel_object->star.ra;
         dec = cel_object->star.dec;
         radvel = cel_object->star.radialvelocity;

         if (cel_object->star.parallax <= 0.0)
         {
            for (i = 0; i < 3; i++)
            {
               v[i] = 0.0;
            }
         }
         break;

      default:    /* Data not used for solar system objects. */
         ra = 0.0;
         dec = 0.0;
         radvel = 0.0;
   }

/*
   Compute length of position vector = distance to object in AU.
*/

   posmag = sqrt (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

/*
   Compute unit vector toward object.
*/

   for (i = 0; i < 3; i++)
   {
      uk[i] = pos[i] / posmag;
   }

/*
   Compute velocity-squared factors.
*/

   v2 =  (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * toms2;
   vo2 = (vel_obs[0] * vel_obs[0] + vel_obs[1] * vel_obs[1] +
      vel_obs[2] * vel_obs[2]) * toms2;

/*
   Compute geopotential at observer, unless observer is geocentric.
*/

      r = d_obs_geo * AU;

      if (r > 1.0e6)
      {
         phigeo = GE / r;
      }
       else
      {
         phigeo = 0.0;
      }


/*
   Compute solar potential at observer.
*/

      r = d_obs_sun * AU;

      if (r > 1.0e8)
      {
         phisun = GS / r;
      }
       else
      {
         phisun = 0.0;
      }

/*
   Compute relativistic potential and velocity factor for observer.
*/

   if ((d_obs_geo != 0.0) || (d_obs_sun != 0.0))
   {

/*
   Lindegren & Dravins eq. (41), second factor in parentheses.
*/
      rel = 1.0 - (phigeo + phisun) / c2 - 0.5 * vo2 / c2;
   }
    else
   {

/*
   Lindegren & Dravins eq. (42), inverse.
*/

      rel = 1.0 - 1.550e-8;
   }

/*
   Complete radial velocity calculation.
*/

   switch (cel_object->type)
   {
      case 2:     /* Objects outside the solar system. */

/*
   For stars, update barycentric radial velocity measure for change
   in view angle.
*/

         rar = ra * 15.0 * DEG2RAD;
         dcr = dec * DEG2RAD;
         cosdec = cos (dcr);
         du[0] = uk[0] - (cosdec * cos (rar));
         du[1] = uk[1] - (cosdec * sin (rar));
         du[2] = uk[2] - sin (dcr);
         zc = radvel * 1.0e3 +
            (v[0] * du[0] + v[1] * du[1] + v[2] * du[2]) * toms;

/*
   Compute observed radial velocity measure of a star (inverse of
   Lindegren & Dravins eq. (41)).
*/

         zb1 = 1.0 + zc / C;
         kvobs = (uk[0] * vel_obs[0] + uk[1] * vel_obs[1] +
            uk[2] * vel_obs[2]) * toms;
         zobs1 = zb1 * rel / (1.0 + kvobs / C);
         break;

      case 0:     /* Objects in the solar system */
      case 1:
      default:

/*
   Compute solar potential at object, if within solar system.
*/

         r = d_obj_sun * AU;

         if ((r > 1.0e8) && (r < 1.0e16))
         {
            phisun = GS / r;
         }
          else
         {
            phisun = 0.0;
         }

/*
   Compute observed radial velocity measure of a planet or other
   object -- including a nearby star -- where kinematic barycentric
   velocity vector is known and gravitational red shift is negligible
   (Lindegren & Dravins eq. (40), applied as per S. Klioner private
   communication (2006)).
*/

         kv = (uk[0] * vel[0] + uk[1] * vel[1] + uk[2] * vel[2]) * toms;
         zb1 = (1.0 + kv / C ) / (1.0 - phisun / c2 - 0.5 * v2 / c2);
         kvobs = (uk[0] * vel_obs[0] + uk[1] * vel_obs[1] +
            uk[2] * vel_obs[2] ) * toms;
         zobs1 = zb1 * rel / (1.0 + kvobs / C);
   }

/*
   Convert observed radial velocity measure to kilometers/second.
*/

   *rv = (zobs1 - 1.0) * C / 1000.0;

   return;
}

/********precession */

short int precession (double jd_tdb1, double *pos1, double jd_tdb2,

                      double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      Precesses equatorial rectangular coordinates from one epoch to
      another.  One of the two epochs must be J2000.0.  The coordinates
      are referred to the mean dynamical equator and equinox of the two
      respective epochs.

   REFERENCES:
      Explanatory Supplement To The Astronomical Almanac, pp. 103-104.
      Capitaine, N. et al. (2003), Astronomy And Astrophysics 412,
         pp. 567-586.
      Hilton, J. L. et al. (2006), IAU WG report, Celest. Mech., 94,
         pp. 351-367.

   INPUT
   ARGUMENTS:
      jd_tdb1 (double)
         TDB Julian date of first epoch.  See Note 1 below.
      pos1[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to mean dynamical equator and equinox of first epoch.
      jd_tdb2 (double)
         TDB Julian date of second epoch.  See Note 1 below.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to mean dynamical equator and equinox of second epoch.

   RETURNED
   VALUE:
      (short int)
         = 0 ... everything OK.
         = 1 ... Precession not to or from J2000.0; 'jd_tdb1' or 'jd_tdb2'
                 not 2451545.0.

   GLOBALS
   USED:
      T0, ASEC2RAD       novascon.c

   FUNCTIONS
   CALLED:
      fabs               math.h
      sin                math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.
      V1.2/03-98/JAB (USNO/AA) Change function type from 'short int' to
                               'void'.
      V1.3/12-99/JAB (USNO/AA) Precompute trig terms for greater
                               efficiency.
      V2.0/10-03/JAB (USNO/AA) Update function for consistency with
                               IERS (2000) Conventions.
      V2.1/01-05/JAB (USNO/AA) Update expressions for the precession
                               angles (extra significant digits).
      V2.2/04-06/JAB (USNO/AA) Update model to 2006 IAU convention.
                               This is model "P03" of second reference.
      V2.3/03-10/JAB (USNO/AA) Implement 'first-time' to fix bug when
                                'jd_tdb2' is 'T0' on first call to
                                function.

   NOTES:
      1. Either 'jd_tdb1' or 'jd_tdb2' must be 2451545.0 (J2000.0) TDB.
      2. This function is the C version of NOVAS Fortran routine
      'preces'.

------------------------------------------------------------------------
*/
{
   static short int first_time = 1;
   short int error = 0;

   static double t_last = 0.0;
   static double xx, yx, zx, xy, yy, zy, xz, yz, zz;
   double eps0 = 84381.406;
   double  t, psia, omegaa, chia, sa, ca, sb, cb, sc, cc, sd, cd;

/*
   Check to be sure that either 'jd_tdb1' or 'jd_tdb2' is equal to T0.
*/

   if ((jd_tdb1 != T0) && (jd_tdb2 != T0))
      return (error = 1);

/*
   't' is time in TDB centuries between the two epochs.
*/

   t = (jd_tdb2 - jd_tdb1) / 36525.0;

   if (jd_tdb2 == T0)
      t = -t;

   if ((fabs (t - t_last) >= 1.0e-15) || (first_time == 1))
   {

/*
   Numerical coefficients of psi_a, omega_a, and chi_a, along with
   epsilon_0, the obliquity at J2000.0, are 4-angle formulation from
   Capitaine et al. (2003), eqs. (4), (37), & (39).
*/

      psia   = ((((-    0.0000000951  * t
                   +    0.000132851 ) * t
                   -    0.00114045  ) * t
                   -    1.0790069   ) * t
                   + 5038.481507    ) * t;

      omegaa = ((((+    0.0000003337  * t
                   -    0.000000467 ) * t
                   -    0.00772503  ) * t
                   +    0.0512623   ) * t
                   -    0.025754    ) * t + eps0;

      chia   = ((((-    0.0000000560  * t
                   +    0.000170663 ) * t
                   -    0.00121197  ) * t
                   -    2.3814292   ) * t
                   +   10.556403    ) * t;

      eps0 = eps0 * ASEC2RAD;
      psia = psia * ASEC2RAD;
      omegaa = omegaa * ASEC2RAD;
      chia = chia * ASEC2RAD;

      sa = sin (eps0);
      ca = cos (eps0);
      sb = sin (-psia);
      cb = cos (-psia);
      sc = sin (-omegaa);
      cc = cos (-omegaa);
      sd = sin (chia);
      cd = cos (chia);
/*
   Compute elements of precession rotation matrix equivalent to
   R3(chi_a) R1(-omega_a) R3(-psi_a) R1(epsilon_0).
*/

      xx =  cd * cb - sb * sd * cc;
      yx =  cd * sb * ca + sd * cc * cb * ca - sa * sd * sc;
      zx =  cd * sb * sa + sd * cc * cb * sa + ca * sd * sc;
      xy = -sd * cb - sb * cd * cc;
      yy = -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc;
      zy = -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc;
      xz =  sb * sc;
      yz = -sc * cb * ca - sa * cc;
      zz = -sc * cb * sa + cc * ca;

      t_last = t;
      first_time = 0;
   }

   if (jd_tdb2 == T0)
   {

/*
   Perform rotation from epoch to J2000.0.
*/
      pos2[0] = xx * pos1[0] + xy * pos1[1] + xz * pos1[2];
      pos2[1] = yx * pos1[0] + yy * pos1[1] + yz * pos1[2];
      pos2[2] = zx * pos1[0] + zy * pos1[1] + zz * pos1[2];
   }
    else
   {

/*
   Perform rotation from J2000.0 to epoch.
*/

      pos2[0] = xx * pos1[0] + yx * pos1[1] + zx * pos1[2];
      pos2[1] = xy * pos1[0] + yy * pos1[1] + zy * pos1[2];
      pos2[2] = xz * pos1[0] + yz * pos1[1] + zz * pos1[2];
   }

   return (error = 0);
}

/********nutation */

void nutation (double jd_tdb, short int direction, short int accuracy,
               double *pos,

               double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:
      Nutates equatorial rectangular coordinates from mean equator and
      equinox of epoch to true equator and equinox of epoch. Inverse
      transformation may be applied by setting flag 'direction'.

   REFERENCES:
      Explanatory Supplement To The Astronomical Almanac, pp. 114-115.

   INPUT
   ARGUMENTS:
      jd_tdb (double)
         TDB Julian date of epoch.
      direction (short int)
         Flag determining 'direction' of transformation;
            direction  = 0 transformation applied, mean to true.
            direction != 0 inverse transformation applied, true to mean.
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy
      pos[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to mean equator and equinox of epoch.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to true equator and equinox of epoch.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      DEG2RAD, ASEC2RAD  novascon.c

   FUNCTIONS
   CALLED:
      e_tilt             novas.c
      cos                math.h
      sin                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.
      V1.2/11-03/JAB (USNO/AA) Remove returned value.
      V1.3/01-06/WKP (USNO/AA) Changed 'mode' to 'accuracy'.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'nutate'.

------------------------------------------------------------------------
*/
{
   double cobm, sobm, cobt, sobt, cpsi, spsi, xx, yx, zx, xy, yy, zy,
      xz, yz, zz, oblm, oblt, eqeq, psi, eps;

/*
   Call 'e_tilt' to get the obliquity and nutation angles.
*/

   e_tilt (jd_tdb,accuracy, &oblm,&oblt,&eqeq,&psi,&eps);

   cobm = cos (oblm * DEG2RAD);
   sobm = sin (oblm * DEG2RAD);
   cobt = cos (oblt * DEG2RAD);
   sobt = sin (oblt * DEG2RAD);
   cpsi = cos (psi * ASEC2RAD);
   spsi = sin (psi * ASEC2RAD);

/*
   Nutation rotation matrix follows.
*/

   xx = cpsi;
   yx = -spsi * cobm;
   zx = -spsi * sobm;
   xy = spsi * cobt;
   yy = cpsi * cobm * cobt + sobm * sobt;
   zy = cpsi * sobm * cobt - cobm * sobt;
   xz = spsi * sobt;
   yz = cpsi * cobm * sobt - sobm * cobt;
   zz = cpsi * sobm * sobt + cobm * cobt;

   if (!direction)
   {

/*
   Perform rotation.
*/

      pos2[0] = xx * pos[0] + yx * pos[1] + zx * pos[2];
      pos2[1] = xy * pos[0] + yy * pos[1] + zy * pos[2];
      pos2[2] = xz * pos[0] + yz * pos[1] + zz * pos[2];
   }
    else
   {

/*
   Perform inverse rotation.
*/

      pos2[0] = xx * pos[0] + xy * pos[1] + xz * pos[2];
      pos2[1] = yx * pos[0] + yy * pos[1] + yz * pos[2];
      pos2[2] = zx * pos[0] + zy * pos[1] + zz * pos[2];
   }

   return;
}

/********nutation_angles */

void nutation_angles (double t, short int accuracy,

                      double *dpsi, double *deps)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function returns the values for nutation in longitude and
      nutation in obliquity for a given TDB Julian date.  The nutation
      model selected depends upon the input value of 'accuracy'.  See
      notes below for important details.

   REFERENCES:
      Kaplan, G. (2005), US Naval Observatory Circular 179.

   INPUT
   ARGUMENTS:
      t (double)
         TDB time in Julian centuries since J2000.0
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      dpsi (double)
         Nutation in longitude in arcseconds.
      deps (double)
         Nutation in obliquity in arcseconds.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      T0, ASEC2RAD       novascon.c

   FUNCTIONS
   CALLED:
      iau2000a           nutation.c
      iau2000b           nutation.c
      nu2000k            nutation.c

   VER./DATE/
   PROGRAMMER:
      V1.0/12-04/JAB (USNO/AA)
      V1.1/01-06/WKP (USNO/AA): Changed 'mode' to 'accuracy'.
      V1.2/02-06/WKP (USNO/AA): Fixed units bug.
      V1.3/01-07/JAB (USNO/AA): Implemented 'low_acc_choice' construct.

   NOTES:
      1. This function selects the nutation model depending first upon
      the input value of 'accuracy'.  If 'accuracy' = 0 (full accuracy),
      the IAU 2000A nutation model is used.  If 'accuracy' = 1 (reduced
      accuracy, the model used depends upon the value of local
      variable 'low_acc_choice', which is set below.
      2. If local variable 'low_acc_choice' = 1 (the default), a
      specially truncated version of IAU 2000A, called 'NU2000K' is
      used.  If 'low_acc_choice' = 2, the IAU 2000B nutation model is
      used.
      3.  See the prologs of the nutation functions in file 'nutation.c'
      for details concerning the models.
      4. This function is the C version of NOVAS Fortran routine
      'nod'.

------------------------------------------------------------------------
*/
{

/*
   Set the value of 'low_acc_choice' according to the rules explained
   under NOTES in the prolog.
*/

   short int low_acc_choice = 1;

   double t1;

   t1 = t * 36525.0;

/*
   High accuracy mode -- use IAU 2000A.
*/

   if (accuracy == 0)
   {
      iau2000a (T0,t1, dpsi,deps);
   }

/*
   Low accuracy mode -- model depends upon value of 'low_acc_choice'.
*/

    else
   {
      switch (low_acc_choice)
      {
         case 2:   /*  use IAU 2000B */
            iau2000b (T0,t1, dpsi,deps);
            break;

         case 1:   /*  use NU2000K  */
         default:
            nu2000k (T0,t1, dpsi,deps);
      }
   }

/*
   Convert output to arcseconds.
*/

   *dpsi /= ASEC2RAD;
   *deps /= ASEC2RAD;

   return;
}

/********fund_args */

void fund_args (double t,

                double a[5])
/*
------------------------------------------------------------------------

   PURPOSE:
      To compute the fundamental arguments (mean elements) of the Sun
      and Moon.

   REFERENCES:
      Simon et al. (1994) Astronomy and Astrophysics 282, 663-683,
         esp. Sections 3.4-3.5.

   INPUT
   ARGUMENTS:
      t (double)
         TDB time in Julian centuries since J2000.0

   OUTPUT
   ARGUMENTS:
      a[5] (double)
         Fundamental arguments, in radians:
          a[0] = l (mean anomaly of the Moon)
          a[1] = l' (mean anomaly of the Sun)
          a[2] = F (mean argument of the latitude of the Moon)
          a[3] = D (mean elongation of the Moon from the Sun)
          a[4] = Omega (mean longitude of the Moon's ascending node);
                 from Simon section 3.4(b.3),
                 precession = 5028.8200 arcsec/cy)

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      ASEC2RAD, ASEC360  novascon.c

   FUNCTIONS
   CALLED:
      fmod               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/10-97/JAB (USNO/AA)
      V1.1/07-98/JAB (USNO/AA): Place arguments in the range 0-TWOPI
                                radians.
      V1.2/09-03/JAB (USNO/AA): Incorporate function 'norm_ang'.
      V1.3/11-03/JAB (USNO/AA): Update with Simon et al. expressions.
      V1.4/01-06/JAB (USNO/AA): Remove function 'norm_ang'; rewrite for
                                consistency with Fortran.
      V1.5/02-11/WKP (USNO/AA): Clarified a[4] description in prolog.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'funarg'.

------------------------------------------------------------------------
*/
{

   a[0] = fmod (485868.249036 +
             t * (1717915923.2178 +
             t * (        31.8792 +
             t * (         0.051635 +
             t * (       - 0.00024470)))), ASEC360) * ASEC2RAD;

   a[1] = fmod (1287104.79305 +
             t * ( 129596581.0481 +
             t * (       - 0.5532 +
             t * (         0.000136 +
             t * (       - 0.00001149)))), ASEC360) * ASEC2RAD;

   a[2] = fmod (335779.526232 +
             t * (1739527262.8478 +
             t * (      - 12.7512 +
             t * (      -  0.001037 +
             t * (         0.00000417)))), ASEC360) * ASEC2RAD;

   a[3] = fmod (1072260.70369 +
             t * (1602961601.2090 +
             t * (       - 6.3706 +
             t * (         0.006593 +
             t * (       - 0.00003169)))), ASEC360) * ASEC2RAD;

   a[4] = fmod (450160.398036 +
             t * ( - 6962890.5431 +
             t * (         7.4722 +
             t * (         0.007702 +
             t * (       - 0.00005939)))), ASEC360) * ASEC2RAD;

   return;
}

/********mean_obliq */

double mean_obliq (double jd_tdb)
/*
------------------------------------------------------------------------

   PURPOSE:
      To compute the mean obliquity of the ecliptic.

   REFERENCES:
      Capitaine et al. (2003), Astronomy and Astrophysics 412, 567-586.

   INPUT
   ARGUMENTS:
      jd_tdb (double)
         TDB Julian Date.

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (double)
         Mean obliquity of the ecliptic in arcseconds.

   GLOBALS
   USED:
      T0                 novascon.c

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/12-04/JAB (USNO/AA)
      V2.0/04-06/JAB (USNO/AA) Update the expression for mean obliquity
                               using data from the reference.

   NOTES:
      None.

------------------------------------------------------------------------
*/
{
   double t, epsilon;

/*
   Compute time in Julian centuries from epoch J2000.0.
*/

   t = (jd_tdb - T0) / 36525.0;

/*
   Compute the mean obliquity in arcseconds.  Use expression from the
   reference's eq. (39) with obliquity at J2000.0 taken from eq. (37)
   or Table 8.
*/

   epsilon = (((( -  0.0000000434   * t
                  -  0.000000576  ) * t
                  +  0.00200340   ) * t
                  -  0.0001831    ) * t
                  - 46.836769     ) * t + 84381.406;

   return (epsilon);
}

/********vector2radec */

short int vector2radec (double *pos,

                        double *ra, double *dec)
/*
------------------------------------------------------------------------

   PURPOSE:
      Converts an vector in equatorial rectangular coordinates to
      equatorial spherical coordinates.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.

   INPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector, equatorial rectangular coordinates.

   OUTPUT
   ARGUMENTS:
      *rightascension (double)
         Right ascension in hours.
      *declination (double)
         Declination in degrees.

   RETURNED
   VALUE:
      (short int)
         = 0 ... Everything OK.
         = 1 ... All vector components are zero; 'ra' and 'dec' are
                 indeterminate.
         = 2 ... Both pos[0] and pos[1] are zero, but pos[2] is nonzero;
                 'ra' is indeterminate.
   GLOBALS
   USED:
      ASEC2RAD           novascon.c

   FUNCTIONS
   CALLED:
      sqrt               math.h
      atan2              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.
      V1.2/01-04/JAB (USNO/AA) Remove calls to 'pow'.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'angles'.

------------------------------------------------------------------------
*/
{
   double xyproj;

   xyproj = sqrt (pos[0] * pos[0] + pos[1] * pos[1]);
   if ((xyproj == 0.0) && (pos[2] == 0))
   {
      *ra = 0.0;
      *dec = 0.0;
      return 1;
   }
    else if (xyproj == 0.0)
   {
      *ra = 0.0;
      if (pos[2] < 0.0)
         *dec = -90.0;
       else
         *dec = 90.0;
      return 2;
   }
    else
   {
      *ra = atan2 (pos[1], pos[0]) / ASEC2RAD / 54000.0;
      *dec = atan2 (pos[2], xyproj) / ASEC2RAD / 3600.0;

      if (*ra < 0.0)
         *ra += 24.0;
   }
   return 0;
}

/********radec2vector */

void radec2vector (double ra, double dec, double dist,

                   double *vector)
/*
------------------------------------------------------------------------

   PURPOSE:
      Converts equatorial spherical coordinates to a vector (equatorial
      rectangular coordinates).

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      ra (double)
         Right ascension (hours).
      dec (double)
         Declination (degrees).
      dist (double)
         Distance (AU)

   OUTPUT
   ARGUMENTS:
      vector[3] (double)
         Position vector, equatorial rectangular coordinates (AU).

   RETURNED
   VALUE:
      (short int)
         = 0 ... Everything OK.

   GLOBALS
   USED:
      DEG2RAD            novascon.c

   FUNCTIONS
   CALLED:
      cos                math.h
      sin                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/05-92/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.
      V1.2/08-09/JLB (USNO/AA) Documented "dist" in prolog

   NOTES:
      None.

------------------------------------------------------------------------
*/
{

   vector[0] = dist * cos (DEG2RAD * dec) * cos (DEG2RAD * 15.0 * ra);
   vector[1] = dist * cos (DEG2RAD * dec) * sin (DEG2RAD * 15.0 * ra);
   vector[2] = dist * sin (DEG2RAD * dec);

   return;
}

/********starvectors */

void starvectors (cat_entry *star,

                  double *pos, double *vel)
/*
------------------------------------------------------------------------

   PURPOSE:
      Converts angular quantities for stars to vectors.

   REFERENCES:
      Kaplan, G. H. et. al. (1989). Astron. Journ. 97, 1197-1210.

   INPUT
   ARGUMENTS:
      *star (struct cat_entry)
         Pointer to catalog entry structure containing ICRS catalog
         data (defined in novas.h).

   OUTPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector, equatorial rectangular coordinates,
         components in AU.
      vel[3] (double)
         Velocity vector, equatorial rectangular coordinates,
         components in AU/Day.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      ASEC2RAD, DEG2RAD  novascon.c
      AU_KM, C           novascon.c

   FUNCTIONS
   CALLED:
      sin                math.h
      cos                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Updated to C programming standards.
      V1.2/11-03/JAB (USNO/AA) Update to account for parallax in units
                               of milliarcseconds and proper motion
                               in milliarcseconds/year.
      V1.3/01-06/WKP (USNO/AA) Changed 'a[3]' to 'd' and added Doppler
                               factor.
      V1.4/02-06/WKP (USNO/AA) Corrected units of C to km/sec.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'vectrs'.

------------------------------------------------------------------------
*/
{
   double paralx, dist, r, d, cra, sra, cdc, sdc, pmr, pmd, rvl, k;

/*
   If parallax is unknown, undetermined, or zero, set it to 1e-6
   milliarcsecond, corresponding to a distance of 1 gigaparsec.
*/

   paralx = star->parallax;

   if (star->parallax <= 0.0)
      paralx = 1.0e-6;

/*
   Convert right ascension, declination, and parallax to position vector
   in equatorial system with units of AU.
*/

   dist = 1.0 / sin (paralx * 1.0e-3 * ASEC2RAD);
   r = (star->ra) * 15.0 * DEG2RAD;
   d = (star->dec) * DEG2RAD;
   cra = cos (r);
   sra = sin (r);
   cdc = cos (d);
   sdc = sin (d);

   pos[0] = dist * cdc * cra;
   pos[1] = dist * cdc * sra;
   pos[2] = dist * sdc;

/*
   Compute Doppler factor, which accounts for change in
   light travel time to star.
*/

   k = 1.0 / (1.0 - star->radialvelocity / C * 1000.0);

/*
   Convert proper motion and radial velocity to orthogonal components of
   motion with units of AU/Day.
*/

   pmr = star->promora  / (paralx * 365.25) * k;
   pmd = star->promodec / (paralx * 365.25) * k;
   rvl = star->radialvelocity * 86400.0 / AU_KM * k;

/*
   Transform motion vector to equatorial system.
*/

   vel[0] = - pmr * sra - pmd * sdc * cra + rvl * cdc * cra;
   vel[1] =   pmr * cra - pmd * sdc * sra + rvl * cdc * sra;
   vel[2] =   pmd * cdc + rvl * sdc;

   return;
}

/********tdb2tt */

void tdb2tt (double tdb_jd,

             double *tt_jd, double *secdiff)
/*
------------------------------------------------------------------------

   PURPOSE:
      Computes the Terrestrial Time (TT) or Terrestrial Dynamical Time
      (TDT) Julian date corresponding to a Barycentric Dynamical Time
      (TDB) Julian date.

   REFERENCES:
      Fairhead, L. & Bretagnon, P. (1990) Astron. & Astrophys. 229, 240.
      Kaplan, G. (2005), US Naval Observatory Circular 179.

   INPUT
   ARGUMENTS:
      tdb_jd (double)
         TDB Julian date.

   OUTPUT
   ARGUMENTS:
      *tt_jd (double)
         TT Julian date.
      *secdiff (double)
         Difference 'tdb_jd'-'tt_jd', in seconds.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      T0                 novascon.c

   FUNCTIONS
   CALLED:
      sin                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/07-92/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.
      V1.2/06-98/JAB (USNO/AA) Adopt new model (Explanatory Supplement
                               to the Astronomical Almanac, pp. 42-44
                               and p. 316.)
      V1.3/11-03/JAB (USNO/AA) Changed variable names of the input
                               Julian dates to make more descriptive.
      V1.4/01-07/JAB (USNO/AA) Adopt Fairhead & Bretagnon expression.

   NOTES:
      1. Expression used in this function is a truncated form of a
      longer and more precise series given in the first reference.  The
      result is good to about 10 microseconds.
      2. This function is the C version of NOVAS Fortran routine
      'times'.
------------------------------------------------------------------------
*/
{
   double t;

   t = (tdb_jd - T0) / 36525.0;

/*
   Expression given in USNO Circular 179, eq. 2.6.
*/

   *secdiff = 0.001657 * sin ( 628.3076 * t + 6.2401)
            + 0.000022 * sin ( 575.3385 * t + 4.2970)
            + 0.000014 * sin (1256.6152 * t + 6.1969)
            + 0.000005 * sin ( 606.9777 * t + 4.0212)
            + 0.000005 * sin (  52.9691 * t + 0.4444)
            + 0.000002 * sin (  21.3299 * t + 5.5431)
            + 0.000010 * t * sin ( 628.3076 * t + 4.2490);

   *tt_jd = tdb_jd - *secdiff / 86400.0;

    return;
}

/********cio_ra */

short int cio_ra (double jd_tt, short int accuracy,

                  double *ra_cio)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function computes the true right ascension of the celestial
      intermediate origin (CIO) at a given TT Julian date.  This is
      -(equation of the origins).

   REFERENCES:
      Kaplan, G. (2005), US Naval Observatory Circular 179.

   INPUT
   ARGUMENTS:
      jd_tt (double)
         TT Julian date.
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra_cio (double)
         Right ascension of the CIO, with respect to the true equinox
         of date, in hours (+ or -).

   RETURNED
   VALUE:
      (short int)
         = 0  ... everything OK.
         = 1  ... invalid value of 'accuracy'.
         > 10 ... 10 + the error code from function 'cio_location'.
         > 20 ... 20 + the error code from function 'cio_basis'.

   GLOBALS
   USED:
      RAD2DEG, T0        novascon.c

   FUNCTIONS
   CALLED:
      tdb2tt             novas.c
      cio_location       novas.c
      cio_basis          novas.c
      nutation           novas.c
      precession         novas.c
      frame_tie          novas.c
      atan2              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-05/JAB (USNO/AA)
      V1.1/12-05/WKP (USNO/AA) Corrected minor syntax problems.
      V1.2/01-06/WKP (USNO/AA) Changed 'mode' to 'accuracy'.
      V1.3/07-06/JAB (USNO/AA) Implement 'cio_location' construct.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'ciora'.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int rs;

   double unitx[3] = {1.0, 0.0, 0.0};
   double jd_tdb, t, secdif, x[3], y[3], z[3], w1[3], w2[3],
       eq[3], az, r_cio;

/*
   Check for valid value of 'accuracy'.
*/

   if ((accuracy < 0) || (accuracy > 1))
      return (error = 1);

/*
   'jd_tdb' is the TDB Julian date.
*/

      jd_tdb = jd_tt;
      tdb2tt (jd_tdb, &t,&secdif);
      jd_tdb = jd_tt + secdif / 86400.0;

/*
   Obtain the basis vectors, in the GCRS, for the celestial intermediate
   system defined by the CIP (in the z direction) and the CIO (in the
   x direction).
*/

   if ((error = cio_location (jd_tdb,accuracy, &r_cio,&rs)) != 0)
   {
      *ra_cio = 0.0;
      return (error += 10);
   }

   if ((error = cio_basis (jd_tdb,r_cio,rs,accuracy, x,y,z)) != 0)
   {
      return (error += 20);
   }

/*
   Compute the direction of the true equinox in the GCRS.
*/

   nutation (jd_tdb,-1,accuracy,unitx, w1);
   precession (jd_tdb,w1,T0, w2);
   frame_tie (w2,-1, eq);

/*
   Compute the RA-like coordinate of the true equinox in the celestial
   intermediate system.
*/

   az = atan2 (eq[0] * y[0] + eq[1] * y[1] + eq[2] * y[2],
      eq[0] * x[0] + eq[1] * x[1] + eq[2] * x[2]) * RAD2DEG;

/*
   The RA of the CIO is minus this coordinate.
*/

   *ra_cio = -az / 15.0;

   return (error);
}

/********cio_location */

short int cio_location (double jd_tdb, short int accuracy,

                        double *ra_cio, short int *ref_sys)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function returns the location of the celestial
      intermediate origin (CIO) for a given Julian date, as a
      right ascension with respect to either the GCRS (geocentric ICRS)
      origin or the true equinox of date.  The CIO is always located on
      the true equator (= intermediate equator) of date.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd_tdb (double)
         TDB Julian date.
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *ra_cio (double)
         Right ascension of the CIO, in hours.
      *ref_sys (short int)
         Reference system in which right ascension is given
            = 1 ... GCRS
            = 2 ... True equator and equinox of date.

   RETURNED
   VALUE:
      (short int)
         = 0  ... everything OK.
         = 1  ... unable to allocate memory for the 'cio' array.
         > 10 ... 10 + the error code from function 'cio_array'.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      cio_array          novas.c
      ira_equinox        novas.c
      fopen              stdio.h
      fclose             stdio.h
      fabs               math.h
      calloc             stdlib.h

   VER./DATE/
   PROGRAMMER:
      V1.0/07-06/JAB (USNO/AA)

   NOTES:
      1. If an external file of CIO right ascensions is available,
      it will be used and 'ref_sys' will be set to 1.  Otherwise an
      internal computation will be used and 'ref_sys' will be set to 2.
      2. The external binary file of CIO right ascensions is assumed
      to be named 'cio_ra.bin'.  Utility program 'cio_file.c', provided
      with the NOVAS-C package, creates this file from a text file also
      provided with NOVAS-C.
      3. This function is the C version of NOVAS Fortran routine
      'cioloc'.

------------------------------------------------------------------------
*/
{
   static short int first_call = 1;
   static short int ref_sys_last = 0;
   static short int use_file = 0;
   short int error = 0;

   long int n_pts = 6;
   long int i, j;

   static double t_last = 0.0;
   static double ra_last;
   double p, eq_origins;

   size_t cio_size;

   static ra_of_cio *cio;

   static FILE *cio_file;

/*
   Check if the input external binary file exists and can be read.
*/

   if (first_call)
   {
      if ((cio_file = fopen ("cio_ra.bin", "rb")) == NULL)
      {
         use_file = 0;
      }
       else
      {
         use_file = 1;
         fclose (cio_file);
      }
   }

/*
   Check if previously computed RA value can be used.
*/

   if ((fabs (jd_tdb - t_last) <= 1.0e-8))
   {
      *ra_cio = ra_last;
      *ref_sys = ref_sys_last;
      return (error = 0);
   }

/*
   Compute the RA of the CIO.
*/

   switch (use_file)
   {

/*
   -----------------------------
   Interpolate values from file.
   -----------------------------
*/

      case 1:

/*
   Allocate memory for the array 'cio'.  This array contains the values
   to be interpolated, extracted from the CIO file.
*/

         if (first_call)
         {
            cio_size = sizeof (ra_of_cio);
            cio = (ra_of_cio *) calloc ((size_t) n_pts, cio_size);
            if (cio == NULL)
               return (error = 1);
             else
               first_call = 0;
         }

/*
   Get array of values to interpolate.
*/

         if ((error = cio_array (jd_tdb,n_pts, cio)) != 0)
         {
            *ra_cio = 0.0;
            return (error += 10);
         }

/*
   Perform Lagrangian interpolation for the RA at 'tdb_jd'.
*/

         *ra_cio = 0.0;
         for (j = 0L; j < n_pts; j++)
         {
            p = 1.0;
            for (i = 0L; i < n_pts; i++)
            {
               if (i != j)
                  p *= ((jd_tdb - cio[i].jd_tdb) /
                       (cio[j].jd_tdb - cio[i].jd_tdb));
            }
            *ra_cio += (p * cio[j].ra_cio);
         }

         *ra_cio /= 54000.0;
         *ref_sys = 1;

         break;

/*
   -------------------------
   Use internal computation.
   -------------------------
*/

      case 0:

/*
   Compute equation of the origins.
*/

         if (first_call)
            first_call = 0;

         eq_origins = ira_equinox (jd_tdb,1,accuracy);

         *ra_cio = -eq_origins;
         *ref_sys = 2;

         break;
   }

   t_last = jd_tdb;
   ra_last = *ra_cio;
   ref_sys_last = *ref_sys;

   return (error);
}

/********cio_basis */

short int cio_basis (double jd_tdb, double ra_cio, short int ref_sys,
                     short int accuracy,

                     double *x, double *y, double *z)
/*
------------------------------------------------------------------------

   PURPOSE:
      To compute the orthonormal basis vectors, with
      respect to the GCRS (geocentric ICRS), of the celestial
      intermediate system defined by the celestial intermediate pole
      (CIP) (in the z direction) and the celestial intermediate origin
      (CIO) (in the x direction).  A TDB Julian date and the right
      ascension of the CIO at that date is required as input.  The
      right ascension of the CIO can be with respect to either the
      GCRS origin or the true equinox of date -- different algorithms
      are used in the two cases.

   REFERENCES:
      Kaplan, G. (2005), US Naval Observatory Circular 179.

   INPUT
   ARGUMENTS:
      jd_tdb (double)
         TDB Julian date of epoch.
      ra_cio (double)
         Right ascension of the CIO at epoch (hours).
      ref_sys (short int)
         Reference system in which right ascension is given (output
         from function 'cio_location')
            = 1 ... GCRS
            = 2 ... True equator and equinox of date.
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      *x (double)
         Unit vector toward the CIO, equatorial rectangular
         coordinates, referred to the GCRS.
      *y (double)
         Unit vector toward the y-direction, equatorial rectangular
         coordinates, referred to the GCRS.
      *z (double)
         Unit vector toward north celestial pole (CIP), equatorial
         rectangular coordinates, referred to the GCRS.

   RETURNED
   VALUE:
      (short int)
         = 0 ... everything OK.
         = 1 ... invalid value of input variable 'ref_sys'.

   GLOBALS
   USED:
      T0, DEG2RAD        novascon.c

   FUNCTIONS
   CALLED:
      nutation           novas.c
      precession         novas.c
      frame_tie          novas.c
      fabs               math.h
      sin                math.h
      cos                math.h
      sqrt               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/08-04/JAB (USNO/AA)
      V1.1/01-06/WKP (USNO/AA) Changed 'mode' to 'accuracy'.
      V1.2/07-06/JAB (USNA/AA) Incorporate code to use 'ref_sys' input.
      V1.3/06-08/WKP (USNO/AA) Changed value of direction argument in
                               calls to 'nutation' from 1 to -1 for
                               consistency.

   NOTES:
      1. This function effectively constructs the matrix C in eq. (3)
      of the reference.
      2. This function is the C version of NOVAS Fortran routine
      'ciobas'.

------------------------------------------------------------------------
*/
{
   static short int ref_sys_last = 0;
   short int error = 0;
   short int i;

   static double t_last = 0.0;
   static double xx[3], yy[3], zz[3];
   double z0[3] = {0.0, 0.0, 1.0};
   double w0[3], w1[3], w2[3], sinra, cosra, xmag;

/*
   Compute unit vector z toward celestial pole.
*/

   if (((fabs (jd_tdb - t_last) > 1.0e-8)) || (ref_sys != ref_sys_last))
   {
      nutation (jd_tdb,-1,accuracy,z0, w1);
      precession (jd_tdb,w1,T0, w2);
      frame_tie (w2,-1, zz);

      t_last = jd_tdb;
      ref_sys_last = ref_sys;
   }
    else
   {
      for (i = 0; i < 3; i++)
      {
         x[i] = xx[i];
         y[i] = yy[i];
         z[i] = zz[i];
      }
      return (error);
   }

/*
   Now compute unit vectors x and y.  Method used depends on the
   reference system in which right ascension of the CIO is given.
*/

   switch (ref_sys)
   {

/*
   ----------------------------
   RA of CIO expressed in GCRS.
   ----------------------------
*/

      case 1:

/*
   Compute vector x toward CIO in GCRS.
*/

         sinra = sin (ra_cio * 15.0 * DEG2RAD);
         cosra = cos (ra_cio * 15.0 * DEG2RAD);
         xx[0] =  zz[2] * cosra;
         xx[1] =  zz[2] * sinra;
         xx[2] = -zz[0] * cosra - zz[1] * sinra;

/*
   Normalize vector x.
*/

         xmag = sqrt (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]);
         xx[0] /= xmag;
         xx[1] /= xmag;
         xx[2] /= xmag;

/*
   Compute unit vector y orthogonal to x and z (y = z cross x).
*/

         yy[0] = zz[1] * xx[2] - zz[2] * xx[1];
         yy[1] = zz[2] * xx[0] - zz[0] * xx[2];
         yy[2] = zz[0] * xx[1] - zz[1] * xx[0];

         break;

/*
   ----------------------------------------------------------
   RA of CIO expressed in equator-and-equinox of date system.
   ----------------------------------------------------------
*/

      case 2:

/*
   Construct unit vector toward CIO in equator-and-equinox-of-date
   system.
*/

          w0[0] = cos (ra_cio * 15.0 * DEG2RAD);
          w0[1] = sin (ra_cio * 15.0 * DEG2RAD);
          w0[2] = 0.0;

/*
   Rotate the vector into the GCRS to form unit vector x.
*/

          nutation (jd_tdb,-1,accuracy,w0, w1);
          precession (jd_tdb,w1,T0, w2);
          frame_tie (w2,-1, xx);

/*
   Compute unit vector y orthogonal to x and z (y = z cross x).
*/

          yy[0] = zz[1] * xx[2] - zz[2] * xx[1];
          yy[1] = zz[2] * xx[0] - zz[0] * xx[2];
          yy[2] = zz[0] * xx[1] - zz[1] * xx[0];

          break;

/*
   ---------------------------
   Invalid value of 'ref_sys'.
   ---------------------------
*/

      default:

         for (i = 0; i < 3; i++)
         {
            xx[i] = 0.0;
            yy[i] = 0.0;
            zz[i] = 0.0;
         }

         error = 1;
         break;
   }

/*
   Load the x, y, and z arrays.
*/

   for (i = 0; i < 3; i++)
   {
      x[i] = xx[i];
      y[i] = yy[i];
      z[i] = zz[i];
   }

   return (error);
}

/********cio_array */

short int cio_array (double jd_tdb, long int n_pts,

                     ra_of_cio *cio)
/*
------------------------------------------------------------------------

   PURPOSE:
      Given an input TDB Julian date and the number of data points
      desired, this function returns a set of Julian dates and
      corresponding values of the GCRS right ascension of the celestial
      intermediate origin (CIO).  The range of dates is centered (at
      least approximately) on the requested date.  The function obtains
      the data from an external data file.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd_tdb (double)
         TDB Julian date.
      n_pts (long int)
         Number of Julian dates and right ascension values requested
         (not less than 2 or more than 20).

   OUTPUT
   ARGUMENTS:
      *cio (struct ra_of_cio)
         An time series (array) of the right ascension of the Celestial
         Intermediate Origin (CIO) with respect to the GCRS (structure
         defined in novas.h).

   RETURNED
   VALUE:
      (short int)
         = 0 ... everything OK
         = 1 ... error opening the 'cio_ra.bin' file.
         = 2 ... 'jd_tdb' not in the range of the CIO file.
         = 3 ... 'n_pts' out of range.
         = 4 ... unable to allocate memory for the internal 't' array.
         = 5 ... unable to allocate memory for the internal 'ra' array.
         = 6 ... 'jd_tdb' is too close to either end of the CIO file;
                 unable to put 'n_pts' data points into the output
                 structure.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      fopen              stdio.h
      fread              stdio.h
      fclose             stdio.h
      abs                math.h
      free               stdlib.h
      calloc             stdlib.h
      fseek              stdio.h

   VER./DATE/
   PROGRAMMER:
      V1.0/09-04/JAB (USNO/AA)
      V1.1/12-05/WKP (USNO/AA) Changed struct type from 'cio_ra' to
                               'ra_of_cio' to avoid conflicts.
      V1.2/02-08/JAB (USNO/AA) Fix file-read strategy "Case 2" and
                               improve documentation.

   NOTES:
      1. This function assumes that binary, random-access file
      'cio_ra.bin' has been created and is in the same directory as
      your executable.  This file is created by program 'cio_file.c',
      included in the NOVAS-C package.  On the first call to this
      function, file 'cio_ra.bin' is opened in read mode.

------------------------------------------------------------------------
*/
{
   static short int first_call = 1;
   short int error = 0;

   static long int last_index_rec = -50L;
   static long int last_n_pts = 0L;
   static long int header_size, record_size, n_recs;
   long int min_pts = 2;
   long int max_pts = 20;
   long int  del_n_pts, index_rec, half_int, lo_limit, hi_limit,
      del_index, abs_del_index, bytes_to_lo, n_swap, n_read, i, j;

   static double jd_beg, jd_end, t_int, *t, *ra;
   double t_temp, ra_temp;

   static size_t double_size, long_size;

   static FILE *cio_file;

/*
   Set the sizes of the file header and data records, open the CIO file,
   and read the file header on the first call to this function.
*/

   if (first_call)
   {
      double_size = sizeof (double);
      long_size = sizeof (long int);
      header_size = (long) ((size_t) 3 * double_size + long_size);
      record_size = (long) ((size_t) 2 * double_size);

/*
   Open the input (binary, random-access) file.
*/

      if ((cio_file = fopen ("cio_ra.bin", "rb")) == NULL)
         return (error = 1);

/*
   Read the file header.
*/

      fread (&jd_beg, double_size, (size_t) 1, cio_file);
      fread (&jd_end, double_size, (size_t) 1, cio_file);
      fread (&t_int, double_size, (size_t) 1, cio_file);
      fread (&n_recs, long_size, (size_t) 1, cio_file);
   }

/*
   Check the input data against limits.
*/

   if ((jd_tdb < jd_beg) || (jd_tdb > jd_end))
      return (error = 2);

   if ((n_pts < min_pts) || (n_pts > max_pts))
      return (error = 3);

/*
   Calculate the difference between the current value of 'n_pts' and
   the last value of 'n_pts'.
*/

   del_n_pts = abs (n_pts - last_n_pts);

/*
   Allocate memory for the 't' and 'ra' arrays.
*/

   if (del_n_pts != 0L)
   {
      if (!first_call)
      {
         free (t);
         free (ra);
      }

      t = (double *) calloc ((size_t) n_pts, double_size);
      if (t == NULL )
      {
         fclose (cio_file);
         return (error = 4);
      }

      ra = (double *) calloc ((size_t) n_pts, double_size);
      if (ra == NULL )
      {
         free (t);
         fclose (cio_file);
         return (error = 5);
      }

      first_call = 0;
   }

/*
   Calculate the record number of the record immediately preceding
   the date of interest: the "index record".
*/

   index_rec = (long int) ((jd_tdb - jd_beg) / t_int) + 1L;

/*
   Test the range of 'n_pts' values centered on 'index_rec' to be sure
   the range of values requested falls within the file limits.
*/

   half_int = (n_pts / 2L) - 1L;
   lo_limit = index_rec - half_int;
   hi_limit = index_rec + (n_pts - half_int - 1L);

   if ((lo_limit < 1L) || (hi_limit > n_recs))
      return (error = 6);

/*
   Compute the number of bytes from the beginning of the file to
   the 'lo_limit'.
*/

   bytes_to_lo = header_size + (lo_limit - 1L) * record_size;

/*
   Compare the current index record with the previous index record.
*/

   del_index = index_rec - last_index_rec;
   abs_del_index = abs (del_index);

/*
   Determine the file read strategy.
*/

/*
   Case 1: The input value of 'n_pts' changed since the last entry,
   or there are no data in the current arrays that can be re-used in the
   new arrays.  In this case, read all new data points into the arrays.
*/

   if ((abs_del_index > n_pts) || (del_n_pts != 0))
   {
      fseek (cio_file, bytes_to_lo, SEEK_SET);

      for (i = 0L; i < n_pts; i++)
      {
         fread (&t[i], double_size, (size_t) 1, cio_file);
         fread (&ra[i], double_size, (size_t) 1, cio_file);
      }
   }

/*
   Case 2: The new index is close enough to the previous index that
   there are some data in the current arrays that can be re-used
   in the new arrays.  The remaining data will be read from the CIO
   file.

   Note that if the new index is the same as the previous index (i.e.,
   'del_index' == 0), neither Case 2a nor 2b is satisfied, and the
   program flow goes directly to load the output arrays with the same
   values as in the current arrays.
*/

    else if ((abs_del_index <= n_pts) && (del_n_pts == 0))
   {
      n_swap = abs (n_pts - abs_del_index);
      n_read = abs_del_index;

/*
   Case 2a: The new index is less than the previous one.  Put the "old"
   data at the end of the new arrays, and read "new" data into the
   beginning of the new arrays.
*/

      if (del_index < 0L)
      {
         for (i = 0L; i < n_swap; i++)
         {
            t_temp = t[i];
            ra_temp = ra[i];

            j = i + abs_del_index;
            t[j] = t_temp;
            ra[j] = ra_temp;
         }

         fseek (cio_file, bytes_to_lo, SEEK_SET);

         for (i = 0L; i < n_read; i++)
         {
            fread (&t[i], double_size, (size_t) 1, cio_file);
            fread (&ra[i], double_size, (size_t) 1, cio_file);
         }
      }

/*
   Case 2b: The new index is greater than the previous one.  Put the
   "old" data at the beginning of the new arrays, and read "new" data
   into the end of the new arrays.
*/

       else if (del_index > 0L)
      {
         for (i = 0L; i < n_swap; i++)
         {
            j = i + abs_del_index;
            t_temp = t[j];
            ra_temp = ra[j];

            t[i] = t_temp;
            ra[i] = ra_temp;
         }

         fseek (cio_file, bytes_to_lo + (n_swap * record_size),
            SEEK_SET);

         j = i++;
         for (i = j; i < n_pts; i++)
         {
            fread (&t[i], double_size, (size_t) 1, cio_file);
            fread (&ra[i], double_size, (size_t) 1, cio_file);
         }
      }
   }

/*
   Load the output 'cio' array with the values in the 't' and 'ra'
   arrays.

   Note that if the input value of 'n_pts' has not changed since the
   last entry, all data in the current arrays can be re-used in
   the new arrays. The if statements above are bypassed and the new
   arrays are the same as the current arrays.
*/

   for (i = 0L; i < n_pts; i++)
   {
      cio[i].jd_tdb = t[i];
      cio[i].ra_cio = ra[i];
   }

/*
   Set values of 'last_index_rec' and 'last_n_pts'.
*/

   last_index_rec = index_rec;
   last_n_pts = n_pts;

   return (error);
}

/********ira_equinox */

double ira_equinox (double jd_tdb, short int equinox,
                    short int accuracy)
/*
------------------------------------------------------------------------

   PURPOSE:
      To compute the intermediate right ascension of the equinox at
      the input Julian date, using an analytical expression for the
      accumulated precession in right ascension.  For the true equinox,
      the result is the equation of the origins.

   REFERENCES:
      Capitaine, N. et al. (2003), Astronomy and Astrophysics 412,
         567-586, eq. (42).

   INPUT
   ARGUMENTS:
      jd_tdb (double)
         TDB Julian date.
      equinox (short int)
         Equinox selection flag:
            = 0 ... mean equinox
            = 1 ... true equinox.
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (double)
         Intermediate right ascension of the equinox, in hours (+ or -).
         If 'equinox' = 1 (i.e true equinox), then the returned value is
         the equation of the origins.

   GLOBALS
   USED:
      T0                 novascon.c

   FUNCTIONS
   CALLED:
      e_tilt             novas.c
      fabs               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/07-06/JAB (USNO/AA)

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'eqxra'.

------------------------------------------------------------------------
*/
{
   static short int acc_last = 99;

   static double t_last = 0.0;
   static double eq_eq = 0.0;
   double t, u, v, w, x, prec_ra, ra_eq;

/*
   Compute time in Julian centuries.
*/

   t = (jd_tdb - T0 ) / 36525.0;

/*
   For the true equinox, obtain the equation of the equinoxes in time
   seconds, which includes the 'complementary terms'.
*/

   if (equinox == 1)
   {
      if (((fabs (jd_tdb - t_last)) > 1.0e-8) || (accuracy != acc_last))
      {
         e_tilt (jd_tdb,accuracy, &u, &v, &eq_eq, &w, &x);
         t_last = jd_tdb;
         acc_last = accuracy;
      }
   }
    else
   {
      eq_eq = 0.0;
   }

/*
   Precession in RA in arcseconds taken from the reference.
*/

   prec_ra = 0.014506 +
      (((( -    0.0000000368   * t
           -    0.000029956  ) * t
           -    0.00000044   ) * t
           +    1.3915817    ) * t
           + 4612.156534     ) * t;

   ra_eq = - (prec_ra / 15.0 + eq_eq) / 3600.0;

   return (ra_eq);
}

/********ephemeris */

short int ephemeris (double jd[2], object *cel_obj, short int origin,
                     short int accuracy,

                     double *pos, double *vel)
/*
------------------------------------------------------------------------

   PURPOSE:
      Retrieves the position and velocity of a solar system body from
      a fundamental ephemeris.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      jd[2] (double)
         TDB Julian date split into two parts, where the sum
         jd[0] + jd[1] is the TDB Julian date.
      *cel_obj (struct object)
         Pointer to structure containing the designation of the body
         of interest (defined in novas.h).
      origin (short int)
         Origin code
            = 0 ... solar system barycenter
            = 1 ... center of mass of the Sun
      accuracy (short int)
         Selection for accuracy
            = 0 ... full accuracy
            = 1 ... reduced accuracy

   OUTPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector of the body at 'jd_tdb'; equatorial rectangular
         coordinates in AU referred to the ICRS.
      vel[3] (double)
         Velocity vector of the body at 'jd_tdb'; equatorial rectangular
         coordinates in AU/day referred to the ICRS.

   RETURNED
   VALUE:
      (short int)
         0    ... Everything OK.
         1    ... Invalid value of 'origin'.
         2    ... Invalid value of 'type' in 'cel_obj'.
         3    ... Unable to allocate memory.
         10+n ... where n is the error code from 'solarsystem'.
         20+n ... where n is the error code from 'readeph'.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      solarsystem         novas.c
      solarsystem_hp      novas.c
      readeph             readeph.c
      strlen              string.h
      strcpy              string.h
      malloc              stdlib.h
      free                stdlib.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-97/JAB (USNO/AA)
      V1.1/10-97/JAB (USNO/AA): Support error code from 'readeph'.
      V1.2/08-98/JAB (USNO/AA): Add computation of barycentric
                                coordinates of a minor planet; support
                                new 'readeph' argument list.
      V1.3/12-99/JAB (USNO/AA): Add error handling to call to
                                'solarsystem' (case 1).
      V1.4/12-02/JAB (USNO/AA): Fix memory leaks with 'posvel' (case 1).
      V1.5/10-08/JAB (USNO/AA): Incorporate higher-precision call to
                                the solar system ephemeris, primarily
                                to support light-time calculations.
      V1.6/02-11/JLB (USNO/AA): Reformatted description of origin for
                                consistency with other documentation.

   NOTES:
      1. It is recommended that the input structure 'cel_obj' be
      created using function 'make_object' in file novas.c.

------------------------------------------------------------------------
*/
{
   char *mp_name;

   int err = 0;
   int mp_number;

   short int error = 0;
   short int ss_number, i;

   double jd_tdb, *posvel, *sun_pos, *sun_vel;

   size_t name_len;

/*
   Check the value of 'origin'.
*/

   if ((origin < 0) || (origin > 1))
      return (error = 1);

/*
   Invoke the appropriate ephemeris access software depending upon the
   type of object.
*/

   switch (cel_obj->type)
   {

/*
   Get the position and velocity of a major planet, Pluto, Sun, or Moon.
   When high accuracy is specified, use function 'solarsystem_hp' rather
   than 'solarsystem'.
*/

      case 0:
         ss_number = cel_obj->number;

         if (accuracy == 0)
         {
            if ((error = solarsystem_hp (jd,ss_number,origin, pos,vel))
               != 0)
               error += 10;
         }
          else
         {
            jd_tdb = jd[0] + jd[1];
            if ((error = solarsystem (jd_tdb,ss_number,origin, pos,vel))
            != 0)
            error += 10;
         }
         break;

/*
   Get the position and velocity of a minor planet.

*/

      case 1:
         mp_number = (int) cel_obj->number;

         name_len = (strlen (cel_obj->name) + 1L) * sizeof (char);
         mp_name = (char *) malloc (name_len);
         if (mp_name == NULL)
            return (error = 3);
         strcpy (mp_name, cel_obj->name);

/*  The USNO minor planet software returns heliocentric positions and
    velocities.
*/

         jd_tdb = jd[0] + jd[1];
         posvel = readeph (mp_number,mp_name,jd_tdb, &err);
         if (posvel == NULL)
         {
            free (mp_name);
            return (error = 3);
         }

         if (err != 0)
         {
            free (mp_name);
            free (posvel);
            return ((short int) (20 + err));
         }

/*  Barycentric coordinates of the minor planet, if desired, are
    computed via the barycentric coordinates of the Sun, obtained
    from the solar system ephemeris.
*/

         if (origin == 0)
         {
            sun_pos = (double *) malloc (3L * sizeof (double));
            if (sun_pos == NULL)
            {
               free (mp_name);
               free (posvel);
               return (error = 3);
            }

            sun_vel = (double *) malloc (3L * sizeof (double));
            if (sun_vel == NULL)
            {
               free (mp_name);
               free (posvel);
               free (sun_pos);
               return (error = 3);
            }

            if ((error = solarsystem (jd_tdb,10,0, sun_pos,sun_vel)) != 0)
            {
               free (mp_name);
               free (posvel);
               free (sun_pos);
               free (sun_vel);
               return (error += 10);
            }

            for (i = 0; i < 3; i++)
            {
               posvel[i] += sun_pos[i];
               posvel[i+3] += sun_vel[i];
            }

            free (sun_pos);
            free (sun_vel);
         }

/*
   Break up 'posvel' into separate position and velocity vectors.
*/

         for (i = 0; i < 3; i++)
         {
            pos[i] = posvel[i];
            vel[i] = posvel[i+3];
         }

         free (mp_name);
         free (posvel);
         break;

/*
   Invalid type of object.
*/

      default:
         error = 2;
         break;

   }

   return (error);
}

/********transform_hip */

void transform_hip (cat_entry *hipparcos,

                    cat_entry *hip_2000)
/*
------------------------------------------------------------------------

   PURPOSE:
      To convert Hipparcos catalog data at epoch J1991.25 to epoch
      J2000.0, for use within NOVAS.  To be used only for Hipparcos or
      Tycho stars with linear space motion.  Both input and output data
      is in the ICRS.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      *hipparcos (struct cat_entry)
         An entry from the Hipparcos catalog, at epoch J1991.25, with
         all members having Hipparcos catalog units.  See Note 1
         below (struct defined in novas.h).

   OUTPUT
   ARGUMENTS:
      *hip_2000 (struct cat_entry)
         The transformed input entry, at epoch J2000.0.  See Note 2
         below (struct defined in novas.h).


   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      T0                 novascon.c

   FUNCTIONS
   CALLED:
      transform_cat      novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/03-98/JAB (USNO/AA)
      V1.1/11-03/JAB (USNO/AA) Update for ICRS

   NOTES:
      1. Input (Hipparcos catalog) epoch and units:
         Epoch: J1991.25
         Right ascension (RA): degrees
         Declination (Dec): degrees
         Proper motion in RA: milliarcseconds per year
         Proper motion in Dec: milliarcseconds per year
         Parallax: milliarcseconds
         Radial velocity: kilometers per second (not in catalog)
      2. Output (modified Hipparcos) epoch and units:
         Epoch: J2000.0
         Right ascension: hours
         Declination: degrees
         Proper motion in RA: milliarcseconds per year
         Proper motion in Dec: milliarcseconds per year
         Parallax: milliarcseconds
         Radial velocity: kilometers per second
      3. This function is the C version of NOVAS Fortran routine
         'gethip'.

------------------------------------------------------------------------
*/
{
   const double epoch_hip = 2448349.0625;

   cat_entry scratch;

/*
   Set up a "scratch" catalog entry containing Hipparcos data in
   "NOVAS units."
*/

   strcpy (scratch.starname, hipparcos->starname);
   scratch.starnumber = hipparcos->starnumber;
   scratch.dec = hipparcos->dec;
   scratch.promora = hipparcos->promora;
   scratch.promodec = hipparcos->promodec;
   scratch.parallax = hipparcos->parallax;
   scratch.radialvelocity = hipparcos->radialvelocity;

   strcpy (scratch.catalog, "SCR");

/*
   Convert right ascension from degrees to hours.
*/

   scratch.ra = hipparcos->ra / 15.0;

/*
   Change the epoch of the Hipparcos data from J1991.25 to J2000.0.
*/

   transform_cat (1,epoch_hip,&scratch,T0,"HP2", hip_2000);

   return;
}

/********transform_cat */

short int transform_cat (short int option, double date_incat,
                         cat_entry *incat, double date_newcat,
                         char newcat_id[SIZE_OF_CAT_NAME],

                         cat_entry *newcat)
/*
------------------------------------------------------------------------

   PURPOSE:
      To transform a star's catalog quantities for a change of epoch
      and/or equator and equinox.  Also used to rotate catalog
      quantities on the dynamical equator and equinox of J2000.0 to the
      ICRS or vice versa.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      option (short int)
         Transformation option
            = 1 ... change epoch; same equator and equinox
            = 2 ... change equator and equinox; same epoch
            = 3 ... change equator and equinox and epoch
            = 4 ... change equator and equinox J2000.0 to ICRS
            = 5 ... change ICRS to equator and equinox of J2000.0
      date_incat (double)
         TT Julian date, or year, of input catalog data.
      *incat (struct cat_entry)
         An entry from the input catalog, with units as given in
         the struct definition (struct defined in novas.h).
      date_newcat (double)
         TT Julian date, or year, of transformed catalog data.
      newcat_id[SIZE_OF_CAT_NAME] (char)
         Catalog identifier ((SIZE_OF_CAT_NAME - 1) characters maximum)
         e.g. HIP = Hipparcos, TY2 = Tycho-2.

   OUTPUT
   ARGUMENTS:
      *newcat (struct cat_entry)
         The transformed catalog entry, with units as given in
         the struct definition (struct defined in novas.h).

   RETURNED
   VALUE:
      = 0 ... Everything OK.
      = 1 ... Invalid value of an input date for option 2 or 3 (see
              Note 1 below).
      = 2 ... length of 'newcat_id' out of bounds.

   GLOBALS
   USED:
      SIZE_OF_CAT_NAME   novas.h
      T0, ASEC2RAD       novascon.c
      AU_KM, C           novascon.c

   FUNCTIONS
   CALLED:
      precession         novas.c
      frame_tie          novas.c
      sin                math.h
      cos                math.h
      sqrt               math.h
      atan2              math.h
      asin               math.h
      strcpy             string.h

   VER./DATE/
   PROGRAMMER:
      V1.0/03-98/JAB (USNO/AA)
      V1.1/11-03/JAB (USNO/AA) Update for ICRS; add options 4 and 5.
      V1.2/06-08/WKP (USNO/AA) Changed 'a[3]' to 'd', added Doppler factor,
                               and corrected 'frame_tie' direction.
      V1.3/06-08/JAB (USNO/AA) Error check on dates for options 2 and 3.
      V1.4/02-11/WKP (USNO/AA) Implement SIZE_OF_CAT_NAME.

   NOTES:
      1. 'date_incat' and 'date_newcat' may be specified either as a
      Julian date (e.g., 2433282.5) or a Julian year and fraction
      (e.g., 1950.0).  Values less than 10000 are assumed to be years.
      For 'option' = 2 or 'option' = 3, either 'date_incat' or
      'date_newcat' must be 2451545.0 or 2000.0 (J2000.0).  For
      'option' = 4 and 'option' = 5, 'date_incat' and 'date_newcat'
      are ignored.
      2. 'option' = 1 updates the star's data to account for the
      star's space motion between the first and second dates, within a
      fixed reference frame.
         'option' = 2 applies a rotation of the reference frame
      corresponding to precession between the first and second dates,
      but leaves the star fixed in space.
         'option' = 3 provides both transformations.
         'option' = 4 and 'option' = 5 provide a a fixed rotation
      about very small angles (<0.1 arcsecond) to take data from the
      dynamical system of J2000.0 to the ICRS ('option' = 4) or vice
      versa ('option' = 5).
      3. For 'option' = 1, input data can be in any fixed reference
      system. for 'option' = 2 or 'option' = 3, this function assumes
      the input data is in the dynamical system and produces output in
      the dynamical system.  for 'option' = 4, the input data must be
      on the dynamical equator and equinox of J2000.0.  for
     'option' = 5, the input data must be in the ICRS.
      4. This function cannot be properly used to bring data from
      old star catalogs into the modern system, because old catalogs
      were compiled using a set of constants that are
      incompatible with modern values.  In particular, it should not
      be used for catalogs whose positions and proper motions were
      derived by assuming a precession constant significantly different
      from the value implicit in function 'precession'.
      5. This function is the C version of NOVAS Fortran routine
      'catran'.
      6. SIZE_OF_CAT_NAME is defined in novas.h.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int j;

   double jd_incat, jd_newcat, paralx, dist, r, d, cra, sra, cdc,
      sdc, k, pos1[3], term1, pmr, pmd, rvl, vel1[3], pos2[3],
      vel2[3], xyproj;

/*
   If necessary, compute Julian dates.

   This function uses TDB Julian dates internally, but no
   distinction between TDB and TT is necessary.
*/

   if (date_incat < 10000.0)
      jd_incat = T0 + (date_incat - 2000.0) * 365.25;
    else
      jd_incat = date_incat;

   if (date_newcat < 10000.0)
      jd_newcat = T0 + (date_newcat - 2000.0) * 365.25;
    else
      jd_newcat = date_newcat;

/*
   Convert input angular components to vectors

   If parallax is unknown, undetermined, or zero, set it to 1.0e-6
   milliarcsecond, corresponding to a distance of 1 gigaparsec.
*/

   paralx = incat->parallax;
   if (paralx <= 0.0)
      paralx = 1.0e-6;

/*
   Convert right ascension, declination, and parallax to position
   vector in equatorial system with units of AU.
*/

   dist = 1.0 / sin (paralx * 1.0e-3 * ASEC2RAD);
   r = incat->ra * 54000.0 * ASEC2RAD;
   d = incat->dec * 3600.0 * ASEC2RAD;
   cra = cos (r);
   sra = sin (r);
   cdc = cos (d);
   sdc = sin (d);
   pos1[0] = dist * cdc * cra;
   pos1[1] = dist * cdc * sra;
   pos1[2] = dist * sdc;

/*
   Compute Doppler factor, which accounts for change in light
   travel time to star.
*/

   k = 1.0 / (1.0 - incat->radialvelocity / C * 1000.0);

/*
   Convert proper motion and radial velocity to orthogonal components
   of motion, in spherical polar system at star's original position,
   with units of AU/day.
*/

   term1 = paralx * 365.25;
   pmr = incat->promora  / term1 * k;
   pmd = incat->promodec / term1 * k;
   rvl = incat->radialvelocity * 86400.0 / AU_KM * k;

/*
   Transform motion vector to equatorial system.
*/

   vel1[0] = - pmr * sra - pmd * sdc * cra + rvl * cdc * cra;
   vel1[1] =   pmr * cra - pmd * sdc * sra + rvl * cdc * sra;
   vel1[2] =               pmd * cdc       + rvl * sdc;

/*
   Update star's position vector for space motion (only if 'option' = 1
   or 'option' = 3).
*/

   if ((option == 1) || (option == 3))
   {
      for (j = 0; j < 3; j++)
      {
         pos2[j] = pos1[j] + vel1[j] * (jd_newcat - jd_incat);
         vel2[j] = vel1[j];
      }
   }
    else
   {
      for (j = 0; j < 3; j++)
      {
           pos2[j] = pos1[j];
           vel2[j] = vel1[j];
      }
   }

/*
   Precess position and velocity vectors (only if 'option' = 2 or
   'option' = 3).
*/

   if ((option == 2) || (option == 3))
   {
      for (j = 0; j < 3; j++)
      {
           pos1[j] = pos2[j];
           vel1[j] = vel2[j];
      }
      if ((error = precession (jd_incat,pos1,jd_newcat, pos2)) != 0)
         return (error);
      precession (jd_incat,vel1,jd_newcat, vel2);
   }

/*
   Rotate dynamical J2000.0 position and velocity vectors to ICRS
   (only if 'option' = 4).
*/

   if (option == 4)
   {
      frame_tie (pos1,-1, pos2);
      frame_tie (vel1,-1, vel2);
   }

/*
   Rotate ICRS position and velocity vectors to dynamical J2000.0
   (only if 'option' = 5).
*/
   if (option == 5)
   {
      frame_tie (pos1,1, pos2);
      frame_tie (vel1,1, vel2);
   }

/*
   Convert vectors back to angular components for output.

   From updated position vector, obtain star's new position expressed
   as angular quantities.
*/

   xyproj = sqrt (pos2[0] * pos2[0] + pos2[1] * pos2[1]);

   if (xyproj > 0.0)
      r = atan2 (pos2[1], pos2[0]);
    else
      r = 0.0;
   newcat->ra = r / ASEC2RAD / 54000.0;
   if (newcat->ra < 0.0)
      newcat->ra += 24.0;
   if (newcat->ra >= 24.0)
      newcat->ra -= 24.0;

   d = atan2 (pos2[2], xyproj);
   newcat->dec = d / ASEC2RAD / 3600.0;

   dist = sqrt (pos2[0] * pos2[0] + pos2[1] * pos2[1] +
      pos2[2] * pos2[2]);

   paralx = asin (1.0 / dist) / ASEC2RAD * 1000.0;
   newcat->parallax = paralx;

/*
   Transform motion vector back to spherical polar system at star's
   new position.
*/

   cra = cos (r);
   sra = sin (r);
   cdc = cos (d);
   sdc = sin (d);
   pmr = - vel2[0] * sra       + vel2[1] * cra;
   pmd = - vel2[0] * cra * sdc - vel2[1] * sra * sdc + vel2[2] * cdc;
   rvl =   vel2[0] * cra * cdc + vel2[1] * sra * cdc + vel2[2] * sdc;

/*
   Convert components of motion to from AU/day to normal catalog units.
*/

   newcat->promora  = pmr * paralx * 365.25 / k;
   newcat->promodec = pmd * paralx * 365.25 / k;
   newcat->radialvelocity = rvl * (AU_KM / 86400.0) / k;

/*
  Take care of zero-parallax case.
*/

   if (newcat->parallax <= 1.01e-6)
   {
      newcat->parallax = 0.0;
      newcat->radialvelocity = incat->radialvelocity;
   }

/*
   Set the catalog identification code for the transformed catalog
   entry.
*/

   if ((short int) strlen (newcat_id) > SIZE_OF_CAT_NAME - 1)
      return (2);
    else
      strcpy (newcat->catalog, newcat_id);

/*
   Copy unchanged quantities from the input catalog entry to the
   transformed catalog entry.
*/

   strcpy (newcat->starname, incat->starname);
   newcat->starnumber = incat->starnumber;

   return (error);
}

/********limb_angle */

void limb_angle (double pos_obj[3], double pos_obs[3],

                 double *limb_ang, double *nadir_ang)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function determines the angle of an object above or below
      the Earth's limb (horizon).  The geometric limb is computed,
      assuming the Earth to be an airless sphere (no refraction or
      oblateness is included).  The observer can be on or above the
      Earth.  For an observer on the surface of the Earth, this
      function returns the approximate unrefracted altitude.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      pos_obj[3] (double)
         Position vector of observed object, with respect to origin at
         geocenter, components in AU.
      pos_obs[3] (double)
         Position vector of observer, with respect to origin at
         geocenter, components in AU.

   OUTPUT
   ARGUMENTS:
      *limb_ang (double)
         Angle of observed object above (+) or below (-) limb
         in degrees.
      *nadir_ang (double)
         Nadir angle of observed object as a fraction of apparent
         radius of limb:
            < 1.0 ... below the limb
            = 1.0 ... on the limb
            > 1.0 ... above the limb


   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      TWOPI, ERAD        novascon.c
      AU, RAD2DEG        novascon.c

   FUNCTIONS
   CALLED:
      sqrt               math.h
      asin               math.h
      acos               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/09-04/JAB (USNO/AA)

   NOTES:
      1.This function is the C version of NOVAS Fortran routine
      'limang'.

------------------------------------------------------------------------
*/
{
   static short int first_entry = 1;

   static double pi, halfpi, rade;
   double disobj, disobs, aprad, zdlim, coszd, zdobj;

   if (first_entry)
   {
      pi = TWOPI / 2.0;
      halfpi = pi / 2.0;
      rade = ERAD / AU;
      first_entry = 0;
   }

/*
   Compute the distance to the object and the distance to the observer.
*/

   disobj = sqrt (pos_obj[0] * pos_obj[0] +
                  pos_obj[1] * pos_obj[1] +
                  pos_obj[2] * pos_obj[2]);

   disobs = sqrt (pos_obs[0] * pos_obs[0] +
                  pos_obs[1] * pos_obs[1] +
                  pos_obs[2] * pos_obs[2]);

/*
   Compute apparent angular radius of Earth's limb.
*/

   if (disobs >= rade)
   {
      aprad = asin (rade / disobs);
   }
    else
   {
      aprad = halfpi;
   }

/*
   Compute zenith distance of Earth's limb.
*/

   zdlim = pi - aprad;

/*
   Compute zenith distance of observed object.
*/

   coszd = (pos_obj[0] * pos_obs[0] + pos_obj[1] * pos_obs[1] +
      pos_obj[2] * pos_obs[2]) / (disobj * disobs);

   if (coszd <= -1.0)
   {
      zdobj = pi;
   }
    else if (coszd >= 1.0)
   {
      zdobj = 0.0;
   }
    else
   {
      zdobj = acos (coszd);
   }

/*
   Angle of object wrt limb is difference in zenith distances.
*/

   *limb_ang = (zdlim - zdobj) * RAD2DEG;

/*
   Nadir angle of object as a fraction of angular radius of limb.
*/

   *nadir_ang = (pi - zdobj) / aprad;

   return;
}

/********refract */

double refract (on_surface *location, short int ref_option,
                double zd_obs)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function computes atmospheric refraction in zenith
      distance.  This version computes approximate refraction for
      optical wavelengths.

   REFERENCES:
      Explanatory Supplement to the Astronomical Almanac, p. 144.
      Bennett, G. (1982), Journal of Navigation (Royal Institute) 35,
         pp. 255-259.

   INPUT
   ARGUMENTS:
      *location (struct on_surface)
         Pointer to structure containing observer's location.  This
         structure also contains weather data (optional) for the
         observer's location (defined in novas.h).
      ref_option (short int)
         = 1 ... Use 'standard' atmospheric conditions.
         = 2 ... Use atmospheric parameters input in the 'location'
                 structure.
      zd_obs (double)
         Observed zenith distance, in degrees.

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (double)
         Atmospheric refraction, in degrees.

   GLOBALS
   USED:
      DEG2RAD            novascon.c

   FUNCTIONS
   CALLED:
      exp                math.h
      tan                math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-98/JAB (USNO/AA)

   NOTES:
      1. This function can be used for planning observations or
      telescope pointing, but should not be used for the reduction
      of precise observations.
      2. This function is the C version of NOVAS Fortran routine
      'refrac'.

------------------------------------------------------------------------
*/
{

/*
   's' is the approximate scale height of atmosphere in meters.
*/

   const double s = 9.1e3;
   double refr, p, t, h, r;

/*
   Compute refraction only for zenith distances between 0.1 and
   91 degrees.
*/

   if ((zd_obs < 0.1) || (zd_obs > 91.0))
      refr = 0.0;
    else
   {

/*
   If observed weather data are available, use them.  Otherwise, use
   crude estimates of average conditions.
*/

      if (ref_option == 2)
      {
         p = location->pressure;
         t = location->temperature;
      }
       else
      {
         p = 1010.0 * exp (-location->height / s);
         t = 10.0;
      }

      h = 90.0 - zd_obs;
      r = 0.016667 / tan ((h + 7.31 / (h + 4.4)) * DEG2RAD);
      refr = r * (0.28 * p / (t + 273.0));
   }

   return (refr);
}

/********julian_date */

double julian_date (short int year, short int month, short int day,
                    double hour)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function will compute the Julian date for a given calendar
      date (year, month, day, hour).

   REFERENCES:
      Fliegel, H. & Van Flandern, T.  Comm. of the ACM, Vol. 11, No. 10,
         October 1968, p. 657.

   INPUT
   ARGUMENTS:
      year (short int)
         Year.
      month (short int)
         Month number.
      day (short int)
         Day-of-month.
      hour (double)
         Hour-of-day.

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (double)
         Julian date.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/06-98/JAB (USNO/AA)
      V1.1/03-08/WKP (USNO/AA) Updated prolog.

   NOTES:
      1. This function is the C version of NOVAS Fortran routine
      'juldat'.
      2. This function makes no checks for a valid input calendar
      date.
      3. Input calendar date must be Gregorian.
      4. Input time value can be based on any UT-like time scale
      (UTC, UT1, TT, etc.) - output Julian date will have the same basis.
------------------------------------------------------------------------
*/
{
   long int jd12h;

   double tjd;

   jd12h = (long) day - 32075L + 1461L * ((long) year + 4800L
      + ((long) month - 14L) / 12L) / 4L
      + 367L * ((long) month - 2L - ((long) month - 14L) / 12L * 12L)
      / 12L - 3L * (((long) year + 4900L + ((long) month - 14L) / 12L)
      / 100L) / 4L;
   tjd = (double) jd12h - 0.5 + hour / 24.0;

   return (tjd);
}

/********cal_date */

void cal_date (double tjd,

               short int *year, short int *month, short int *day,
               double *hour)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function will compute a date on the Gregorian calendar given
      the Julian date.

   REFERENCES:
      Fliegel, H. & Van Flandern, T.  Comm. of the ACM, Vol. 11, No. 10,
         October 1968, p. 657.

   INPUT
   ARGUMENTS:
      tjd (double)
         Julian date.

   OUTPUT
   ARGUMENTS:
      *year (short int)
         Year.
      *month (short int)
         Month number.
      *day (short int)
         Day-of-month.
      *hour (double)
         Hour-of-day.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      fmod               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-98/JAB (USNO/AA)

   NOTES:
      1. This routine valid for any 'jd' greater than zero.
      2. Input Julian date can be based on any UT-like time scale
      (UTC, UT1, TT, etc.) - output time value will have same basis.
      3. This function is the C version of NOVAS Fortran routine
      'caldat'.


------------------------------------------------------------------------
*/
{
   long int jd, k, m, n;

   double djd;

   djd = tjd + 0.5;
   jd = (long int) djd;

   *hour = fmod (djd,1.0) * 24.0;

   k     = jd + 68569L;
   n     = 4L * k / 146097L;

   k     = k - (146097L * n + 3L) / 4L;
   m     = 4000L * (k + 1L) / 1461001L;
   k     = k - 1461L * m / 4L + 31L;

   *month = (short int) (80L * k / 2447L);
   *day   = (short int) (k - 2447L * (long int) *month / 80L);
   k      = (long int) *month / 11L;

   *month = (short int) ((long int) *month + 2L - 12L * k);
   *year  = (short int) (100L * (n - 49L) + m + k);

   return;
}

/********norm_ang */

double norm_ang (double angle)
/*
------------------------------------------------------------------------

   PURPOSE:
      Normalize angle into the range 0 <= angle < (2 * pi).

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      angle (double)
         Input angle (radians).

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (double)
          The input angle, normalized as described above (radians).

   GLOBALS
   USED:
      TWOPI              novascon.c

   FUNCTIONS
   CALLED:
      fmod               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/09-03/JAB (USNO/AA)

   NOTES:
      None.

------------------------------------------------------------------------
*/
{
   double a;

   a = fmod (angle,TWOPI);
   if (a < 0.0)
         a += TWOPI;

   return (a);
}

/********make_cat_entry */

short int make_cat_entry (char star_name[SIZE_OF_OBJ_NAME],
                          char catalog[SIZE_OF_CAT_NAME],
                          long int star_num, double ra, double dec,
                          double pm_ra, double pm_dec, double parallax,
                          double rad_vel,

                          cat_entry *star)
/*
------------------------------------------------------------------------

   PURPOSE:
      To create a structure of type 'cat_entry' containing catalog
      data for a star or "star-like" object.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      star_name[SIZE_OF_OBJ_NAME] (char)
         Object name ((SIZE_OF_OBJ_NAME - 1) characters maximum).
      catalog[SIZE_OF_CAT_NAME] (char)
         Catalog identifier ((SIZE_OF_CAT_NAME - 1) characters maximum)
         e.g. HIP = Hipparcos, TY2 = Tycho-2.
      star_num (long int)
         Object number in the catalog.
      ra (double)
         Right ascension of the object (hours).
      dec (double)
         Declination of the object (degrees).
      pm_ra (double)
         Proper motion in right ascension (milliarcseconds/year).
      pm_dec (double)
         Proper motion in declination (milliarcseconds/year).
      parallax (double)
         Parallax (milliarcseconds).
      rad_vel (double)
         Radial velocity (kilometers/second).

   OUTPUT
   ARGUMENTS:
      *star (struct cat_entry)
         Structure containing the input data (structure defined in
         novas.h).

   RETURNED
   VALUE:
      (short int)
         = 0 ... no problems.
         = 1 ... length of 'star_name' out of bounds.
         = 2 ... length of 'catalog' out of bounds.

   GLOBALS
   USED:
      SIZE_OF_OBJ_NAME   novas.h
      SIZE_OF_CAT_NAME   novas.h

   FUNCTIONS
   CALLED:
      strlen             string.h
      strcpy             string.h

   VER./DATE/
   PROGRAMMER:
      V1.0/03-98/JAB (USNO/AA)
      V1.1/08-04/JAB (USNO/AA)
      V1.3/05-10/JAB (USNO/AA): Fix bugs in set-up of 'star->starname'
                                and 'star->catalog'.

   NOTES:
      1. This utility function creates a single data structure
      encapsulating the input data.
      2. SIZE_OF_OBJ_NAME and SIZE_OF_CAT_NAME are defined in novas.h.

------------------------------------------------------------------------
*/
{

/*
   Set up the 'star' structure.
*/

   if ((short int) strlen (star_name) > SIZE_OF_OBJ_NAME - 1)
      return (1);
    else
      strcpy (star->starname, star_name);

   if ((short int) strlen (catalog) > SIZE_OF_CAT_NAME - 1)
      return (2);
    else
      strcpy (star->catalog, catalog);

   star->starnumber = star_num;
   star->ra = ra;
   star->dec = dec;
   star->promora = pm_ra;
   star->promodec = pm_dec;
   star->parallax = parallax;
   star->radialvelocity = rad_vel;

   return (0);
}

/********make_object */

short int make_object (short int type, short int number,
                       char name[SIZE_OF_OBJ_NAME], cat_entry *star_data,

                       object *cel_obj)
/*
------------------------------------------------------------------------

   PURPOSE:
      Makes a structure of type 'object' - specifying a celestial
      object - based on the input parameters.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      type (short int)
         Type of object
            = 0 ... major planet, Pluto, Sun, or Moon
            = 1 ... minor planet
            = 2 ... object located outside the solar system
                    (e.g. star, galaxy, nebula, etc.)
      number (short int)
         Body number
            For 'type' = 0: Mercury = 1,...,Pluto = 9, Sun = 10,
                            Moon = 11
            For 'type' = 1: minor planet number
            For 'type' = 2: set to 0 (zero)
      name[SIZE_OF_OBJ_NAME] (char)
         Name of the object ((SIZE_OF_OBJ_NAME - 1) characters maximum).
      *star_data (struct cat_entry)
         Structure containing basic astrometric data for any celestial
         object located outside the solar system; the catalog
         data for a star (defined in novas.h).

   OUTPUT
   ARGUMENTS:
      struct object *cel_obj
         Structure containing the object definition (defined in novas.h).

   RETURNED
   VALUE:
      (short int)
         = 0 ... everything OK
         = 1 ... invalid value of 'type'
         = 2 ... 'number' out of range
         = 3 ... Initialization of 'cel_obj' failed (object name).
         = 4 ... Initialization of 'cel_obj' failed (catalog name).
         = 5 ... 'name' is out of string bounds.

   GLOBALS
   USED:
      SIZE_OF_OBJ_NAME   novas.h
      SIZE_OF_CAT_NAME   novas.h

   FUNCTIONS
   CALLED:
      strlen             string.h
      strcpy             string.h
      toupper            ctype.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-97/JAB (USNO/AA)
      V1.1/10-98/JAB (USNO/AA): Change object name to mixed case.
      V1.2/12-02/JAB (USNO/AA): Fix bug by adding brackets in the "Set
                                the object number" 'if' block.
      V1.3/07-04/JAB (USNO/AA): Change name from 'set_body'; expand to
                                include objects outside solar system.
      V1.4/05-10/JAB (USNO/AA): Implement SIZE_OF_OBJ_NAME and
                                SIZE_OF_CAT_NAME; check length of
                                'name'.
      V1.5/09-10/WKP (USNO/AA): Explicitly cast 'toupper' return value
                                to type 'char' to silence compiler.

   NOTES:
      1. SIZE_OF_OBJ_NAME and SIZE_OF_CAT_NAME are defined in novas.h.

------------------------------------------------------------------------
*/
{
   short int error = 0;
   short int i;

/*
   Initialize the structure.
*/

   cel_obj->type = 0;
   cel_obj->number = 0;
   strcpy (cel_obj->name, "  ");

   if ((short int) strlen ("  ") > SIZE_OF_OBJ_NAME - 1)
      return (error = 3);
    else
   {
      strcpy (cel_obj->name, "  ");
      strcpy (cel_obj->star.starname, "  ");
   }

   if ((short int) strlen ("  ") > SIZE_OF_CAT_NAME - 1)
      return (error = 4);
    else
      strcpy (cel_obj->star.catalog, "  ");

   cel_obj->star.starnumber = 0L;
   cel_obj->star.ra = 0.0;
   cel_obj->star.dec = 0.0;
   cel_obj->star.promora = 0.0;
   cel_obj->star.promodec = 0.0;
   cel_obj->star.parallax = 0.0;
   cel_obj->star.radialvelocity = 0.0;

/*
   Set the object type.
*/

   if ((type < 0) || (type > 2))
      return (error = 1);
    else
      cel_obj->type = type;

/*
   Set the object number.
*/

   if (type == 0)
   {
      if ((number <= 0) || (number > 11))
         return (error = 2);
   }
    else if (type == 1)
   {
      if (number <= 0)
         return (error = 2);
   }
    else
   {
      number = 0;
   }

   cel_obj->number = number;

/*
   Check length of 'name'.  Set the object name in upper case.
*/

   if ((short int) strlen (name) > SIZE_OF_OBJ_NAME - 1)
      return (error = 5);

   for (i = 0; i < SIZE_OF_OBJ_NAME - 1; i++)
   {
      cel_obj->name[i] = (char) toupper (name[i]);
      if (name[i] == '\0')
         break;
   }
   cel_obj->name[i] = '\0';

/*
   Populate the astrometric-data structure if the object is 'type' = 2.
*/

   if (type == 2)
   {
      strcpy (cel_obj->star.starname, star_data->starname);
      strcpy (cel_obj->star.catalog, star_data->catalog);
      cel_obj->star.starnumber = star_data->starnumber;
      cel_obj->star.ra = star_data->ra;
      cel_obj->star.dec = star_data->dec;
      cel_obj->star.promora = star_data->promora;
      cel_obj->star.promodec = star_data->promodec;
      cel_obj->star.parallax = star_data->parallax;
      cel_obj->star.radialvelocity = star_data->radialvelocity;
   }

    return (error);
}

/********make_observer */

short int make_observer (short int where, on_surface *obs_surface,
                         in_space *obs_space,

                         observer *obs)
/*
------------------------------------------------------------------------

   PURPOSE:
      Makes a structure of type 'observer' - specifying the location of
      the observer.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      where (short int)
         Integer code specifying location of observer.
            = 0: observer at geocenter
            = 1: observer on surface of earth
            = 2: observer on near-earth spacecraft
      *obs_surface (struct on_surface)
         Structure containing data for an observer's location on the
         surface of the Earth; used when 'where' = 1 (defined in
         novas.h).
      *obs_space (struct in_space)
         Data for an observer's location on a near-Earth spacecraft;
         used when 'where' = 2  (defined in novas.h).

   OUTPUT
   ARGUMENTS:
      struct observer *obs
         Structure specifying the location of the observer (defined in
         novas.h).

   RETURNED
   VALUE:
      (short int)
         = 0 ... everything OK.
         = 1 ... input value of 'where' is out-of-range.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/08-04/JAB (USNO/AA)

   NOTES:
      None.

------------------------------------------------------------------------
*/
{
   short int error = 0;

/*
   Initialize the output structure.
*/

   obs->where = where;
   obs->on_surf.latitude = 0.0;
   obs->on_surf.longitude = 0.0;
   obs->on_surf.height = 0.0;
   obs->on_surf.temperature = 0.0;
   obs->on_surf.pressure = 0.0;
   obs->near_earth.sc_pos[0] = 0.0;
   obs->near_earth.sc_pos[1] = 0.0;
   obs->near_earth.sc_pos[2] = 0.0;
   obs->near_earth.sc_vel[0] = 0.0;
   obs->near_earth.sc_vel[1] = 0.0;
   obs->near_earth.sc_vel[2] = 0.0;

/*
   Populate the output structure based on the value of 'where'.
*/

   switch (where)
   {
      case (0):  /* Geocentric */
         break;

      case (1):  /* On surface of Earth */
         obs->on_surf.latitude = obs_surface->latitude;
         obs->on_surf.longitude = obs_surface->longitude;
         obs->on_surf.height = obs_surface->height;
         obs->on_surf.temperature = obs_surface->temperature;
         obs->on_surf.pressure = obs_surface->pressure;
         break;

      case (2):  /* In near-Earth spacecraft */
         obs->near_earth.sc_pos[0] = obs_space->sc_pos[0];
         obs->near_earth.sc_pos[1] = obs_space->sc_pos[1];
         obs->near_earth.sc_pos[2] = obs_space->sc_pos[2];
         obs->near_earth.sc_vel[0] = obs_space->sc_vel[0];
         obs->near_earth.sc_vel[1] = obs_space->sc_vel[1];
         obs->near_earth.sc_vel[2] = obs_space->sc_vel[2];
         break;

      default:
         error = 1;
         break;
   }

   return (error);
}

/********make_observer_at_geocenter */

void make_observer_at_geocenter (

                                 observer *obs_at_geocenter)
/*
------------------------------------------------------------------------

   PURPOSE:
      Makes a structure of type 'observer' specifying an observer at the
      geocenter.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      None.

   OUTPUT
   ARGUMENTS:
      struct observer *obs_at_geocenter
         Structure specifying the location of the observer at the geocenter
         (defined in novas.h).

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      make_in_space      novas.c
      make_on_surface    novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/03-07/JAB (USNO/AA)

   NOTES:
      1. The output data structure is one of the inputs to NOVAS-C
      function 'place'.

------------------------------------------------------------------------
*/
{
   double latitude = 0.0;
   double longitude = 0.0;
   double height = 0.0;
   double temperature = 0.0;
   double pressure = 0.0;
   double sat_pos[3] = {0.0, 0.0, 0.0};
   double sat_vel[3] = {0.0, 0.0, 0.0};

   in_space sat_state;

   on_surface surface_loc;

   make_in_space (sat_pos,sat_vel, &sat_state);
   make_on_surface (latitude,longitude,height,temperature,pressure,
      &surface_loc);

   obs_at_geocenter->where = 0;
   obs_at_geocenter->on_surf = surface_loc;
   obs_at_geocenter->near_earth = sat_state;
}

/********make_observer_on_surface */

void make_observer_on_surface (double latitude, double longitude,
                               double height, double temperature,
                               double pressure,

                               observer *obs_on_surface)
/*
------------------------------------------------------------------------

   PURPOSE:
      Makes a structure of type 'observer' specifying the location
      of and weather for an observer on the surface of the Earth.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      latitude (double)
         Geodetic (ITRS) latitude in degrees; north positive.
      longitude (double)
         Geodetic (ITRS) longitude in degrees; east positive.
      height (double)
         Height of the observer (meters).
      temperature (double)
         Temperature (degrees Celsius).
      pressure (double)
         Atmospheric pressure (millibars).

   OUTPUT
   ARGUMENTS:
      struct observer *obs_on_surface
         Structure containing the location of and weather for an
         observer on the surface of the Earth (defined in novas.h).

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      make_in_space      novas.c
      make_on_surface    novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/03-07/JAB (USNO/AA)

   NOTES:
      1. The output data structure is one of the inputs to NOVAS-C
      function 'place'.

------------------------------------------------------------------------
*/
{
   double sat_pos[3] = {0.0, 0.0, 0.0};
   double sat_vel[3] = {0.0, 0.0, 0.0};

   in_space sat_state;

   on_surface surface_loc;

   make_in_space (sat_pos,sat_vel, &sat_state);
   make_on_surface (latitude,longitude,height,temperature,pressure,
      &surface_loc);

   obs_on_surface->where = 1;
   obs_on_surface->on_surf = surface_loc;
   obs_on_surface->near_earth = sat_state;
}

/********make_observer_in_space */

void make_observer_in_space (double sc_pos[3], double sc_vel[3],

                             observer *obs_in_space)
/*
------------------------------------------------------------------------

   PURPOSE:
      Makes a structure of type 'observer' specifying the position and
      velocity of an observer situated on a near-Earth spacecraft.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      sc_pos[3] (double)
         Geocentric position vector (x, y, z) in km.
      sc_vel[3] (double)
         Geocentric velocity vector (x_dot, y_dot, z_dot) in km/s.

   OUTPUT
   ARGUMENTS:
      struct observer *obs_in_space
         Structure containing the position and velocity of an observer
         situated on a near-Earth spacecraft (defined in novas.h).

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      make_in_space      novas.c
      make_on_surface    novas.c

   VER./DATE/
   PROGRAMMER:
      V1.0/03-07/JAB (USNO/AA)

   NOTES:
      1. Both input vectors are with respect to true equator and equinox
      of date.
      2. The output data structure is one of the inputs to NOVAS-C
      function 'place'.

------------------------------------------------------------------------
*/
{
   const double latitude = 0.0;
   const double longitude = 0.0;
   const double height = 0.0;
   const double temperature = 0.0;
   const double pressure = 0.0;

   in_space sat_state;

   on_surface surface_loc;

   make_in_space (sc_pos,sc_vel, &sat_state);
   make_on_surface (latitude,longitude,height,temperature,pressure,
      &surface_loc);

   obs_in_space->where = 2;
   obs_in_space->on_surf = surface_loc;
   obs_in_space->near_earth = sat_state;
}

/********make_on_surface */

void make_on_surface (double latitude, double longitude, double height,
                      double temperature, double pressure,

                      on_surface *obs_surface)
/*
------------------------------------------------------------------------

   PURPOSE:
      Makes a structure of type 'on_surface' - specifying the location
      of and weather for an observer on the surface of the Earth.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      latitude (double)
         Geodetic (ITRS) latitude in degrees; north positive.
      longitude (double)
         Geodetic (ITRS) longitude in degrees; east positive.
      height (double)
         Height of the observer (meters).
      temperature (double)
         Temperature (degrees Celsius).
      pressure (double)
         Atmospheric pressure (millibars).

   OUTPUT
   ARGUMENTS:
      struct on_surface *obs_surface
         Structure containing the location of and weather for an
         observer on the surface of the Earth (defined in novas.h).

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/08-04/JAB (USNO/AA)

   NOTES:
      None.

------------------------------------------------------------------------
*/
{
   obs_surface->latitude = latitude;
   obs_surface->longitude = longitude;
   obs_surface->height = height;
   obs_surface->temperature = temperature;
   obs_surface->pressure = pressure;
}

/********make_in_space */

void make_in_space (double sc_pos[3], double sc_vel[3],

                    in_space *obs_space)
/*
------------------------------------------------------------------------

   PURPOSE:
      Makes a structure of type 'in_space' - specifying the position and
      velocity of an observer situated on a near-Earth spacecraft.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      sc_pos[3] (double)
         Geocentric position vector (x, y, z) in km.
      sc_vel[3] (double)
         Geocentric velocity vector (x_dot, y_dot, z_dot) in km/s.

   OUTPUT
   ARGUMENTS:
      struct in_space *obs_space
         Structure containing the position and velocity of an observer
         situated on a near-Earth spacecraft (defined in novas.h).

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/08-04/JAB (USNO/AA)

   NOTES:
      1. Both vectors with respect to true equator and equinox of date.

------------------------------------------------------------------------
*/
{
   obs_space->sc_pos[0] = sc_pos[0];
   obs_space->sc_pos[1] = sc_pos[1];
   obs_space->sc_pos[2] = sc_pos[2];

   obs_space->sc_vel[0] = sc_vel[0];
   obs_space->sc_vel[1] = sc_vel[1];
   obs_space->sc_vel[2] = sc_vel[2];
}
