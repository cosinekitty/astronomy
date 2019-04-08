/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1

  solsys1.c: Interface to JPL ephemerides for use with eph_manager

  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC
  http://www.usno.navy.mil/USNO/astronomical-applications
*/


#ifndef _NOVAS_
   #include "novas.h"
#endif

#ifndef _EPHMAN_
   #include "eph_manager.h"
#endif


/********solarsystem */

short int solarsystem (double tjd, short int body, short int origin,

                       double *position, double *velocity)
/*
------------------------------------------------------------------------

   PURPOSE:
      Provides an interface between the JPL direct-access solar system
      ephemerides and NOVAS-C.

   REFERENCES:
      JPL. 2007, "JPL Planetary and Lunar Ephemerides: Export Information,"
        (Pasadena, CA: JPL) http://ssd.jpl.nasa.gov/?planet_eph_export.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      tjd (double)
         Julian date of the desired time, on the TDB time scale.
      body (short int)
         Body identification number for the solar system object of
         interest;Mercury = 1, ..., Pluto= 9, Sun= 10, Moon = 11.
      origin (short int)
         Origin code
            = 0 ... solar system barycenter
            = 1 ... center of mass of the Sun
            = 2 ... center of Earth

   OUTPUT
   ARGUMENTS:
      position[3] (double)
         Position vector of 'body' at tjd; equatorial rectangular
         coordinates in AU referred to the ICRS.
      velocity[3] (double)
         Velocity vector of 'body' at tjd; equatorial rectangular
         system referred to the ICRS, in AU/day.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      planet_ephemeris          eph_manager.h

   VER./DATE/
   PROGRAMMER:
      V1.0/03-93/WTH (USNO/AA): Convert FORTRAN to C.
      V1.1/07-93/WTH (USNO/AA): Update to C standards.
      V2.0/07-93/WTH (USNO/AA): Update to C standards.
      V2.1/06-99/JAB (USNO/AA): Minor style and documentation mods.
      V2.2/11-06/JAB (USNO/AA): Update to handle split-JD input
                                now supported in 'Planet_Ephemeris'.
      V2.3/09-10/WKP (USNO/AA): Initialized local variable 'center'
                                to silence compiler warning.
      V2.4/10-10/WKP (USNO/AA): Changed 'Planet_Ephmeris' function
                                name to lower case to comply with
                                C coding standards.
      V2.5/02-11/JLB (USNO/AA): Reformatted description of origin for
                                consistency with other documentation.
      V2.6/02-11/WKP (USNO/AA): More minor prolog changes for
                                consistency among all solsysn.c files.


   NOTES:
      1. This function and function 'planet_ephemeris' were designed
         to work with the 1997 version of the JPL ephemerides, as
         noted in the references.
      2. The user must create the binary ephemeris files using
         software from JPL, and open the file using function
         'ephem_open' in eph_manager.h, prior to calling this
         function.
      3. This function places the entire Julian date in the first
         element of the input time to 'planet_ephemeris'. This is
         adequate for all but the highest precision applications.  For
         highest precision, use function 'solarsystem_hp' in file
         'solsys1.c'.
      4. Function 'planet_ephemeris' is a C rewrite of the JPL Fortran
         subroutine 'pleph'.

------------------------------------------------------------------------
*/
{
   short int target, center = 0;

   double jd[2];

/*
   Perform sanity checks on the input body and origin.
*/

   if ((body < 1) || (body > 11))
      return 1;
    else if ((origin < 0) || (origin > 2))
      return 2;

/*
   Select 'target' according to value of 'body'.
*/

   switch (body)
  {
      case 10:
         target = 10;
         break;
      case 11:
         target = 9;
         break;
      default:
         target = body - 1;
   }

/*
   Select 'center' according to the value of 'origin'.
*/

   if (!origin)
      center = 11;
    else if (origin == 1)
      center = 10;
    else if (origin == 2)
      center = 2;

/*
   Obtain position and velocity vectors.  The entire Julian date is
   contained in the first element of 'jd'.  This is adequate for all
   but the highest precision applications.
*/

   jd[0] = tjd;
   jd[1] = 0.0;

   planet_ephemeris (jd,target,center, position,velocity);

   return 0;
}

/********solarsystem_hp */

short int solarsystem_hp (double tjd[2], short int body,
                          short int origin,

                          double *position, double *velocity)
/*
------------------------------------------------------------------------

   PURPOSE:
      Provides an interface between the JPL direct-access solar system
      ephemerides and NOVAS-C for highest precision applications.

   REFERENCES:
      JPL. 2007, "JPL Planetary and Lunar Ephemerides: Export Information,"
        (Pasadena, CA: JPL) http://ssd.jpl.nasa.gov/?planet_eph_export.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      tjd[2] (double)
         Two-element array containing the Julian date, which may be
         split any way (although the first element is usually the
         "integer" part, and the second element is the "fractional"
         part).  Julian date is on the TDB or "T_eph" time scale.
      body (short int)
         Body identification number for the solar system object of
         interest;Mercury = 1, ..., Pluto= 9, Sun= 10, Moon = 11.
      origin (short int)
         Origin code
            = 0 ... solar system barycenter
            = 1 ... center of mass of the Sun
            = 2 ... center of Earth

   OUTPUT
   ARGUMENTS:
      position[3] (double)
         Position vector of 'body' at tjd; equatorial rectangular
         coordinates in AU referred to the ICRS.
      velocity[3] (double)
         Velocity vector of 'body' at tjd; equatorial rectangular
         system referred to the ICRS, in AU/day.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      planet_ephemeris          eph_manager.h

   VER./DATE/
   PROGRAMMER:
      V1.0/11-06/JAB (USNO/AA): Update to handle split-JD input
                                now supported in 'Planet_Ephemeris'.
      V1.1/09-10/WKP (USNO/AA): Initialized local variable 'center'
                                to silence compiler warning.
      V1.2/10-10/WKP (USNO/AA): Changed 'Planet_Ephmeris' function
                                name to lower case to comply with
                                C coding standards.
      V1.3/02-11/JLB (USNO/AA): Reformatted description of origin for
                                consistency with other documentation.
      V1.4/02-11/WKP (USNO/AA): More minor prolog changes for
                                consistency among all solsysn.c files.


   NOTES:
      1. This function and function 'planet_ephemeris' were designed
         to work with the 1997 version of the JPL ephemerides, as
         noted in the references.
      2. The user must create the binary ephemeris files using
         software from JPL, and open the file using function
         'ephem_open' in eph_manager.h, prior to calling this
         function.
      3. This function supports the "split" Julian date feature of
         function 'planet_ephemeris' for highest precision.  For
         usual applications, use function 'solarsystem' in file
         'solsys1.c'.
      4. Function 'planet_ephemeris' is a C rewrite of the JPL Fortran
         subroutine 'pleph'.

------------------------------------------------------------------------
*/
{
   short int target, center = 0;

/*
   Perform sanity checks on the input body and origin.
*/

   if ((body < 1) || (body > 11))
      return 1;
    else if ((origin < 0) || (origin > 2))
      return 2;

/*
   Select 'target' according to value of 'body'.
*/

   switch (body)
  {
      case 10:
         target = 10;
         break;
      case 11:
         target = 9;
         break;
      default:
         target = body - 1;
   }

/*
   Select 'center' according to the value of 'origin'.
*/

   if (!origin)
      center = 11;
    else if (origin == 1)
      center = 10;
    else if (origin == 2)
      center = 2;

/*
   Obtain position and velocity vectors.  The Julian date is split
   between two double-precision elements for highest precision.
*/

   planet_ephemeris (tjd,target,center, position,velocity);

   return 0;
}
