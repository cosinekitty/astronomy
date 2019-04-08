/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1

  solsys2.c: Interface to JPL ephemeris-access software for use with jplint.f

  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC
  http://www.usno.navy.mil/USNO/astronomical-applications
*/


#ifndef _NOVAS_
   #include "novas.h"
#endif

/*
   Dummy function prototype for Fortran subroutine 'jplint'.
*/

void jplint_ (double *tjd, long int *targ, long int *cent,

            double *posvel, long int *err_flg);

void jplihp_ (double *jed, long int *targ, long int *cent,

            double *posvel, long int *err_flg);




/********solarsystem */

short int solarsystem (double tjd, short int body, short int origin,

                       double *position, double *velocity)
/*
------------------------------------------------------------------------

   PURPOSE:
      Provides an interface between the JPL direct-access solar system
      ephemerides and NOVAS-C.

   REFERENCES:
      JPL. 2007, JPL Planetary and Lunar Ephemerides: Export Information,
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
         interest;  Mercury = 1,...,Pluto = 9, Sun = 10, Moon = 11.
      origin (short int)
         Origin code
            = 0 ... solar system barycenter
            = 1 ... center of mass of the Sun

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
      (short int)
         0...Everything OK.
         1...Invalid value of body or origin.
         2...Error detected by JPL software.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      jplint_ (A Fortran subroutine that serves as the interface to
              JPL's Fortran code that accesses the solar system
              ephemerides)

   VER./DATE/
   PROGRAMMER:
      V1.0/10-97/JAB (USNO/AA)
      V1.1/06-08/WKP (USNO/AA) Added underscore to end of jplint_
                               function.
      V1.2/11-09/WKP (USNO/AA) Updated output argument names to
                               'position' and 'velocity'.
      V1.3/02-11/JLB (USNO/AA): Reformatted description of origin for
                                consistency with other documentation.
      V1.4/02-11/WKP (USNO/AA): More minor prolog changes for
                                consistency among all solsysn.c files.


   NOTES:
      1. This function is based on 'solsys2d.c' from version 1 of
         NOVAS-C.  It generalizes access to the JPL software by calling
         a Fortran interface subroutine, 'jplint', instead of making a
         direct call to the JPL subroutine 'pleph', whose arguments
         have changed several times throughout the years.  This way,
         any future change to the arguments can be accommodated in
         'jplint' rather than in this function.
------------------------------------------------------------------------
*/
{

   short int i, error = 0;

   long int targ, cent, err_flg = 0;

   double posvel[6] = {0.0,0.0,0.0,0.0,0.0,0.0};

/*
   Perform sanity checks on the input body and origin.
*/

   if ((body < 1) || (body > 11))
      return 1;
    else if ((origin < 0) || (origin > 1))
      return 1;

/*
   Select 'targ' according to the value of 'body'.
*/

   if (body == 10)
      targ = 11L;
     else if (body == 11)
      targ = 10L;
     else
      targ = (long int) body;

/*
   Select 'cent' according to the value of 'origin'.
*/

   if (origin == 0)
      cent = 12L;
    else if (origin == 1)
      cent = 11L;
    else
      return 1;

/*
   Call Fortran subroutine 'jplint' to obtain position and velocity
   array 'posvel'.  This is the only point in the NOVAS-C package
   where the Fortran/C interface occurs.
   Note that arguments must be sent to Fortran by reference, not by
   value.
*/

   jplint_ (&tjd, &targ, &cent, posvel, &err_flg);
   if (err_flg)
      return (error = 2);

/*
   Decompose 'posvel' into 'position' and 'velocity'.
*/

   for (i = 0; i < 3; i++)
   {
      position[i] = posvel[i];
      velocity[i] = posvel[i+3];
   }

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
         interest;  Mercury = 1,...,Pluto = 9, Sun = 10, Moon = 11.
      origin (short int)
         Origin code
            = 0 ... solar system barycenter
            = 1 ... center of mass of the Sun

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
      (short int)
         0...Everything OK.
         1...Invalid value of body or origin.
         2...Error detected by JPL software.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      jplihp_ (A Fortran subroutine that serves as the interface to
              JPL's Fortran code that accesses the solar system
              ephemerides)

   VER./DATE/
   PROGRAMMER:
      V1.0/12-06/JAB (USNO/AA)
      V1.1/06-08/WKP (USNO/AA) Added underscore to end of jplihp_
                               function.
      V1.2/11-09/WKP (USNO/AA) Updated output argument names to
                               'position' and 'velocity'.
      V1.3/02-11/JLB (USNO/AA): Reformatted description of origin for
                                consistency with other documentation.
      V1.4/02-11/WKP (USNO/AA): More minor prolog changes for
                                consistency among all solsysn.c files.


   NOTES:
      1. This function is based on 'solsys2d.c' from version 1 of
         NOVAS-C.  It generalizes access to the JPL software by calling
         a Fortran interface subroutine, 'jplihp', instead of making a
         direct call to the JPL subroutine 'dpleph', whose arguments
         have changed several times throughout the years.  This way,
         any future change to the arguments can be accommodated in
         'jplihp' rather than in this function.
      2. This function supports the "split" Julian date feature of
         JPL subroutine 'dpleph' for highest precision.  For
         usual applications, use function 'solarsystem' in file
         'solsys2.c'.

------------------------------------------------------------------------
*/
{

   short int i, error = 0;

   long int targ, cent, err_flg = 0;

   double posvel[6] = {0.0,0.0,0.0,0.0,0.0,0.0};

/*
   Perform sanity checks on the input body and origin.
*/

   if ((body < 1) || (body > 11))
      return 1;
    else if ((origin < 0) || (origin > 1))
      return 1;

/*
   Select 'targ' according to the value of 'body'.
*/

   if (body == 10)
      targ = 11L;
     else if (body == 11)
      targ = 10L;
     else
      targ = (long int) body;

/*
   Select 'cent' according to the value of 'origin'.
*/

   if (origin == 0)
      cent = 12L;
    else if (origin == 1)
      cent = 11L;
    else
      return 1;

/*
   Call Fortran subroutine 'jplihp' to obtain position and velocity
   array 'posvel'.  This is the only point in the NOVAS-C package
   where the Fortran/C interface occurs.
   Note that arguments must be sent to Fortran by reference, not by
   value.
*/

   jplihp_ (tjd, &targ, &cent, posvel, &err_flg);
   if (err_flg)
      return (error = 2);

/*
   Decompose 'posvel' into 'position' and 'velocity'.
*/

   for (i = 0; i < 3; i++)
   {
      position[i] = posvel[i];
      velocity[i] = posvel[i+3];
   }

   return 0;
}

