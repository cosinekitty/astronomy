/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1
 
  novascon.h: Header file for novascon.c 
 
  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC 
  http://www.usno.navy.mil/USNO/astronomical-applications
*/

#ifndef _CONSTS_
   #define _CONSTS_

   extern const short int FN1;
   extern const short int FN0;

/*
   TDB Julian date of epoch J2000.0.
*/

   extern const double T0;

/*
   Speed of light in meters/second is a defining physical constant.
*/

   extern const double C;

/*
   Light-time for one astronomical unit (AU) in seconds, from DE-405.
*/

   extern const double AU_SEC;

/*
   Speed of light in AU/day.  Value is 86400 / AU_SEC.
*/

   extern const double C_AUDAY;

/*
   Astronomical unit in meters.  Value is AU_SEC * C.
*/

   extern const double AU;

/*
   Astronomical Unit in kilometers.
*/

   extern const double AU_KM;

/*
   Heliocentric gravitational constant in meters^3 / second^2, from
   DE-405.
*/

   extern const double GS;

/*
   Geocentric gravitational constant in meters^3 / second^2, from
   DE-405.
*/

   extern const double GE;

/*
   Radius of Earth in kilometers from IERS Conventions (2003).
*/

   extern const double ERAD;

/*
   Earth ellipsoid flattening from IERS Conventions (2003).
   Value is 1 / 298.25642.
*/

   extern const double F;

/*
   Rotational angular velocity of Earth in radians/sec from IERS
   Conventions (2003).
*/

   extern const double ANGVEL;

/*
   Reciprocal masses of solar system bodies, from DE-405
   (Sun mass / body mass).
   MASS[0] = Earth/Moon barycenter, MASS[1] = Mercury, ...,
   MASS[9] = Pluto, MASS[10] = Sun, MASS[11] = Moon.
*/

   extern const double RMASS[12];

/*
   Value of 2 * pi in radians.
*/

   extern const double TWOPI;

/*
   Number of arcseconds in 360 degrees.
*/

   extern const double ASEC360;

/*
   Angle conversion constants.
*/

   extern const double ASEC2RAD;
   extern const double DEG2RAD;
   extern const double RAD2DEG;

#endif
