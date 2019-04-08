/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1

  solarsystem.h: Header file for solsys1.c, solsys2.c, & solsys3.c

  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC
  http://www.usno.navy.mil/USNO/astronomical-applications
*/

#ifndef _SOLSYS_
   #define _SOLSYS_

/*
   Function prototypes
*/

   short int solarsystem (double tjd, short int body, short int origin,

                          double *position, double *velocity);

   short int solarsystem_hp (double tjd[2], short body, short origin,

                             double *position, double *velocity);


#endif
