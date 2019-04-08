/*
   Naval Observatory Vector Astrometry Software (NOVAS)
   C Edition, Version 3.1
   
   nutation.h: Header file for nutation models

   U. S. Naval Observatory
   Astronomical Applications Dept.
   Washington, DC 
   http://www.usno.navy.mil/USNO/astronomical-applications
*/

#ifndef _NUTATION_
   #define _NUTATION_

/*
   Function prototypes
*/

   void iau2000a (double jd_high, double jd_low,

                  double *dpsi, double *deps);

   void iau2000b (double jd_high, double jd_low,

                  double *dpsi, double *deps);

   void nu2000k (double jd_high, double jd_low,

                 double *dpsi, double *deps);

#endif
