/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1

  eph_manager.c: C version of JPL Ephemeris Manager for use with solsys1

  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC
  http://www.usno.navy.mil/USNO/astronomical-applications
*/

#ifndef _EPHMAN_
   #include "eph_manager.h"
#endif

/*
   Define global variables
*/

short int KM;

/*
   IPT and LPT defined as int to support 64 bit systems.
*/

int IPT[3][12], LPT[3];

long int NRL, NP, NV;
long int RECORD_LENGTH;

double SS[3], JPLAU, PC[18], VC[18], TWOT, EM_RATIO;
double *BUFFER;

FILE *EPHFILE = NULL;

/********ephem_open */

short int ephem_open (char *ephem_name,

                      double *jd_begin, double *jd_end,
                      short int *de_number)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function opens a JPL planetary ephemeris file and
      sets initial values.  This function must be called
      prior to calls to the other JPL ephemeris functions.

   REFERENCES:
      Standish, E.M. and Newhall, X X (1988). "The JPL Export
         Planetary Ephemeris"; JPL document dated 17 June 1988.

   INPUT
   ARGUMENTS:
      *ephem_name (char)
         Name of the direct-access ephemeris file.

   OUTPUT
   ARGUMENTS:
      *jd_begin (double)
         Beginning Julian date of the ephemeris file.
      *jd_end (double)
         Ending Julian date of the ephemeris file.
      *de_number (short int)
         DE number of the ephemeris file opened.

   RETURNED
   VALUE:
      (short int)
          0   ...file exists and is opened correctly.
          1   ...file does not exist/not found.
          2-10...error reading from file header.
          11  ...unable to set record length; ephemeris (DE number)
                 not in look-up table.

   GLOBALS
   USED:
      SS                eph_manager.h
      JPLAU             eph_manager.h
      PC                eph_manager.h
      VC                eph_manager.h
      TWOT              eph_manager.h
      EM_RATIO          eph_manager.h
      BUFFER            eph_manager.h
      IPT               eph_manager.h
      LPT               eph_manager.h
      NRL               eph_manager.h
      KM                eph_manager.h
      NP                eph_manager.h
      NV                eph_manager.h
      RECORD_LENGTH     eph_manager.h
      EPHFILE           eph_manager.h

   FUNCTIONS
   CALLED:
      fclose            stdio.h
      free              stdlib.h
      fopen             stdio.h
      fread             stdio.h
      calloc            stdlib.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-90/JAB (USNO/NA)
      V1.1/06-92/JAB (USNO/AA): Restructure and add initializations.
      V1.2/07-98/WTH (USNO/AA): Modified to open files for different
                                ephemeris types. (200,403,404,405,406)
      V1.3/11-07/WKP (USNO/AA): Updated prolog.
      V1.4/09-10/WKP (USNO/AA): Changed ncon and denum variables and
                                sizeof ipt array to type 'int' for
                                64-bit system compatibility.
      V1.5/09-10/WTH (USNO/AA): Added support for DE421, default case
                                for switch, close file on error.
      V1.6/10-10/WKP (USNO/AA): Renamed function to lowercase to
                                comply with coding standards.

   NOTES:
      KM...flag defining physical units of the output states.
         = 1, km and km/sec
         = 0, AU and AU/day
      Default value is 0 (KM determines time unit for nutations.
                          Angle unit is always radians.)

------------------------------------------------------------------------
*/
{
   char ttl[252], cnam[2400];

   short int i, j;

   int ncon, denum;

   ephem_close();

/*
   Open file ephem_name.
*/

   if ((EPHFILE = fopen (ephem_name, "rb")) == NULL)
   {
      return 1;
   }
    else
   {

/*
   File found. Set initializations and default values.
*/

      KM = 0;

      NRL = 0;

      NP = 2;
      NV = 3;
      TWOT = 0.0;

      for (i = 0; i < 18; i++)
      {
         PC[i] = 0.0;
         VC[i] = 0.0;
      }

      PC[0] = 1.0;
      VC[1] = 1.0;

/*
   Read in values from the first record, aka the header.
*/

      if (fread (ttl, sizeof ttl, 1, EPHFILE) != 1)
      {
         ephem_close();
         return 2;
      }
      if (fread (cnam, sizeof cnam, 1, EPHFILE) != 1)
      {
         ephem_close();
         return 3;
      }
      if (fread (SS, sizeof SS, 1, EPHFILE) != 1)
      {
         ephem_close();
         return 4;
      }
      if (fread (&ncon, sizeof ncon, 1, EPHFILE) != 1)
      {
         ephem_close();
         return 5;
      }
      if (fread (&JPLAU, sizeof JPLAU, 1, EPHFILE) != 1)
      {
         ephem_close();
         return 6;
      }
      if (fread (&EM_RATIO, sizeof EM_RATIO, 1, EPHFILE) != 1)
      {
         ephem_close();
         return 7;
      }
      for (i = 0; i < 12; i++)
         for (j = 0; j < 3; j++)
            if (fread (&IPT[j][i], sizeof(int), 1, EPHFILE) != 1)
            {
               ephem_close();
               return 8;
            }
      if (fread (&denum, sizeof denum, 1, EPHFILE) != 1)
      {
         ephem_close();
         return 9;
      }
      if (fread (LPT, sizeof LPT, 1, EPHFILE) != 1)
      {
         ephem_close();
         return 10;
      }

/*
   Set the value of the record length according to what JPL ephemeris is
   being opened.
*/

      switch (denum)
      {
         case 200:
            RECORD_LENGTH = 6608;
            break;
         case 403: case 405:
         case 421:
            RECORD_LENGTH = 8144;
            break;
         case 404: case 406:
            RECORD_LENGTH = 5824;
            break;

/*
   An unknown DE file was opened. Close the file and return an error
   code.
*/

         default:
            *jd_begin = 0.0;
            *jd_end = 0.0;
            *de_number = 0;
            ephem_close();
            return 11;
      }

      BUFFER = (double *) calloc (RECORD_LENGTH / 8, sizeof(double));

      *de_number = (short int) denum;
      *jd_begin = SS[0];
      *jd_end = SS[1];
   }

   return 0;
}

/********ephem_close */

short int ephem_close (void)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function closes a JPL planetary ephemeris file and frees the
      memory.

   REFERENCES:
      Standish, E.M. and Newhall, X X (1988). "The JPL Export
         Planetary Ephemeris"; JPL document dated 17 June 1988.

   INPUT
   ARGUMENTS:
      None.

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (short int)
          0  ...file was already closed or closed correctly.
          EOF...error closing file.

   GLOBALS
   USED:
      BUFFER            eph_manager.h
      EPHFILE           eph_manager.h

   FUNCTIONS
   CALLED:
      fclose            stdio.h
      free              stdlib.h

   VER./DATE/
   PROGRAMMER:
      V1.0/11-07/WKP (USNO/AA)
      V1.1/09-10/WKP (USNO/AA): Explicitly cast fclose return value to
                                type 'short int'.
      V1.2/10-10/WKP (USNO/AA): Renamed function to lowercase to
                                comply with coding standards.

   NOTES:
      None.

------------------------------------------------------------------------
*/
{
   short int error = 0;

   if (EPHFILE)
   {
      error =  (short int) fclose (EPHFILE);
      EPHFILE = NULL;   /* modified per errata in https://aa.usno.navy.mil/software/novas/novas_faq.php */
   }

   if (BUFFER)
   {
      free (BUFFER);
      BUFFER = NULL;
   }

   return error;
}

/********planet_ephemeris */

short int planet_ephemeris (double tjd[2], short int target,
                            short int center,

                            double *position, double *velocity)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function accesses the JPL planetary ephemeris to give the
      position and velocity of the target object with respect to the
      center object.

   REFERENCES:
      Standish, E.M. and Newhall, X X (1988). "The JPL Export
         Planetary Ephemeris"; JPL document dated 17 June 1988.

   INPUT
   ARGUMENTS:
      tjd[2] (double)
         Two-element array containing the Julian date, which may be
         split any way (although the first element is usually the
         "integer" part, and the second element is the "fractional"
         part).  Julian date is in the TDB or "T_eph" time scale.
      target (short int)
         Number of 'target' point.
      center (short int)
         Number of 'center' (origin) point.
         The numbering convention for 'target' and'center' is:
            0 = Mercury           7 = Neptune
            1 = Venus             8 = Pluto
            2 = Earth             9 = Moon
            3 = Mars             10 = Sun
            4 = Jupiter          11 = Solar system bary.
            5 = Saturn           12 = Earth-Moon bary.
            6 = Uranus           13 = Nutations (long int. and obliq.)
            (If nutations are desired, set 'target' = 13;
             'center' will be ignored on that call.)

   OUTPUT
   ARGUMENTS:
      *position (double)
         Position vector array of target relative to center, measured
         in AU.
      *velocity (double)
         Velocity vector array of target relative to center, measured
         in AU/day.

   RETURNED
   VALUE:
      (short int)
         0  ...everything OK.
         1,2...error returned from State.

   GLOBALS
   USED:
      EM_RATIO          eph_manager.h

   FUNCTIONS
   CALLED:
      state             eph_manager.h

   VER./DATE/
   PROGRAMMER:
      V1.0/03-93/WTH (USNO/AA): Convert FORTRAN to C.
      V1.1/07-93/WTH (USNO/AA): Update to C standards.
      V2.0/07-98/WTH (USNO/AA): Modified for ease of use and linearity.
      V3.0/11-06/JAB (USNO/AA): Allowed for use of input 'split' Julian
                                date for higher precision.
      V3.1/11-07/WKP (USNO/AA): Updated prolog and error codes.
      V3.1/12-07/WKP (USNO/AA): Removed unreferenced variables.
      V3.2/10-10/WKP (USNO/AA): Renamed function to lowercase to
                                comply with coding standards.

   NOTES:
      None.

------------------------------------------------------------------------
*/
{
   short int i, error = 0, earth = 2, moon = 9;
   short int do_earth = 0, do_moon = 0;

   double jed[2];
   double pos_moon[3] = {0.0,0.0,0.0}, vel_moon[3] = {0.0,0.0,0.0},
          pos_earth[3] = {0.0,0.0,0.0}, vel_earth[3] = {0.0,0.0,0.0};
   double target_pos[3] = {0.0,0.0,0.0}, target_vel[3] = {0.0,0.0,0.0},
          center_pos[3] = {0.0,0.0,0.0}, center_vel[3] = {0.0,0.0,0.0};

/*
   Initialize 'jed' for 'state' and set up component count.
*/

   jed[0] = tjd[0];
   jed[1] = tjd[1];

/*
   Check for target point = center point.
*/

   if (target == center)
   {
      for (i = 0; i < 3; i++)
      {
         position[i] = 0.0;
         velocity[i] = 0.0;
      }
      return 0;
   }

/*
   Check for instances of target or center being Earth or Moon,
   and for target or center being the Earth-Moon barycenter.
*/

   if ((target == earth) || (center == earth))
      do_moon = 1;
   if ((target == moon) || (center == moon))
      do_earth = 1;
   if ((target == 12) || (center == 12))
      do_earth = 1;

   if (do_earth)
   {
      error = state (jed,2, pos_earth,vel_earth);
      if (error)
         return error;
   }

   if (do_moon)
   {
      error = state (jed,9, pos_moon,vel_moon);
      if (error)
         return error;
   }

/*
   Make call to State for target object.
*/

   if (target == 11)
   {
      for (i = 0; i < 3; i++)
      {
         target_pos[i] = 0.0;
         target_vel[i] = 0.0;
      }
   }
    else if (target == 12)
   {
      for (i = 0; i < 3; i++)
      {
         target_pos[i] = pos_earth[i];
         target_vel[i] = vel_earth[i];
      }
   }
    else
      error = state (jed,target, target_pos,target_vel);

   if (error)
      return error;

/*
   Make call to State for center object.
*/

/*
   If the requested center is the Solar System barycenter,
   then don't bother with a second call to State.
*/

   if (center == 11)
   {
      for (i = 0; i < 3; i++)
      {
         center_pos[i] = 0.0;
         center_vel[i] = 0.0;
      }
   }

/*
   Center is Earth-Moon barycenter, which was already computed above.
*/

    else if (center == 12)
   {
      for (i = 0; i < 3; i++)
      {
         center_pos[i] = pos_earth[i];
         center_vel[i] = vel_earth[i];
      }
   }
    else
      error = state (jed,center, center_pos,center_vel);

   if (error)
      return error;

/*
   Check for cases of Earth as target and Moon as center or vice versa.
*/

   if ((target == earth) && (center == moon))
   {
      for (i = 0; i < 3; i++)
      {
         position[i] = -center_pos[i];
         velocity[i] = -center_vel[i];
      }
      return 0;
   }
    else if ((target == moon) && (center == earth))
   {
      for (i = 0; i < 3; i++)
      {
         position[i] = target_pos[i];
         velocity[i] = target_vel[i];
      }
      return 0;
   }

/*
   Check for Earth as target, or as center.
*/

    else if (target == earth)
   {
      for (i = 0; i < 3; i++)
      {
         target_pos[i] = target_pos[i] - (pos_moon[i] /
            (1.0 + EM_RATIO));
         target_vel[i] = target_vel[i] - (vel_moon[i] /
            (1.0 + EM_RATIO));
      }
   }
    else if (center == earth)
   {
      for (i = 0; i < 3; i++)
      {
         center_pos[i] = center_pos[i] - (pos_moon[i] /
            (1.0 + EM_RATIO));
         center_vel[i] = center_vel[i] - (vel_moon[i] /
            (1.0 + EM_RATIO));
      }
   }

/*
   Check for Moon as target, or as center.
*/

    else if (target == moon)
   {
      for (i = 0; i < 3; i++)
      {
         target_pos[i] = (pos_earth[i] - (target_pos[i] /
            (1.0 + EM_RATIO))) + target_pos[i];
         target_vel[i] = (vel_earth[i] - (target_vel[i] /
            (1.0 + EM_RATIO))) + target_vel[i];
      }
   }
    else if (center == moon)
   {
      for (i = 0; i < 3; i++)
      {
         center_pos[i] = (pos_earth[i] - (center_pos[i] /
            (1.0 + EM_RATIO))) + center_pos[i];
         center_vel[i] = (vel_earth[i] - (center_vel[i] /
            (1.0 + EM_RATIO))) + center_vel[i];
      }
   }

/*
   Compute position and velocity vectors.
*/

   for (i = 0; i < 3; i++)
   {
      position[i] = target_pos[i] - center_pos[i];
      velocity[i] = target_vel[i] - center_vel[i];
   }

   return 0;
}


/********state */

short int state (double *jed, short int target,

                 double *target_pos, double *target_vel)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function reads and interpolates the JPL planetary
      ephemeris file.

   REFERENCES:
      Standish, E.M. and Newhall, X X (1988). "The JPL Export
         Planetary Ephemeris"; JPL document dated 17 June 1988.

   INPUT
   ARGUMENTS:
      *jed (double)
         2-element Julian date (TDB) at which interpolation is wanted.
         Any combination of jed[0]+jed[1] which falls within the time
         span on the file is a permissible epoch.  See Note 1 below.
      target (short int)
         The requested body to get data for from the ephemeris file.
         The designation of the astronomical bodies is:
                 0 = Mercury                    6 = Uranus
                 1 = Venus                      7 = Neptune
                 2 = Earth-Moon barycenter      8 = Pluto
                 3 = Mars                       9 = geocentric Moon
                 4 = Jupiter                   10 = Sun
                 5 = Saturn

   OUTPUT
   ARGUMENTS:
      *target_pos (double)
         The barycentric position vector array of the requested object,
         in AU.
         (If target object is the Moon, then the vector is geocentric.)
      *target_vel (double)
         The barycentric velocity vector array of the requested object,
         in AU/Day.

         Both vectors are referenced to the Earth mean equator and
         equinox of epoch.

   RETURNED
   VALUE:
      (short int)
         0...everything OK.
         1...error reading ephemeris file.
         2...epoch out of range.
         3...target out of range.

   GLOBALS
   USED:
      KM                eph_manager.h
      EPHFILE           eph_manager.h
      IPT               eph_manager.h
      BUFFER            eph_manager.h
      NRL               eph_manager.h
      RECORD_LENGTH     eph_manager.h
      SS                eph_manager.h
      JPLAU             eph_manager.h

   FUNCTIONS
   CALLED:
      split             eph_manager.h
      fseek             stdio.h
      fread             stdio.h
      interpolate       eph_manager.h
      ephem_close       eph_manager.h

   VER./DATE/
   PROGRAMMER:
      V1.0/03-93/WTH (USNO/AA): Convert FORTRAN to C.
      V1.1/07-93/WTH (USNO/AA): Update to C standards.
      V2.0/07-98/WTH (USNO/AA): Modify to make position and velocity
                                two distinct vector arrays.  Routine set
                                to compute one state per call.
      V2.1/11-07/WKP (USNO/AA): Updated prolog.
      V2.2/10-10/WKP (USNO/AA): Renamed function to lowercase to
                                comply with coding standards.

   NOTES:
      1. For ease in programming, the user may put the entire epoch in
         jed[0] and set jed[1] = 0. For maximum interpolation accuracy,
         set jed[0] = the most recent midnight at or before
         interpolation epoch, and set jed[1] = fractional part of a day
         elapsed between jed[0] and epoch. As an alternative, it may
         prove convenient to set jed[0] = some fixed epoch, such as
         start of the integration and jed[1] = elapsed interval between
         then and epoch.

------------------------------------------------------------------------
*/
{
   short int i;

   long int nr, rec;

   double t[2], aufac = 1.0, jd[4], s;

/*
   Check for invalid targets.
*/
   if (target < 0 || target > 10)
      return 3;

/*
   Set units based on value of the 'KM' flag.
*/

   if (KM)
      t[1] = SS[2] * 86400.0;
    else
   {
      t[1] = SS[2];
      aufac = 1.0 / JPLAU;
   }

/*
   Check epoch.
*/

   s = jed[0] - 0.5;
   split (s, &jd[0]);
   split (jed[1], &jd[2]);
   jd[0] += jd[2] + 0.5;
   jd[1] += jd[3];
   split (jd[1], &jd[2]);
   jd[0] += jd[2];

/*
   Return error code if date is out of range.
*/

   if ((jd[0] < SS[0]) || ((jd[0] + jd[3]) > SS[1]))
      return 2;

/*
   Calculate record number and relative time interval.
*/

   nr = (long int) ((jd[0] - SS[0]) / SS[2]) + 3;
   if (jd[0] == SS[1])
      nr -= 2;
   t[0] = ((jd[0] - ((double) (nr-3) * SS[2] + SS[0])) + jd[3]) / SS[2];

/*
   Read correct record if it is not already in memory.
*/

   if (nr != NRL)
   {
      NRL = nr;
      rec = (nr - 1) * RECORD_LENGTH;
      fseek (EPHFILE, rec, SEEK_SET);
      if (!fread (BUFFER, RECORD_LENGTH, 1, EPHFILE))
      {
         ephem_close ();
         return 1;
      }
   }

/*
   Check and interpolate for requested body.
*/

   interpolate (&BUFFER[IPT[0][target]-1],t,IPT[1][target],
      IPT[2][target], target_pos,target_vel);

   for (i = 0; i < 3; i++)
   {
      target_pos[i] *= aufac;
      target_vel[i] *= aufac;
   }

   return 0;
}

/********interpolate */

void interpolate (double *buf, double *t, long int ncf, long int na,

                  double *position, double *velocity)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function differentiates and interpolates a set of
      Chebyshev coefficients to give position and velocity.

   REFERENCES:
      Standish, E.M. and Newhall, X X (1988). "The JPL Export
         Planetary Ephemeris"; JPL document dated 17 June 1988.

   INPUT
   ARGUMENTS:
      *buf (double)
         Array of Chebyshev coefficients of position.
      *t (double)
         t[0] is fractional time interval covered by coefficients at
         which interpolation is desired (0 <= t[0] <= 1).
         t[1] is length of whole interval in input time units.
      ncf (long int)
         Number of coefficients per component.
      na (long int)
         Number of sets of coefficients in full array
         (i.e., number of sub-intervals in full interval).

   OUTPUT
   ARGUMENTS:
      *position (double)
         Position array of requested object.
      *velocity (double)
         Velocity array of requested object.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      NP                eph_manager.h
      NV                eph_manager.h
      PC                eph_manager.h
      VC                eph_manager.h
      TWOT              eph_manager.h

   FUNCTIONS
   CALLED:
      fmod              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/03-93/WTH (USNO/AA): Convert FORTRAN to C.
      V1.1/07-93/WTH (USNO/AA): Update to C standards.
      V1.2/07-98/WTH (USNO/AA): Modify to make position and velocity
                                two distinct vector arrays.
      V1.3/11-07/WKP (USNO/AA): Updated prolog.
      V1.4/12-07/WKP (USNO/AA): Changed ncf and na arguments from short
                                int to long int.
      V1.5/10-10/WKP (USNO/AA): Renamed function to lowercase to
                                comply with coding standards.

   NOTES:
      None.

------------------------------------------------------------------------
*/
{
   long int i, j, k, l;

   double dna, dt1, temp, tc, vfac;

/*
   Get correct sub-interval number for this set of coefficients and
   then get normalized Chebyshev time within that subinterval.
*/

   dna = (double) na;
   dt1 = (double) ((long int) t[0]);
   temp = dna * t[0];
   l = (long int) (temp - dt1);

/*
   'tc' is the normalized Chebyshev time (-1 <= tc <= 1).
*/

   tc = 2.0 * (fmod (temp, 1.0) + dt1) - 1.0;

/*
   Check to see whether Chebyshev time has changed, and compute new
   polynomial values if it has.  (The element PC[1] is the value of
   t1[tc] and hence contains the value of 'tc' on the previous call.)
*/

   if (tc != PC[1])
   {
      NP = 2;
      NV = 3;
      PC[1] = tc;
      TWOT = tc + tc;
   }

/*
   Be sure that at least 'ncf' polynomials have been evaluated and
   are stored in the array 'PC'.
*/

   if (NP < ncf)
   {
      for (i = NP; i < ncf; i++)
         PC[i] = TWOT * PC[i-1] - PC[i-2];
      NP = ncf;
   }

/*
   Interpolate to get position for each component.
*/

   for (i = 0; i < 3; i++)
   {
      position[i] = 0.0;
      for (j = ncf-1; j >= 0; j--)
      {
         k = j + (i * ncf) + (l * (3 * ncf));
         position[i] += PC[j] * buf[k];
      }
   }

/*
   If velocity interpolation is desired, be sure enough derivative
   polynomials have been generated and stored.
*/

   vfac = (2.0 * dna) / t[1];
   VC[2] = 2.0 * TWOT;
   if (NV < ncf)
   {
      for (i = NV; i < ncf; i++)
         VC[i] = TWOT * VC[i-1] + PC[i-1] + PC[i-1] - VC[i-2];
      NV = ncf;
   }

/*
   Interpolate to get velocity for each component.
*/

   for (i = 0; i < 3; i++)
   {
      velocity[i] = 0.0;
      for (j = ncf-1; j > 0; j--)
      {
         k = j + (i * ncf) + (l * (3 * ncf));
         velocity[i] += VC[j] * buf[k];
      }
      velocity[i] *= vfac;
   }

   return;
}

/********split */

void split (double tt,

            double *fr)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function breaks up a double number into a double integer
      part and a fractional part.

   REFERENCES:
      Standish, E.M. and Newhall, X X (1988). "The JPL Export
         Planetary Ephemeris"; JPL document dated 17 June 1988.

   INPUT
   ARGUMENTS:
      tt (double)
         Input number.

   OUTPUT
   ARGUMENTS:
      *fr (double)
         2-element output array;
            fr[0] contains integer part,
            fr[1] contains fractional part.
         For negative input numbers,
            fr[0] contains the next more negative integer;
            fr[1] contains a positive fraction.

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
      V1.0/06-90/JAB (USNO/NA): CA coding standards
      V1.1/03-93/WTH (USNO/AA): Convert to C.
      V1.2/07-93/WTH (USNO/AA): Update to C standards.
      V1.3/10-10/WKP (USNO/AA): Renamed function to lowercase to
                                comply with coding standards.

   NOTES:
      None.

------------------------------------------------------------------------
*/
{

/*
   Get integer and fractional parts.
*/

   fr[0] = (double)((long int) tt);
   fr[1] = tt - fr[0];

/*
   Make adjustments for negative input number.
*/

   if ((tt >= 0.0) || (fr[1] == 0.0))
      return;
    else
   {
      fr[0] = fr[0] - 1.0;
      fr[1] = fr[1] + 1.0;
   }

   return;
}
