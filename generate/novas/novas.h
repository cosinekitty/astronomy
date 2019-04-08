/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1

  novas.h: Header file for novas.c

  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC
  http://www.usno.navy.mil/USNO/astronomical-applications
*/

#ifndef _NOVAS_
   #define _NOVAS_

   #ifndef __STDIO__
      #include <stdio.h>
   #endif

   #ifndef __MATH__
      #include <math.h>
   #endif

   #ifndef __STRING__
      #include <string.h>
   #endif

   #ifndef __STDLIB__
      #include <stdlib.h>
   #endif

   #ifndef __CTYPE__
      #include <ctype.h>
   #endif

   #ifndef _CONSTS_
      #include "novascon.h"
   #endif

   #ifndef _SOLSYS_
      #include "solarsystem.h"
   #endif

   #ifndef _NUTATION_
      #include "nutation.h"
   #endif


/*
   Structures
*/

/*
   struct cat_entry:  basic astrometric data for any celestial object
                      located outside the solar system; the catalog
                      data for a star

   starname[SIZE_OF_OBJ_NAME] = name of celestial object
   catalog[SIZE_OF_CAT_NAME]  = catalog designator (e.g., HIP)
   starnumber                 = integer identifier assigned to object
   ra                         = ICRS right ascension (hours)
   dec                        = ICRS declination (degrees)
   promora                    = ICRS proper motion in right ascension
                                (milliarcseconds/year)
   promodec                   = ICRS proper motion in declination
                                (milliarcseconds/year)
   parallax                   = parallax (milliarcseconds)
   radialvelocity             = radial velocity (km/s)

   SIZE_OF_OBJ_NAME and SIZE_OF_CAT_NAME are defined below.  Each is the
   number of characters in the string (the string length) plus the null
   terminator.
*/

   #define SIZE_OF_OBJ_NAME 51
   #define SIZE_OF_CAT_NAME 4

   typedef struct
   {
      char starname[SIZE_OF_OBJ_NAME];
      char catalog[SIZE_OF_CAT_NAME];
      long int starnumber;
      double ra;
      double dec;
      double promora;
      double promodec;
      double parallax;
      double radialvelocity;
   } cat_entry;

/*
   struct object:    specifies the celestial object of interest

   type              = type of object
                     = 0 ... major planet, Pluto, Sun, or Moon
                     = 1 ... minor planet
                     = 2 ... object located outside the solar system
                             (star, nebula, galaxy, etc.)
   number            = object number
                       For 'type' = 0: Mercury = 1, ..., Pluto = 9,
                                       Sun = 10, Moon = 11
                       For 'type' = 1: minor planet number
                       For 'type' = 2: set to 0 (object is
                       fully specified in 'struct cat_entry')
   name              = name of the object (limited to
                       (SIZE_OF_OBJ_NAME - 1) characters)
   star              = basic astrometric data for any celestial object
                       located outside the solar system; the catalog
                       data for a star
*/

   typedef struct
   {
      short int type;
      short int number;
      char name[SIZE_OF_OBJ_NAME];
      cat_entry star;
   } object;

/*
   struct on_surface: data for an observer's location on the surface of
                      the Earth.  The atmospheric parameters are used
                      only by the refraction function called from
                      function 'equ2hor'. Additional parameters can be
                      added to this structure if a more sophisticated
                      refraction model is employed.

   latitude           = geodetic (ITRS) latitude; north positive (degrees)
   longitude          = geodetic (ITRS) longitude; east positive (degrees)
   height             = height of the observer (meters)
   temperature        = temperature (degrees Celsius)
   pressure           = atmospheric pressure (millibars)
*/

   typedef struct
   {
      double latitude;
      double longitude;
      double height;
      double temperature;
      double pressure;
   } on_surface;

/*
   struct in_space:   data for an observer's location on a near-Earth
                      spacecraft

   sc_pos[3]          = geocentric position vector (x, y, z), components
                        in km
   sc_vel[3]          = geocentric velocity vector (x_dot, y_dot,
                        z_dot), components in km/s

                        Both vectors with respect to true equator and
                        equinox of date
*/

   typedef struct
   {
      double sc_pos[3];
      double sc_vel[3];
   } in_space;

/*
   struct observer:   data specifying the location of the observer

   where              = integer code specifying location of observer
                        = 0: observer at geocenter
                        = 1: observer on surface of earth
                        = 2: observer on near-earth spacecraft
   on_surface         = structure containing data for an observer's
                        location on the surface of the Earth (where = 1)
   near_earth         = data for an observer's location on a near-Earth
                        spacecraft (where = 2)
*/

   typedef struct
   {
      short int where;
      on_surface on_surf;
      in_space near_earth;
   } observer;

/*
   struct sky_pos:    data specifying a celestial object's place on the
                      sky; contains the output from function 'place'

   r_hat[3]           = unit vector toward object (dimensionless)
   ra                 = apparent, topocentric, or astrometric
                        right ascension (hours)
   dec                = apparent, topocentric, or astrometric
                        declination (degrees)
   dis                = true (geometric, Euclidian) distance to solar
                        system body or 0.0 for star (AU)
   rv                 = radial velocity (km/s)
*/

   typedef struct
   {
      double r_hat[3];
      double ra;
      double dec;
      double dis;
      double rv;
   } sky_pos;

/*
   struct ra_of_cio:  right ascension of the Celestial Intermediate
                      Origin (CIO) with respect to the GCRS

   jd_tdb             = TDB Julian date
   ra_cio             = right ascension of the CIO with respect
                        to the GCRS (arcseconds)
*/

   typedef struct
   {
      double jd_tdb;
      double ra_cio;
   } ra_of_cio;


/*
   Define "origin" constants.
*/

   #define BARYC  0
   #define HELIOC 1

/*
   Function prototypes
*/

   double *readeph (int mp, char *name, double jd,

                    int *err);

   short int app_star (double jd_tt, cat_entry *star,
                       short int accuracy,

                       double *ra, double *dec);

   short int virtual_star (double jd_tt, cat_entry *star,
                           short int accuracy,

                           double *ra, double *dec);

   short int astro_star (double jd_tt, cat_entry *star,
                         short int accuracy,

                         double *ra, double *dec);

   short int app_planet (double jd_tt, object *ss_body,
                         short int accuracy,

                         double *ra, double *dec, double *dis);

   short int virtual_planet (double jd_tt, object *ss_body,
                             short int accuracy,

                             double *ra, double *dec, double *dis);

   short int astro_planet (double jd_tt, object *ss_body,
                           short int accuracy,

                           double *ra, double *dec, double *dis);

   short int topo_star (double jd_tt, double delta_t, cat_entry *star,
                        on_surface *position, short int accuracy,

                        double *ra, double *dec);

   short int local_star (double jd_tt, double delta_t, cat_entry *star,
                         on_surface *position, short int accuracy,

                         double *ra, double *dec);

   short int topo_planet (double jd_tt, object *ss_body, double delta_t,
                          on_surface *position, short int accuracy,

                          double *ra, double *dec, double *dis);

   short int local_planet (double jd_tt, object *ss_body,
                           double delta_t, on_surface *position,
                           short int accuracy,

                           double *ra, double *dec, double *dis);

   short int mean_star (double jd_tt, double ra, double dec,
                        short int accuracy,

                        double *ira, double *idec);

   short int place (double jd_tt, object *cel_object,
                    observer *location, double delta_t,
                    short int coord_sys, short int accuracy,

                    sky_pos *output);

   void equ2gal (double rai, double deci,

                 double *glon, double *glat);

   short int equ2ecl (double jd_tt, short int coord_sys,
                      short int accuracy, double ra, double dec,

                      double *elon, double *elat);

   short int equ2ecl_vec (double jd_tt, short int coord_sys,
                          short int accuracy, double *pos1,

                          double *pos2);

   short int ecl2equ_vec (double jd_tt, short int coord_sys,
                          short int accuracy, double *pos1,

                          double *pos2);

   void equ2hor (double jd_ut1, double delta_t, short int accuracy,
                 double xp, double yp, on_surface *location, double ra,
                 double dec, short int ref_option,

                 double *zd, double *az, double *rar, double *decr);

   short int gcrs2equ (double jd_tt, short int coord_sys,
                       short int accuracy, double rag, double decg,

                       double *ra, double *dec);

   short int sidereal_time (double jd_high, double jd_low,
                            double delta_t, short int gst_type,
                            short int method, short int accuracy,

                            double *gst);

   double era (double jd_high, double jd_low);

   short int ter2cel (double jd_ut_high, double jd_ut_low,
                      double delta_t, short int method,
                      short int accuracy, short int option, double xp,
                      double yp, double *vec1,

                      double *vec2);

   short int cel2ter (double jd_ut_high, double jd_ut_low,
                      double delta_t, short int method,
                      short int accuracy, short int option,
                      double xp, double yp, double *vec1,

                      double *vec2);

   void spin (double angle, double *pos1,

              double *pos2);

   void wobble (double tjd, short int direction, double xp, double yp,
                double *pos1,

                double *pos2);

   void terra (on_surface *location, double st,

               double *pos, double *vel);

   void e_tilt (double jd_tdb, short int accuracy,

                double *mobl, double *tobl, double *ee, double *dpsi,
                double *deps);

   short int cel_pole (double tjd, short int type, double dpole1,
                       double dpole2);

   double ee_ct (double jd_high, double jd_low, short int accuracy);

   void frame_tie (double *pos1, short int direction,

                   double *pos2);

   void proper_motion (double jd_tdb1, double *pos, double *vel,
                       double jd_tdb2,

                       double *pos2);

   void bary2obs (double *pos, double *pos_obs,

                  double *pos2, double *lighttime);

   short int geo_posvel (double jd_tt, double delta_t,
                         short int accuracy, observer *obs,

                         double *pos, double *vel);

   short int light_time (double jd_tdb, object *ss_object,
                         double pos_obs[3], double tlight0,
                         short int accuracy,

                         double pos[3], double *tlight);

   double d_light (double *pos1, double *pos_obs);

   short int grav_def (double jd_tdb, short int loc_code,
                       short int accuracy, double *pos1, double *pos_obs,

                       double *pos2);

   void grav_vec (double *pos1, double *pos_obs, double *pos_body,
                  double rmass,

                  double *pos2);

   void aberration (double *pos, double *ve, double lighttime,

                    double *pos2);

   void rad_vel (object *cel_object, double *pos, double *vel,
                 double *vel_obs, double d_obs_geo, double d_obs_sun,
                 double d_obj_sun,

                 double *rv);

   short int precession (double jd_tdb1, double *pos1, double jd_tdb2,

                         double *pos2);

   void nutation (double jd_tdb, short int direction, short int accuracy,
                  double *pos,

                  double *pos2);

   void nutation_angles (double t, short int accuracy,

                         double *dpsi, double *deps);

   void fund_args (double t,

                   double a[5]);

   double mean_obliq (double jd_tdb);

   short int vector2radec (double *pos,

                           double *ra, double *dec);

   void radec2vector (double ra, double dec, double dist,

                      double *vector);

   void starvectors (cat_entry *star,

                     double *pos, double *vel);

   void tdb2tt (double tdb_jd,

                double *tt_jd, double *secdiff);

   short int cio_ra (double jd_tt, short int accuracy,

                     double *ra_cio);

   short int cio_location (double jd_tdb, short int accuracy,

                           double *ra_cio, short int *ref_sys);

   short int cio_basis (double jd_tdb, double ra_cio, short int ref_sys,
                        short int accuracy,

                        double *x, double *y, double *z);

   short int cio_array (double jd_tdb, long int n_pts,

                        ra_of_cio *cio);

   double ira_equinox (double jd_tdb, short int equinox,
                       short int accuracy);

   short int ephemeris (double jd[2], object *cel_obj, short int origin,
                        short int accuracy,

                        double *pos, double *vel);

   void transform_hip (cat_entry *hipparcos,

                       cat_entry *hip_2000);

   short int transform_cat (short int option, double date_incat,
                            cat_entry *incat, double date_newcat,
                            char newcat_id[SIZE_OF_CAT_NAME],

                            cat_entry *newcat);

   void limb_angle (double pos_obj[3], double pos_obs[3],

                    double *limb_ang, double *nadir_ang);

   double refract (on_surface *location, short int ref_option,
                   double zd_obs);

   double julian_date (short int year, short int month, short int day,
                       double hour);

   void cal_date (double tjd,

                  short int *year, short int *month, short int *day,
                  double *hour);

   double norm_ang (double angle);

   short int make_cat_entry (char star_name[SIZE_OF_OBJ_NAME],
                             char catalog[SIZE_OF_CAT_NAME],
                        long int star_num, double ra, double dec,
                        double pm_ra, double pm_dec, double parallax,
                        double rad_vel,

                        cat_entry *star);

   short int make_object (short int type, short int number,
                          char name[SIZE_OF_OBJ_NAME],
                          cat_entry *star_data,

                          object *cel_obj);

   short int make_observer (short int where, on_surface *obs_surface,
                            in_space *obs_space,

                            observer *obs);

   void make_observer_at_geocenter (

                                 observer *obs_at_geocenter);

   void make_observer_on_surface (double latitude, double longitude,
                                  double height, double temperature,
                                  double pressure,

                                  observer *obs_on_surface);

   void make_observer_in_space (double sc_pos[3], double sc_vel[3],

                                observer *obs_in_space);

   void make_on_surface (double latitude, double longitude,
                         double height,
                         double temperature, double pressure,

                         on_surface *obs_surface);

   void make_in_space (double sc_pos[3], double sc_vel[3],

                       in_space *obs_space);


#endif
