/*
    galactic.c  -  Don Cross  -  2021-06-10

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This program illustrates how to convert a location
    in the sky expressed in IAU 1958 galactic coordinates
    into the local altitude and azimuth of someone wanting
    to aim a radio dish at it.
*/

#include <stdio.h>
#include <math.h>
#include "astro_demo_common.h"


int GalaticToHorizontal(
  astro_time_t time,
  astro_observer_t observer,
  double glat,
  double glon,
  double *altitude,
  double *azimuth) {
  astro_rotation_t rot, adjust_rot;
  astro_spherical_t gsphere, hsphere;
  astro_vector_t gvec, hvec;

  /*
        Calculate a rotation matrix that converts
        galactic coordinates to J2000 equatorial coordinates.
    */
  rot = Astronomy_Rotation_GAL_EQJ();

  /*
        Adjust the rotation matrix to convert galatic to horizontal (HOR).
    */
  adjust_rot = Astronomy_Rotation_EQJ_HOR(&time, observer);
  rot = Astronomy_CombineRotation(rot, adjust_rot);

  /*
        Convert the galactic coordinates from angles to a unit vector.
    */
  gsphere.status = ASTRO_SUCCESS;
  gsphere.lat = glat;
  gsphere.lon = glon;
  gsphere.dist = 1.0;
  gvec = Astronomy_VectorFromSphere(gsphere, time);
  if (gvec.status != ASTRO_SUCCESS) {
    Serial.printf("Astronomy_VectorFromSphere returned error %d\n", gvec.status);
    return 1;
  }

  /*
        Use the rotation matrix to convert the galactic vector to a horizontal vector.
    */
  hvec = Astronomy_RotateVector(rot, gvec);
  if (hvec.status != ASTRO_SUCCESS) {
    Serial.printf("Astronomy_RotateVector returned error %d\n", hvec.status);
    return 1;
  }

  /*
        Convert the horizontal vector back to angular coordinates: altitude and azimuth.
        Assuming this is a radio source (not optical), do not correct for refraction.
    */
  hsphere = Astronomy_HorizonFromVector(hvec, REFRACTION_NONE);
  if (hsphere.status != ASTRO_SUCCESS) {
    Serial.printf("Astronomy_HorizonFromVector returned error %d\n", hsphere.status);
    return 1;
  }

  *altitude = hsphere.lat;
  *azimuth = hsphere.lon;
  return 0;
}

void setup() {
  Serial.begin(9600);
}

void loop() {
  astro_observer_t observer;
  astro_time_t time;
  double glat, glon;
  double azimuth, altitude;


  const char *argvs[] = { "22.2", "33.8", "2023-10-07T00:38:57+03:30" };
  glat = 23.3;  //-90.0 to 90.0
  glon = 17.6;  //-90.0 to 90.0
  int argcs = 3;
  int error;
  error = ParseArgs(argcs, argvs, &observer, &time);


  if (glat < -90.0 || glat > +90.0) {
    Serial.printf("ERROR: Invalid galatic latitude '%s' on command line\n", glat);
    error = 1;
  }


  if (glon <= -360.0 || glon >= +360.0) {
    Serial.printf("ERROR: Invalid galatic longitude '%s' on command line\n", glon);
    error = 1;
  }

  if (GalaticToHorizontal(time, observer, glat, glon, &altitude, &azimuth)) {
    Serial.println("error.......");
  }
  if (error) {
    Serial.println("error.......");
  }
  Serial.printf("altitude = %10.3lf, azimuth = %10.3lf\n", altitude, azimuth);

  delay(1000);
}
