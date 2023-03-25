/*
    camera.c  -  Don Cross  -  2021-03-21

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Suppose you want to photograph the Moon,
    and you want to know what it will look like in the photo.
    Given a location on the Earth, and a date/time,
    this program calculates the orientation of the sunlit
    side of the Moon with respect to the top of your
    photo image. It assumes the camera faces directly
    toward the Moon's azimuth and tilts upward to its
    altitude angle above the horizon.
*/

#include <stdio.h>
#include <math.h>
#include "astro_demo_common.h"

#define FAIL(...)   do{fprintf(stderr, __VA_ARGS__); error = 1; goto fail;}while(0)

static int CameraImage(astro_observer_t observer, astro_time_t time);

int main(int argc, const char *argv[])
{
    int error;
    astro_observer_t observer;
    astro_time_t time;

    error = ParseArgs(argc, argv, &observer, &time);
    if (error)
        return error;

    return CameraImage(observer, time);
}


static int CameraImage(astro_observer_t observer, astro_time_t time)
{
    int error;
    astro_equatorial_t moon_equ;
    astro_equatorial_t sun_equ;
    astro_horizon_t moon_hor;
    astro_rotation_t rot;
    astro_vector_t vec;
    astro_illum_t illum;
    astro_angle_result_t angle;
    double radius;
    double tilt;
    const double tolerance = 1.0e-15;

    /* Calculate the topocentric equatorial coordinates of date for the Moon. */
    /* Assume aberration does not matter because the Moon is so close and has such a small relative velocity. */
    moon_equ = Astronomy_Equator(BODY_MOON, &time, observer, EQUATOR_OF_DATE, NO_ABERRATION);
    if (moon_equ.status != ASTRO_SUCCESS)
        FAIL("Error %d calculating Moon position.\n", moon_equ.status);

    /* Also calculate the Sun's topocentric position in the same coordinate system. */
    sun_equ = Astronomy_Equator(BODY_SUN, &time, observer, EQUATOR_OF_DATE, NO_ABERRATION);
    if (sun_equ.status != ASTRO_SUCCESS)
        FAIL("Error %d calculating Sun position.\n", sun_equ.status);

    /* Get the Moon's horizontal coordinates, so we know how much to pivot azimuth and altitude. */
    moon_hor = Astronomy_Horizon(&time, observer, moon_equ.ra, moon_equ.dec, REFRACTION_NONE);
    printf("Moon horizontal position: azimuth = %0.3lf, altitude = %0.3lf\n", moon_hor.azimuth, moon_hor.altitude);

    /* Get the rotation matrix that converts equatorial to horizontal coordintes for this place and time. */
    rot = Astronomy_Rotation_EQD_HOR(&time, observer);

    /*
        Modify the rotation matrix in two steps:
        First, rotate the orientation so we are facing the Moon's azimuth.
        We do this by pivoting around the zenith axis.
        Horizontal axes are: 0 = north, 1 = west, 2 = zenith.
        Tricky: because the pivot angle increases counterclockwise, and azimuth
        increases clockwise, we undo the azimuth by adding the positive value.
    */
    rot = Astronomy_Pivot(rot, 2, moon_hor.azimuth);
    if (rot.status != ASTRO_SUCCESS)
        FAIL("Error %d in Astronomy_Pivot(azimuth)\n", rot.status);

    /*
        Second, pivot around the leftward axis to bring the Moon to the camera's altitude level.
        From the point of view of the leftward axis, looking toward the camera,
        adding the angle is the correct sense for subtracting the altitude.
    */
    rot = Astronomy_Pivot(rot, 1, moon_hor.altitude);
    if (rot.status != ASTRO_SUCCESS)
        FAIL("Error %d in Astronomy_Pivot(altitude)\n", rot.status);

    /* As a sanity check, apply this rotation to the Moon's equatorial (EQD) coordinates and verify x=0, y=0. */
    vec = Astronomy_RotateVector(rot, moon_equ.vec);
    if (vec.status != ASTRO_SUCCESS)
        FAIL("Error %d in Astronomy_RotateVector(moon)\n", vec.status);

    /* Convert to unit vector. */
    radius = Astronomy_VectorLength(vec);
    vec.x /= radius;
    vec.y /= radius;
    vec.z /= radius;
    printf("Moon check: x = %0.6lf, y = %0.6lf, z = %0.6lf\n", vec.x, fabs(vec.y), fabs(vec.z));
    if (!isfinite(vec.x) || fabs(vec.x - 1.0) > tolerance)
        FAIL("Excessive error in moon check (x)\n");

    if (!isfinite(vec.y) || fabs(vec.y) > tolerance)
        FAIL("Excessive error in moon check (y)\n");

    if (!isfinite(vec.z) || fabs(vec.z) > tolerance)
        FAIL("Excessive error in moon check (z)\n");

    /* Apply the same rotation to the Sun's equatorial vector. */
    /* The x- and y-coordinates now tell us which side appears sunlit in the camera! */

    vec = Astronomy_RotateVector(rot, sun_equ.vec);
    if (vec.status != ASTRO_SUCCESS)
        FAIL("Error %d in Astronomy_RotateVector(sun)\n", vec.status);

    /* Don't bother normalizing the Sun vector, because in AU it will be close to unit anyway. */
    printf("Sun vector: x = %0.6lf, y = %0.6lf, z = %0.6lf\n", vec.x, vec.y, vec.z);

    /* Calculate the tilt angle of the sunlit side, as seen by the camera. */
    /* The x-axis is now pointing directly at the object, z is up in the camera image, y is to the left. */
    tilt = RAD2DEG * atan2(vec.y, vec.z);
    printf("Tilt angle of sunlit side of the Moon = %0.3lf degrees counterclockwise from up.\n", tilt);

    illum = Astronomy_Illumination(BODY_MOON, time);
    if (illum.status != ASTRO_SUCCESS)
        FAIL("Error %d trying to calculate Moon illumination.\n", illum.status);

    printf("Moon magnitude = %0.2lf, phase angle = %0.2lf degrees.\n", illum.mag, illum.phase_angle);

    angle = Astronomy_AngleFromSun(BODY_MOON, time);
    if (angle.status != ASTRO_SUCCESS)
        FAIL("Error %d trying to calculate angle between Moon and Sun\n", angle.status);

    printf("Angle between Moon and Sun as seen from Earth = %0.2lf degrees.\n", angle.angle);

    error = 0;
fail:
    return error;
}
