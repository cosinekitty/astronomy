/*
    triangulate.c  -  Don Cross  -  2021-06-22

    Example C program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astronomy.h"

char printBuffer[80];
static const char * UsageText = "\n"
"USAGE:  triangulate  lat1 lon1 elv1 az1 alt1  lat2 lon2 elv2 az2 alt2\n"
"\n"
"Calculate the best-fit location of a point as observed\n"
"from two different locations on or near the Earth's surface.\n"
"\n"
"lat1, lat2 = Geographic latitudes in degrees north of the equator.\n"
"lon1, lon2 = Geographic longitudes in degrees east of the prime meridian.\n"
"elv1, elv2 = Elevations above sea level in meters.\n"
"az1,  az2  = Azimuths toward observed object in degrees clockwise from north.\n"
"alt1, alt2 = Altitude angles toward observed object in degrees above horizon.\n"
"\n"
"This program extrapolates lines in the given directions from the two\n"
"geographic locations and finds the location in space where they\n"
"come closest to intersecting. It then prints out the coordinates\n"
"of that triangulation point, along with the error radius in meters.\n"
"\n";


static double DotProduct(astro_vector_t a, astro_vector_t b)
{
    if (a.status != ASTRO_SUCCESS || b.status != ASTRO_SUCCESS)
    {
        snprintf(printBuffer, sizeof printBuffer, "DotProduct(FATAL): a.status=%d, b.status=%d\n", a.status, b.status);
        Serial.println(printBuffer);
        exit(1);
    }

    return a.x*b.x + a.y*b.y + a.z*b.z;
}


static astro_vector_t AddScale(double sa, astro_vector_t va, double sb, astro_vector_t vb)
{
    astro_vector_t v;

    if (va.status != ASTRO_SUCCESS || vb.status != ASTRO_SUCCESS)
    {
        snprintf(printBuffer, sizeof printBuffer, "AddScale(FATAL): va.status=%d, vb.status=%d\n", va.status, vb.status);
        Serial.println(printBuffer);
        exit(1);
    }

    v.x = sa*va.x + sb*vb.x;
    v.y = sa*va.y + sb*vb.y;
    v.z = sa*va.z + sb*vb.z;
    v.t = va.t;
    v.status = ASTRO_SUCCESS;

    return v;
}


static double ScanArg(const char *name, const char *text)
{
    double value;
    if (sscanf(text, "%lf", &value) != 1 || !isfinite(value))
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Invalid value for %s: '%s'\n", name, text);
        Serial.println(printBuffer);
        exit(1);
    }
    return value;
}


static int DirectionVector(
    astro_time_t time,
    astro_observer_t observer,
    double altitude,
    double azimuth,
    astro_vector_t *equ_vector)
{
    astro_spherical_t hor;
    astro_vector_t hor_vector;
    astro_rotation_t rot;

    /* Convert horizontal angles to a horizontal unit vector. */

    hor.lat = altitude;
    hor.lon = azimuth;
    hor.dist = 1.0;
    hor.status = ASTRO_SUCCESS;

    hor_vector = Astronomy_VectorFromHorizon(hor, time, REFRACTION_NONE);
    if (hor_vector.status != ASTRO_SUCCESS)
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Astronomy_VectorFromHorizon returned %d\n", hor_vector.status);
        Serial.println(printBuffer);
        return 1;
    }

    /* Find the rotation matrix that converts horizontal vectors to equatorial vectors. */
    rot = Astronomy_Rotation_HOR_EQD(&time, observer);
    if (rot.status != ASTRO_SUCCESS)
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Astronomy_Rotation_HOR_EQD returned %d\n", rot.status);
        Serial.println(printBuffer);
        return 1;
    }

    /* Rotate the horizontal (HOR) vector to an equator-of-date (EQD) vector. */
    *equ_vector = Astronomy_RotateVector(rot, hor_vector);
    if (equ_vector->status != ASTRO_SUCCESS)
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Astronomy_RotateVector returned %d\n", equ_vector->status);
        Serial.println(printBuffer);
        return 1;
    }

    return 0;
}


static int Intersect(
    astro_vector_t pos1,
    astro_vector_t dir1,
    astro_vector_t pos2,
    astro_vector_t dir2)
{
    double E, F, G, denom, u, v, dist;
    astro_vector_t amb, a, b, c, miss;
    astro_observer_t obs;

    F = DotProduct(dir1, dir2);
    amb = AddScale(+1.0, pos1, -1.0, pos2);       /* amb = pos1 - pos2 */
    E = DotProduct(dir1, amb);
    G = DotProduct(dir2, amb);
    denom = 1.0 - F*F;
    if (denom == 0.0)
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Cannot solve because directions are parallel.\n");
        Serial.println(printBuffer);
        return 1;
    }

    u = (F*G - E) / denom;
    v = G + F*u;
    if (u < 0.0 || v < 0.0)
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Lines of sight do not converge.\n");
        Serial.println(printBuffer);
        return 1;
    }

    a = AddScale(1.0, pos1, u, dir1);   /*  a = pos1 + u*dir1   */
    b = AddScale(1.0, pos2, v, dir2);   /*  b = pos2 + v*dir2   */
    c = AddScale(0.5, a, 0.5, b);       /*  c = (a+b)/2         */
    miss = AddScale(1.0, a, -1.0, b);   /*  miss = a-b          */

    dist = (KM_PER_AU * 1000.0 / 2.0) * Astronomy_VectorLength(miss);   /* error radius in meters */
    obs = Astronomy_VectorObserver(&c, EQUATOR_OF_DATE);

    snprintf(printBuffer, sizeof printBuffer,"Solution: lat = %0.6lf, lon = %0.6lf, elv = %0.3lf meters; error = %0.3lf meters\n",
        obs.latitude, obs.longitude, obs.height, dist);
    Serial.println(printBuffer);

    return 0;
}


void setup() {
  Serial.begin(9600);
}

void loop() {
    
    double lat1, lon1, elv1, az1, alt1;
    double lat2, lon2, elv2, az2, alt2;
    astro_time_t time;
    astro_observer_t obs1, obs2;
    astro_vector_t pos1, pos2;
    astro_vector_t dir1, dir2;

    /* Validate and parse command line arguments. */
    const char *argvs[] = {"10","10","10","10","10","10","10","10","10","10","10"};
    int args = sizeof(argvs) / sizeof(argvs[1]);
    if (args != 11)
    {
        char printBufferError[200];
        snprintf(printBufferError, sizeof printBufferError, "%s", UsageText);
        Serial.println(printBufferError);
    }

    lat1 = ScanArg("lat1", argvs[ 1]);
    lon1 = ScanArg("lon1", argvs[ 2]);
    elv1 = ScanArg("elv1", argvs[ 3]);
    az1  = ScanArg("az1",  argvs[ 4]);
    alt1 = ScanArg("alt1", argvs[ 5]);
    lat2 = ScanArg("lat2", argvs[ 6]);
    lon2 = ScanArg("lon2", argvs[ 7]);
    elv2 = ScanArg("elv2", argvs[ 8]);
    az2  = ScanArg("az2",  argvs[ 9]);
    alt2 = ScanArg("alt2", argvs[10]);

    obs1 = Astronomy_MakeObserver(lat1, lon1, elv1);
    obs2 = Astronomy_MakeObserver(lat2, lon2, elv2);

    /* Convert geographic coordinates into 3D vectors. */
    /* We can use an arbitrary time because we don't care about */
    /* the Earth's rotation with respect to extraterrestrial bodies. */
    time = Astronomy_TimeFromDays(0.0);

    /* Convert first observer's geographic coordinates to a geocentric vector. */
    pos1 = Astronomy_ObserverVector(&time, obs1, EQUATOR_OF_DATE);
    if (pos1.status != ASTRO_SUCCESS)
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Astronomy_ObserverVector returned error %d for first observer.\n", pos1.status);
        Serial.println(printBuffer);

    }

    /* Convert second observer's geographic coordinates to a geocentric vector. */
    pos2 = Astronomy_ObserverVector(&time, obs2, EQUATOR_OF_DATE);

    if (pos2.status != ASTRO_SUCCESS)
    {
        snprintf(printBuffer, sizeof printBuffer, "ERROR: Astronomy_ObserverVector returned error %d for second observer.\n", pos2.status);
        Serial.println(printBuffer);

    }
    /* Convert observers' angular lines of sight into equatorial unit vectors. */

    if (DirectionVector(time, obs1, alt1, az1, &dir1))
        Serial.println("DirectionVector...A..");
    if (DirectionVector(time, obs2, alt2, az2, &dir2))
        Serial.println("DirectionVector...B..");

    /* Solve for the target point. */
    Intersect(pos1, dir1, pos2, dir2);
    delay(2000);
}
