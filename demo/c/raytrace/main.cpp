/*
    main.cpp  -  Entry point for Jupiter raytracer.
    by Don Cross <cosinekitty.com>
    https://github.com/cosinekitty/astronomy
*/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include "algebra.h"
#include "imager.h"
#include "astro_demo_common.h"

bool Verbose = false;

static const char UsageText[] =
"\n"
"USAGE:\n"
"\n"
"raytrace outfile.png width height planet yyyy-mm-ddThh:mm:ssZ [-f]\n"
"\n"
"where\n"
"    outfile.png = name of output PNG image\n"
"    width  = number of horizontal pixels in the output image\n"
"    height = number of vertical pixels in the output image\n"
"    planet = Jupiter\n"
"\n"
"options:\n"
"    -f       =  flip the image (match inverted telescope view)\n"
"    -s<ang>  =  spin the image by the specified angle in degrees\n"
"    -s       =  auto-spin the image to make the planet's north pole upward\n"
"    -z<fac>  =  zoom in by the given multiplication factor\n"
"\n";


const double AUTO_SPIN = 999.0;

const astro_time_t dummyTime = Astronomy_TimeFromDays(0.0);

astro_vector_t MakeAstroVector(double x, double y, double z)
{
    astro_vector_t a;

    a.x = x;
    a.y = y;
    a.z = z;
    a.t = dummyTime;
    a.status = ASTRO_SUCCESS;

    return a;
}


class RotationMatrixAimer : public Imager::Aimer
{
private:
    astro_rotation_t rotation;
    int flip;
    double xspin, yspin;

public:
    RotationMatrixAimer(astro_rotation_t _rotation, int _flip, double _spinAngleDegrees)
        : rotation(_rotation)
        , flip(_flip)
        , xspin(cos(_spinAngleDegrees * DEG2RAD))
        , yspin(sin(_spinAngleDegrees * DEG2RAD))
    {
    }

    virtual Imager::Vector Aim(const Imager::Vector& raw) const
    {
        astro_vector_t v;
        double x = flip ? -raw.x : raw.x;
        double y = raw.y;
        v.x = xspin*x - yspin*y;
        v.y = yspin*x + xspin*y;
        v.z = raw.z;
        v.t = dummyTime;
        v.status = ASTRO_SUCCESS;
        astro_vector_t rv = Astronomy_RotateVector(rotation, v);
        return Imager::Vector(rv.x, rv.y, rv.z);
    }
};


int JupiterImage(const char *filename, int width, int height, astro_time_t time, int flip, double spin, double zoom)
{
    using namespace Imager;

    // Calculate the geocentric position of Jupiter, corrected for light travel time.
    astro_vector_t jupiter = Astronomy_GeoVector(BODY_JUPITER, time, ABERRATION);
    if (jupiter.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Error %d calculating Jupiter geocentric position\n", jupiter.status);
        return 1;
    }

    // Calculate the geocentric position of the Sun, as our light source.
    astro_vector_t sun = Astronomy_GeoVector(BODY_SUN, time, ABERRATION);
    if (sun.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Error %d calculating Sun geocentric position\n", sun.status);
        return 1;
    }

    // Calculate the time light left the Jupiter system to be seen on Earth.
    double light_travel_time = Astronomy_VectorLength(jupiter) / C_AUDAY;
    astro_time_t depart = Astronomy_AddDays(time, -light_travel_time);

    // Calculate the orientation of Jupiter's rotation axis.
    astro_axis_t axis = Astronomy_RotationAxis(BODY_JUPITER, depart);
    if (axis.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Error %d calculating Jupiter's rotation axis.\n", axis.status);
        return 1;
    }

    // Calculate the position of Jupiter's moons at the backdated time.
    astro_jupiter_moons_t jm = Astronomy_JupiterMoons(depart);

    Scene scene(Color(0.0, 0.0, 0.0));

    const double km_scale = 10000.0;   // makes Jupiter's radius about 7 units.
    const double au_scale = km_scale / KM_PER_AU;

    const double equ_radius = JUPITER_EQUATORIAL_RADIUS_KM / km_scale;
    const double pol_radius = JUPITER_POLAR_RADIUS_KM / km_scale;
    Spheroid *planet = new Spheroid(equ_radius, equ_radius, pol_radius);
    scene.AddSolidObject(planet);
    planet->SetFullMatte(Color(1.0, 0.95, 0.85));
    planet->Move(Vector(jupiter.x, jupiter.y, jupiter.z) / au_scale);
    planet->SetTag("Jupiter");

    // Reorient Jupiter's rotation axis to match the calculated orientation.
    planet->RotateY(90.0 - axis.dec);       // ?? not sure about direction
    planet->RotateZ(15.0 * axis.ra);        // ?? not sure about direction

    // Add Jupiter's moons to the scene.

    // Io
    Vector io_center(
        jm.moon[JM_IO].x + jupiter.x,
        jm.moon[JM_IO].y + jupiter.y,
        jm.moon[JM_IO].z + jupiter.z
    );
    Sphere *io = new Sphere(io_center/au_scale, IO_RADIUS_KM/km_scale);
    scene.AddSolidObject(io);
    io->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    io->SetTag("Io");

    // Europa
    Vector eu_center(
        jm.moon[JM_EUROPA].x + jupiter.x,
        jm.moon[JM_EUROPA].y + jupiter.y,
        jm.moon[JM_EUROPA].z + jupiter.z
    );
    Sphere *europa = new Sphere(eu_center/au_scale, EUROPA_RADIUS_KM/km_scale);
    scene.AddSolidObject(europa);
    europa->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    europa->SetTag("Europa");

    // Ganymede
    Vector gan_center(
        jm.moon[JM_GANYMEDE].x + jupiter.x,
        jm.moon[JM_GANYMEDE].y + jupiter.y,
        jm.moon[JM_GANYMEDE].z + jupiter.z
    );
    Sphere *ganymede = new Sphere(gan_center/au_scale, GANYMEDE_RADIUS_KM/km_scale);
    scene.AddSolidObject(ganymede);
    ganymede->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???

    // Callisto
    Vector cal_center(
        jm.moon[JM_CALLISTO].x + jupiter.x,
        jm.moon[JM_CALLISTO].y + jupiter.y,
        jm.moon[JM_CALLISTO].z + jupiter.z
    );
    Sphere *callisto = new Sphere(cal_center/au_scale, CALLISTO_RADIUS_KM/km_scale);
    scene.AddSolidObject(callisto);
    callisto->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    callisto->SetTag("Callisto");

    // Add the Sun as the point light source.
    scene.AddLightSource(LightSource(Vector(sun.x, sun.y, sun.z) / au_scale, Color(1.0, 1.0, 1.0)));

    // Aim the camera at the center of Jupiter.
    // Start with an identity matrix, which leaves the camera pointing in the -z direction,
    // i.e. <0, 0, -1>.
    astro_rotation_t rotation = Astronomy_IdentityMatrix();

    // Convert Jupiter's rectangular coordinates to angular spherical coordinates.
    astro_spherical_t sph = Astronomy_SphereFromVector(jupiter);

    // The camera starts aimed at the direction <0, 0, -1>,
    // i.e. pointing away from the z-axis.
    // Perform a series of rotations that aims the camera toward
    // the center of Jupiter, but keeping the left/right direction
    // in the image parallel to the equatorial reference plane (the x-y plane).
    //
    // Start by rotating 90 degrees around the x-axis, to bring the camera center
    // to point in the y-direction. This leaves the image oriented with
    // the x-y plane horizontal, and the x-axis to the right.
    // Add the extra angle needed to bring the camera to Jupiter's declination.
    rotation = Astronomy_Pivot(rotation, 0, sph.lat + 90.0);

    // Now rotate around the z-axis to bring the camera to the
    // same right ascension as Jupiter. Subtract 90 degrees because
    // we are aimed at the y-axis, not the x-axis.
    rotation = Astronomy_Pivot(rotation, 2, sph.lon - 90.0);

    Vector planetUnitVector(jupiter.x/sph.dist, jupiter.y/sph.dist, jupiter.z/sph.dist);

    // Check for the sentinel value that indicates the caller
    // wants us to calculate the spin angle that shows the planet's
    // north pole in the upward-facing direction.
    if (spin == AUTO_SPIN)
    {
        // Earth-equatorial north is aligned with the camera's vertical direction.
        // We want to spin the camera around its aim point so that the planet's
        // north pole is aligned with the vertical direction instead.
        // So we calculate two angles:
        // (1) The angle of the Earth's north pole projected onto the camera's focal plane.
        // (2) The angle of the planet's north pole projected onto the camera's focal plane.
        // Subtract these two angles to obtain the desired spin angle.
        // We know angle (1) trivially as 90 degrees counterclockwise from the image's right.

        // The rotation matrix in 'rotation' converts from camera coordinates to EQJ.
        // Calculate the inverse matrix, which converts from EQJ to camera coordinates.
        astro_rotation_t inv = Astronomy_InverseRotation(rotation);

        // As a sanity check, reorient the Earth's north pole axis and verify
        // it is still pointing upward in the camera's focal plane.
        astro_vector_t earthNorthPole = MakeAstroVector(0, 0, 1);
        astro_vector_t earthCheck = Astronomy_RotateVector(inv, earthNorthPole);
        if (Verbose) printf("Earth north pole check: (%lf, %lf, %lf)\n", earthCheck.x, earthCheck.y, earthCheck.z);
        if (fabs(earthCheck.x) > 1.0e-6 || earthCheck.y <= 0.0)
        {
            fprintf(stderr, "FAIL Earth north pole check.\n");
            return 1;
        }

        // Now do the same thing with the planet's north pole axis.
        astro_vector_t axisShadow = Astronomy_RotateVector(inv, axis.north);
        if (Verbose) printf("Planet north pole shadow: (%lf, %lf, %lf)\n", axisShadow.x, axisShadow.y, axisShadow.z);
        spin = (RAD2DEG * atan2(axisShadow.y, axisShadow.x)) - 90.0;
        if (Verbose) printf("Auto-spin angle = %0.3lf degrees\n", spin);
    }

    RotationMatrixAimer aimer(rotation, flip, spin);
    // Verify that the aimer redirects the vector <0, 0, -1> directly
    // toward the center of Jupiter.
    Vector aimTest = aimer.Aim(Vector(0.0, 0.0, -1.0));
    if (Verbose)
    {
        printf("Aim Test: x = %12.8lf, y = %12.8lf, z = %12.8lf\n", aimTest.x, aimTest.y, aimTest.z);

        printf("Jupiter : x = %12.8lf, y = %12.8lf, z = %12.8lf\n",
            planetUnitVector.x,
            planetUnitVector.y,
            planetUnitVector.z);
    }

    const double diff = (aimTest - planetUnitVector).Magnitude();
    if (diff > 1.0e-15)
    {
        fprintf(stderr, "FAIL aim test: diff = %le\n", diff);
        return 1;
    }

    scene.SetAimer(&aimer);

    scene.SaveImage(filename, (size_t)width, (size_t)height, zoom, 4);

    return 0;
}


int main(int argc, const char *argv[])
{
    using namespace std;

    if (argc >= 6)
    {
        int flip = 0;
        double spin = 0.0;
        double zoom = 200.0;

        for (int i = 6; i < argc; ++i)
        {
            if (!strcmp(argv[i], "-f"))
            {
                flip = 1;
            }
            else if (!strcmp(argv[i], "-v"))
            {
                Verbose = true;
            }
            else if (argv[i][0] == '-' && argv[i][1] == 's')
            {
                if (argv[i][2] == '\0')
                {
                    // The "auto-spin" option. Set a sentinel value for 'spin' to indicate
                    // that the imager should calculate the ideal spin angle to make the planet's
                    // north pole appear toward the top of the generated image.
                    spin = AUTO_SPIN;
                }
                else if (1 != sscanf(&argv[i][2], "%lf", &spin) || !isfinite(spin) || spin < -360 || spin > +360)
                {
                    fprintf(stderr, "ERROR: invalid spin angle after '-s'\n");
                    return 1;
                }
            }
            else if (argv[i][0] == '-' && argv[i][1] == 'z')
            {
                if (1 != sscanf(&argv[i][2], "%lf", &zoom) || !isfinite(zoom) || zoom < 1 || zoom > 1.0e+6)
                {
                    fprintf(stderr, "ERROR: invalid zoom factor after '-z'\n");
                    return 1;
                }
            }
            else
            {
                fprintf(stderr, "ERROR: Unknown option: %s\n", argv[i]);
                return 1;
            }
        }

        const char *filename = argv[1];
        int width = atoi(argv[2]);
        if (width < 100 || width > 60000)
        {
            fprintf(stderr, "ERROR: Pixel width must be in the range 100..60000\n");
            return 1;
        }
        int height = atoi(argv[3]);
        if (height < 100 || height > 60000)
        {
            fprintf(stderr, "ERROR: Pixel height must be in the range 100..60000\n");
            return 1;
        }
        const char *planet = argv[4];
        astro_time_t time;
        if (Verbose) printf("time = [%s]\n", argv[5]);
        if (ParseTime(argv[5], &time))
            return 1;

        if (!strcmp(planet, "Jupiter"))
            return JupiterImage(filename, width, height, time, flip, spin, zoom);

        fprintf(stderr, "ERROR: Unknown planet '%s'\n", planet);
        return 1;
    }

    fprintf(stderr, "%s", UsageText);
    return 1;
}
