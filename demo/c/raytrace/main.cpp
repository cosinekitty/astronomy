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

static const char UsageText[] =
"\n"
"USAGE:\n"
"\n"
"raytrace outfile.png width planet yyyy-mm-ddThh:mm:ssZ [-f]\n"
"\n"
"where\n"
"    outfile.png = name of output PNG image\n"
"    width = number of pixels wide to make the image\n"
"    planet = Jupiter\n"
"\n"
"options:\n"
"    -f   =  flip the image (match inverted telescope view)\n"
"\n";


class RotationMatrixAimer : public Imager::Aimer
{
private:
    astro_rotation_t rotation;
    astro_time_t dummyTime;
    int flip;
    double xspin, yspin;

public:
    RotationMatrixAimer(astro_rotation_t _rotation, int _flip, double _spinAngleDegrees)
        : rotation(_rotation)
        , dummyTime(Astronomy_TimeFromDays(0.0))
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


int JupiterImage(const char *filename, int width, astro_time_t time, int flip, double spin, double zoom)
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
    const double io_radius = IO_RADIUS_KM / km_scale;
    Spheroid *io = new Spheroid(io_radius, io_radius, io_radius);
    scene.AddSolidObject(io);
    io->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    io->Move(Vector(
        jm.moon[JM_IO].x + jupiter.x,
        jm.moon[JM_IO].y + jupiter.y,
        jm.moon[JM_IO].z + jupiter.z) / au_scale
    );
    io->SetTag("Io");

    // Europa
    const double eu_radius = EUROPA_RADIUS_KM / km_scale;
    Spheroid *europa = new Spheroid(eu_radius, eu_radius, eu_radius);
    scene.AddSolidObject(europa);
    europa->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    europa->Move(Vector(
        jm.moon[JM_EUROPA].x + jupiter.x,
        jm.moon[JM_EUROPA].y + jupiter.y,
        jm.moon[JM_EUROPA].z + jupiter.z) / au_scale
    );
    europa->SetTag("Europa");

    // Ganymede
    const double gan_radius = GANYMEDE_RADIUS_KM / km_scale;
    Spheroid *ganymede = new Spheroid(gan_radius, gan_radius, gan_radius);
    scene.AddSolidObject(ganymede);
    ganymede->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    ganymede->Move(Vector(
        jm.moon[JM_GANYMEDE].x + jupiter.x,
        jm.moon[JM_GANYMEDE].y + jupiter.y,
        jm.moon[JM_GANYMEDE].z + jupiter.z) / au_scale
    );

    // Callisto
    const double cal_radius = CALLISTO_RADIUS_KM / km_scale;
    Spheroid *callisto = new Spheroid(cal_radius, cal_radius, cal_radius);
    scene.AddSolidObject(callisto);
    callisto->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    callisto->Move(Vector(
        jm.moon[JM_CALLISTO].x + jupiter.x,
        jm.moon[JM_CALLISTO].y + jupiter.y,
        jm.moon[JM_CALLISTO].z + jupiter.z) / au_scale
    );
    callisto->SetTag("Callisto");

    // Add the Sun as the point light source.
    scene.AddLightSource(LightSource(Vector(sun.x, sun.y, sun.z) / au_scale, Color(1.0, 1.0, 1.0)));

    // Aim the camera at the center of Jupiter.
    // Start with an identity matrix, which leaves the camera pointing in the -z direction,
    // i.e. <0, 0, -1>.
    astro_rotation_t rotation = Astronomy_IdentityMatrix();

    // Convert Jupiter's rectangular coordinates to angular spherical coordinates.
    astro_spherical_t sph = Astronomy_SphereFromVector(jupiter);

    // Rotate 90 degrees around the y-axis, plus the declination of Jupiter,
    // to bring the camera up to the declination level of Jupiter.
    rotation = Astronomy_Pivot(rotation, 1, -(90.0 + sph.lat));

    // Rotate around the z-axis to aim the camera at Jupiter's right ascension.
    rotation = Astronomy_Pivot(rotation, 2, sph.lon);

    RotationMatrixAimer aimer(rotation, flip, spin);
    // Verify that the aimer redirects the vector <0, 0, -1> directly
    // toward the center of Jupiter.
    Vector aimTest = aimer.Aim(Vector(0.0, 0.0, -1.0));
    printf("Aim Test: x = %12.8lf, y = %12.8lf, z = %12.8lf\n", aimTest.x, aimTest.y, aimTest.z);

    printf("Jupiter : x = %12.8lf, y = %12.8lf, z = %12.8lf\n",
        jupiter.x / sph.dist,
        jupiter.y / sph.dist,
        jupiter.z / sph.dist);

    scene.SetAimer(&aimer);

    scene.SaveImage(filename, (size_t)width, (size_t)width, zoom, 4);

    return 0;
}


int main(int argc, const char *argv[])
{
    using namespace std;

    if (argc >= 5)
    {
        int flip = 0;
        double spin = 0.0;
        double zoom = 200.0;

        for (int i = 5; i < argc; ++i)
        {
            if (!strcmp(argv[i], "-f"))
            {
                flip = 1;
            }
            else if (argv[i][0] == '-' && argv[i][1] == 's')
            {
                if (1 != sscanf(&argv[i][2], "%lf", &spin) || !isfinite(spin) || spin < -360 || spin > +360)
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
        const char *planet = argv[3];
        astro_time_t time;
        printf("time = [%s]\n", argv[4]);
        if (ParseTime(argv[4], &time))
            return 1;

        if (!strcmp(planet, "Jupiter"))
            return JupiterImage(filename, width, time, flip, spin, zoom);

        fprintf(stderr, "ERROR: Unknown planet '%s'\n", planet);
        return 1;
    }

    fprintf(stderr, "%s", UsageText);
    return 1;
}
