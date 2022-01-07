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
"    -z<fac>  =  zoom in by the given multiplication factor\n"
"\n";

const double KM_SCALE = 1000.0;
const double AU_SCALE = KM_SCALE / KM_PER_AU;
const double AUTO_ZOOM = 0.0;
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


double BodyEquatorialRadiusKm(astro_body_t body)
{
    switch (body)
    {
    case BODY_MERCURY:
        return MERCURY_EQUATORIAL_RADIUS_KM;

    case BODY_VENUS:
        return VENUS_RADIUS_KM;

    case BODY_MOON:
        return MOON_EQUATORIAL_RADIUS_KM;

    case BODY_MARS:
        return MARS_EQUATORIAL_RADIUS_KM;

    case BODY_JUPITER:
        return JUPITER_EQUATORIAL_RADIUS_KM;

    case BODY_SATURN:
        return SATURN_EQUATORIAL_RADIUS_KM;

    case BODY_URANUS:
        return URANUS_EQUATORIAL_RADIUS_KM;

    case BODY_NEPTUNE:
        return NEPTUNE_EQUATORIAL_RADIUS_KM;

    case BODY_PLUTO:
        return PLUTO_RADIUS_KM;

    default:
        return -1.0;    // error code: unsupported body
    }
}


double BodyPolarRadiusKm(astro_body_t body)
{
    switch (body)
    {
    case BODY_MERCURY:
        return MERCURY_POLAR_RADIUS_KM;

    case BODY_VENUS:
        return VENUS_RADIUS_KM;

    case BODY_MOON:
        return MOON_POLAR_RADIUS_KM;

    case BODY_MARS:
        return MARS_POLAR_RADIUS_KM;

    case BODY_JUPITER:
        return JUPITER_POLAR_RADIUS_KM;

    case BODY_SATURN:
        return SATURN_POLAR_RADIUS_KM;

    case BODY_URANUS:
        return URANUS_POLAR_RADIUS_KM;

    case BODY_NEPTUNE:
        return NEPTUNE_POLAR_RADIUS_KM;

    case BODY_PLUTO:
        return PLUTO_RADIUS_KM;

    default:
        return -1.0;    // error code: unsupported body
    }
}


Imager::Color BodyColor(astro_body_t body)
{
    switch (body)
    {
    case BODY_MERCURY:
        return Imager::Color(168.0 / 255.0, 175.0 / 255.0, 168.0 / 255.0);

    case BODY_VENUS:
        return Imager::Color(1.0, 1.0, 1.0);

    case BODY_MOON:
        return Imager::Color(112.0 / 255.0, 108.0 / 255.0, 107.0 / 255.0);

    case BODY_MARS:
        return Imager::Color(255.0 / 255.0, 133.0 / 255.0, 94.0 / 255.0);

    case BODY_JUPITER:
        return Imager::Color(210.0 / 255.0, 198.0 / 255.0, 174.0 / 255.0);

    case BODY_SATURN:
        return Imager::Color(169.0 / 255.0, 149.0 / 255.0,  99.0 / 255.0, 0.45);

    case BODY_URANUS:
        return Imager::Color(142.0 / 255.0, 164.0 / 255.0, 177.0 / 255.0);

    case BODY_NEPTUNE:
        return Imager::Color(119.0 / 255.0, 154.0 / 255.0, 192.0 / 255.0);

    case BODY_PLUTO:
        return Imager::Color(216.0 / 255.0, 201.0 / 255.0, 180.0 / 255.0);

    default:
        return Imager::Color(1.0, 1.0, 1.0);
    }
}


static Imager::SolidObject* CreateSaturnRings()
{
    using namespace Imager;

    struct RingData
    {
        double  innerRadiusKm;
        double  outerRadiusKm;
        double  red;
        double  green;
        double  blue;
    };

    static const RingData ringData[] =
    {
        {  92154.0,  117733.0,  125.0,  118.0,   99.0 },
        { 122405.0,  133501.0,  118.0,  112.0,  100.0 },
        { 134085.0,  136888.0,  100.0,   94.0,   82.0 }
    };

    SolidObject* ringSystem = nullptr;

    static const size_t NUM_RINGS = sizeof(ringData) / sizeof(ringData[0]);
    for (size_t i=0; i < NUM_RINGS; ++i)
    {
        const Color color(
            ringData[i].red / 255.0,
            ringData[i].green / 255.0,
            ringData[i].blue / 255.0
        );

        ThinRing* ringSolid = new ThinRing(
            ringData[i].innerRadiusKm / KM_SCALE,
            ringData[i].outerRadiusKm / KM_SCALE
        );

        ringSolid->SetFullMatte(color);

        if (ringSystem != nullptr)
            ringSystem = new SetUnion(Vector(), ringSolid, ringSystem);
        else
            ringSystem = ringSolid;
    }

    return ringSystem;
}


static void AddJupiterMoons(
    Imager::Scene &scene,
    astro_vector_t geo_planet,
    astro_time_t depart)
{
    using namespace Imager;

    // Calculate the position of Jupiter's moons at the backdated time.
    astro_jupiter_moons_t jm = Astronomy_JupiterMoons(depart);

    // Add Jupiter's moons to the scene.

    // Io
    Vector io_center(
        jm.moon[JM_IO].x + geo_planet.x,
        jm.moon[JM_IO].y + geo_planet.y,
        jm.moon[JM_IO].z + geo_planet.z
    );
    Sphere *io = new Sphere(io_center/AU_SCALE, IO_RADIUS_KM/KM_SCALE);
    scene.AddSolidObject(io);
    io->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    io->SetTag("Io");

    // Europa
    Vector eu_center(
        jm.moon[JM_EUROPA].x + geo_planet.x,
        jm.moon[JM_EUROPA].y + geo_planet.y,
        jm.moon[JM_EUROPA].z + geo_planet.z
    );
    Sphere *europa = new Sphere(eu_center/AU_SCALE, EUROPA_RADIUS_KM/KM_SCALE);
    scene.AddSolidObject(europa);
    europa->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    europa->SetTag("Europa");

    // Ganymede
    Vector gan_center(
        jm.moon[JM_GANYMEDE].x + geo_planet.x,
        jm.moon[JM_GANYMEDE].y + geo_planet.y,
        jm.moon[JM_GANYMEDE].z + geo_planet.z
    );
    Sphere *ganymede = new Sphere(gan_center/AU_SCALE, GANYMEDE_RADIUS_KM/KM_SCALE);
    scene.AddSolidObject(ganymede);
    ganymede->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???

    // Callisto
    Vector cal_center(
        jm.moon[JM_CALLISTO].x + geo_planet.x,
        jm.moon[JM_CALLISTO].y + geo_planet.y,
        jm.moon[JM_CALLISTO].z + geo_planet.z
    );
    Sphere *callisto = new Sphere(cal_center/AU_SCALE, CALLISTO_RADIUS_KM/KM_SCALE);
    scene.AddSolidObject(callisto);
    callisto->SetFullMatte(Color(1.0, 1.0, 1.0));     // ??? Actual moon colors ???
    callisto->SetTag("Callisto");
}


int PlanetImage(
    astro_body_t body,
    const char *filename,
    int width,
    int height,
    astro_time_t time,
    int flip,
    double spin,
    double zoom)
{
    using namespace Imager;

    // Calculate the geocentric position of the planet, corrected for light travel time.
    astro_vector_t geo_planet = Astronomy_GeoVector(body, time, ABERRATION);
    if (geo_planet.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Error %d calculating planet geocentric position\n", geo_planet.status);
        return 1;
    }

    // Calculate the geocentric position of the Sun, as our light source.
    astro_vector_t sun = Astronomy_GeoVector(BODY_SUN, time, ABERRATION);
    if (sun.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Error %d calculating Sun geocentric position\n", sun.status);
        return 1;
    }

    // Calculate the time light left the planet to be seen on Earth.
    double planet_distance_au = Astronomy_VectorLength(geo_planet);
    double light_travel_time = planet_distance_au / C_AUDAY;
    astro_time_t depart = Astronomy_AddDays(time, -light_travel_time);

    // Calculate the orientation of the planet's rotation axis.
    astro_axis_t axis = Astronomy_RotationAxis(body, depart);
    if (axis.status != ASTRO_SUCCESS)
    {
        fprintf(stderr, "Error %d calculating planet's rotation axis.\n", axis.status);
        return 1;
    }

    Scene scene(Color(0.0, 0.0, 0.0));

    const double equ_radius = BodyEquatorialRadiusKm(body) / KM_SCALE;
    const double pol_radius = BodyPolarRadiusKm(body) / KM_SCALE;
    if (equ_radius < 0.0 || pol_radius < 0.0)
    {
        fprintf(stderr, "Error: cannot find radius data for requested body.\n");
        return 1;
    }
    SolidObject *planet = new Spheroid(equ_radius, equ_radius, pol_radius);
    planet->SetFullMatte(BodyColor(body));

    switch (body)
    {
    case BODY_SATURN:
        // Replace the spheroid with a union of the spheroid and its system of rings.
        planet = new SetUnion(Vector(), planet, CreateSaturnRings());
        break;

    case BODY_JUPITER:
        AddJupiterMoons(scene, geo_planet, depart);
        break;

    default:
        break;
    }

    scene.AddSolidObject(planet);
    planet->Move(Vector(geo_planet.x, geo_planet.y, geo_planet.z) / AU_SCALE);
    planet->SetTag(Astronomy_BodyName(body));

    // Reorient the planet's rotation axis to match the calculated orientation.
    planet->RotateY(90.0 - axis.dec);
    planet->RotateZ(15.0 * axis.ra);

    // Add the Sun as the point light source.
    scene.AddLightSource(LightSource(Vector(sun.x, sun.y, sun.z) / AU_SCALE, Color(1.0, 1.0, 1.0)));

    // Aim the camera at the planet's center.
    // Start with an identity matrix, which leaves the camera pointing in the -z direction,
    // i.e. <0, 0, -1>.
    astro_rotation_t rotation = Astronomy_IdentityMatrix();

    // Convert the planet's rectangular coordinates to angular spherical coordinates.
    astro_spherical_t sph = Astronomy_SphereFromVector(geo_planet);

    // The camera starts aimed at the direction <0, 0, -1>,
    // i.e. pointing away from the z-axis.
    // Perform a series of rotations that aims the camera toward
    // the center of the planet, but keeping the left/right direction
    // in the image parallel to the equatorial reference plane (the x-y plane).
    //
    // Start by rotating 90 degrees around the x-axis, to bring the camera center
    // to point in the y-direction. This leaves the image oriented with
    // the x-y plane horizontal, and the x-axis to the right.
    // Add the extra angle needed to bring the camera to the planet's declination.
    rotation = Astronomy_Pivot(rotation, 0, sph.lat + 90.0);

    // Now rotate around the z-axis to bring the camera to the
    // same right ascension as the planet. Subtract 90 degrees because
    // we are aimed at the y-axis, not the x-axis.
    rotation = Astronomy_Pivot(rotation, 2, sph.lon - 90.0);

    Vector planetUnitVector(geo_planet.x/sph.dist, geo_planet.y/sph.dist, geo_planet.z/sph.dist);

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
    // toward the center of the planet.
    Vector aimTest = aimer.Aim(Vector(0.0, 0.0, -1.0));
    if (Verbose)
    {
        printf("Aim Test: x = %12.8lf, y = %12.8lf, z = %12.8lf\n", aimTest.x, aimTest.y, aimTest.z);

        printf("Planet : x = %12.8lf, y = %12.8lf, z = %12.8lf\n",
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

    // Calculate the zoom factor that will fit the planet into
    // the smaller pixel dimension. To be more precise, calculate
    // the equatorial angular diameter with a 10% cushion for the
    // smaller of the two pixel dimensions (width or height).
    double diameter_radians = (2.0 * equ_radius) / (planet_distance_au / AU_SCALE);
    double factor = 0.9 / diameter_radians;
    if (zoom == AUTO_ZOOM)
        zoom = factor;
    else
        zoom *= factor;

    scene.SaveImage(filename, (size_t)width, (size_t)height, zoom, 4);

    return 0;
}


int main(int argc, const char *argv[])
{
    using namespace std;

    if (argc >= 6)
    {
        int flip = 0;
        double spin = AUTO_SPIN;
        double zoom = AUTO_ZOOM;

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
                if (1 != sscanf(&argv[i][2], "%lf", &spin) || !isfinite(spin) || spin < -360 || spin > +360)
                {
                    fprintf(stderr, "ERROR: invalid spin angle after '-s'\n");
                    return 1;
                }
            }
            else if (argv[i][0] == '-' && argv[i][1] == 'z')
            {
                if (1 != sscanf(&argv[i][2], "%lf", &zoom) || !isfinite(zoom) || zoom < 0.001 || zoom > 1000.0)
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

        astro_body_t body = Astronomy_BodyCode(planet);
        if (body == BODY_INVALID)
        {
            fprintf(stderr, "ERROR: Unknown body '%s'", planet);
            return 1;
        }

        return PlanetImage(body, filename, width, height, time, flip, spin, zoom);
    }

    fprintf(stderr, "%s", UsageText);
    return 1;
}
