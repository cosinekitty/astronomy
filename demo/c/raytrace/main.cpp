/*
    main.cpp  -  Entry point for Jupiter raytracer.
    by Don Cross <cosinekitty.com>
    https://github.com/cosinekitty/astronomy
*/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include "algebra.h"
#include "planet.h"
#include "astro_demo_common.h"

static const char UsageText[] =
"\n"
"USAGE:\n"
"\n"
"raytrace outfile.png width planet yyyy-mm-ddThh:mm:ssZ\n"
"\n"
"where\n"
"    outfile.png = name of output PNG image\n"
"    width = number of pixels wide to make the image\n"
"    planet = Jupiter\n"
"\n";

void SpheroidTest()
{
    using namespace Imager;

    Scene scene(Color(0.0, 0.0, 0.0));

    Spheroid* spheroid = new Spheroid(4.0, 2.0, 1.0);
    spheroid->Move(0.0, 0.0, -50.0);
    spheroid->RotateX(-12.0);
    spheroid->RotateY(-60.0);

    scene.AddSolidObject(spheroid);
    scene.AddLightSource(LightSource(Vector(+35.0, +50.0, +20.0), Color(0.2, 1.0, 1.0)));
    scene.AddLightSource(LightSource(Vector(-47.0, -37.0, +12.0), Color(1.0, 0.2, 0.2)));

    const char *filename = "spheroid.png";
    scene.SaveImage(filename, 300, 300, 8.0, 2);
    std::cout << "Wrote " << filename << std::endl;
}

void SaturnTest()
{
    using namespace Imager;

    Scene scene(Color(0.0, 0.0, 0.0));
    Saturn* saturn = new Saturn();
    saturn->Move(0.0, 0.0, -100.0);
    saturn->RotateY(-15.0);
    saturn->RotateX(-83.0);
    scene.AddSolidObject(saturn);
    scene.AddLightSource(LightSource(Vector(+30.0, +26.0, +20.0), Color(1.0, 1.0, 1.0)));

    const char *filename = "saturn.png";
    scene.SaveImage(filename, 500, 250, 4.0, 4);
    std::cout << "Wrote " << filename << std::endl;
}

int JupiterImage(const char *filename, int width, astro_time_t time)
{
    return 0;
}

int main(int argc, const char *argv[])
{
    using namespace std;

    if (argc == 5)
    {
        const char *filename = argv[1];
        int width = atoi(argv[2]);
        if (width < 100 || width > 60000)
        {
            fprintf(stderr, "ERROR: Pixel width must be in the range 100..60000\n");
            return 1;
        }
        const char *planet = argv[3];
        astro_time_t time;
        if (ParseTime(argv[4], &time))
            return 1;

        if (!strcmp(planet, "Jupiter"))
            return JupiterImage(filename, width, time);

        fprintf(stderr, "ERROR: Unknown planet '%s'\n", planet);
        return 1;
    }

    fprintf(stderr, "%s", UsageText);
    return 1;
}
