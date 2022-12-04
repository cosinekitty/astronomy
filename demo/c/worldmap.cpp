/*
    worldmap.cpp  -  Don Cross  -  2022-04-09

    Example C++ program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Efficiently calculate a Mercator projection map
    of the Earth, showing areas where the Sun and Moon
    are visible.
*/

#include <stdio.h>
#include <math.h>
#include <vector>
#include "astro_demo_common.h"
#include "lodepng.h"

static const char UsageText[] =
"\n"
"USAGE:\n"
"\n"
"worldmap outfile.png [yyyy-mm-ddThh:mm:ssZ]\n"
"\n"
"Draws a Mercator projection of the Earth showing areas\n"
"where the Sun or Moon are visible. Writes the result to\n"
"the output specified PNG image file.\n"
"\n"
"If the observation time is not specified on the command line,\n"
"this program uses the computer's current date and time.\n"
"\n";

typedef unsigned char   byte;

struct Pixel
{
    byte    red;
    byte    green;
    byte    blue;
    byte    alpha;

    Pixel(): red(0x00), green(0x00), blue(0x00), alpha(0xff) {}
};

class Image
{
private:
    std::vector<Pixel> buffer;

public:
    const int width;
    const int height;

    Image(int pixelsWide, int pixelsHigh)
        : buffer(pixelsWide * pixelsHigh)
        , width(pixelsWide)
        , height(pixelsHigh)
        {}

    Pixel& pixel(int x, int y)
    {
        if (x < 0 || x >= width)
            throw "Invalid x coordinate";
        if (y < 0 || y >= height)
            throw "Invalid y coordinate";
        return buffer.at(x + y*width);
    }
};

const int PixelsPerDegree = 5;
const int MinLatitude = -80;
const int MaxLatitude = +80;
const int MinLongitude = -180;
const int MaxLongitude = +180;
const int PixelsWide = (MaxLongitude - MinLongitude) * PixelsPerDegree;
const int PixelsHigh = (MaxLatitude - MinLatitude) * PixelsPerDegree;

int Render(Image &image, astro_time_t time);


int main(int argc, const char *argv[])
{
    const char *outFileName;
    const char *timeString;
    astro_time_t time;

    // Parse the command line parameters.
    switch (argc)
    {
    case 2:
        outFileName = argv[1];
        time = Astronomy_CurrentTime();
        break;

    case 3:
        outFileName = argv[1];
        timeString = argv[2];
        if (!strcmp(timeString, "default"))
            timeString = "2022-04-09T16:05:35Z";    // hack to get around vscode unhelpful modifying of my command-line arguments when debugging!
        if (0 != ParseTime(timeString, &time))
            return 1;
        break;

    default:
        printf("%s", UsageText);
        return 1;
    }

    // Create a world map image in memory for the given time.
    Image image(PixelsWide, PixelsHigh);
    if (0 != Render(image, time))
        return 1;

    // Write the memory image to a PNG file.
    const unsigned char *firstBytePtr = &image.pixel(0, 0).red;
    unsigned error = lodepng::encode(outFileName, firstBytePtr, PixelsWide, PixelsHigh);
    if (error != 0)
    {
        fprintf(stderr, "FATAL: lodepng::encode returned %u\n", error);
        return 1;
    }

    return 0;
}


double VerticalComponent(
    const astro_rotation_t &rot,
    astro_vector_t &ovec,   // vector from center of Earth to the observer
    astro_vector_t &bvec)   // vector from center of Earth to the celestial body (Sun/Moon)
{
    // Subtract the vectors to get a vector from the observer to the Sun/Moon.
    // This is called a topocentric vector.
    astro_vector_t topo;
    topo.status = ASTRO_SUCCESS;
    topo.t = bvec.t;
    topo.x = bvec.x - ovec.x;
    topo.y = bvec.y - ovec.y;
    topo.z = bvec.z - ovec.z;

    // Rotate the topocentric vector into horizontal orientation.
    astro_vector_t hor = Astronomy_RotateVector(rot, topo);

    // Return the normalized vertical component of the object.
    // This is a number from -1 (straight down) to +1 (straight up).
    // If the vertical component is 0, it means the object is on the observer's horizon.
    return hor.z / sqrt(hor.x*hor.x + hor.y*hor.y + hor.z*hor.z);
}


void ColorPixel(Pixel &pixel, double vert, double red, double green, double blue)
{
    if (vert > 0.0)     // can the observer see this object at all?
    {
        double sumRed   = pixel.red   + (0xff * red   * vert);
        double sumGreen = pixel.green + (0xff * green * vert);
        double sumBlue  = pixel.blue  + (0xff * blue  * vert);

        if (sumRed > 0xff)      sumRed = 0xff;
        if (sumGreen > 0xff)    sumGreen = 0xff;
        if (sumBlue > 0xff)     sumBlue = 0xff;

        pixel.red   = (byte) sumRed;
        pixel.green = (byte) sumGreen;
        pixel.blue  = (byte) sumBlue;
    }
}


void PrintZenithPoint(const char *name, astro_vector_t geovec)
{
    // Convert the geocentric vector into an observer location.
    astro_observer_t observer = Astronomy_VectorObserver(&geovec, EQUATOR_OF_DATE);
    printf("%-4s is at zenith for: lat = %7.3lf, lon = %8.3lf\n", name, observer.latitude, observer.longitude);
    fflush(stdout);
}


int Render(Image &image, astro_time_t time)
{
    // To minimize the work for each pixel, we calculate the
    // geocentric positions of the Sun and Moon only once.
    // We convert them to equator-of-date coordinates (EQD)
    // to make it very easy to find an altitude angle for each pixel.

    // We do not need aberration correction for the Sun,
    // because 22 arcseconds is too small to notice on a map.
    // For the Moon, the aberration parameter has no effect.
    astro_vector_t geo_sun_eqj  = Astronomy_GeoVector(BODY_SUN,  time, NO_ABERRATION);
    astro_vector_t geo_moon_eqj = Astronomy_GeoVector(BODY_MOON, time, NO_ABERRATION);

    // Align the Sun/Moon vectors with the Earth's axis at the time of observation.
    astro_rotation_t rot = Astronomy_Rotation_EQJ_EQD(&time);
    astro_vector_t geo_sun_eqd  = Astronomy_RotateVector(rot, geo_sun_eqj);
    astro_vector_t geo_moon_eqd = Astronomy_RotateVector(rot, geo_moon_eqj);

    // Just for fun, find the geographic locations where the Sun/Moon
    // appear directly overhead (at the zenith = straight up).
    PrintZenithPoint("Sun",  geo_sun_eqd);
    PrintZenithPoint("Moon", geo_moon_eqd);

    astro_observer_t observer;
    observer.height = 0.0;

    for (int i = 0; i < PixelsWide; ++i)
    {
        observer.longitude = (i / (double)PixelsPerDegree) + MinLongitude;
        for (int j = 0; j < PixelsHigh; ++j)
        {
            observer.latitude = ((PixelsHigh - (j + 1)) / (double)PixelsPerDegree) + MinLatitude;

            // Calculate a vector from the center of the Earth to the observer's
            // geographic location that corresponds to this pixel.
            astro_vector_t ovec = Astronomy_ObserverVector(&time, observer, EQUATOR_OF_DATE);

            // Create a rotation matrix that converts equator-of-date vectors
            // into horizontal vectors for this pixel's location on the Earth.
            astro_rotation_t rot = Astronomy_Rotation_EQD_HOR(&time, observer);

            // Use the rotation matrix and observer vector to calculate
            // a vertical component value in the range -1.0 .. +1.0.
            // This number tells whether the Sun/Moon are above that observer's
            // horizon, and if so, how high it appears in the sky.
            double sunVert  = VerticalComponent(rot, ovec, geo_sun_eqd);
            double moonVert = VerticalComponent(rot, ovec, geo_moon_eqd);

            // Add yellow for sunlight intensity for this pixel.
            ColorPixel(image.pixel(i, j), sunVert,  1.0, 1.0, 0.0);

            // Add blue for moonlight intensity for this pixel.
            ColorPixel(image.pixel(i, j), moonVert, 0.0, 0.0, 0.5);
        }
    }

    return 0;
}
