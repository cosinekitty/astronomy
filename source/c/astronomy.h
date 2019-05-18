/*
    Astronomy library for C/C++.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef __ASTRONOMY_H
#define __ASTRONOMY_H

#ifdef __cplusplus
extern "C" {
#endif

/*---------- types ----------*/

typedef enum
{
    ASTRO_SUCCESS,
    ASTRO_INVALID_BODY,
    ASTRO_NO_CONVERGE
}
astro_status_t;

typedef struct
{
    double ut;
    double tt;
}
astro_time_t;

typedef struct
{
    astro_status_t status;
    double x;
    double y;
    double z;
    astro_time_t t;
}
astro_vector_t;

typedef enum
{
    BODY_MERCURY,
    BODY_VENUS,
    BODY_EARTH,
    BODY_MARS,
    BODY_JUPITER,
    BODY_SATURN,
    BODY_URANUS,
    BODY_NEPTUNE,
    BODY_PLUTO,
    BODY_SUN,
    BODY_MOON
}
astro_body_t;

#define MIN_BODY    BODY_MERCURY
#define MAX_BODY    BODY_MOON

typedef struct
{
    double latitude;
    double longitude;
    double elevation;
}
astro_observer_t;


/*---------- functions ----------*/

const char *Astronomy_BodyName(astro_body_t body);
astro_observer_t Astronomy_MakeObserver(double latitude, double longitude, double elevation);
astro_time_t Astronomy_MakeTime(int year, int month, int day, int hour, int minute, double second);
astro_time_t Astronomy_AddDays(astro_time_t time, double days);
astro_vector_t Astronomy_HelioVector(astro_body_t body, astro_time_t time);
astro_vector_t Astronomy_GeoVector(astro_body_t body, astro_time_t time);
astro_vector_t Astronomy_GeoMoon(astro_time_t time);
double Astronomy_VectorLength(astro_vector_t vector);

#ifdef __cplusplus
}
#endif

#endif  /* ifndef __ASTRONOMY_H */
