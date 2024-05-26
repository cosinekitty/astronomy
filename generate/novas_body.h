/*
    novas_body.h

    Defines body/target integer values as used by NOVAS C 3.1,
    along with a few extensions (Earth, Moon) not directly represented
    in the NOVAS ephemerides.

    MIT License

    Copyright (c) 2019-2024 Don Cross <cosinekitty@gmail.com>

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
#ifndef __DDC_NOVAS_H
#define __DDC_NOVAS_H

#define BODY_INVALID    (-1)

/* The following are supported by the NOVAS state() function... */
#define BODY_MERCURY     0
#define BODY_VENUS       1
#define BODY_EMB         2      /* Earth/Moon Barycenter */
#define BODY_MARS        3
#define BODY_JUPITER     4
#define BODY_SATURN      5
#define BODY_URANUS      6
#define BODY_NEPTUNE     7
#define BODY_PLUTO       8
#define BODY_GM          9      /* Geocentric Moon */
#define BODY_SUN        10

/* The following are extensions supported by my own code... */
#define BODY_EARTH      11
#define BODY_MOON       12
#define BODY_SSB        13     /* Solar System Barycenter */

#endif /* __DDC_NOVAS_H */
