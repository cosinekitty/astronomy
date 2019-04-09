/*
    novas_body.h  -  Don Cross  -  2019-03-28

    Defines body/target integer values as used by NOVAS C 3.1,
    along with a few extensions (Earth, Moon) not directly represented
    in the NOVAS ephemerides.
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

#endif /* __DDC_NOVAS_H */
