/*
    earth.h  -  Don Cross  -  2019-03-23

    Special-purpose logic for calculating Earth-Moon Barycenter.
    Based on sun_eph() from NOVAS C 3.1 file solsys3.c.
*/
#ifndef __DDC_EARTH_H
#define __DDC_EARTH_H

extern const double EarthMoonMassRatio;

void CalcEarth(double jd, double pos[3]);
int NovasEarth(double jd, double pos[3]);

#endif /* __DDC_EARTH_H */
