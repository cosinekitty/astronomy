### Data for testing EQJ/EQD conversion.

I generated this data using [JPL Horzions](https://ssd.jpl.nasa.gov/horizons.cgi)
to test Astronomy Engine's conversion between J2000 equatorial coordinates (EQJ)
and equator-of-date coordinates (EQD).

The first pair of RA/DEC are in ICRF, which is very close to Astronomy Engine EQJ.
The second pair of RA/DEC are apparent coordinates, which is close to
Astronomy Engine EQD, but also with aberration correction.

So my test algorithm will need to correct for aberration also.

Here are the batch settings to reproduce the data:

```
!$$SOF
COMMAND= '499'
CENTER= 'coord'
COORD_TYPE= 'GEODETIC'
SITE_COORD= '-77.06580,+38.92056,0'
MAKE_EPHEM= 'YES'
TABLE_TYPE= 'OBSERVER'
START_TIME= '2021-06-06'
STOP_TIME= '2031-06-06'
STEP_SIZE= '10 d'
CAL_FORMAT= 'JD'
TIME_DIGITS= 'MINUTES'
ANG_FORMAT= 'DEG'
OUT_UNITS= 'KM-S'
RANGE_UNITS= 'AU'
APPARENT= 'AIRLESS'
SUPPRESS_RANGE_RATE= 'NO'
SKIP_DAYLT= 'NO'
EXTRA_PREC= 'YES'
R_T_S_ONLY= 'NO'
REF_SYSTEM= 'J2000'
CSV_FORMAT= 'NO'
OBJ_DATA= 'YES'
QUANTITIES= '1,2'
!$$EOF
```
