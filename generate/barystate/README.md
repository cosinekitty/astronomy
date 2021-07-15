### Data for testing barycentric state calculation.

This [JPL Horizons](https://ssd.jpl.nasa.gov/horizons.cgi) data is for testing calculation of the barycentric state
(position and velocity vectors) using the Astronomy Engine `BaryState` function.
There are test files for all the supported bodies.
The files use varying start time, stop time, and time step values,
but they are otherwise configured the same.

As an example, here is the batch data for generating `Sun.txt`:

```
!$$SOF
COMMAND= '10'
CENTER= '500@0'
MAKE_EPHEM= 'YES'
TABLE_TYPE= 'VECTORS'
START_TIME= '1900-01-01'
STOP_TIME= '2100-01-01'
STEP_SIZE= '20 d'
OUT_UNITS= 'AU-D'
REF_PLANE= 'FRAME'
REF_SYSTEM= 'J2000'
VECT_CORR= 'NONE'
VEC_LABELS= 'YES'
VEC_DELTA_T= 'NO'
CSV_FORMAT= 'NO'
OBJ_DATA= 'YES'
VEC_TABLE= '2'
!$$EOF
```
