/*
    moon_north_south.js  -  by Don Cross  -  2023-08-22

    Calculates when the Moon reaches extreme declination
    (most north or south with respect to the plane of the Earth's equator)
    and when the Moon reaches the most extreme ecliptic latitude
    (most north or south with respect to the plane of the Earth's orbit around the Sun).
*/

const Astronomy = require('./astronomy.js');

function Ecliptic(time) {
    // Return the Moon's ecliptic latitude at the given time.
    const vec = Astronomy.GeoMoon(time);
    const ecl = Astronomy.Ecliptic(vec);
    return ecl.elat;
}

function Equatorial(time) {
    // Return the Moon's declination angle at the given time.
    // Start with the Moon's position vector in J2000 coordinates.
    const eqj = Astronomy.GeoMoon(time);
    // Find rotation matrix to convert J2000 coordinates to equator-of-date.
    const rot = Astronomy.Rotation_EQJ_EQD(time);
    // Transform coordinates into equator-of-date.
    const eqd = Astronomy.RotateVector(rot, eqj);
    // Convert to angular coordinates to find declination angle.
    const equ = Astronomy.EquatorFromVector(eqd);
    return equ.dec;
}

function Search(time1, direction, func) {
    // Create a callback function that reports the rate of change of the desired variable.
    function f(t) {
        const dt = 1 / 86400;       // one second, expressed in days
        const x1 = func(t.AddDays(-dt));
        const x2 = func(t.AddDays(+dt));
        return direction*(x2 - x1);
    }

    return Astronomy.Search(f, time1, time1.AddDays(30));
}

function Solve(time1, direction, func, comment) {
    const time = Search(time1, direction, func);
    const angle = func(time);
    console.log(`${time}  Moon next reaches ${comment} = ${angle.toFixed(7)}.`);
}


function ParseDate(text) {
    const d = new Date(text);
    if (!Number.isFinite(d.getTime())) {
        console.error(`ERROR: Not a valid date: "${text}"`);
        process.exit(1);
    }
    return d;
}

function Demo() {
    if (process.argv.length === 3) {
        const time1 = Astronomy.MakeTime(ParseDate(process.argv[2]));
        Solve(time1, -1, Ecliptic, 'maximum ecliptic latitude');
        Solve(time1, +1, Ecliptic, 'minimum ecliptic latitude');
        Solve(time1, -1, Equatorial, 'maximum declination');
        Solve(time1, +1, Equatorial, 'minimum declination');
        process.exit(0);
    } else {
        console.log('USAGE: node moon_north_south yyyy-mm-dd[Thh:mm:ssZ]');
        process.exit(1);
    }
}

Demo();

/*
    ---------------------------------------------------------------------------------------

    JPL Horizons data for maximum ecliptic latitude:

        time                  longitude     latitude
        2023-Sep-10 07:35     115.5157583   5.1662606
        2023-Sep-10 07:40     115.5572561   5.1662659
        2023-Sep-10 07:45     115.5987527   5.1662686   <== max
        2023-Sep-10 07:50     115.6402480   5.1662685
        2023-Sep-10 07:55     115.6817421   5.1662657

    This program says:
        2023-09-10T07:47:00.420Z  Moon next reaches maximum ecliptic latitude = 5.1662436.

    ---------------------------------------------------------------------------------------

    JPL Horizons data for minimum ecliptic latitude:

        time                  longitude     latitude
        2023-Aug-28 08:15     296.0963021  -5.1108070
        2023-Aug-28 08:20     296.1479665  -5.1108143
        2023-Aug-28 08:25     296.1996351  -5.1108175  <== min
        2023-Aug-28 08:30     296.2513079  -5.1108166
        2023-Aug-28 08:35     296.3029851  -5.1108115

    This program says:
        2023-08-28T08:26:14.233Z  Moon next reaches minimum ecliptic latitude = -5.1107593.

    ---------------------------------------------------------------------------------------

    JPL Horizons data for maximum declination:

        time                   RA(deg)   DEC(deg)
        2023-Sep-08 13:00      94.62655  28.17823
        2023-Sep-08 13:05      94.67459  28.17827
        2023-Sep-08 13:10      94.72262  28.17829
        2023-Sep-08 13:15      94.77064  28.17829   <== max
        2023-Sep-08 13:20      94.81867  28.17828
        2023-Sep-08 13:25      94.86669  28.17824

    This program says:
        2023-09-08T13:13:01.641Z  Moon next reaches maximum declination = 28.1783302.

    ---------------------------------------------------------------------------------------

    JPL Horizons data for maximum declination:

        time                   RA(deg)   DEC(deg)
        2023-Aug-26 20:05     274.66967 -28.10661
        2023-Aug-26 20:10     274.72558 -28.10665
        2023-Aug-26 20:15     274.78150 -28.10668
        2023-Aug-26 20:20     274.83743 -28.10668   <==  min
        2023-Aug-26 20:25     274.89336 -28.10667
        2023-Aug-26 20:30     274.94930 -28.10662

    This program says:
        2023-08-26T20:18:11.883Z  Moon next reaches minimum declination = -28.1066748.

    ---------------------------------------------------------------------------------------

    JPL Batch Data for calculating ecliptic latitudes:

!$$SOF
MAKE_EPHEM=YES
COMMAND=301
EPHEM_TYPE=OBSERVER
CENTER='500@399'
START_TIME='2023-08-28'
STOP_TIME='2023-08-29'
STEP_SIZE='5 MINUTES'
QUANTITIES='31'
REF_SYSTEM='ICRF'
CAL_FORMAT='CAL'
CAL_TYPE='M'
TIME_DIGITS='MINUTES'
ANG_FORMAT='HMS'
APPARENT='AIRLESS'
RANGE_UNITS='AU'
SUPPRESS_RANGE_RATE='NO'
SKIP_DAYLT='NO'
SOLAR_ELONG='0,180'
EXTRA_PREC='NO'
R_T_S_ONLY='NO'
CSV_FORMAT='NO'
OBJ_DATA='YES'

    ---------------------------------------------------------------------------------------

    JPL Batch data for calculating equatorial coordinates:

!$$SOF
MAKE_EPHEM=YES
COMMAND=301
EPHEM_TYPE=OBSERVER
CENTER='500@399'
START_TIME='2023-08-26'
STOP_TIME='2023-08-27'
STEP_SIZE='5 MINUTES'
QUANTITIES='2'
REF_SYSTEM='ICRF'
CAL_FORMAT='CAL'
CAL_TYPE='M'
TIME_DIGITS='MINUTES'
ANG_FORMAT='DEG'
APPARENT='AIRLESS'
RANGE_UNITS='AU'
SUPPRESS_RANGE_RATE='NO'
SKIP_DAYLT='NO'
SOLAR_ELONG='0,180'
EXTRA_PREC='NO'
R_T_S_ONLY='NO'
CSV_FORMAT='NO'
OBJ_DATA='YES'

*/
