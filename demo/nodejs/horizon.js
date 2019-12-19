/*
    horizon.js  -  Don Cross  -  2019-12-14

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    This is a more advanced example. It shows how to use coordinate
    transforms and a binary search to find the two azimuths where the
    ecliptic intersects with an observer's horizon at a given date and time.

    node horizon.js latitude longitude [date]
*/
'use strict';

const Astronomy = require('../../source/js/astronomy.js');

const NUM_SAMPLES = 4;

function ECLIPLON(i) {
    return (360 * i) / NUM_SAMPLES;
}

function HorizontalCoords(ecliptic_longitude, time, rot_ecl_hor) {
    const eclip = Astronomy.MakeSpherical(
        0.0,                    /* being "on the ecliptic plane" means ecliptic latitude is zero. */
        ecliptic_longitude,
        1.0);                   /* any positive distance value will work fine. */

    /* Convert ecliptic angular coordinates to ecliptic vector. */
    const ecl_vec = Astronomy.VectorFromSphere(eclip, time);

    /* Use the rotation matrix to convert ecliptic vector to horizontal vector. */
    const hor_vec = Astronomy.RotateVector(rot_ecl_hor, ecl_vec);

    /* Find horizontal angular coordinates, correcting for atmospheric refraction. */
    return Astronomy.HorizonFromVector(hor_vec, 'normal');
}


function Search(time, rot_ecl_hor, e1, e2)
{
    const tolerance = 1.0e-6;        /* one-millionth of a degree is close enough! */

    /*
        Binary search: find the ecliptic longitude such that the horizontal altitude
        ascends through a zero value. The caller must pass e1, e2 such that the altitudes
        bound zero in ascending order.
    */

    for(;;)
    {
        const e3 = (e1 + e2) / 2.0;
        const h3 = HorizontalCoords(e3, time, rot_ecl_hor);

        if (Math.abs(e2-e1) < tolerance)
        {
            /* We have found the horizon crossing within tolerable limits. */
            return { ex:e3, h:h3 };
        }

        if (h3.lat < 0.0)
            e1 = e3;
        else
            e2 = e3;
    }
}


function FindEclipticCrossings(observer, time) {
    /*
        The ecliptic is a celestial circle that describes the mean plane of
        the Earth's orbit around the Sun. We use J2000 ecliptic coordinates,
        meaning the x-axis is defined to where the plane of the Earth's
        equator on January 1, 2000 at noon UTC intersects the ecliptic plane.
        The positive x-axis points toward the March equinox.
        Calculate a rotation matrix that converts J2000 ecliptic vectors
        to horizontal vectors for this observer and time.
    */
    const rot = Astronomy.Rotation_ECL_HOR(time, observer);

    /*
        Sample several points around the ecliptic.
        Remember the horizontal coordinates for each sample.
    */
    const hor = [];
    let i;

    for (i=0; i < NUM_SAMPLES; ++i) {
        hor.push(HorizontalCoords(ECLIPLON(i), time, rot));
    }

    for (i=0; i < NUM_SAMPLES; ++i) {
        const a1 = hor[i].lat;
        const a2 = hor[(i+1) % NUM_SAMPLES].lat;
        const e1 = ECLIPLON(i);
        const e2 = ECLIPLON(i+1);
        if (a1 * a2 <= 0) {
            let s, direction;

            if (a2 > a1)
                s = Search(time, rot, e1, e2);
            else
                s = Search(time, rot, e2, e1);

            if (s.h.lon > 0 && s.h.lon < 180)
                direction = 'ascends';
            else
                direction = 'descends';

            console.log(`Ecliptic longitude ${s.ex.toFixed(4)} ${direction} through horizon az ${s.h.lon.toFixed(4)}, alt ${s.h.lat.toExponential(4)}`);
        }
    }
}

function Demo() {
    if (process.argv.length === 4 || process.argv.length === 5) {
        const latitude = parseFloat(process.argv[2]);
        const longitude = parseFloat(process.argv[3]);
        const observer = Astronomy.MakeObserver(latitude, longitude, 0);
        const time = Astronomy.MakeTime((process.argv.length === 5) ? new Date(process.argv[4]) : new Date());
        FindEclipticCrossings(observer, time);
        process.exit(0);
    } else {
        console.log('USAGE: node horizon.js latitude longitude [date]');
        process.exit(1);
    }
}

Demo();
