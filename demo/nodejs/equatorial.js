/*
    equatorial.js  -  by Don Cross - 2021-03-27

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy

    Given an observer's location on the Earth, a
    date/time, and horizontal coordinates (azimuth, altitude)
    for that observer, this program works backwards
    to figure out the equatorial coordinates for that
    location in the sky. It provides two solutions:
    one that includes atmospheric refraction, another
    that ignores atmospheric refraction.

    To execute, run the command:
    node equatorial latitude longitude azimuth altitude [date]
*/

const Astronomy = require('./astronomy.js');

function ParseNumber(text, name) {
    const x = Number(text);
    if (!Number.isFinite(x)) {
        console.error(`ERROR: Not a valid numeric value for ${name}: "${text}"`);
        process.exit(1);
    }
    return x;
}

function ParseDate(text) {
    const d = new Date(text);
    if (!Number.isFinite(d.getTime())) {
        console.error(`ERROR: Not a valid date: "${text}"`);
        process.exit(1);
    }
    return d;
}


function Format(x, length, digits) {
    let s = x.toFixed(digits);
    while (s.length < length)
        s = ' ' + s;
    return s;
}


function Solve(refract, observer, time, azimuth, altitude) {
    // Convert the angular horizontal coordinates (azimuth, altitude)
    // to a horizontal vector (north, west, zenith).
    const hor_sphere = new Astronomy.Spherical(altitude, azimuth, 1);
    const refraction_option = refract ? 'normal' : null;
    const hor_vec = Astronomy.VectorFromHorizon(hor_sphere, time, refraction_option);

    // Make a rotation matrix for this observer and date/time that converts
    // horizontal coordinates (HOR) to equatorial coordinates in the J2000 epoch (EQJ).
    const rot_hor_eqj = Astronomy.Rotation_HOR_EQJ(time, observer);

    // Use the rotation matrix to convert the horizontal vector to an equatorial vector.
    const eqj_vec = Astronomy.RotateVector(rot_hor_eqj, hor_vec);

    // Convert the equatorial vector to equatorial angular coordinates (RA, DEC).
    const eqj = Astronomy.EquatorFromVector(eqj_vec);

    // Self-check the answers by converting back to horizontal coordinates,
    // using a different algorithm that has been tested to work.

    // First we need to convert J2000 equatorial (EQJ) to equator-of-date (EQD),
    // because the Horizon function expects EQD.
    const rot_eqj_eqd = Astronomy.Rotation_EQJ_EQD(time);
    const eqd_vec = Astronomy.RotateVector(rot_eqj_eqd, eqj_vec);
    const eqd = Astronomy.EquatorFromVector(eqd_vec);

    const check_hor = Astronomy.Horizon(time, observer, eqd.ra, eqd.dec, refraction_option);
    const alt_error = Math.abs(check_hor.altitude - altitude);
    const az_error = Math.abs(check_hor.azimuth - azimuth);

    let line = (refract ? '   yes' : '   no ');
    line += '    ' + Format(eqj.ra, 10, 4);
    line += ' ' + Format(eqj.dec, 10, 4);
    line += ' ' + Format(eqd.ra, 10, 4);
    line += ' ' + Format(eqd.dec, 10, 4);
    line += ' ' + Format(alt_error, 10, 6);
    line += ' ' + Format(az_error, 10, 6);
    console.log(line);
}


function Demo() {
    if (process.argv.length === 6 || process.argv.length === 7) {
        const latitude  = ParseNumber(process.argv[2]);
        const longitude = ParseNumber(process.argv[3]);
        const observer = new Astronomy.Observer(latitude, longitude, 0);
        const azimuth = ParseNumber(process.argv[4]);
        const altitude = ParseNumber(process.argv[5]);
        const time = Astronomy.MakeTime((process.argv.length === 7) ? ParseDate(process.argv[6]) : new Date());

        // Print a common header for both solutions.
        console.log('Refract?    J2000_RA  J2000_DEC  OFDATE_RA OFDATE_DEC  ALT_error   AZ_error');

        // Solve once ignoring atmospheric refraction.
        Solve(false, observer, time, azimuth, altitude);

        // Solve again considering atmospheric refraction.
        Solve(true, observer, time, azimuth, altitude);

        process.exit(0);
    } else {
        console.log('USAGE: node equatorial latitude longitude azimuth altitude [date]');
        process.exit(1);
    }
}

Demo();
