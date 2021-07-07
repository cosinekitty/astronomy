/*
    equator_of_date.js  -  by Don Cross  -  2021-07-06

    Example Node.js program for Astronomy Engine:
    https://github.com/cosinekitty/astronomy
*/

const UsageText = `
USAGE:  node equator_of_date.js ra dec [yyyy-mm-ddThh:mm:ssZ]

Converts J2000 equatorial coordinates to
equator-of-date coordinates.

ra  = J2000 Right ascension in sidereal hours (0 .. 24).
dec = J2000 Declination in degrees (-90 .. +90).
yyyy-mm-ddThh:mm:ssZ = Optional date and time in UTC.
(If omitted, the current date and time are used.)

This program prints out the right ascension and declination
of the same point in the sky, but expressed in the Earth's
equator at the given date and time.
`;

const Astronomy = require('./astronomy.js');

function ParseNumber(name, text, minValue, maxValue) {
    const x = Number(text);
    if (!Number.isFinite(x) || (x < minValue) || (x > maxValue)) {
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

function Demo() {
    if (process.argv.length < 4 || process.argv.length > 5) {
        console.log(UsageText);
        process.exit(1);
    } else {
        // Parse the command line arguments.
        const ra = ParseNumber("RA", process.argv[2], 0, 24);
        const dec = ParseNumber("DEC", process.argv[3], -90, +90);
        const date = (process.argv.length > 4) ? ParseDate(process.argv[4]) : new Date();
        const time = Astronomy.MakeTime(date);
        console.log(`time = ${time}`);

        // Create a rotation matrix that converts J2000 equatorial (EQJ)
        // orientation to equator-of-date (EQD) orientation.
        const rot = Astronomy.Rotation_EQJ_EQD(time);

        // Convert the spherical angular EQJ coordinates to a unit vector.
        // Multiply ra by 15 to convert sidereal hours to degrees
        const eqj_sphere = new Astronomy.Spherical(dec, 15*ra, 1);
        const eqj_vec = Astronomy.VectorFromSphere(eqj_sphere, time);

        // Use the rotation matrix to re-orient the EQJ vector to a EQD vector.
        const eqd_vec = Astronomy.RotateVector(rot, eqj_vec);

        // Convert the EQJ vector back to spherical angular coordinates.
        const eqd_sphere = Astronomy.SphereFromVector(eqd_vec);

        // Print out the converted angular coordinates.
        const eqd_ra = eqd_sphere.lon / 15;     // convert degrees to sidereal hours
        const eqd_dec = eqd_sphere.lat;
        console.log(`Equator-of-date coordinates: RA=${eqd_ra.toFixed(4)}, DEC=${eqd_dec.toFixed(4)}`);

        // Success!
        process.exit(0);
    }
}

Demo();
